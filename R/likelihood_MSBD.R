#' Likelihood calculation for randomly sampled trees
#' 
#' Calculates the negative log likelihood of a multi-states model given a tree. 
#' This function is designed to work with constant extant and/or extinct sampling.
#' 
#' @param tree Phylogenetic tree (in ape format) to calculate the likelihood on.
#' @param shifts Matrix describing the positions (edges and times) of shifts. See 'Details'.
#' @param gamma Rate of state change.
#' @param lambdas Birth rates of all states.
#' @param mus Death rates of all states.
#' @param lambda_rates Rates of decay of birth rate for all states. To use exponential decay, \code{stepsize} should also be provided.
#' @param stepsize Size of the step to use for time discretization with exponential decay, default NULL. To use exponential decay, \code{lambda_rates} should also be provided.
#' @param uniform_weights Whether all states are weighted uniformly in shifts, default TRUE. If FALSE, the weights of states are calculated from the distributions \code{p_lambda} and \code{p_mu}. See 'Details'.
#' @param p_lambda Prior probability distribution on lambdas, used if \code{uniform_weights = FALSE}.
#' @param p_mu Prior probability distribution on mus, used if \code{uniform_weights = FALSE}.
#' @param rho Sampling proportion on extant tips, default 1.
#' @param sigma Sampling probability on extinct tips (tips are sampled upon extinction), default 0.
#' @param rho_sampling Whether the most recent tips should be considered extant tips, sampled with sampling proportion \code{rho}. If FALSE, all tips will be considered extinct tips, sampled with sampling probability \code{sigma}. Should be TRUE for most macroevolution datasets and FALSE for most epidemiology datasets.
#' @param unresolved Whether this tree is the backbone of an unresolved tree. This is an internal variable used in calculations for unresolved trees.
#' @param add_time The time between the most recent tip and the end of the process (>=0). This is an internal variable used in calculations for unresolved trees.
#' 
#' @return The value of the negative log likelihood of the model given the tree.
#' 
#' @details It is to be noted that all times are counted backwards, with the most recent tip positioned at 0. \cr\cr
#' The 'shifts' matrix is composed of 3 columns and a number of rows. Each row describes a shift: the first column is the index of the edge on which the shift happens, 
#' the second column is the time of the shift, and the third column is the index of the new state. For example the row vector (3,0.5,2) specifies a shift on edge number 3, at time 0.5, 
#' towards the state that has parameters \code{lambdas[2]}, \code{lambda_rates[2]} and \code{mus[2]}. \cr\cr
#' The weights w are used for calculating the transition rates q from each state i to j: \eqn{q_{i,j}=\gamma*w_{i,j}}{q(i,j)=\gamma*w(i,j)}. 
#' If \code{uniform_weights = TRUE}, \eqn{w_{i,j} = \frac{1}{N-1}}{w(i,j)=1/(N-1)} for all i,j, where N is the total number of states. 
#' If \code{uniform_weights = FALSE}, \eqn{w_{i,j} = \frac{p_\lambda(\lambda_j)p_\mu(\mu_j)}{sum_{k \ne i}p_\lambda(\lambda_k)p_\mu(\mu_k)}}{w(i,j)=p\lambda(\lambdaj)p\mu(\muj)/sum(p\lambda(\lambdak)p\mu(\muk)) for all k!=i}
#' where the distributions \eqn{p_\lambda}{p\lambda} and \eqn{p_\mu}{p\mu} are provided by the inputs \code{p_lambda} and \code{p_mu}.
#' 
#' @examples
#' # Input a phylogeny
#' tree <- ape::read.tree(text = "(((t4:0.7293960718,(t1:0.450904974,t3:0.09259337652)
#'         :0.04068535892):0.4769176776,t8:0.1541864066):0.7282000314,((t7:0.07264320855,
#'         (((t5:0.8231869878,t6:0.3492440532):0.2380232813,t10:0.2367582193):0.5329497182,
#'         t9:0.1016243151):0.5929288475):0.3003101915,t2:0.8320755605):0.2918686506);")
#' 
#' # Calculate the log likelihood under a constant birth-death model (i.e, no shifts) 
#' # with full extant & extinct sampling
#' likelihood_MSBD(tree, shifts = c(), gamma = 0, lambdas = 10, mus = 1, sigma = 1)
#' # Calculate the log likelihood under a multi-states model with 2 states 
#' # and full extant & extinct sampling
#' likelihood_MSBD(tree, shifts = matrix(c(2,1.8,2), nrow = 1), 
#'                 gamma = 0.05, lambdas = c(10, 6), mus = c(1, 0.5), sigma = 1)
#' # Calculate the log likelihood under a multi-states model with 2 states and exponential decay 
#' # with full extant & extinct sampling
#' likelihood_MSBD(tree, shifts = matrix(c(2,1.8,2), nrow = 1), 
#'                 gamma = 0.05, lambdas = c(10, 6), mus = c(1, 0.5), 
#'                 sigma = 1, stepsize = 0.01, lambda_rates = c(0.1, 0.1))
#' 
#' @export

likelihood_MSBD = function(tree,shifts,gamma,lambdas,mus,
                           lambda_rates = NULL,stepsize = NULL,
                           uniform_weights = TRUE,p_lambda=0,p_mu=0,
                           rho = 1, sigma = 0, rho_sampling = TRUE,
                           add_time = 0, unresolved = FALSE) {
  
  if(length(lambdas)!=length(mus)) {
    stop("Numbers of lambdas and mus are not consistent")
  }
  if(!is.null(lambda_rates) && length(lambdas)!=length(lambda_rates)) {
    stop("Numbers of lambdas and lambda rates are not consistent")
  }
  if(rho>1 || rho<0 || sigma>1 || sigma<0) stop("Invalid sampling proportions")
  if(!rho_sampling && rho != 0) rho=0
  
  if(!is.null(stepsize) && is.null(lambda_rates)) {
    warning("Stepsize provided but no lambda rates, will default to no decay")
    stepsize = NULL
  }
  if(is.null(stepsize) && !is.null(lambda_rates)) {
    warning("Lambda rates provided but no stepsize, will default to stepsize = min(tree$edge.length)")
    stepsize = min(tree$edge.length)
  }
  
  ntips = length(tree$tip.label)
  depths = ape::node.depth.edgelength(tree)
  totalnodes = ntips + tree$Nnode
  tor = max(depths)
  
  states = rep(1,totalnodes)
  inittimestates = rep(0,length(lambdas))
  
  desc = .get_desc(tree)
  org_nodes = totalnodes
  
  #process tree to add nodes for rate shifts (and update states accordingly)
  if(length(shifts[,2])>1) {
    shifts = shifts[order(shifts[,2],decreasing=TRUE),] #only works if shifts contains at least two lines
  }
  for(i in seq_along(shifts[,1])) {
    e = shifts[i,1] #edge concerned
    t = shifts[i,2] #time of the event
    j = shifts[i,3] #new state
    inittimestates[j]=t
    if(j>length(lambdas)) {
      stop("New state has no parameters associated")
    }
    edgee=tree$edge[e,]
    tree$edge[e,2]=totalnodes+1
    tree$edge.length[e] = tor-t-depths[edgee[1]]
    tree$edge = rbind(tree$edge,c(totalnodes+1,edgee[2]))
    tree$edge.length = c(tree$edge.length,t-tor+depths[edgee[2]])
    depths[totalnodes+1] = tor-t
    totalnodes=totalnodes+1
    tree$Nnode=tree$Nnode+1
    desc[[totalnodes-ntips]]=edgee[2]
    
    for(i in seq_along(shifts[,1])) {
      if(shifts[i,1]==e && shifts[i,2]<t) {
        shifts[i,1] = length(tree$edge.length)
      }
    }
    prop = totalnodes
    while(length(prop)>0) {
      states[prop[1]]=j
      if(prop[1]>ntips)  prop=c(prop,desc[[prop[1]-ntips]])
      prop=prop[-1]
    }
  }
  if(!is.null(tree$root.edge)) {
    ca = tree$edge[1,1]
    tree$edge=rbind(c(totalnodes+1,ca),tree$edge)
    tree$edge.length = c(tree$root.edge,tree$edge.length)
    states=c(states,1)
    tor = tor + tree$root.edge
    depths = c(depths+tree$root.edge,0)
  }
  
  inittimestates[1] = tor
  times = tor-depths
  
  if(!is.null(stepsize)) {
    pprecalc = c()
    anc_states = states
    for(i in seq_along(shifts[,1])) {
      anc_states[org_nodes+i] = states[tree$edge[which(tree$edge[,2]==org_nodes+i),1]]
    }
    for(i in 1:length(lambdas)) {
      nodes = which(anc_states==i)
      pprecalc[nodes] = .initial_calc_p_decay(lambdas[i],lambda_rates[i],mus[i],gamma,times[nodes]+add_time,inittimestates[i],rho,sigma,stepsize)
    }
  }
  
  #process tree for division in steps
  likelihood=0
  
  for(e in 1:length(tree$edge[,1])) {
    node = tree[[1]][e,1]
    edge_ti = times[node]
    node2 = tree[[1]][e,2]
    edge_te = times[node2]
    
    if(tree$edge.length[e]<0) stop(paste0("Tree edge ",e," has negative length"))
    s = states[node]
    
    if(is.null(stepsize)) { #no decay, constant parameters
      likelihood = likelihood -log(.fN_ratio(lambdas[s],mus[s],gamma,edge_ti+add_time,edge_te+add_time,rho,sigma))
    }
    else {
      t0 = inittimestates[s]
      z = lambda_rates[s]
      nsteps = ceiling(tree$edge.length[e]/stepsize)
      
      if(nsteps == 0) { #i.e ts==te
        likelihood = likelihood -log(.fN_ratio_decay(lambdas[s],z,mus[s],gamma,edge_ti+add_time,edge_te+add_time,t0,rho,sigma))
      }
      else {
        #pre-calculate all p values
        grid = seq(edge_ti,edge_te,length.out = nsteps+1)
        pval = .calc_p_decay(lambdas[s],z,mus[s],gamma,grid,t0,rho,sigma,pprecalc[node2])
        
        for(i in 1:nsteps) {
          step_ti = grid[i]
          step_te = grid[i+1]
          ival = pval[i+1]
          
          likelihood = likelihood -log(.fN_ratio_decay(lambdas[s],z,mus[s],gamma,step_ti+add_time,step_te+add_time,t0,rho,sigma,ival))
        }
      }
    }
    
    #what is the event at the end of the edge
    if(node2>ntips) {
      d = desc[[node2-ntips]]
      desc_count = length(d)
    }
    else desc_count=0
    
    if(desc_count==2) { #speciation event
      if(is.null(stepsize)) likelihood = likelihood - log(lambdas[s])
      else likelihood = likelihood - log(lambdas[s]*exp(lambda_rates[s]*(edge_te - inittimestates[s])))
    }
    else if(desc_count==1) { #rate shift event
      ns = states[node2]
      if(uniform_weights) {wj = 1/(length(lambdas)-1)}
      else {
        ws = p_lambda(lambdas)*p_mu(mus)
        ws[s]=0
        if(sum(ws)==0) {wj = 0}
        else {wj=ws[ns]/sum(ws)}
      }
      likelihood = likelihood - log(gamma*wj)
    }
    else if (desc_count==0) { #leaf
      if(unresolved) {} #unresolved tree, no sampling
      else {
        if(edge_te < 0.0001 && rho_sampling) { #Extant sampling
          likelihood = likelihood - log(rho)
        }
        else { #Extinct sampling
          likelihood = likelihood - log(sigma*mus[s])
        }
      }
    }
    else { #oops
      print("A node was found with more than two children")
    }
  }
  
  likelihood
}