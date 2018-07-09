#' Likelihood calculation for unresolved trees
#' 
#' Calculates the negative log likelihood of a multi-states model given a tree. 
#' This function is designed to work with unresolved trees, where tips represent collapsed clades. This sampling scheme is not recommended for epidemiology datasets.
#' The MRCA times of collapsed clades and the number of collapsed lineages need to be provided for all tips. If neither is provided the function will default to random sampling. 
#' Extinct tips can be present outside of the unresolved parts, but not below the time(s) set for \code{tcut}.
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
#' @param lineage_counts Number of lineages collapsed on each tip. Should be set to 1 for extinct tips.
#' @param tcut Times of clade collapsing for each tip (i.e time of the MRCA of all collapsed lineages). Can be a single number or a vector of length the number of tips.
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
#' # Simulate a random phylogeny
#' set.seed(24)
#' tree <- ape::rcoal(10)
#' 
#' # Calculate the log likelihood under a constant birth-death model (i.e, no shifts) 
#' # with unresolved tips
#' likelihood_MSBD_unresolved(tree, shifts = c(), gamma = 0, lambdas = 10, mus = 1, 
#'                            lineage_counts = c(2,5,1,3,1,1,1,1,2,6), tcut = 0.05)
#' # Calculate the log likelihood under a multi-states model with 2 states and unresolved tips
#' likelihood_MSBD_unresolved(tree, shifts = matrix(c(2,0.7,2), nrow = 1), 
#'                            gamma = 0.05, lambdas = c(10, 5), mus = c(1, 1), 
#'                            lineage_counts = c(2,5,1,3,1,1,1,1,2,6), tcut = 0.05)
#' 
#' @export

likelihood_MSBD_unresolved = function(tree,shifts,gamma,lambdas,mus,
                                      lambda_rates = NULL,stepsize = NULL,
                                      uniform_weights = TRUE,p_lambda=0,p_mu=0,
                                      rho = 1, sigma = 0, rho_sampling = TRUE,
                                      lineage_counts = c(), tcut = NULL) {
  
  if(rho>1 || rho<0 || sigma>1 || sigma<0) stop("Invalid sampling proportions")
  
  unresolved = (length(lineage_counts) > 0 && !is.null(tcut))
  
  lik = 0
  if(unresolved) {
    ntips = length(tree$tip.label)
    if(length(lineage_counts) != ntips) stop("The vector of number of collapsed species doesn't match with the number of tips")
    
    if(is.null(tcut)) stop("Time(s) of clade collapsing need to be provided")
    if(length(tcut) ==1) tcut = rep(tcut, ntips)
    if(length(tcut) != ntips) stop("The vector of times of clade collapsing doesn't match with the number of tips")
    
    depths = ape::node.depth.edgelength(tree)
    tor = max(depths)
    times = tor-depths
    if(!is.null(tree$root.edge)) tor = tor + tree$root.edge
    
    ntips = length(tree$tip.label)
    alldesc = .get_desc(tree)
    states = rep(1,ntips)
    for(i in seq_along(shifts[,1])) {
      e = shifts[i,1]
      desc = tree$edge[e,2]
      while(length(desc)>0) {
        if(desc[1] <= ntips) states[desc[1]]=shifts[i,3]
        else desc=c(desc,alldesc[[desc[1] - ntips]])
        desc=desc[-1]
      }
    }
    
    tip_edges = which(tree$edge[,2] %in% 1:ntips)
    for(e in tip_edges) {
      tip = tree$edge[e,2]
      if(lineage_counts[tip] > 1) { #cutoff
        if(tree$edge.length[e] < tcut[tip]+times[tip]) stop(paste0("Invalid cut-off time for tip ",tree$tip.label[tip]))
        tree$edge.length[e] = tree$edge.length[e]-tcut[tip]+times[tip]
        lik = lik - log(.pn(lambdas[states[tip]],mus[states[tip]], lineage_counts[tip], tcut[tip]-times[tip], rho=1))
        #rho set to 1 in the unresolved parts <=> assumption of all species taken into account
        #Warning: sigma set to 0 in the unresolved parts <=> assumption of no fossil data there
        #Warning: gamma set to 0 in the unresolved parts <=> assumption of no rate shifts there
      }
      else { #fully resolved tip - need to account for sampling
        tcut[tip] = times[tip] #to be used in min(tcut)
        if(times[tip]<0.00001 && rho_sampling) lik = lik - log(rho)
        else lik = lik - log(sigma*mus[states[tip]])
      }
    }
    mint = min(tcut)
    lik = lik + likelihood_MSBD(tree,shifts,gamma,lambdas,mus,lambda_rates,stepsize,
                                uniform_weights,p_lambda,p_mu,
                                rho = rho, sigma = sigma, rho_sampling = TRUE, add_time = mint, unresolved = TRUE)
  }
  
  else lik = likelihood_MSBD(tree,shifts,gamma,lambdas,mus,lambda_rates,stepsize,
                             uniform_weights,p_lambda,p_mu,
                             rho = rho, sigma = sigma, rho_sampling = rho_sampling)
  return(lik) 
}