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