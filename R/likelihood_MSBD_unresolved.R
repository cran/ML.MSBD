likelihood_MSBD_unresolved = function(tree,shifts,gamma,lambdas,mus,
                                      lambda_rates = NULL,stepsize = NULL,
                                      uniform_weights = TRUE,p_lambda=0,p_mu=0,
                                      rho = 1, sigma = 0, rho_sampling = TRUE,
                                      unresolved = FALSE, lineage_counts = c(), tcut = NULL) {
  
  if(rho>1 || rho<0 || sigma>1 || sigma<0) stop("Invalid sampling proportions")
  
  lik = 0
  if(unresolved) {
    ntips = length(tree$tip.label)
    if(length(lineage_counts) != ntips) stop("The vector of number of collapsed species doesn't match with the number of tips")
    
    if(is.null(tcut)) stop("Times of cutoff not provided")
    if(length(tcut) ==1) tcut = rep(tcut, ntips)
    if(length(tcut) != ntips) stop("The vector of times of cutoff doesn't match with the number of tips")
    
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
        x = tree$edge[desc[1],2]
        if(x<=ntips) states[x]=shifts[i,3]
        else desc=c(desc,alldesc[[x - ntips]])
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
        #Warning: sigma set to 0 in the unresolved parts <=> assumption of no fossile data in the unresolved
        #Warning: gamma set to 0 in the unresolved parts <=> assumption of no rate shifts in the unresolved
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