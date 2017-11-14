ML_MSBD = function(tree,initial_values,
                   uniform_weights=TRUE,p_lambda=0,p_mu=0,
                   rho = 1, sigma=0, rho_sampling = TRUE,
                   unresolved = FALSE, lineage_counts = c(), tcut = 0,
                   optim_control = list(),attempt_remove=TRUE,max_nshifts=Inf,
                   stepsize=NULL, no_extinction=FALSE, fixed_gamma=NULL,
                   unique_lambda = FALSE, unique_mu = FALSE,
                   saved_state = NULL, save_path = NULL,
                   time_mode = c("real","3pos","tip","mid","root"),
                   fast_optim = FALSE) {
  
  if(time_mode %in% c("real","tip","mid","root")) time_positions = time_mode
  else if (time_mode == "3pos") time_positions = c("tip","mid","root")
  else stop("Invalid time positions mode, available are tip, mid, root, 3pos and real")
  
  if(rho>1 || rho<0 || sigma>1 || sigma<0) stop("Invalid sampling proportions")
  if(!rho_sampling && rho != 0) rho=0
  
  if(is.null(saved_state)) {
    initial_values = .get_initial_values_from_vector(initial_values,stepsize, no_extinction, fixed_gamma)
    
    #first test without shifts
    temp = .ML_optim(tree,c(),initial_values, c(),
                     uniform_weights,p_lambda,p_mu,
                     rho,sigma,rho_sampling,
                     unresolved,lineage_counts,tcut,
                     optim_control, "fixed",
                     stepsize,no_extinction,fixed_gamma,
                     unique_lambda,unique_mu, fast_optim)
    
    if(is.na(temp$l)) {
      stop("First optimization failure")
    }
    else {
      saved_state = list(initial_screen = TRUE, bestl = temp$l, pars = temp$p, initial_values = initial_values)
      if(!is.null(save_path)) save(saved_state, file=save_path)
    }
  }
  
  #initial conditions screen
  if(!is.null(saved_state$initial_screen)) {
    while(TRUE) {
      done = TRUE
      for(i in 1:length(saved_state$initial_values)) {
        #no screen on gamma, not useful
        if(names(saved_state$initial_values)[i]=="gamma") next
        
        iv = saved_state$initial_values
        iv[[i]] = iv[[i]]*10
        temp = .ML_optim(tree,c(),iv, c(),uniform_weights,p_lambda,p_mu,rho,sigma,rho_sampling,unresolved,lineage_counts,tcut,
                         optim_control, "fixed",stepsize,no_extinction,fixed_gamma,unique_lambda,unique_mu, fast_optim)
        if(!is.na(temp$l) && temp$l<saved_state$bestl) {
          done = FALSE
          saved_state$bestl = temp$l
          saved_state$initial_values = iv
          saved_state$pars = temp$p
          if(!is.null(save_path)) save(saved_state, file=save_path)
        }
        iv[[i]] = iv[[i]]/100
        temp = .ML_optim(tree,c(),iv, c(),uniform_weights,p_lambda,p_mu,rho,sigma,rho_sampling,unresolved,lineage_counts,tcut,
                         optim_control, "fixed",stepsize,no_extinction,fixed_gamma,unique_lambda,unique_mu, fast_optim)
        if(!is.na(temp$l) && temp$l<saved_state$bestl) {
          done = FALSE
          saved_state$bestl = temp$l
          saved_state$initial_values = iv
          saved_state$pars = temp$p
          if(!is.null(save_path)) save(saved_state, file=save_path)
        }
      }
      if(done) break
    } 
    
    newsaved_state = list()
    newsaved_state$params = saved_state$pars
    newsaved_state$likelihood = saved_state$bestl
    newsaved_state$best_models = c(saved_state$bestl)
    newsaved_state$edges = c()
    newsaved_state$times = c()
    newsaved_state$initial_values = saved_state$pars
    if(is.null(fixed_gamma)) newsaved_state$initial_values$gamma = saved_state$initial_values$gamma
    
    saved_state = newsaved_state
    if(!is.null(save_path)) save(saved_state, file=save_path)
  }
  
  while(TRUE) {
    #if no more edges free, stop
    if(length(saved_state$edges)==length(tree$edge) || length(saved_state$edges) == max_nshifts) break
    
    if(is.null(saved_state$partial)) {
      saved_state$partial = list(edge_min = 0, time_min = NULL, p_min = NULL, min_lik = Inf, tested_edges = c())
    }
    
    #test max estimates for all edges
    for(i in 1:length(tree$edge[,1])) {
      
      if(i %in% saved_state$partial$tested_edges) next #already tested edge
      if(is.element(i,saved_state$edges)) next #disallow multiple shifts on the same edge
      if(tree$edge.length[i] == 0) { #zero-length edge
        saved_state$partial$tested_edges = c(saved_state$partial$tested_edges,i)
        next
      }
      
      if(fast_optim) anc_state = .find_ancestral_state(tree, i, saved_state$edges)
      
      for(time_pos in time_positions) {
        temp = .ML_optim(tree,c(saved_state$edges,i),
                         saved_state$initial_values,saved_state$times,
                         uniform_weights,p_lambda,p_mu,
                         rho,sigma,rho_sampling,
                         unresolved,lineage_counts,tcut,
                         optim_control, time_pos,
                         stepsize,no_extinction,fixed_gamma,
                         unique_lambda,unique_mu, fast_optim, saved_state$params, anc_state)
        
        if(!is.na(temp$l) && temp$l<saved_state$partial$min_lik) {
          saved_state$partial$edge_min = i
          saved_state$partial$time_min = temp$t
          saved_state$partial$p_min = temp$p
          saved_state$partial$min_lik = temp$l
        }
      }
      
      saved_state$partial$tested_edges = c(saved_state$partial$tested_edges,i)
      if(!is.null(save_path)) save(saved_state, file=save_path)
    }
    
    saved_state$best_models = c(saved_state$best_models, saved_state$partial$min_lik)
    
    if(saved_state$likelihood>saved_state$partial$min_lik) {
      saved_state$edges = c(saved_state$edges,saved_state$partial$edge_min)
      saved_state$times = c(saved_state$times,saved_state$partial$time_min)
      saved_state$params = saved_state$partial$p_min
      saved_state$likelihood = saved_state$partial$min_lik
      saved_state$initial_values$gamma = saved_state$params$gamma #gamma needs to increase with the number of clusters
    }
    else {
      break
    }
    
    saved_state$partial = NULL
    if(!is.null(save_path)) save(saved_state, file=save_path)
  }
  
  #optional : try removing stuff from the shifts and test
  if(attempt_remove && length(saved_state$edges)>0) {
    removal = TRUE
    while(removal==TRUE) {
      removal = FALSE
      for(i in length(saved_state$edges):1) {
        edgetemp = saved_state$edges[-i]
        time_tmp = saved_state$times[-i]
        
        if(fast_optim) anc_state = .find_ancestral_state(tree, saved_state$edges[i], edgetemp)
        values_temp = .remove_state(saved_state$params, i+1, stepsize, no_extinction, unique_lambda, unique_mu)
        
        temp = .ML_optim(tree, edgetemp,
                         saved_state$initial_values, time_tmp,
                         uniform_weights,p_lambda,p_mu,
                         rho,sigma,rho_sampling,
                         unresolved,lineage_counts,tcut,
                         optim_control, "fixed",
                         stepsize,no_extinction,fixed_gamma,
                         unique_lambda,unique_mu,
                         fast_optim,values_temp, anc_state)
        
        if(!is.na(temp$l) && temp$l<saved_state$likelihood) {
          print("shift removed")
          saved_state$params = temp$p
          saved_state$likelihood = temp$l
          saved_state$best_models[length(edgetemp)+1] = temp$l
          saved_state$edges = edgetemp
          saved_state$times = time_tmp
          removal = TRUE
        }
        else {
          if(!is.na(temp$l) && saved_state$best_models[length(edgetemp)+1]> temp$l) saved_state$best_models[length(edgetemp)+1] = temp$l
        }
      }
      
      if(!is.null(save_path)) save(saved_state, file=save_path)
    }
  }
  #decombine parameter vector in exploitable form
  result = saved_state$params
  
  result$likelihood = saved_state$likelihood
  result$shifts.edge = c(0,saved_state$edges)
  result$shifts.time = c(0,saved_state$times)
  result$best_models = saved_state$best_models
  
  return(result)
}

.ML_optim = function(tree,edges,initial_values, times,
                     uniform_weights=TRUE,p_lambda=0,p_mu=0,
                     rho=1, sigma=0, rho_sampling = TRUE,
                     unresolved = FALSE, lineage_counts = c(), tcut = 0,
                     optim_control=list(), time_mode = c("real","tip","mid","root","fixed"),
                     stepsize=NULL, no_extinction=FALSE, fixed_gamma=NULL,
                     unique_lambda = FALSE, unique_mu = FALSE,
                     fast_optim = FALSE, current_values = NULL, anc_state=0) {
  
  n = length(edges)
  
  if(time_mode != "fixed") { #fixed means all times are given as input
    e = edges[n]
    depths = ape::node.depth.edgelength(tree)
    tor = max(depths)
    tmin = tor-depths[tree$edge[e,2]] #times are backward
    tmax = tor-depths[tree$edge[e,1]]
    
    #if applicable, fix time of last shift
    if(time_mode == "tip") times[n] = tmin*0.01+tmax*0.99
    if(time_mode == "mid") times[n] = (tmin+tmax)/2
    if(time_mode == "root") times[n] = tmin*0.99+tmax*0.01
  }
  
  if(no_extinction && !uniform_weights) p_mu = function(x) {
    cond = function(c) { if(c==0){1}else{0}}
    sapply(x,cond)
  }
  
  #create function of parameters
  auxfun = function(p) {
    shifts = NULL
    
    if(!fast_optim) values = .get_params_from_vector(p, n+1, stepsize, no_extinction, fixed_gamma, unique_lambda, unique_mu)
    else values = .get_params_from_vector_fast(p, n+1, current_values, anc_state,
                                               stepsize, no_extinction, fixed_gamma, unique_lambda, unique_mu)
    
    #test validity of parameters
    if(values$gamma<0) return(NA)
    if(sum(values$lambdas<0)>0) return(NA)
    if(!is.null(values$lambda_rates) && sum(values$lambda_rates<0)>0) return(NA)
    if(sum(values$mus<0)>0) return(NA)
    
    if(time_mode == "real") { #time of last shift is inferred
      times[n] = p[length(p)]
      if(times[n]<tmin || times[n]>tmax) return(NA)
    }
    
    for(j in seq(along=edges)) {
      e = edges[j]
      t = times[j]
      shifts = rbind(shifts,c(e,t,1+j))
    }
    
    res = likelihood_MSBD_unresolved(tree,shifts,values$gamma,values$lambdas,values$mus,
                                     values$lambda_rates,stepsize,
                                     uniform_weights,p_lambda,p_mu,
                                     rho,sigma,rho_sampling,
                                     unresolved,lineage_counts,tcut)
    return(res)
  }
  
  #initial parameters values
  if(!fast_optim) initp = .get_vector_from_initial_values(initial_values, n+1, stepsize, no_extinction, fixed_gamma, unique_lambda, unique_mu)
  else initp = .get_vector_from_initial_values_fast(initial_values, n+1, stepsize, no_extinction, fixed_gamma, unique_lambda, unique_mu)
  if(time_mode == "real") initp = c(initp, (tmin+tmax)/2)
  
  #optim on auxfun
  out=try(stats::optim(initp, auxfun,control=optim_control))
  
  if (class(out) == "try-error") {
    print("Optimization failed")
    result = list(l=NA,p=out)
  }
  else {
    while(out$convergence==1) {
      if(is.null(optim_control$maxit)) optim_control$maxit = 1000
      else optim_control$maxit = 2*optim_control$maxit
      out=stats::optim(out$par, auxfun,control=optim_control)
    }
    
    if(time_mode == "real") {
      times[n] = out$par[length(out$par)]
      out$par = out$par[-length(out$par)]
    }
    result = list(l=out$value, t=times[n])
    if(!fast_optim) result$p=.get_params_from_vector(out$par, n+1, stepsize, no_extinction, fixed_gamma, unique_lambda, unique_mu)
    else result$p=.get_params_from_vector_fast(out$par, n+1, current_values, anc_state,
                                               stepsize, no_extinction, fixed_gamma, unique_lambda, unique_mu)
  }
  return(result)
}