#' Full Maximum Likelihood inference of birth and death rates together with their changes along a phylogeny under a multi-type birth-death model.
#' 
#' Infers a complete MSBD model from a phylogeny, including the most likely number of states, positions and times of state changes, and parameters associated with each state. 
#' Uses a greedy approach to add states and Maximum Likelihood inference for the other parameters.
#' 
#' @param tree Phylogenetic tree (in ape format) to calculate the likelihood on.
#' @param initial_values Initial values for the optimizer, to be provided as a vector in this order: gamma (optional), lambda, lambda decay rate (optional), mu (optional). See 'Details'.
#' @param uniform_weights Whether all states are weighted uniformly in shifts, default TRUE. If FALSE, the weights of states are calculated from the distributions \code{p_lambda} and \code{p_mu}. See 'Details'.
#' @param p_lambda Prior probability distribution on lambdas, used if \code{uniform_weights = FALSE}.
#' @param p_mu Prior probability distribution on mus, used if \code{uniform_weights = FALSE}.
#' @param rho Sampling proportion on extant tips, default 1.
#' @param sigma Sampling probability on extinct tips (tips are sampled upon extinction), default 0.
#' @param rho_sampling Whether the most recent tips should be considered extant tips, sampled with sampling proportion \code{rho}. If FALSE, all tips will be considered extinct tips, sampled with sampling probability \code{sigma}. Should be TRUE for most macroevolution datasets and FALSE for most epidemiology datasets.
#' @param lineage_counts For trees with clade collapsing. Number of lineages collapsed on each tip. Should be set to 1 for extinct tips. 
#' @param tcut For trees with clade collapsing. Times of clade collapsing for each tip (i.e time of the MRCA of all collapsed lineages). Can be a single number or a vector of length the number of tips.
#' @param stepsize Size of the step to use for time discretization with exponential decay, default NULL. To use exponential decay, an initial value for \code{lambda_rates} should also be provided.
#' @param no_extinction Whether to use the Yule process (\code{mu=0}) for all states, default FALSE. If TRUE no initial value for \code{mu} is needed.
#' @param fixed_gamma Value to which \code{gamma} should be fixed, default NULL. If provided no initial value for \code{gamma} is needed.
#' @param unique_lambda Whether to use the same value of \code{lambda} for all states, default FALSE. If TRUE and exponential decay is active all states will also share the same value for \code{lambda_rate}.
#' @param unique_mu Whether to use the same value of \code{mu} for all states, default FALSE.
#' 
#' @param optim_control Control list for the optimizer, corresponds to control input in optim function, see \code{?optim} for details.
#' @param attempt_remove Whether to attempt to remove shifts at the end of the inference, default TRUE. If FALSE, use a pure greedy algorithm.
#' @param max_nshifts Maximum number of shifts to test for, default \code{Inf}.
#' @param saved_state If provided, the inference will be restarted from this state.
#' @param save_path If provided, the progress of the inference will be saved to this path after each optimization step.
#' @param time_mode String controlling the time positions of inferred shifts. See 'Details'.
#' @param fast_optim Whether to use the faster mode of optimization, default FALSE. If TRUE only rates associated with the state currently being added to the tree and its ancestor will be optimized at each step, otherwise all rates are optimized.
#' @param parallel Whether the computation should be run in parallel, default FALSE. Will use a user-defined cluster if one is found, otherwise will define its own.
#' @param ncores Number of cores to use for a parallel computation.
#' 
#' @return Returns a list describing the most likely model found, with the following components:
#' \item{\code{likelihood}}{the negative log likelihood of the model}
#' \item{\code{shifts.edge}}{the indexes of the edges where shifts happen, 0 indicates the root state}
#' \item{\code{shifts.time}}{the time positions of shifts}
#' \item{\code{gamma}}{the rate of state change}
#' \item{\code{lambdas}}{the birth rates of all states}
#' \item{\code{lambda_rates}}{if exponential decay was activated, the rates of decay of birth rate for all states}
#' \item{\code{mus}}{the death rates of all states}
#' \item{\code{best_models}}{a vector containing the negative log likelihood of the best model found for each number of states tested (\code{best_models[i]} corresponds to i states, i.e i-1 shifts)}
#' All vectors are indexed in the same way, so that the state with parameters \code{lambdas[i]}, \code{lambda_rates[i]} and \code{mus[i]} starts on edge \code{shifts.edge[i]} at time \code{shifts.time[i]}.
#' 
#' @details It is to be noted that all times are counted backwards, with the most recent tip positioned at 0. \cr\cr
#' 
#' Five time modes are possible for the input \code{time_mode}. 
#' In \code{tip} mode, the shifts will be placed at 10\% of the length of the edge. 
#' In \code{mid} mode, the shifts will be placed at 50\% of the length of the edge. 
#' In \code{root} mode, the shifts will be placed at 90\% of the length of the edge. 
#' In \code{3pos} mode, the three "tip", "mid" and "root" positions will be tested.\cr\cr
#' 
#' The weights w are used for calculating the transition rates q from each state i to j: \eqn{q_{i,j}=\gamma*w_{i,j}}{q(i,j)=\gamma*w(i,j)}. 
#' If \code{uniform_weights = TRUE}, \eqn{w_{i,j} = \frac{1}{N-1}}{w(i,j)=1/(N-1)} for all i,j, where N is the total number of states. 
#' If \code{uniform_weights = FALSE}, \eqn{w_{i,j} = \frac{p_\lambda(\lambda_j)p_\mu(\mu_j)}{sum_{k \ne i}p_\lambda(\lambda_k)p_\mu(\mu_k)}}{w(i,j)=p\lambda(\lambdaj)p\mu(\muj)/sum(p\lambda(\lambdak)p\mu(\muk)) for all k!=i}
#' where the distributions \eqn{p_\lambda}{p\lambda} and \eqn{p_\mu}{p\mu} are provided by the inputs \code{p_lambda} and \code{p_mu}.\cr\cr
#' 
#' Initial values for the optimization need to be provided as a vector and contain the following elements (in order): 
#' an initial value for gamma, which is required unless \code{fixed_gamma} is provided, 
#' an initial value for lambda which is always required, 
#' an initial value for lambda decay rate, which is required if \code{stepsize} is provided, 
#' and an initial value for mu, which is required unless \code{no_extinction = TRUE}. 
#' An error will be raised if the number of initial values provided does not match the one expected from the rest of the settings, 
#' and the function will fail if the likelihood cannot be calculated at the initial values.
#' 
#' @examples
#' # Input a phylogeny
#' tree <- ape::read.tree(text = "(((t4:0.7293960718,(t1:0.450904974,t3:0.09259337652)
#'         :0.04068535892):0.4769176776,t8:0.1541864066):0.7282000314,((t7:0.07264320855,
#'         (((t5:0.8231869878,t6:0.3492440532):0.2380232813,t10:0.2367582193):0.5329497182,
#'         t9:0.1016243151):0.5929288475):0.3003101915,t2:0.8320755605):0.2918686506);")
#'
#' # Infer the most likely multi-states birth-death model 
#' # with full extant & extinct sampling
#' \dontrun{ML_MSBD(tree, initial_values = c(0.1, 10, 1), sigma = 1, time_mode = "mid") }
#' # Infer the most likely multi-states birth-death model with exponential decay
#' # and full extant & extinct sampling
#' \dontrun{ML_MSBD(tree, initial_values = c(0.1, 10, 0.5, 1), sigma = 1, 
#'                  stepsize = 0.1, time_mode = "mid")}
#' 
#' # Input a phylogeny with extant samples
#' tree2 <- ape::read.tree(text = "(t3:0.9703302342,((t4:0.1999577823,(t2:0.1287530271,
#'         (t7:0.08853561159,(t8:0.07930237712,t9:0.07930237712):0.009233234474):0.04021741549):
#'         0.07120475526):0.4269919425,(((t10:0.0191876225,t5:0.0191876225):0.04849906822,
#'         t6:0.06768669072):0.1672340445,t1:0.2349207353):0.3920289896):0.3433805094);")
#' 
#' # Infer the most likely multi-states Yule model with partial extant sampling
#' \dontrun{ML_MSBD(tree2, initial_values = c(0.1, 10), no_extinction = TRUE, 
#'                   rho = 0.5, time_mode = "mid")}
#' # Infer the most likely multi-states birth-death model with full extant sampling 
#' # and unresolved extant tips
#' \dontrun{ML_MSBD(tree2, initial_values = c(0.1, 10, 1), 
#'                   lineage_counts = c(2,5,1,3,1,1,1,1,2,6), tcut = 0.05, time_mode = "mid")}
#' 
#' @export

ML_MSBD = function(tree,initial_values,
                   uniform_weights=TRUE,p_lambda=0,p_mu=0,
                   rho = 1, sigma=0, rho_sampling = TRUE,
                   lineage_counts = c(), tcut = 0,
                   stepsize=NULL, no_extinction=FALSE, fixed_gamma=NULL,
                   unique_lambda = FALSE, unique_mu = FALSE,
                   optim_control = list(), attempt_remove=TRUE, max_nshifts=Inf,
                   saved_state = NULL, save_path = NULL,
                   time_mode = c("3pos","tip","mid","root"),
                   fast_optim = FALSE, 
                   parallel = FALSE, ncores = getOption('mc.cores', 2L)) {
  
  if(time_mode %in% c("tip","mid","root")) time_positions = time_mode
  else if (time_mode == "3pos") time_positions = c("tip","mid","root")
  else stop("Invalid time positions mode, available are tip, mid, root and 3pos")
  
  if(rho>1 || rho<0 || sigma>1 || sigma<0) stop("Invalid sampling proportions")
  if(!rho_sampling && rho != 0) rho=0
  
  if(length(lineage_counts) > 0 && !is.null(tcut)) {
    print("Clade collapsing detected")
    
    ntips = length(tree$tip.label)
    if(length(lineage_counts) != ntips) stop("The vector of number of collapsed species doesn't match with the number of tips")
    if(is.null(tcut)) stop("Time(s) of clade collapsing need to be provided")
    if(length(tcut) ==1) tcut = rep(tcut, ntips)
    if(length(tcut) != ntips) stop("The vector of times of clade collapsing doesn't match with the number of tips")
  }
  
  if(parallel) {
    if (! requireNamespace("doParallel", quietly = TRUE)) stop("Parallel computation requires the doParallel package.")
    if (! foreach::getDoParRegistered() ) {
      doParallel::registerDoParallel(ncores)
      message('Registered parallel computation with ', ncores, ' workers')
      on.exit(doParallel::stopImplicitCluster())
    } else {
      message('Using parallel computation with existing', foreach::getDoParName(), 
            ' with ', foreach::getDoParWorkers(), ' workers')
    }
    `%d%` <- foreach::`%dopar%`
  } else {
    message('Executing sequential computation')
    `%d%` <- foreach::`%do%`
  }
  
  ptm = proc.time()[3]
  
  if(is.null(saved_state)) {
    initial_values = .get_initial_values_from_vector(initial_values,stepsize, no_extinction, fixed_gamma)
    
    #first test without shifts
    temp = .ML_optim(tree,c(),initial_values, c(),
                     uniform_weights,p_lambda,p_mu,
                     rho,sigma,rho_sampling,
                     lineage_counts,tcut,
                     optim_control, "fixed",
                     stepsize,no_extinction,fixed_gamma,
                     unique_lambda,unique_mu, fast_optim)
    
    if(is.na(temp$l)) {
      stop("First optimization failure")
    }
    else {
      saved_state = list(initial_screen = TRUE, bestl = temp$l, pars = temp$p, initial_values = initial_values)
      if(!is.null(save_path) && proc.time()[3] - ptm > 600) {
        save(saved_state, file=save_path)
        ptm = proc.time()[3]
      }
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
        temp = .ML_optim(tree,c(),iv, c(),uniform_weights,p_lambda,p_mu,rho,sigma,rho_sampling,lineage_counts,tcut,
                         optim_control, "fixed",stepsize,no_extinction,fixed_gamma,unique_lambda,unique_mu, fast_optim)
        if(!is.na(temp$l) && temp$l<saved_state$bestl) {
          done = FALSE
          saved_state$bestl = temp$l
          saved_state$initial_values = iv
          saved_state$pars = temp$p
          if(!is.null(save_path) && proc.time()[3] - ptm > 600) {
            save(saved_state, file=save_path)
            ptm = proc.time()[3]
          }
        }
        iv[[i]] = iv[[i]]/100
        temp = .ML_optim(tree,c(),iv, c(),uniform_weights,p_lambda,p_mu,rho,sigma,rho_sampling,lineage_counts,tcut,
                         optim_control, "fixed",stepsize,no_extinction,fixed_gamma,unique_lambda,unique_mu, fast_optim)
        if(!is.na(temp$l) && temp$l<saved_state$bestl) {
          done = FALSE
          saved_state$bestl = temp$l
          saved_state$initial_values = iv
          saved_state$pars = temp$p
          if(!is.null(save_path) && proc.time()[3] - ptm > 600) {
            save(saved_state, file=save_path)
            ptm = proc.time()[3]
          }
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
    if(!is.null(save_path) && proc.time()[3] - ptm > 600) {
      save(saved_state, file=save_path)
      ptm = proc.time()[3]
    }
  }
  
  while(TRUE) {
    #if no more edges free, stop
    if(length(saved_state$edges)==length(tree$edge) || length(saved_state$edges) == max_nshifts) break
    
    if(!parallel && is.null(saved_state$partial)) {
      saved_state$partial = list(edge_min = 0, time_min = NULL, p_min = NULL, min_lik = Inf, tested_edges = c())
    }
    
    #test max estimates for all edges
    all_edges = foreach::foreach (i = 1:length(tree$edge[,1]), .packages = "ML.MSBD") %d% {

      if(!parallel && i %in% saved_state$partial$tested_edges) return(list(edge = i, lik = Inf)) #already tested edge
      if(is.element(i,saved_state$edges)) return(list(edge = i, lik = Inf)) #disallow multiple shifts on the same edge

      if(tree$edge.length[i] == 0) { #zero-length edge
        if(!parallel) saved_state$partial$tested_edges = c(saved_state$partial$tested_edges,i)
        return(list(edge = i, lik = Inf))
      }
      
      if(fast_optim) anc_state = .find_ancestral_state(tree, i, saved_state$edges)
      
      edge_results = list(edge = i, lik = Inf, pars = NULL, time = NULL)
      for(time_pos in time_positions) {
        temp = .ML_optim(tree,c(saved_state$edges,i),
                         saved_state$initial_values,saved_state$times,
                         uniform_weights,p_lambda,p_mu,
                         rho,sigma,rho_sampling,
                         lineage_counts,tcut,
                         optim_control, time_pos,
                         stepsize,no_extinction,fixed_gamma,
                         unique_lambda,unique_mu, fast_optim, saved_state$params, anc_state)
        
        if(!is.na(temp$l) && temp$l < edge_results$lik) {
          edge_results$lik = temp$l
          edge_results$pars = temp$p
          edge_results$time = temp$t
        }
      }
      
      if(!parallel) {
        if(edge_results$lik < saved_state$partial$min_lik) {
          saved_state$partial$edge_min = i
          saved_state$partial$time_min = temp$t
          saved_state$partial$p_min = temp$p
          saved_state$partial$min_lik = temp$l
        }
        saved_state$partial$tested_edges = c(saved_state$partial$tested_edges,i)
        if(!is.null(save_path) && proc.time()[3] - ptm > 600) {
          save(saved_state, file=save_path)
          ptm = proc.time()[3]
        }
      }
      return(edge_results)
    }
    
    if(parallel) {
      liks = sapply(all_edges, function(x) x$lik)
      best = which(liks == min(liks))
      saved_state$partial = list(edge_min = all_edges[[best]]$edge, time_min = all_edges[[best]]$time, 
                                 p_min = all_edges[[best]]$pars, min_lik = all_edges[[best]]$lik)
    }
    
    saved_state$best_models = c(saved_state$best_models, saved_state$partial$min_lik)
    
    if(saved_state$likelihood > saved_state$partial$min_lik) {
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
                         lineage_counts,tcut,
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
      
      if(!is.null(save_path) && proc.time()[3] - ptm > 600) {
        save(saved_state, file=save_path)
        ptm = proc.time()[3]
      }
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
                     lineage_counts = c(), tcut = 0,
                     optim_control=list(), time_mode = c("tip","mid","root","fixed"),
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
    
    for(j in seq(along=edges)) {
      e = edges[j]
      t = times[j]
      shifts = rbind(shifts,c(e,t,1+j))
    }
    
    res = likelihood_MSBD_unresolved(tree,shifts,values$gamma,values$lambdas,values$mus,
                                     values$lambda_rates,stepsize,
                                     uniform_weights,p_lambda,p_mu,
                                     rho,sigma,rho_sampling,
                                     lineage_counts,tcut)
    return(res)
  }
  
  #initial parameters values
  if(!fast_optim) initp = .get_vector_from_initial_values(initial_values, n+1, stepsize, no_extinction, fixed_gamma, unique_lambda, unique_mu)
  else initp = .get_vector_from_initial_values_fast(initial_values, n+1, stepsize, no_extinction, fixed_gamma, unique_lambda, unique_mu)
  
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
    
    result = list(l=out$value, t=times[n])
    if(!fast_optim) result$p=.get_params_from_vector(out$par, n+1, stepsize, no_extinction, fixed_gamma, unique_lambda, unique_mu)
    else result$p=.get_params_from_vector_fast(out$par, n+1, current_values, anc_state,
                                               stepsize, no_extinction, fixed_gamma, unique_lambda, unique_mu)
  }
  return(result)
}