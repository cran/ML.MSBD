.get_vector_from_initial_values = function(values,nstates,
                                           stepsize=NULL, no_extinction=FALSE, fixed_gamma=NULL, 
                                           unique_lambda = FALSE, unique_mu = FALSE) {
  vector = c()
  
  if(is.null(fixed_gamma)) vector = values$gamma
  
  if(!unique_lambda) vector = c(vector, rep(values$lambdas[1],nstates))
  else vector = c(vector, values$lambdas[1])
  
  if(!is.null(stepsize)) {
    if(!unique_lambda) {
      vector = c(vector, rep(values$lambda_rates[1],nstates)) 
    }
    else { #TODO - lambda rate controlled differently than lambda ?
      vector = c(vector, values$lambda_rates[1])
      warning("Lambda rates will be shared between states as well")
    }
  }
  
  if(!no_extinction) {
    if(!unique_mu) vector = c(vector, rep(values$mus[1],nstates))
    else vector = c(vector, values$mus[1])
  }
  vector
}

.get_vector_from_params = function(values,
                                   stepsize=NULL, no_extinction=FALSE, fixed_gamma=NULL) {
  vector = c()
  
  if(is.null(fixed_gamma)) vector = values$gamma
  vector = c(vector, values$lambdas)
  if(!is.null(stepsize)) vector = c(vector, values$lambda_rates)
  if(!no_extinction) vector = c(vector, values$mus)
  
  vector
}

.get_params_from_vector = function(vector, nstates,
                                   stepsize=NULL, no_extinction=FALSE, fixed_gamma=NULL, 
                                   unique_lambda = FALSE, unique_mu = FALSE) {
  values = list()
  k = 1
  if(is.null(fixed_gamma)) {
    values$gamma = vector[k]
    k = k+1
  }
  else values$gamma = fixed_gamma
  
  if(!unique_lambda) {
    values$lambdas = vector[k:(k+nstates-1)]
    k = k+nstates
  }
  else {
    values$lambdas = rep(vector[k],nstates)
    k = k+1
  }
  
  if(!is.null(stepsize)) {
    if(!unique_lambda) {
      values$lambda_rates = vector[k:(k+nstates-1)]
      k = k+nstates
    }
    else {
      values$lambda_rates = rep(vector[k],nstates)
      k = k+1
    }
  }
  
  if(!no_extinction) {
    if(!unique_mu) values$mus = vector[k:(k+nstates-1)]
    else values$mus = rep(vector[k],nstates)
  }
  else values$mus = rep(0,nstates)
  
  values
}

.get_initial_values_from_vector = function(vector,
                                           stepsize=NULL, no_extinction=FALSE, fixed_gamma=NULL) {
  #length check
  k = 1
  if(is.null(fixed_gamma)) k = k+1
  if(!is.null(stepsize)) k = k+1
  if(!no_extinction) k = k+1
  
  if(length(vector) != k) stop(paste0("Expected ", k, " initial values but got ", length(vector)))
  
  values = list()
  k = 1
  if(is.null(fixed_gamma)) {
    values$gamma = vector[k]
    k = k+1
  }
  
  values$lambdas = vector[k]
  k = k+1
  
  if(!is.null(stepsize)) {
    values$lambda_rates = vector[k]
    k = k+1
  }
  
  if(!no_extinction) values$mus = vector[k]
  
  values
}

.get_params_from_vector_fast = function(vector, nstates, values, anc_state,
                                        stepsize=NULL, no_extinction=FALSE, fixed_gamma=NULL, 
                                        unique_lambda = FALSE, unique_mu = FALSE) {
  k = 1
  if(is.null(fixed_gamma)) {
    values$gamma = vector[k]
    k = k+1
  }
  else values$gamma = fixed_gamma
  
  if(!unique_lambda && nstates>1) {
    values$lambdas[c(anc_state,nstates)] = vector[k:(k+1)]
    k = k+2
  }
  else {
    values$lambdas = rep(vector[k],nstates)
    k = k+1
  }
  
  if(!is.null(stepsize)) {
    if(!unique_lambda && nstates>1) {
      values$lambda_rates[c(anc_state,nstates)] = vector[k:(k+1)]
      k = k+2
    }
    else {
      values$lambda_rates = rep(vector[k],nstates)
      k = k+1
    }
  }
  
  if(!no_extinction) {
    if(!unique_mu && nstates>1) values$mus[c(anc_state,nstates)] = vector[k:(k+1)]
    else {
      values$mus = rep(vector[k],nstates)
    }
  }
  else values$mus = rep(0,nstates)
  
  values
}

.get_vector_from_initial_values_fast = function(values, nstates,
                                                stepsize=NULL, no_extinction=FALSE, fixed_gamma=NULL, 
                                                unique_lambda = FALSE, unique_mu = FALSE) {
  vector = c()
  
  if(is.null(fixed_gamma)) vector = values$gamma
  
  if(!unique_lambda && nstates>1) vector = c(vector, rep(values$lambdas[1],2))
  else vector = c(vector, values$lambdas[1])
  
  if(!is.null(stepsize)) {
    if(!unique_lambda && nstates>1) {
      vector = c(vector, rep(values$lambda_rates[1],2)) 
    }
    else { #TODO - lambda rate controlled differently than lambda ?
      vector = c(vector, values$lambda_rates[1])
      if(unique_lambda) warning("Lambda rates will be shared between states as well")
    }
  }
  
  if(!no_extinction) {
    if(!unique_mu && nstates>1) vector = c(vector, rep(values$mus[1],2))
    else vector = c(vector, values$mus[1])
  }
  vector
}