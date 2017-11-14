.find_ancestral_state = function(tree, edgei, shifts) {
  if(length(shifts)==0) return(1)
  
  x = tree$edge[edgei,1]
  while(x > length(tree$tip.label)+1){
    anc = which(tree$edge[,2]==x)
    if(anc %in% shifts) {
      state = which(shifts==anc)+1
      return(state)
    }
    x = tree$edge[anc,1]
  }
  return(1)
}

.remove_state = function(values, i, 
                         stepsize=NULL, no_extinction=FALSE, unique_lambda = FALSE, unique_mu = FALSE) {
  values$lambdas = values$lambdas[-i]
  if(!is.null(stepsize)) values$lambda_rates = values$lambda_rates[-i]
  values$mus = values$mus[-i]
  
  values
}