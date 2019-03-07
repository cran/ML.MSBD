## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(ML.MSBD)

## ------------------------------------------------------------------------
tree <- ape::read.tree(text = "((((t3:0.04098955599,(t8:0.03016935301,t10:0.03016935301):0.01082020298):0.06041650538,t2:0.1014060614):0.7530620333,(t1:0.5805635547,(t5:0.2225489503,t7:0.2225489503):0.3580146044):0.2739045399):1.005016052,((t9:0.09730025079,t6:0.09730025079):0.0209186338,t4:0.1182188846):1.741265262);")

## ---- eval=FALSE---------------------------------------------------------
#  ML_MSBD(tree, initial_values = c(0.1, 10, 1), time_mode = "mid")

## ------------------------------------------------------------------------
tree_stt <- ape::read.tree(text = "((t6:0.1203204831,t3:0.1251815392):0.527233894,(((t8:0.1153702512,t5:0.2362936393):0.1328287013,((t1:0.6202914508,t4:0.839935567):0.6779029006,t10:0.472941763):0.3049683478):0.1837210727,((t9:0.24675539,t2:0.8737900169):0.6451289295,t7:0.1376576095):0.8951371317):0.5465055171);")

## ---- eval=FALSE---------------------------------------------------------
#  ML_MSBD(tree_stt, initial_values = c(0.1, 10, 1), rho = 0.5, sigma = 0.1, time_mode = "mid")
#  
#  ML_MSBD(tree_stt, initial_values = c(0.1, 10, 1), rho_sampling = FALSE, sigma = 0.1, time_mode = "mid")

## ---- eval=FALSE---------------------------------------------------------
#  ML_MSBD(tree, initial_values = c(0.1, 10, 1), c(0.1, 10, 0.5, 1), sigma = 1, stepsize = 0.1, time_mode = "mid")

## ------------------------------------------------------------------------
likelihood_MSBD(tree, shifts = matrix(c(2,0.8,2), nrow = 1), gamma = 0.05, lambdas = c(10, 6), mus = c(1, 0.5))

## ------------------------------------------------------------------------
likelihood_MSBD(tree_stt, shifts = c(), gamma = 0.05, lambdas = 10, mus = 0.5, sigma = 0.5, rho_sampling = FALSE)
likelihood_MSBD(tree, shifts = c(), gamma = 0.05, lambdas = 10, mus = 0.5, lambda_rates = 0.1, stepsize = 0.05)

## ------------------------------------------------------------------------
tree_collapsed = ape::read.tree(text = "((t1:0.379876463,t3:0.379876463):1.668231124,(t2:0.5653793315,t4:0.5653793315):1.482728255);")
likelihood_MSBD_unresolved(tree_collapsed, shifts = matrix(c(2,0.25,2), nrow = 1), gamma = 0.05, lambdas = c(10, 6), mus = c(1, 0.5), lineage_counts = c(5,1,3,6), tcut = 0.1)
likelihood_MSBD_unresolved(tree_collapsed, shifts = matrix(c(2,0.25,2), nrow = 1), gamma = 0.05, lambdas = c(10, 6), mus = c(1, 0.5), lineage_counts = c(5,1,3,6), tcut = c(0.1,0.0,0.15,0.4))

