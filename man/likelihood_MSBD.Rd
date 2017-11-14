\name{likelihood_MSBD}
\alias{likelihood_MSBD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Likelihood calculation for constantly sampled trees
}
\description{
Calculates the negative log likelihood of a multi-states model given a tree. This function is designed to work with constant extant and/or extinct sampling.
}
\usage{
likelihood_MSBD(tree, shifts, gamma, lambdas, mus, 
                lambda_rates = NULL, stepsize = NULL, 
                uniform_weights = TRUE, p_lambda = 0, p_mu = 0, 
                rho = 1, sigma = 0, rho_sampling = TRUE,
                add_time = 0, unresolved = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
The phylogenetic tree (in ape format) to calculate the likelihood on
}
  \item{shifts}{
Matrix describing the positions (edges and times) of shifts. See 'Details'.
}
  \item{gamma}{
Rate of state change
}
  \item{lambdas}{
Birth rates of all states
}
  \item{mus}{
Death rates of all states
}
  \item{lambda_rates}{
Rates of decay of birth rate for all states. If provided, stepsize should also be provided.
}
  \item{stepsize}{
Size of the step to use for time discretization with exponential decay, default NULL. If provided, lambda_rates should also be provided.
}
  \item{uniform_weights}{
Whether all states are weighted uniformly in shifts, default true. If false, the weights of states are calculated from the provided distributions on lambda and mu. See 'Details'.
}
  \item{p_lambda}{
Prior probability distribution on lambda, used if uniform_weights is false.
}
  \item{p_mu}{
Prior probability distribution on mu, used if uniform_weights is false.
}
  \item{rho}{
Sampling proportion on extant tips, default 1.
}
  \item{sigma}{
Sampling probability on extinct tips (tips are sampled upon extinction), default 0.
}
  \item{rho_sampling}{
If true, the most recent tips will be considered extant tips and use the sampling proportion rho. If false, all tips will be considered extinct tips and use the sampling probability sigma. Should be true for most macroevolution datasets and false for most epidemiology datasets.
}
  \item{add_time}{
The time between the most recent tip and the end of the process (>=0). This is an internal variable used in calculations for unresolved trees.
}
    \item{unresolved}{
Whether to use the 'cutoff' sampling scheme, default false. If true, the times at cutoff and the number of lineages need to be provided for all tips. This sampling scheme is not recommended for epidemiology datasets. This sampling scheme supports extinct tips outside of the unresolved parts.
}
}
\details{
It is to be noted that all times are counted backwards, with the most recent tip positioned at 0. \cr
\cr
The 'shifts' matrix is composed of 3 columns and a number of rows. Each row describes a shift: the first column is the index of the edge on which the shift happens, the second column is the time of the shift, and the third column is the index of the new state. For example the row vector (3,0.5,2) specifies a shift on edge number 3, at time 0.5, towards the state that has parameters lambdas[2], lambda_rates[2] and mus[2]. \cr
\cr
The weights w are used for calculating the transition rates q from each state i to j: \eqn{q_{i,j}=\gamma*w_{i,j}}{q(i,j)=\gamma*w(i,j)}. If \code{uniform_weights} is true, \eqn{w_{i,j} = \frac{1}{N-1}}{w(i,j)=1/(N-1)} for all i,j, where N is the total number of states. If \code{uniform_weights} is false, \eqn{w_{i,j} = \frac{p_\lambda(\lambda_j)p_\mu(\mu_j)}{sum_{k \ne i}p_\lambda(\lambda_k)p_\mu(\mu_k)}}{w(i,j)=p\lambda(\lambdaj)p\mu(\muj)/sum(p\lambda(\lambdak)p\mu(\muk)) for all k!=i} where the distributions \eqn{p_\lambda}{p\lambda} and \eqn{p_\mu}{p\mu} are provided by the inputs \code{plambda} and \code{pmu}.
}
\value{
The value of the negative log likelihood of the model given the tree.
}
\author{
\packageAuthor{ML.MSBD}
}
\examples{
# Simulate a random phylogeny
set.seed(25)
tree <- ape::rtree(10)

# Calculate the log likelihood under a constant birth-death model (i.e, no shifts) 
# with full extant & extinct sampling
likelihood_MSBD(tree, shifts = c(), gamma = 0, lambdas = 10, mus = 1, sigma = 1)
# Calculate the log likelihood under a multi-states model with 2 states 
# and full extant & extinct sampling
likelihood_MSBD(tree, shifts = matrix(c(2,1.8,2), nrow = 1), 
  gamma = 0.05, lambdas = c(10, 6), mus = c(1, 0.5), sigma = 1)
# Calculate the log likelihood under a multi-states model with 2 states and exponential decay 
# with full extant & extinct sampling
likelihood_MSBD(tree, shifts = matrix(c(2,1.8,2), nrow = 1), 
  gamma = 0.05, lambdas = c(10, 6), mus = c(1, 0.5), 
  sigma = 1, stepsize = 0.01, lambda_rates = c(0.1, 0.1))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }% use one of  RShowDoc("KEYWORDS")
