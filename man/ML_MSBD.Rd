\name{ML_MSBD}
\alias{ML_MSBD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Full Maximum Likelihood inference of birth and death rates together with their changes along a phylogeny under a multi-type birth-death model.
}
\description{
Infers complete model from a phylogeny, including number of states, positions and timing of state changes, and parameters associated with each state. Uses a greedy approach to add states and maximum likelihood for the rest of the model.
}
\usage{
ML_MSBD(tree, initial_values,
        uniform_weights=TRUE, p_lambda=0, p_mu=0,
        rho = 1, sigma=0, rho_sampling = TRUE,
        unresolved = FALSE, lineage_counts = c(), tcut = 0,
        optim_control = list(), attempt_remove=TRUE, max_nshifts=Inf,
        stepsize=NULL, no_extinction=FALSE, fixed_gamma=NULL,
        unique_lambda = FALSE, unique_mu = FALSE,
        saved_state = NULL, save_path = NULL,
        time_mode = c("real","3pos","tip","mid","root"),
        fast_optim = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
The phylogenetic tree (in ape format) to run inference on
}
  \item{initial_values}{
Initial values for the optimizer, to be provided as a vector in this order: gamma (optional), lambda, lambda decay rate (optional), mu (optional). See 'Details'.
}
  \item{uniform_weights}{
Whether all states are weighted uniformly in shifts, default true. If false, the weights of states are calculated from the provided distributions on lambda and mu. See 'Details'.
}
  \item{p_lambda}{
Probability distribution on lambda, used if uniform_weights is false. 
}
  \item{p_mu}{
Probability distribution on mu, used if uniform_weights is false.
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
  \item{unresolved}{
Whether to use the 'cutoff' sampling scheme, default false. If true, the times at cutoff and the number of lineages need to be provided for all tips. This sampling scheme is not recommended for epidemiology datasets. This sampling scheme supports extinct tips outside of the unresolved parts.
}
  \item{lineage_counts}{
Number of lineages collapsed on each tip, used if unresolved is true. Should be set to 1 for extinct tips.
}
  \item{tcut}{
Times of collapsing for each tip (i.e time of the mrca of all collapsed lineages), used if unresolved is true. Can be a single number or a vector of length the number of tips (recommended).
}
  \item{optim_control}{
Control list for the optimizer, corresponds to control input in optim function, see ?optim for details.
}
  \item{attempt_remove}{
Whether to attempt to remove shifts at the end of the inference, default true. Full greedy approach if false.
}
  \item{max_nshifts}{
Maximum number of shifts to test for, default Inf.
}
  \item{stepsize}{
Size of the step to use for time discretization with exponential decay. Default NULL, will use exponential decay if a value is provided.
}
  \item{no_extinction}{
Whether to use the Yule process (mu=0) for all states, default false. If true no initial value for mu is needed.
}
  \item{fixed_gamma}{
Value to which gamma should be fixed, default NULL. If provided no initial value for gamma is needed.
}
  \item{unique_lambda}{
Whether to use the same value of lambda for all states, default false. If true and exponential decay is active all states will also share the value of lambda rate.
}
  \item{unique_mu}{
Whether to use the same value of mu for all states, default false.
}
  \item{saved_state}{
If provided, the inference will be restarted from this state.
}
  \item{save_path}{
If provided, the progress of the inference will be saved to this path after each optimization step.
}
  \item{time_mode}{
String controlling how time positions of the shifts will be inferred. See 'Details'.
}
  \item{fast_optim}{
Whether to use the faster mode of optimization, default FALSE. If true only rates associated with the state currently being added to the tree and its ancestor will be optimized at each step, otherwise all rates are optimized.
  }
}
\details{
It is to be noted that all times are counted backwards, with the most recent tip positioned at 0.\cr
\cr
Five time modes are possible for the input \code{time_mode}. 
In \code{real} mode, the time positions of shifts will be estimated in addition to the other parameters. 
In \code{tip} mode, the shifts will be placed at 10\% of the length of the edge. 
In \code{mid} mode, the shifts will be placed at 50\% of the length of the edge. 
In \code{root} mode, the shifts will be placed at 90\% of the length of the edge. 
In \code{3pos} mode, the three "tip", "mid" and "root" positions will be tested.\cr
\cr
The weights w are used for calculating the transition rates q from each state i to j: \eqn{q_{i,j}=\gamma*w_{i,j}}{q(i,j)=\gamma*w(i,j)}. If \code{uniform_weights} is true, \eqn{w_{i,j} = \frac{1}{N-1}}{w(i,j)=1/(N-1)} for all i,j, where N is the total number of states. If \code{uniform_weights} is false, \eqn{w_{i,j} = \frac{p_\lambda(\lambda_j)p_\mu(\mu_j)}{sum_{k \ne i}p_\lambda(\lambda_k)p_\mu(\mu_k)}}{w(i,j)=p\lambda(\lambdaj)p\mu(\muj)/sum(p\lambda(\lambdak)p\mu(\muk)) for all k!=i} where the distributions \eqn{p_\lambda}{p\lambda} and \eqn{p_\mu}{p\mu} are provided by the inputs \code{p_lambda} and \code{p_mu}.\cr
\cr
Initial values for the optimization need to be provided as a vector and contain the following elements (in order): an initial value for gamma, which is required unless \code{fixed_gamma} is provided, an initial value for lambda which is always required, an initial value for lambda decay rate, which is required if \code{stepsize} is provided, and an initial value for mu, which is required unless \code{no_extinction} is true. 
An error will be raised if the number of initial values provided does not match the one expected from the rest of the settings, and the function will fail if the likelihood cannot be calculated at the initial values.
}
\value{
Returns a list describing the most likely model found, with the following components:
\item{likelihood}{the negative log likelihood of the model}
\item{shifts.edge}{the indexes of the edges where shifts happen, 0 indicates the root state}
\item{shifts.time}{the time positions of shifts}
\item{gamma}{the rate of state change}
\item{lambdas}{the birth rates of all states}
\item{lambda_rates}{if exponential decay was activated, the rates of decay of birth rate for all states}
\item{mus}{the death rates of all states}
\item{best_models}{a vector containing the negative log likelihood of the best model found for each number of states tested (best_models[i] corresponds to i states, i.e i-1 shifts)}
All vectors are indexed in the same way, so that the state with parameters lambdas[i], lambda_rates[i] and mus[i] starts on edge shifts.edge[i] at time shifts.time[i].
}
\author{
\packageAuthor{ML.MSBD}
}
\examples{
# Simulate a random phylogeny with extinct samples
set.seed(25)
tree <- ape::rtree(10)

# Infer the most likely multi-states birth-death model 
# with full extant & extinct sampling
\donttest{ML_MSBD(tree, initial_values = c(0.1, 10, 1), sigma = 1, time_mode = "mid") }
# Infer the most likely multi-states birth-death model with exponential decay
# and full extant & extinct sampling
\donttest{ML_MSBD(tree, initial_values = c(0.1, 10, 0.5, 1), sigma = 1, 
  stepsize = 0.1, time_mode = "mid")}

# Simulate a random phylogeny with extant samples
set.seed(24)
tree2 <- ape::rcoal(10)

# Infer the most likely multi-states Yule model with partial extant sampling
\donttest{ML_MSBD(tree2, initial_values = c(0.1, 10), no_extinction = TRUE, 
  rho = 0.5, time_mode = "mid")}
# Infer the most likely multi-states birth-death model with full extant sampling 
# and unresolved extant tips
\donttest{ML_MSBD(tree2, initial_values = c(0.1, 10, 1), unresolved = TRUE, 
  lineage_counts = c(2,5,1,3,1,1,1,1,2,6), tcut = 0.05, time_mode = "mid")}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }% use one of  RShowDoc("KEYWORDS")