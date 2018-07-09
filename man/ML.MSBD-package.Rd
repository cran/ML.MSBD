\name{ML.MSBD-package}
\alias{ML.MSBD-package}
\alias{ML.MSBD}
\docType{package}
\title{
\packageTitle{ML.MSBD}
}
\description{
\packageDescription{ML.MSBD}
}
\details{
\packageDESCRIPTION{ML.MSBD}
\packageIndices{ML.MSBD}
}
\author{
\packageAuthor{ML.MSBD}

Maintainer: \packageMaintainer{ML.MSBD}
}
\references{
J. Barido-Sottani and T. Stadler. Accurate detection of HIV transmission clusters from phylogenetic trees using a multi-state birth-death model, BioRXiv 2017. (https://www.biorxiv.org/content/early/2017/11/10/215491)
}

\keyword{ package }
\seealso{
\code{\link[ape:ape-package]{ape}}
}
\examples{
# Simulate a random phylogeny
set.seed(25)
tree <- ape::rtree(10)

# Calculate the log likelihood under a multi-states model with 2 states 
# and full extant & extinct sampling
likelihood_MSBD(tree, shifts = matrix(c(2,1.8,2), nrow = 1), 
  gamma = 0.05, lambdas = c(10, 6), mus = c(1, 0.5), sigma = 1)

# Infer the most likely multi-states birth-death model with full extant & extinct sampling
\dontrun{ML_MSBD(tree, initial_values = c(0.1, 10, 1), sigma = 1, time_mode = "mid") }
# Infer the most likely multi-states birth-death model with exponential decay 
# and full extant & extinct sampling
\dontrun{ML_MSBD(tree, initial_values = c(0.1, 10, 0.5, 1), sigma = 1, 
  stepsize = 0.1, time_mode = "mid")}
}
