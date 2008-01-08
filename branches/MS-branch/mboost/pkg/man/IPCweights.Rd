\name{IPCweights}
\alias{IPCweights}
\title{ Inverse Probability of Censoring Weights }
\description{
  Compute weights for censored regression via the inverted 
  probability of censoring principle.
}
\usage{
IPCweights(x, maxweight = 5)
}
\arguments{
  \item{x}{ an object of class \code{Surv}.}
  \item{maxweight}{ the maximal value of the returned weights. }
}
\details{

  Inverse probability of censoring weights are one possibility to fit
  models formulated in the \emph{full data world} in the presence of censoring,
  i.e., the \emph{observed data world}, see van der Laan and Robins (2003) for
  the underlying theory and Hothorn et al. (2006) for an application to
  survival analysis.

}
\value{
  A vector of numeric weights. 
}
\references{ 

    Mark J. van der Laan and James M. Robins (2003), 
    \emph{Unified Methods for Censored Longitudinal Data and Causality}, 
    Springer, New York.

    Torsten Hothorn, Peter Buhlmann, Sandrine Dudoit,
    Annette Molinaro and Mark J. van der Laan (2006), Survival ensembles.
    \emph{Biostatistics} \bold{7}(3), 355--373.

    Peter Buhlmann and Torsten Hothorn (2007),
    Boosting algorithms: regularization, prediction and model fitting.
    \emph{Statistical Science}, accepted.
    \url{http://www.imstat.org/sts/future_papers.html}

}
\keyword{survival}
