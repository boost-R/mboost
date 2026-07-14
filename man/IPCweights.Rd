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
  i.e., the \emph{observed data world} \bibcitep{mboost::Laan_Robins_2003}, see 
  \bibcitet{Hothorn:2006:Biostatistics:16344280} for an application to
  survival analysis.

}
\value{
  A vector of numeric weights. 
}
\references{\bibshow{*}}
\keyword{survival}
