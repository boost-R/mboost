\name{boost_control}
\alias{boost_control}
\title{ Control Hyper Parameters for Boosting Algorithms }
\description{
  Definition of the initial number of boosting iterations, step size 
  and other hyper parameters for boosting algorithms.
}
\usage{
boost_control(mstop = 100, nu = 0.1, constraint = FALSE, 
              risk = c("inbag", "oobag", "none"), 
              savedata = TRUE, center = FALSE, trace = FALSE)
}
\arguments{
  \item{mstop}{ an integer giving the number of initial boosting iterations. }
  \item{nu}{ a double (between 0 and 1) defining the step size or shrinkage parameter. } 
  \item{constraint}{ a logical indicating whether the working responses
                     should be restricted to \eqn{(-1, +1)}. }
  \item{risk}{ a character indicating how the empirical risk should be 
               computed for each boosting iteration. \code{inbag} leads
               risks computed for the learning sample (i.e., all non-zero weights),
               \code{oobag} to risks based on the out-of-bag (all observations with 
               zero weights) and \code{none} to no risk computations at all.}
  \item{savedata}{ a logical, should the data be saved in the 
                   returned object?}
  \item{center}{ a logical indicating if the numerical covariates should be mean 
                 centered before fitting. Only implemented for \code{\link{glmboost}} and
                 \code{\link{gamboost}}.}
  \item{trace}{ a logical triggering printout of status information during 
                the fitting process.}
}
\details{

  Objects returned by this function specify hyper parameters of the
  boosting algorithms implemented in \code{\link{glmboost}}, 
  \code{\link{gamboost}} and \code{\link{blackboost}} 
  (via the \code{control} argument).

}
\value{
  An object of class \code{boost_control}, a list.
}
\keyword{misc}
