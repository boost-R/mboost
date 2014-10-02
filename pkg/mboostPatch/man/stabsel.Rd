\name{stabsel}
\alias{stabsel}
\alias{stabsel.mboost}
\alias{stabsel_parameters.mboost}
\title{
    Stability Selection
}
\description{
    Selection of influential variables or model components with error control.
}
\usage{
## a method to compute stability selection paths for fitted mboost models
\method{stabsel}{mboost}(x, cutoff, q, PFER,
        folds = subsample(model.weights(x), B = B),
        B = ifelse(sampling.type == "MB", 100, 50),
        assumption = c("unimodal", "r-concave", "none"),
        sampling.type = c("SS", "MB"),
        papply = mclapply, verbose = TRUE, FWER, eval = TRUE, ...)

## just a wrapper to stabsel(p, ..., eval = FALSE)
\method{stabsel_parameters}{mboost}(p, ...)
}
\arguments{
  \item{x, p}{an fitted model of class \code{"mboost"}.}
  \item{cutoff}{cutoff between 0.5 and 1. Preferably a value between 0.6
    and 0.9 should be used.}
  \item{q}{number of (unique) selected variables (or groups of variables
    depending on the model) that are selected on each subsample.}
  \item{PFER}{upper bound for the per-family error rate. This
    specifies the amount of falsely selected base-learners, which is
    tolerated. See details.}
  \item{folds}{ a weight matrix with number of rows equal to the number
    of observations, see \code{\link{cvrisk}} and
    \code{\link[stabs]{subsample}}. Usually one should not
    change the default here as subsampling with a fraction of \eqn{1/2}
    is needed for the error bounds to hold. One usage scenario where
    specifying the folds by hand might be the case when one has
    dependent data (e.g. clusters) and thus wants to draw clusters
    (i.e., multiple rows together) not individuals.}
  \item{assumption}{ Defines the type of assumptions on the
    distributions of the selection probabilities and simultaneous
    selection probabilities. Only applicable for
    \code{sampling.type = "SS"}. For \code{sampling.type = "MB"} we
    always use code{"none"}.}
  \item{sampling.type}{ use sampling scheme of of Shah & Samworth
    (2013), i.e., with complementarty pairs (\code{sampling.type = "SS"}),
    or the original sampling scheme of Meinshausen & Buehlmann (2010).}
  \item{B}{ number of subsampling replicates. Per default, we use 50
    complementary pairs for the error bounds of Shah & Samworth (2013)
    and 100 for the error bound derived in  Meinshausen & Buehlmann
    (2010). As we use \eqn{B} complementray pairs in the former case
    this leads to \eqn{2B} subsamples.}
  \item{papply}{ (parallel) apply function, defaults to
    \code{\link[parallel]{mclapply}}. Alternatively, \code{parLapply}
    can be used. In the latter case, usually more setup is needed (see
    example of \code{\link{cvrisk}} for some details).}
  \item{verbose}{ logical (default: \code{TRUE}) that determines wether
    \code{warnings} should be issued. }
  \item{FWER}{ deprecated. Only for compatibility with older versions,
    use PFER instead.}
  \item{eval}{ logical. Determines whether stability selection is
    evaluated (\code{eval = TRUE}; default) or if only the parameter
    combination is returned.}
  \item{\dots}{ additional arguments to parallel apply methods such as
    \code{\link{mclapply}} and to \code{\link{cvrisk}}.}
}
\details{

  For details see \code{\link[stabs]{stabsel}} in package \pkg{stabs}.

}
\value{
  An object of class \code{stabsel} with a special \code{print} method.
  The object has the following elements:
  \item{phat}{selection probabilities.}
  \item{selected}{elements with maximal selection probability greater
    \code{cutoff}.}
  \item{max}{maximum of selection probabilities.}
  \item{cutoff}{cutoff used.}
  \item{q}{average number of selected variables used.}
  \item{PFER}{per-family error rate.}
  \item{sampling.type}{the sampling type used for stability selection.}
  \item{assumption}{the assumptions made on the selection
    probabilities.}
  \item{call}{the call.}
}
\references{

  N. Meinshausen and P. Buehlmann (2010), Stability selection.
  \emph{Journal of the Royal Statistical Society, Series B},
  \bold{72}, 417--473.

  R.D. Shah and R.J. Samworth (2013), Variable selection with error
  control: another look at stability selection. \emph{Journal of the Royal
  Statistical Society, Series B}, \bold{75}, 55--80.

}
\seealso{
  \code{\link[stabs]{stabsel}} and
  \code{\link[stabs]{stabsel_parameters}}
}
\examples{
  ## make data set available
  data("bodyfat", package = "TH.data")
  ## set seed
  set.seed(1234)

  ### low-dimensional example
  mod <- glmboost(DEXfat ~ ., data = bodyfat)

  ## compute cutoff ahead of running stabsel to see if it is a sensible
  ## parameter choice.
  ##   p = ncol(bodyfat) - 1 (= Outcome) + 1 ( = Intercept)
  stabsel_parameters(q = 3, PFER = 1, p = ncol(bodyfat) - 1 + 1,
                     sampling.type = "MB")
  ## the same:
  stabsel(mod, q = 3, PFER = 1, sampling.type = "MB", eval = FALSE)

  ## now run stability selection
  (sbody <- stabsel(mod, q = 3, PFER = 1, sampling.type = "MB"))
  opar <- par(mai = par("mai") * c(1, 1, 1, 2.7))
  plot(sbody)
  par(opar)

  plot(sbody, type = "maxsel", ymargin = 6)
}
\keyword{nonparametric}
