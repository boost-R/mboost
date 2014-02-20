\name{stabsel}
\alias{stabsel}
\alias{stabsel_parameters}
\title{
    Stability Selection
}
\description{
    Selection of influential variables or model components with error control.
}
\usage{
stabsel(object, cutoff, q, PFER,
        folds = cv(model.weights(object), type = "subsampling",
                   B = ifelse(sampling.type == "MB", 100, 50)),
        assumption = c("unimodal", "r-concave", "none"),
        sampling.type = c("SS", "MB"),
        papply = mclapply, verbose = TRUE, FWER, ...)

## function to compute missing parameter from the other two parameters
## (internally used within stabsel)
stabsel_parameters(cutoff, q, PFER, p,
                   B = ifelse(sampling.type == "MB", 100, 50),
                   assumption = c("unimodal", "r-concave", "none"),
                   sampling.type = c("SS", "MB"),
                   verbose = FALSE, FWER)
}
\arguments{
  \item{object}{an \code{mboost} object.}
  \item{cutoff}{cutoff between 0.5 and 1. Preferably a value between 0.6
    and 0.9 should be used.}
  \item{q}{number of (unique) selected base-learners per boosting run.}
  \item{PFER}{upper bound for the per-family error rate. This
    specifies the amount of falsely selected base-learners, which is
    tolerated. See details.}
  \item{folds}{ a weight matrix with number of rows equal to the number
    of observations, see \code{\link{cvrisk}}.}
  \item{assumption}{ Defines the type of assumptions on the
    distributions of the selection probabilities and simultaneous
    selection probabilities. Only applicable for
    \code{sampling.type = "SS"}. For \code{sampling.type = "MB"} we
    always use code{"none"}.}
  \item{sampling.type}{ use sampling scheme of of Shah & Samworth
    (2013), i.e., with complementarty pairs (\code{sampling.type = "SS"}),
    or the original sampling scheme of Meinshausen & Buehlmann (2010).}
  \item{p}{ number of possible predictors (including intercept if
    applicable).}
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
  \item{\dots}{additional arguments to \code{\link{cvrisk}} and
    further arguments to parallel apply methods such as
    \code{\link{mclapply}}.}
}
\details{

  This function implements the stability selection procedure
  by Meinshausen and Buehlmann (2010) and the improved error bounds by
  Shah and Samworth (2013).

  Two of the three arguments \code{cutoff}, \code{q} and \code{PFER}
  \emph{must} be specified. The per-family error rate (PFER), i.e., the
  expected number of false positives \eqn{E(V)}, where \eqn{V} is the
  number of false positives, is bounded by the argument \code{PFER}.

  As controlling the PFER is more conservative as controlling the
  family-wise error rate (FWER), the procedure also controlls the FWER,
  i.e., the probability of selecting at least one non-influential
  variable (or model component) is less than \code{PFER}.

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
  \bold{72}:417--473.

  R.D. Shah and R.J. Samworth (2013), Variable selection with error
  control: another look at stability selection. \emph{Journal of the Royal
  Statistical Society, Series B}, \bold{75}:55--80.

}
\examples{

  data("bodyfat", package = "TH.data")

  ### low-dimensional example
  mod <- glmboost(DEXfat ~ ., data = bodyfat)

  ## compute cutoff ahead of running stabsel to see if it is a sensible
  ## parameter choice.
  ##   p = ncol(bodyfat) - 1 (= Outcome) + 1 ( = Intercept)
  stabsel_parameters(q = 3, PFER = 1, p = ncol(bodyfat) - 1 + 1,
                     sampling.type = "MB")

  ## now run stability selection; to make results reproducible
  set.seed(1234)
  (sbody <- stabsel(mod, q = 3, PFER = 1, sampling.type = "MB"))
  opar <- par(mai = par("mai") * c(1, 1, 1, 2.7))
  plot(sbody)
  par(opar)

}
\keyword{nonparametric}
