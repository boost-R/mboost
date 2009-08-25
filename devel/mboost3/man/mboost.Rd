\name{mboost}
\alias{mboost}
\alias{mboost_fit}
\title{ Model-based Gradient Boosting }
\description{
  Gradient boosting for optimizing arbitrary loss functions, where component-wise
  models are utilized as base learners.
}
\usage{
mboost(formula, data = list(), baselearner = c("bbs", "bols", "btree", "bss", "bns"), ...) 
mboost_fit(blg, response, weights = NULL, offset = NULL,
           family = GaussReg(), control = boost_control())
}
\arguments{
  \item{formula}{ a symbolic description of the model to be fit. }
  \item{data}{ a data frame containing the variables in the model. }
  \item{baselearner}{ a character specifying the component-wise base learner to be used:
      \code{\link{bbs}} means P-splines with a B-spline basis (see Schmid and Hothorn
      2007), 
      \code{\link{bols}} linear models, \code{\link{bspatial}} bivariate tensor
      product penalized splines, and \code{\link{brandom}} random effects.
      In addition, \code{\link{btree}} boosts stumps.
      Component-wise smoothing splines
      have been considered in Buhlmann and Yu (2003) and Schmid and Hothorn (2007) investigate P-splines with
      a B-spline basis. Kneib, Hothorn and Tutz (2007) also utilise P-splines with a B-spline basis,
      supplement them with their bivariate tensor product version to estimate interaction surfaces and spatial effects
      and also consider random effects base learners.}
  \item{blg}{ a list of objects of class \code{blg}, as returned by all base-learners.}
  \item{response}{ the response variable.}
  \item{weights}{ an optionally numeric vector of weights.}
  \item{offset}{ an optionally numeric vector to be used as offset.}
  \item{family}{a \code{Family} object.}
  \item{control}{ a list of parameters controlling the algorithm.}
  \item{\dots}{ additional arguments passed to \code{mboost}.}
}
\details{
to be done
}
\value{
  An object of class \code{mboost} with \code{\link{print}},
  \code{\link{AIC}}, \code{\link{plot}} and \code{\link{predict}} methods being available.
}
\references{

  Peter Buhlmann and Bin Yu (2003),
  Boosting with the L2 loss: regression and classification.
  \emph{Journal of the American Statistical Association}, \bold{98},
  324--339.

  Peter Buhlmann and Torsten Hothorn (2007),
  Boosting algorithms: regularization, prediction and model fitting.
  \emph{Statistical Science}, \bold{22}(4), 477--505.

}
\seealso{\code{\link{glmboost}} for boosted linear models and
  \code{\link{blackboost}} for boosted trees. See e.g. \code{\link{bbs}}
  for possible base-learners. See \code{\link{cvrisk}} for
  cross-validated stoping iteration. Furthermore see
  \code{\link{boost_control}}, \code{\link{Family}} and  \code{\link[mboost]{methods}}}
\examples{

  data("bodyfat", package = "mboost3")

  mod <- mboost(DEXfat ~ btree(age) + bols(waistcirc) + bbs(hipcirc), data = bodyfat)

  layout(matrix(1:3, nc = 3))
  plot(mod, ask = FALSE)

}
\keyword{models}
\keyword{nonlinear}
