\name{mboost}
\alias{mboost}
\alias{mboost_fit}
\title{ Model-based Gradient Boosting }
\description{
  Gradient boosting for optimizing arbitrary loss functions, where component-wise
  models are utilized as base-learners.
}
\usage{
mboost(formula, data = list(), 
       baselearner = c("bbs", "bols", "btree", "bss", "bns"), ...) 
mboost_fit(blg, response, weights = NULL, offset = NULL,
           family = Gaussian(), control = boost_control())
}
\arguments{
  \item{formula}{ a symbolic description of the model to be fit. }
  \item{data}{ a data frame containing the variables in the model. }
  \item{baselearner}{ a character specifying the component-wise base  
    learner to be used: \code{\link{bbs}} means P-splines with a
    B-spline basis (see Schmid and Hothorn 2008), \code{\link{bols}}
    linear models abd \code{\link{btree}} boosts stumps.
    \code{bss} and \code{bns} are deprecated.
    Component-wise smoothing splines have been considered in Buehlmann
    and Yu (2003) and Schmid and Hothorn (2008) investigate P-splines
    with a B-spline basis. Kneib, Hothorn and Tutz (2009) also utilize
    P-splines with a B-spline basis, supplement them with their
    bivariate tensor product version to estimate interaction surfaces
    and spatial effects and also consider random effects base
    learners.} 
  \item{blg}{ a list of objects of class \code{blg}, as returned by all base-learners.}
  \item{response}{ the response variable.}
  \item{weights}{ an optionally numeric vector of weights.}
  \item{offset}{ an optionally numeric vector to be used as offset.}
  \item{family}{a \code{Family} object.}
  \item{control}{ a list of parameters controlling the algorithm.}
  \item{\dots}{ additional arguments passed to \code{mboost}.}
}
\details{

  The function implements component-wise functional gradient boosting
  in a generic way. Basically,
  the algorithm is initialized with a function for computing
  the negative gradient of the loss function (via its \code{family} argument)
  and one or more base-learners (given as \code{blg}).

  The algorithm minimized the in-sample empirical risk defined as 
  the weighted sum (by \code{weights}) of the loss function (corresponding
  to the negative gradient) evaluated at the data.

  The structure of the model is determined by the structure 
  of the base-learners. If more than one base-learner is given,
  the model is additive in these components.

  Base-learners can be specified via a formula interface 
  (function \code{mboost}) or as a list of objects of class \code{bl},
  see \code{\link{baselearners}}.

  Note that the more convenient modeling interfaces \code{\link{gamboost}},
  \code{\link{glmboost}} and \code{\link{blackboost}} all call 
  \code{mboost} directly.

}
\value{
  An object of class \code{mboost} with \code{\link{print}},
  \code{\link{AIC}}, \code{\link{plot}} and \code{\link{predict}} methods being available.
}
\references{

  Peter Buehlmann and Bin Yu (2003),
  Boosting with the L2 loss: regression and classification.
  \emph{Journal of the American Statistical Association}, \bold{98},
  324--339.

  Peter Buehlmann and Torsten Hothorn (2007),
  Boosting algorithms: regularization, prediction and model fitting.
  \emph{Statistical Science}, \bold{22}(4), 477--505.

  Yoav Freund and Robert E. Schapire (1996),
  Experiments with a new boosting algorithm.
  In \emph{Machine Learning: Proc. Thirteenth International Conference},
  148--156.

  Jerome H. Friedman (2001),
  Greedy function approximation: A gradient boosting machine.
  \emph{The Annals of Statistics}, \bold{29}, 1189--1232.

}
\seealso{\code{\link{glmboost}} for boosted linear models and
  \code{\link{blackboost}} for boosted trees. See e.g. \code{\link{bbs}}
  for possible base-learners. See \code{\link{cvrisk}} for
  cross-validated stopping iteration. Furthermore see
  \code{\link{boost_control}}, \code{\link{Family}} and  \code{\link[mboost]{methods}}}
\examples{

  data("bodyfat", package = "mboost")

  ### formula interface: additive Gaussian model with
  ### a non-linear step-function in `age', a linear function in `waistcirc'
  ### and a smooth non-linear smooth function in `hipcirc'
  mod <- mboost(DEXfat ~ btree(age) + bols(waistcirc) + bbs(hipcirc), 
                data = bodyfat)
  layout(matrix(1:6, nc = 3, byrow = TRUE))
  plot(mod, ask = FALSE, main = "formula")

  ### the same
  with(bodyfat, 
       mod <- mboost_fit(list(btree(age), bols(waistcirc), bbs(hipcirc)), 
                         response = DEXfat))
  plot(mod, ask = FALSE, main = "base-learner")
}
\keyword{models}
\keyword{nonlinear}
