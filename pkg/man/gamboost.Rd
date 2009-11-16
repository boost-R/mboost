\name{gamboost}
\alias{gamboost}
\title{ Gradient Boosting with Smooth Components }
\description{
  Gradient boosting for optimizing arbitrary loss functions, where component-wise
  smoothing procedures are utilized as base learners.
}
\usage{
gamboost(formula, data = list(), baselearner = c("bbs", "bols", "btree", "bss", "bns"),
         dfbase = 4, ...)
}
\arguments{
  \item{formula}{ a symbolic description of the model to be fit. }
  \item{data}{ a data frame containing the variables in the model. }
  \item{baselearner}{ a character specifying the component-wise base
    learner to be used: \code{\link{bbs}} means P-splines with a
    B-spline basis (see Schmid and Hothorn 2007), \code{\link{bols}}
    linear models abd \code{\link{btree}} boosts stumps.
    \code{bss} and \code{bns} are only used for backward
    compatibility. You should consider to use \code{\link{bbs}} instead.
    Component-wise smoothing splines have been considered in Buhlmann
    and Yu (2003) and Schmid and Hothorn (2007) investigate P-splines
    with a B-spline basis. Kneib, Hothorn and Tutz (2007) also utilize
    P-splines with a B-spline basis, supplement them with their
    bivariate tensor product version to estimate interaction surfaces
    and spatial effects and also consider random effects base
    learners.} 
  \item{dfbase}{ an integer vector giving the degrees of freedom for the smoothing
    spline, either globally for all variables (when its length is one)
    or separately for each single covariate. }
  \item{\dots}{ additional arguments passed to \code{mboost}.}
}
\details{

  A (generalized) additive model is fitted using a boosting algorithm based on
  component-wise univariate base learners. The base-learners can either be
  specified via the \code{formula} object or via the \code{baselearner} argument
  (see \code{\link{bbs}} for an example). If the base learners specified in
  \code{formula} differ from \code{baselearner}, the latter argument will be
  ignored. Furthermore, two additional base-learners can be specified in
  \code{formula}: \code{\link{bspatial}} for bivariate tensor product
  penalized splines and \code{\link{brandom}} for random effects. 
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

  Thomas Kneib, Torsten Hothorn and Gerhard Tutz (2009), Variable selection and
  model choice in geoadditive regression models, \emph{Biometrics}, \bold{65}(2),
  626--634.

  Matthias Schmid and Torsten Hothorn (2008),
  Boosting additive models using component-wise P-splines as
  base-learners. \emph{Computational Statistics \& Data Analysis}, 
  \bold{53}(2), 298--311.

}
\seealso{\code{\link{mboost}} for the generic boosting function and
  \code{\link{glmboost}} for boosted linear models and
  \code{\link{blackboost}} for boosted trees. See e.g. \code{\link{bbs}}
  for possible base-learners. See \code{\link{cvrisk}} for
  cross-validated stopping iteration. Furthermore see
  \code{\link{boost_control}}, \code{\link{Family}} and  \code{\link[mboost]{methods}}}
\examples{

    ### a simple two-dimensional example: cars data
    cars.gb <- gamboost(dist ~ speed, data = cars, dfbase = 4,
                        control = boost_control(mstop = 50))
    cars.gb
    AIC(cars.gb, method = "corrected")

    ### plot fit for mstop = 1, ..., 50
    plot(dist ~ speed, data = cars)
    tmp <- sapply(1:mstop(AIC(cars.gb)), function(i)
        lines(cars$speed, predict(cars.gb[i]), col = "red"))
    lines(cars$speed, predict(smooth.spline(cars$speed, cars$dist),
                              cars$speed)$y, col = "green")

    ### artificial example: sinus transformation
    x <- sort(runif(100)) * 10
    y <- sin(x) + rnorm(length(x), sd = 0.25)
    plot(x, y)
    ### linear model
    lines(x, fitted(lm(y ~ sin(x) - 1)), col = "red")
    ### GAM
    lines(x, fitted(gamboost(y ~ x,
                    control = boost_control(mstop = 500))),
          col = "green")

}
\keyword{models}
\keyword{nonlinear}
