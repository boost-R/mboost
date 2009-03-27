\name{gamboost}
\alias{gamboost}
\alias{gamboost.formula}
\alias{gamboost.matrix}
\alias{gamboost_fit}
\alias{plot.gamboost}
\title{ Gradient Boosting with Smooth Components }
\description{
  Gradient boosting for optimizing arbitrary loss functions, where component-wise
  smoothing procedures are utilized as base learners.
}
\usage{
\method{gamboost}{formula}(formula, data = list(), weights = NULL, 
        na.action = na.omit, ...)
\method{gamboost}{matrix}(x, y, weights = NULL, ...)
gamboost_fit(object, baselearner = c("bbs", "bss", "bols", 
             "bns", "btree"), dfbase = 4, family = GaussReg(), 
             control = boost_control(), weights = NULL)
\method{plot}{gamboost}(x, which = NULL, ask = TRUE && dev.interactive(), 
    type = "b", ylab = expression(f[partial]), add_rug = TRUE, ...)
}
\arguments{
  \item{formula}{ a symbolic description of the model to be fit. }
  \item{data}{ a data frame containing the variables in the model. }
  \item{weights}{ an optional vector of weights to be used in the fitting
                  process. }
  \item{na.action}{a function which indicates what should happen when the data
          contain \code{NA}s.}
  \item{x}{ design matrix (for \code{gamboost.matrix}) or an object returned
            by \code{gamboost} to be plotted via \code{plot}.}
  \item{y}{ vector of responses. }
  \item{object}{ an object of class \code{boost_data}, see \code{\link{boost_dpp}}. }
  \item{baselearner}{ a character specifying the component-wise base learner to be used: 
      \code{\link{bss}} means smoothing splines (see Buhlmann and Yu 2003),
      \code{\link{bbs}} P-splines with a B-spline basis (see Schmid and Hothorn
      2007), \code{\link{bns}} P-splines with a natural spline basis, 
      \code{\link{bols}} linear models, \code{\link{bspatial}} bivariate tensor
      product penalized splines, and \code{\link{brandom}} random effects. 
      In addition, \code{\link{btree}} boosts stumps.
      Component-wise smoothing splines
      have been considered in Buhlmann and Yu (2003) and Schmid and Hothorn (2007) investigate P-splines with 
      a B-spline basis. Kneib, Hothorn and Tutz (2007) also utilise P-splines with a B-spline basis,
      supplement them with their bivariate tensor product version to estimate interaction surfaces and spatial effects
      and also consider random effects base learners.}
  \item{dfbase}{ an integer vector giving the degrees of freedom for the smoothing 
                 spline, either globally for all variables (when its length is one)
                 or separately for each single covariate. }
  \item{family}{ an object of class \code{\link{boost_family-class}},
                implementing the negative gradient corresponding
                to the loss function to be optimized, by default, squared error loss
                for continuous responses is used.}
  \item{control}{ an object of class \code{\link{boost_control}}.}
  \item{which}{ if a subset of the plots is required, specify a subset of the
                variables. Only selected variables are plotted by default.}
  \item{ask}{ logical; if \code{TRUE}, the user is \emph{ask}ed before each plot, see
          \code{par(ask=.)}.}
  \item{type}{what type of plot should be drawn: see \code{\link{plot}}.}
  \item{ylab}{a title for the y axis: see \code{\link{title}}.}
  \item{add_rug}{ logical; if \code{TRUE}, \code{\link{rug}}s are added.}
  \item{\dots}{ additional arguments passed to callies. }
}
\details{

  A (generalized) additive model is fitted using a boosting algorithm based on
  component-wise univariate base learners. The base learner can either be
  specified via the \code{formula} object or via the \code{baselearner} argument
  (see \code{\link{bbs}} for an example). If the base learners specified in
  \code{formula} differ from \code{baselearner}, the latter argument will be
  ignored.

  The function \code{gamboost_fit} provides access to the fitting
  procedure without data pre-processing, e.g. for cross-validation.

  Note that penalized B-splines instead of smoothing splines
  are used as default baselearners as of version 1.1-0.

}
\value{
  An object of class \code{gamboost} with \code{\link{print}},
  \code{\link{AIC}} and \code{\link{predict}} methods being available.
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
  model choice in geoadditive regression models. 
  \emph{Biometrics}, accepted.
  \url{http://epub.ub.uni-muenchen.de/2063/}

  Matthias Schmid and Torsten Hothorn (2009),
  Boosting additive models using component-wise P-splines as
  base-learners. \emph{Computational Statistics \& Data Analysis}, accepted.
  \url{http://epub.ub.uni-muenchen.de/2057/}

}
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
    lines(x, fitted(gamboost(y ~ x - 1, 
                    control = boost_control(mstop = 500))), 
          col = "green")

}
\keyword{models}
\keyword{nonlinear}
