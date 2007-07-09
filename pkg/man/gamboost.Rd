\name{gamboost}
\alias{gamboost}
\alias{gamboost.formula}
\alias{gamboost.matrix}
\alias{gamboost_fit}
\title{ Gradient Boosting with Component-wise Smoothing Splines }
\description{
  Gradient boosting for optimizing arbitrary loss functions where component-wise
  smoothing splines are utilized as base learners.
}
\usage{
\method{gamboost}{formula}(formula, data = list(), weights = NULL, ...)
\method{gamboost}{matrix}(x, y, weights = NULL, ...)
gamboost_fit(object, baselearner = c("ssp", "bsp", "ols"), 
             dfbase = 4, family = GaussReg(), 
             control = boost_control(), weights = NULL)
}
\arguments{
  \item{formula}{ a symbolic description of the model to be fit. }
  \item{data}{ a data frame containing the variables in the model. }
  \item{weights}{ an optional vector of weights to be used in the fitting
                  process. }
  \item{x}{ design matrix.}
  \item{y}{ vector of responses. }
  \item{object}{ an object of class \code{boost_data}, see \code{\link{boost_dpp}}. }
  \item{baselearner}{ an character specifying the component-wise base-learner to be used: 
      \code{ssp} means smoothing splines, \code{bsp} B-splines (see \code{\link[splines]{bs}}
      and \code{ols} means linear models. Please note that only 
      the characteristics of component-wise smoothing
      splines have been investigated theoretically and practically until now.}
  \item{dfbase}{ an integer vector giving the degrees of freedom for the smoothing 
                 spline, either globally for all variables (when its length is one)
                 or separately for each single covariate. }
  \item{family}{ an object of class \code{\link{boost_family-class}},
                implementing the negative gradient corresponding
                to the loss function to be optimized, by default, squared error loss
                for continuous responses is used.}
  \item{control}{ an object of class \code{\link{boost_control}}.}
  \item{\dots}{ additional arguments passed to callies. }
}
\details{

  A (generalized) additive model is fitted using a boosting algorithm based on component-wise
  univariate smoothing splines. The methodology is described in
  Buhlmann and Yu (2003). If \code{dfbase = 1}, a univariate 
  linear model is used as base learner (resulting in a linear partial fit
  for this variable).

  The function \code{gamboost_fit} provides access to the fitting
  procedure without data pre-processing, e.g. for cross-validation.

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
  Statistical Science, accepted.
  \url{ftp://ftp.stat.math.ethz.ch/Research-Reports/Other-Manuscripts/buhlmann/BuehlmannHothorn_Boosting-rev.pdf}

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
