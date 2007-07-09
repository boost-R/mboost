\name{methods}
\alias{print.glmboost}
\alias{coef.glmboost}
\alias{print.gamboost}
\alias{AIC.gb}
\alias{predict.gb}
\alias{mstop}
\alias{mstop.gbAIC}
\alias{mstop.gb}
\alias{mstop.cvrisk}
\alias{mstop.blackboost}
\alias{fitted.gb}
\alias{logLik.gb}
\title{ Methods for Gradient Boosting Objects }
\description{
  Methods for models fitted by boosting algorithms.
}
\usage{
\method{print}{glmboost}(x, ...)
\method{print}{gamboost}(x, ...)
\method{coef}{glmboost}(object, ...)
\method{AIC}{gb}(object, method = c("corrected", "classical"), ...)
\method{mstop}{gbAIC}(object, ...)
\method{mstop}{gb}(object, ...)
\method{mstop}{cvrisk}(object, ...)
\method{mstop}{blackboost}(object, ...)
\method{predict}{gb}(object, newdata = NULL, type = c("lp", "response"), ...)
\method{fitted}{gb}(object, type = c("lp", "response"), ...)
\method{logLik}{gb}(object, ...)
}
\arguments{
  \item{object}{ objects of class \code{glmboost}, \code{gamboost} or \code{gbAIC}. }
  \item{x}{ objects of class \code{glmboost} or \code{gamboost}. }
  \item{newdata}{ optionally, a data frame in which to look for variables with
          which to predict. }
  \item{type}{ a character indicating whether the fit or the response
               (classes) should be predicted in case of classification 
               problems. }
  \item{method}{ a character specifying if the corrected AIC criterion or
                 a classical (-2 logLik + 2 * df) should be computed.}
  \item{\dots}{ additional arguments passed to callies. }
}
\details{

  These functions can be used to extract details from fitted models. \code{print}
  shows a dense representation of the model fit and \code{coef} extracts the
  regression coefficients of a linear model fitted using the \code{\link{glmboost}}
  function.
  
  The \code{predict} function can be used to predict the status of the response variable
  for new observations whereas \code{fitted} extracts the regression fit for the observations
  in the learning sample.

  For (generalized) linear and additive models, the \code{AIC} function can be used
  to compute both the classical and corrected AIC (Hurvich et al., 1998, only available
  when \code{family = GaussReg()} was used), 
  which is useful for the determination
  of the optimal number of boosting iterations to be applied (which can be extracted via
  \code{mstop}).

  Note that \code{logLik} and \code{AIC} only make sense when the corresponding
  \code{\link{Family}} implements the appropriate loss function.  

}
\references{

  Clifford M. Hurvich, Jeffrey S. Simonoff and Chih-Ling Tsai (1998),
  Smoothing parameter selection in nonparametric regression using
  an improved Akaike information criterion.
  \emph{Journal of the Royal Statistical Society, Series B},
  \bold{20}(2), 271--293.

  Peter Buhlmann and Torsten Hothorn (2007),
  Boosting algorithms: regularization, prediction and model fitting.
  Statistical Science, accepted.
  \url{ftp://ftp.stat.math.ethz.ch/Research-Reports/Other-Manuscripts/buhlmann/BuehlmannHothorn_Boosting-rev.pdf}

}
\examples{

    ### a simple two-dimensional example: cars data
    cars.gb <- glmboost(dist ~ speed, data = cars, 
                        control = boost_control(mstop = 2000))
    cars.gb

    ### initial number of boosting iterations
    mstop(cars.gb)

    ### AIC criterion
    aic <- AIC(cars.gb, method = "corrected")
    aic

    ### coefficients for optimal number of boosting iterations
    coef(cars.gb[mstop(aic)])
    plot(cars$dist, predict(cars.gb[mstop(aic)]), 
         ylim = range(cars$dist))
    abline(a = 0, b = 1)

}
\keyword{methods}
