\name{methods}
\alias{print.glmboost}
\alias{coef.glmboost}
\alias{coef.gamboost}
\alias{print.gamboost}
\alias{AIC.gamboost}
\alias{AIC.glmboost}
\alias{predict.gb}
\alias{predict.blackboost}
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
\method{coef}{gamboost}(object, ...)
\method{AIC}{gamboost}(object, method = c("corrected", "classical", "gMDL"), ..., k = 2)
\method{AIC}{glmboost}(object, method = c("corrected", "classical", "gMDL"),
    df = c("trace", "actset"), ..., k = 2)
\method{mstop}{gbAIC}(object, ...)
\method{mstop}{gb}(object, ...)
\method{mstop}{cvrisk}(object, ...)
\method{mstop}{blackboost}(object, ...)
\method{predict}{gb}(object, newdata = NULL, type = c("lp", "response"),
       allIterations = FALSE, ...)
\method{predict}{blackboost}(object, newdata = NULL, type = c("lp", "response"),
       allIterations = FALSE, ...)
\method{fitted}{gb}(object, type = c("lp", "response"), ...)
\method{logLik}{gb}(object, ...)
}
\arguments{
  \item{object}{ objects of class \code{glmboost}, \code{gamboost},
    \code{blackboost} or \code{gbAIC}. }
  \item{x}{ objects of class \code{glmboost} or \code{gamboost}. }
  \item{newdata}{ optionally, a data frame in which to look for variables with
          which to predict. }
  \item{type}{ a character indicating whether the fit or the response
               (classes) should be predicted in case of classification
               problems. }
  \item{allIterations}{ computes the (linear) predictor for all boosting
               iterations simultaneously and returns a matrix.}
  \item{method}{ a character specifying if the corrected AIC criterion or
                 a classical (-2 logLik + k * df) should be computed.}
  \item{df}{ a character specifying how degrees of freedom should be computed:
             \code{trace} defines degrees of freedom by the trace of the
             boosting hat matrix and \code{actset} uses the number of
             non-zero coefficients for each boosting iteration.}
  \item{k}{  numeric, the \emph{penalty} per parameter to be used; the default
             \code{k = 2} is the classical AIC. Only used when \code{method = "classical"}.}
  \item{\dots}{ additional arguments passed to callies. }
}
\details{

  These functions can be used to extract details from fitted models. \code{print}
  shows a dense representation of the model fit and \code{coef} extracts the
  regression coefficients of a linear model fitted using the \code{\link{glmboost}}
  function (for \code{\link{gamboost}} objects, \code{coef} only works
  when \code{\link{bols}} is used as baselearner).

  The \code{predict} function can be used to predict the status of the response variable
  for new observations whereas \code{fitted} extracts the regression fit for the observations
  in the learning sample. When \code{allIterations = TRUE}, the matrix
  of all (linear) predictors for boosting iterations 1 to \code{mstop} is returned.

  For (generalized) linear and additive models, the \code{AIC} function can be used
  to compute both the classical and corrected AIC (Hurvich et al., 1998, only available
  when \code{family = GaussReg()} was used), which is useful for the determination
  of the optimal number of boosting iterations to be applied (which can be extracted via
  \code{mstop}). The degrees of freedom are either computed via the trace of the
  boosting hat matrix (which is rather slow even for moderate sample sizes)
  or the number of variables (non-zero coefficients) that entered the model so far
  (faster but only meaningful for linear models fitted via \code{\link{gamboost}}
  (see Hastie, 2007).

  In addition, the general Minimum Description Length criterion (Buhlmann and Yu, 2006)
  can be computed using function \code{AIC}.

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
  \emph{Statistical Science}, \bold{22}(4), 477--505.

  Travor Hastie (2007), Discussion of ``Boosting Algorithms: Regularization, Prediction and Model Fitting''
  by Peter Buhlmann and Torsten Hothorn. \emph{Statistical Science}, \bold{22}(4), 505.

  Peter Buhlmann and Bin Yu (2006),
  Sparse Boosting. \emph{Journal of Machine Learning Research}, \bold{7}, 1001--1024.

}
\seealso{ \code{\link{gamboost}}, \code{\link{glmboost}} and
  \code{\link{blackboost}} for model fitting. See \code{\link{cvrisk}} for
  cross-validated stoping iteration.}
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
