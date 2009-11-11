\name{methods}
\alias{print.glmboost}
\alias{print.mboost}

\alias{summary.mboost}

\alias{coef.glmboost}
\alias{coef.mboost}

\alias{[.mboost}

\alias{AIC.mboost}
\alias{predict.mboost}
\alias{predict.glmboost}

\alias{mstop}
\alias{mstop.gbAIC}
\alias{mstop.mboost}
\alias{mstop.cvrisk}

\alias{fitted.mboost}
\alias{logLik.mboost}
\alias{hatvalues.mboost}
\alias{hatvalues.glmboost}

\title{ Methods for Gradient Boosting Objects }
\description{
  Methods for models fitted by boosting algorithms.
}
\usage{
\method{print}{glmboost}(x, ...)
\method{print}{mboost}(x, ...)

\method{summary}{mboost}(object, ...)

\method{coef}{glmboost}(object, which = NULL,
    aggregate = c("sum", "cumsum", "none"), ...)
\method{coef}{mboost}(object, which = NULL,
    aggregate = c("sum", "cumsum", "none"), ...)

\method{[}{mboost}(x, i, return = TRUE, ...)

\method{AIC}{mboost}(object, method = c("corrected", "classical", "gMDL"),
    df = c("trace", "actset"), ..., k = 2)

\method{mstop}{mboost}(object, ...)
\method{mstop}{gbAIC}(object, ...)
\method{mstop}{cvrisk}(object, ...)

\method{predict}{mboost}(object, newdata = NULL,
    type = c("link", "response", "class"), which = NULL,
    aggregate = c("sum", "cumsum", "none"), ...)
\method{predict}{glmboost}(object, newdata = NULL,
    type = c("link", "response", "class"), which = NULL,
    aggregate = c("sum", "cumsum", "none"), ...)

\method{fitted}{mboost}(object, ...)
\method{logLik}{mboost}(object, ...)
\method{hatvalues}{mboost}(model, ...)
\method{hatvalues}{glmboost}(model, ...)
}
\arguments{
  \item{object}{ objects of class \code{glmboost}, \code{gamboost},
    \code{blackboost} or \code{gbAIC}. }
  \item{x}{ objects of class \code{glmboost} or \code{gamboost}. }
  \item{model}{objects of class mboost}
  \item{newdata}{ optionally, a data frame in which to look for variables with
          which to predict. In case the model was fitted using the \code{matrix}
          interface to \code{\link{glmboost}}, \code{newdata} must be a \code{matrix}
          as well (an error is given otherwise).}
  \item{which}{ a subset of base-learners to take into account for computing
                predictions or coefficients. If \code{which} is given
                (as an integer vector or characters corresponding
                 to base-learners) a list is returned.}
  \item{type}{ the type of prediction required.  The default is on the scale
          of the predictors; the alternative \code{"response"} is on
          the scale of the response variable.  Thus for a
          binomial model the default predictions are of log-odds
          (probabilities on logit scale) and \code{type = "response"} gives
          the predicted probabilities.  The \code{"class"} option returns
          predicted classes.}
  \item{aggregate}{ a character specifying how to aggregate predictions
                    of single base-learners. The default returns the prediction
                    for the final number of boosting iterations. \code{"cumsum"}
                    returns a matrix with the predictions for all iterations
                    simultaneously (in columns). \code{"none"} returns
                    a matrix where the jth columns contains the predictions
                    of the base-learner of the jth boosting iteration.}
  \item{i}{ integer. Index specifying the model to extract. If \code{i}
            is smaller than the initial \code{mstop}, a subset is used.
            If \code{i} is larger than the initial \code{mstop},
            additional boosting steps are performed until step \code{i}
            is reached. See details for more information. }
  \item{return}{ a logical indicating whether the changed object is
                 returned. }
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
  shows a dense representation of the model fit and \code{summary} gives
  a more detailed representation. The function \code{coef} extracts the
  regression coefficients of a linear model fitted using the \code{\link{glmboost}}
  function or an additive model fitted using the \code{\link{gamboost}}.

  The \code{predict} function can be used to predict the status of the response variable
  for new observations whereas \code{fitted} extracts the regression fit for the observations
  in the learning sample.

  The \code{[.mstop} function can be used to enhance or restrict a given
  boosting model to the specified boosting iteration \code{i}. Note that
  in both cases the original \code{x} will be changed to reduce the
  memory footage. If the boosting model is enhanced by specifying an
  index that is larger than the initial \code{mstop}, only the missing
  \code{i - mstop} steps are fitted. If the model is restricted,
  the spare steps are not dropped, i.e. if we increase mstop again,
  these boosting steps are immediately available.

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
\note{
  The \code{[.mstop} function changes the original object, i.e.
  \code{gbmodel[10]} changes \code{gbmodel} directly!
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
  cross-validated stopping iteration.}
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

    ### enhance or restrict model
    cars.gb <- gamboost(dist ~ speed, data = cars,
                        control = boost_control(mstop = 100, trace=TRUE))
    cars.gb[10]
    cars.gb[100, return=FALSE] # no refitting required
    cars.gb[150, return=FALSE] # only iterations 101 to 150 are newly fitted

    ### coefficients for optimal number of boosting iterations
    coef(cars.gb[mstop(aic)])
    plot(cars$dist, predict(cars.gb[mstop(aic)]),
         ylim = range(cars$dist))
    abline(a = 0, b = 1)
}
\keyword{methods}
