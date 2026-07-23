\name{survFit}
\alias{survFit}
\alias{survFit.mboost}
\alias{plot.survFit}
\title{ Survival Curves for a Cox Proportional Hazards Model }
\description{
  Computes the predicted survivor function for a Cox proportional hazards model.
}
\usage{
\method{survFit}{mboost}(object, newdata = NULL, ...)
\method{plot}{survFit}(x, xlab = "Time", ylab = "Probability", \dots)
}
\arguments{
  \item{object}{ an object of class \code{mboost} which is assumed to have a \code{\link{CoxPH}}
                 family component. }
  \item{newdata}{ an optional data frame in which to look for variables with
                  which to predict the survivor function. }
  \item{x}{ an object of class \code{survFit} for plotting. }
  \item{xlab}{ the label of the x axis. }
  \item{ylab}{ the label of the y axis. }
  \item{...}{ additional arguments passed to callies.}
}
\details{

  If \code{newdata = NULL}, the survivor function of the Cox proportional
  hazards model is computed for a centered linear predictor based on the
  \code{\link{blackboost}}, \code{\link{gamboost}}, or \code{\link{glmboost}}
  call. The Breslow estimator is used for computing the baseline survivor
  function. If \code{newdata} is a data frame, the \code{\link{predict}} method
  of \code{object}, along with the Breslow estimator, is used for computing the
  predicted survivor function for each row in \code{newdata}.

  If the model was fitted with case weights, the weighted Breslow
  estimator is used for computing the baseline survivor function, i.e.,
  each observation contributes to the number of events and to the risk
  sets according to its weight, and the survivor function for
  \code{newdata = NULL} is computed at the weighted mean of the linear
  predictor.

}
\value{
 An object of class \code{survFit} containing the following components:
  \item{surv}{ the estimated survival probabilities at the time points
               given in \code{time}. }
  \item{time}{ the time points at which the survivor functions are
               evaluated. }
  \item{n.event}{ the number of events observed at each time point given
                  in \code{time}. For models fitted with case weights,
                  this is the weighted number of events.}
}
\seealso{ \code{\link{gamboost}}, \code{\link{glmboost}} and
  \code{\link{blackboost}} for model fitting.}
\examples{

library("survival")
data("cancer", package = "survival")


fm <- Surv(futime,fustat) ~ age + resid.ds + rx + ecog.ps
fit <- glmboost(fm, data = ovarian, family = CoxPH(),
    control=boost_control(mstop = 500))

S1 <- survFit(fit)
S1
newdata <- ovarian[c(1,3,12),]
S2 <- survFit(fit, newdata = newdata)
S2

plot(S1)
}