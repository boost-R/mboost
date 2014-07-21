\name{confint.mboost}

\alias{confint.mboost}
\alias{confint.glmboost}

\alias{plot.mboost.ci}
\alias{lines.mboost.ci}
\alias{print.glmboost.ci}

\title{
  Pointwise Bootstrap Confidence Intervals
}
\description{
  Compute and display pointwise confidence intervals
}
\usage{
\method{confint}{mboost}(object, parm = NULL, level = 0.95, B = 1000,
        B.mstop = 25, newdata = NULL, which = parm, ...)
\method{plot}{mboost.ci}(x, which, level = x$level, ylim = NULL, type = "l", col = "black",
     ci.col = rgb(170, 170, 170, alpha = 85, maxColorValue = 255),
     raw = FALSE, ...)
\method{lines}{mboost.ci}(x, which, level = x$level,
     col = rgb(170, 170, 170, alpha = 85, maxColorValue = 255),
     raw = FALSE, ...)


\method{confint}{glmboost}(object, parm = NULL, level = 0.95,
        B = 1000, B.mstop = 25, which = parm, ...)
\method{print}{glmboost.ci}(x, which = NULL, level = x$level, pe = FALSE, ...)
}

\arguments{
  \item{object}{
    a fitted model object of class \code{glmboost}, \code{gamboost} or
    \code{mboost} for which the confidence intervals should be computed.
  }
  \item{parm, which}{
    a subset of base-learners to take into account for computing
    confidence intervals. See \code{\link{mboost_methods}} for details.
    \code{parm} is just a synonyme for \code{which} to be in line with
    the generic \code{confint} function. Preferably use \code{which}.
  }
  \item{level}{
    the confidence level required.
  }
  \item{B}{
    number of outer bootstrap replicates used to compute the empirical
    bootstrap confidence intervals.
  }
  \item{newdata}{
    optionally, a data frame on which to compute the predictions for the
    confidence intervals.
  }
  \item{B.mstop}{
    number of inner bootstrap replicates used to determine the optimal
    mstop on each of the \code{B} bootstrap samples.
  }
  \item{x}{
    a confidence interval object.
  }
  \item{ylim}{
    limits of the y scale. Per default computed from the data to plot.
  }
  \item{type}{
    type of graphic for the point estimate, i.e., the predicted function.
methods  }
  \item{col}{
    color of the point estimate, i.e., the predicted function.
  }
  \item{ci.col}{
    color of the confidence interval.
  }
  \item{raw}{
    logical, should the raw function estimates or the derived confidence
    estimates be plotted?
  }
  \item{pe}{
    logical, should the point estimtate (PE) be also returned?
  }
  \item{\dots}{
    additional arguments
  }
}
\details{
  Use a nested boostrap approach to compute pointwise confidence
  intervals for the predicted partial functions or regression
  parameters.
}
\value{
  An object of class \code{glmboost.ci} or \code{mboost.ci} with special
  \code{print} and/or \code{plot} functions.
}
\references{
  %% ~put references to the literature/web site here ~
}
\author{
  Benjamin Hofner <benjamin.hofner@fau.de>
}
\seealso{
  \code{\link{cvrisk}} for crossvalidation approaches and
  \code{\link{mboost_methods}} for other methods.
}
\examples{
### a simple linear example
set.seed(1907)
data <- data.frame(x1 = rnorm(100), x2 = rnorm(100),
                   z = factor(sample(1:3, 100, replace = TRUE)))
data$y <- rnorm(100, mean = data$x1 - data$x2 - 1 * (data$z == 2) +
                            1 * (data$z == 3), sd = 0.1)
linmod <- glmboost(y ~ x1 + x2 + z, data = data,
                   control = boost_control(mstop = 200))

## compute confidence interval from 10 samples. Usually one should use
## at least 1000 samples.
CI <- confint(linmod, B = 10, level = 0.9)
CI

## to compute a confidence interval for another level simply change the
## level in the print function:
print(CI, level = 0.8)
## or print a subset (with point estimates):
print(CI, level = 0.8, pe = TRUE, which = "z")


### a simple smooth example
set.seed(1907)
data <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
data$y <- rnorm(100, mean = data$x1^2 - sin(data$x2), sd = 0.1)
gam <- gamboost(y ~ x1 + x2, data = data,
                control = boost_control(mstop = 200))

## compute confidence interval from 10 samples. Usually one should use
## at least 1000 samples.
CI_gam <- confint(gam, B = 10, level = 0.9)

par(mfrow = c(1, 2))
plot(CI_gam, which = 1)
plot(CI_gam, which = 2)
## to compute a confidence interval for another level simply change the
## level in the plot or lines function:
lines(CI_gam, which = 2, level = 0.8)
}
\keyword{methods}
