\name{plot}
\alias{plot}

\alias{plot.glmboost}
\alias{plot.mboost}
\alias{lines.mboost}

\title{
  Plot effect estimates of boosting models
}
\description{
  Plot coefficient plots for \code{glmboost} models and partial effect
  plots for all other \code{mboost} models.
}
\usage{

\method{plot}{glmboost}(x, main = deparse(x$call), col = NULL,
     off2int = FALSE, ...)

\method{plot}{mboost}(x, which = NULL, newdata = NULL,
     type = "b", rug = TRUE, rugcol = "black",
     ylim = NULL, xlab = NULL, ylab = expression(f[partial]),
     add = FALSE, ...)

\method{lines}{mboost}(x, which = NULL, type = "l", rug = FALSE, ...)

}

\arguments{
  \item{x}{
    object of class \code{glmboost} or an object inheriting from
    \code{mboost} for plotting.
  }
  \item{main}{
    a title for the plot.
  }
  \item{col}{
    (a vector of) colors for plotting the lines representing the
    coefficient paths.
  }
  \item{off2int}{
    logical indicating whether the offset should be added to the
    intercept (if there is any) or if the offset is neglected for
    plotting (default).
  }
  \item{which}{
    a subset of base-learners used for plotting. If \code{which} is
    given (as an integer vector or characters corresponding
    to base-learners) only the corresponding partial effect plots are
    depicted. Per default all selected base-learners are plotted.
  }
  \item{newdata}{
    optionally, a data frame in which to look for variables with
    which to make predictions that are then plotted. This is especially
    useful if the data that was used to fit the model shows some larger
    gaps as effect plots are linearly interpolated between observations.
    For an example using \code{newdata} see below.
  }
  \item{type}{
    character string giving the type of plot desired. Per default,
    points and lines are plotted (\code{"b"}). Other useful options are
    points (\code{"p"}) or lines (\code{"l"}). See
    \code{\link{plot.default}} for details.
  }
  \item{rug}{
    logical. Should a rug be added to the x-axis?
  }
  \item{rugcol}{
    color for the rug.
  }
  \item{ylim}{
    the y limits of the plot.
  }
  \item{xlab}{
    a label for the x axis.
  }
  \item{ylab}{
    a label for the y axis.
  }
  \item{add}{
    logical. Should the plot be added to the previous plot?
  }
  \item{\dots}{
    Additional arguments to the \code{plot} functions. E.g. one can
    specify the x limits \code{xlim} or the color of the plot using
    \code{col}.
  }

}
\details{
  The coefficient paths for \code{glmboost} models show how the
  coefficient estimates evolve with increasing \code{mstop}. Each line
  represents one parameter estimate. Parameter estimates are only
  depicted when they they are selected at any time in the boosting
  model. Parameters that are not selected are droped from the figure
  (see example).

  Models specified with \code{gamboost} or \code{mboost} are plotted as
  partial effects. Only the effect of the current bossting iteration is
  depicted instead of the coefficient paths as for linear models. The
  function \code{lines} is just a wrapper to \code{plot(\ldots , add =
    TRUE)} where per default the effect is plotted as line and the
  \code{rug} is set to \code{FALSE}.

  Spatial effects can be also plotted using the function \code{plot}
  for mboost models (using \code{lattice} graphics). More complex
  effects reuquire manual plotting: One needs to predict the effects on
  a disired grid and plot the effect estimates.
}
\value{
  A plot of the fitted model.
}
\references{
  Benjamin Hofner, Andreas Mayr, Nikolay Robinzonov and Matthias Schmid
  (2014). Model-based Boosting in R: A Hands-on Tutorial Using the R
  Package mboost. \emph{Computational Statistics}, \bold{29}, 3--35.\cr
  \url{http://dx.doi.org/10.1007/s00180-012-0382-5}
}
\seealso{
  \code{\link{mboost_methods}} for further methods.
}
\examples{

### a simple example: cars data with one random variable
set.seed(1234)
cars$z <- rnorm(50)

########################################
## Plot linear models
########################################

## fit a linear model
cars.lm <- glmboost(dist ~ speed + z, data = cars)

## plot coefficient paths of glmboost
par(mfrow = c(3, 1), mar = c(4, 4, 4, 8))
plot(cars.lm,
     main = "Coefficient paths (offset not included)")
plot(cars.lm, off2int = TRUE,
     main = "Coefficient paths (offset included in intercept)")

## plot coefficient paths only for the first 15 steps,
## i.e., bevore z is selected
mstop(cars.lm) <- 15
plot(cars.lm, off2int = TRUE, main = "z is not yet selected")


########################################
## Plot additive models; basics
########################################

## fit an additive model
cars.gam <- gamboost(dist ~ speed + z, data = cars)

## plot effects
par(mfrow = c(1, 2), mar = c(4, 4, 0.1, 0.1))
plot(cars.gam)

## use same y-lims
plot(cars.gam, ylim = c(-50, 50))

## plot only the effect of speed
plot(cars.gam, which = "speed")
## as partial matching is used we could also use
plot(cars.gam, which = "sp")


########################################
## More complex plots
########################################

## Let us use more boosting iterations and compare the effects.

## We change the plot type and plot both effects in one figure:
par(mfrow = c(1, 1), mar = c(4, 4, 4, 0.1))
mstop(cars.gam) <- 100
plot(cars.gam, which = 1, col = "red", type = "l", rug = FALSE,
     main = "Compare effect for various models")

## Now the same model with 1000 iterations
mstop(cars.gam) <- 1000
lines(cars.gam, which = 1, col = "grey", lty = "dotted")

## There are some gaps in the data. Use newdata to get a smoother curve:
newdata <- data.frame(speed = seq(min(cars$speed), max(cars$speed),
                                  length = 200))
lines(cars.gam, which = 1, col = "grey", lty = "dashed",
      newdata = newdata)

## The model with 1000 steps seems to overfit the data.
## Usually one should use e.g. cross-validation to tune the model.

## Finally we refit the model using linear effects as comparison
cars.glm <- gamboost(dist ~ speed + z, baselearner = bols, data = cars)
lines(cars.glm, which = 1, col = "black")
## We see that all effects are more or less linear.

## Add a legend
legend("topleft", title = "Model",
       legend = c("... with mstop = 100", "... with mstop = 1000",
         "... with linear effects"),
       lty = c("solid", "dashed", "solid"),
       col = c("red", "grey", "black"))

}
\keyword{ methods }
