\name{glmboost}
\alias{glmboost}
\alias{glmboost.formula}
\alias{glmboost.matrix}
\alias{glmboost.default}
\alias{plot.glmboost}
\title{ Gradient Boosting with Component-wise Linear Models }
\description{
  Gradient boosting for optimizing arbitrary loss functions where component-wise
  linear models are utilized as base-learners.
}
\usage{
\method{glmboost}{formula}(formula, data = list(), weights = NULL,
          na.action = na.pass, contrasts.arg = NULL,
          center = FALSE, control = boost_control(), ...)
\method{glmboost}{matrix}(x, y, center = FALSE, control = boost_control(), ...)
\method{glmboost}{default}(x,  ...)
\method{plot}{glmboost}(x, main = deparse(x$call), col = NULL,
                        off2int = FALSE, ...)
}
\arguments{
  \item{formula}{ a symbolic description of the model to be fit. }
  \item{data}{ a data frame containing the variables in the model. }
  \item{weights}{ an optional vector of weights to be used in the fitting
                  process. }
  \item{contrasts.arg}{a list, whose entries are contrasts suitable for input
          to the \code{contrasts} replacement function and whose names are
          the names of columns of \code{data} containing factors.
          See \code{\link{model.matrix.default}}.}
  \item{na.action}{a function which indicates what should happen when the data
          contain \code{NA}s.}
  \item{center}{logical indicating of the predictor variables are centered before fitting.}
  \item{control}{ a list of parameters controlling the algorithm.}
  \item{x}{ design matrix or an object of class \code{glmboost} for plotting.
            Sparse matrices of class \code{Matrix} can be used as well.}
  \item{y}{ vector of responses. }
  \item{main}{ a title for the plot.}
  \item{col}{ (a vector of) colors for plotting the lines representing
             the coefficient paths.}
  \item{off2int}{ logical indicating whether the offset should be
    added to the intercept (if there is any) or if the offset is
    neglected for plotting (default).}
  \item{\dots}{ additional arguments passed to \code{\link{mboost_fit}},
    including \code{weights}, \code{offset}, \code{family} and
    \code{control}. For default values see \code{\link{mboost_fit}}.}
}
\details{

  A (generalized) linear model is fitted using a boosting algorithm based on component-wise
  univariate linear models. The fit, i.e., the regression coefficients, can be
  interpreted in the usual way. The methodology is described in
  Buehlmann and Yu (2003), Buehlmann (2006), and Buehlmann and Hothorn (2007).

}
\value{
  An object of class \code{glmboost} with \code{\link{print}}, \code{\link{coef}},
  \code{\link{AIC}} and \code{\link{predict}} methods being available.
  For inputs with longer variable names, you might want to change
  \code{par("mai")} before calling the \code{plot} method of \code{glmboost}
  objects visualizing the coefficients path.
}
\references{

  Peter Buehlmann and Bin Yu (2003),
  Boosting with the L2 loss: regression and classification.
  \emph{Journal of the American Statistical Association}, \bold{98},
  324--339.

  Peter Buehlmann (2006), Boosting for high-dimensional linear models.
  \emph{The Annals of Statistics}, \bold{34}(2), 559--583.

  Peter Buehlmann and Torsten Hothorn (2007),
  Boosting algorithms: regularization, prediction and model fitting.
  \emph{Statistical Science}, \bold{22}(4), 477--505.

}
\seealso{\code{\link{mboost}} for the generic boosting function and
  \code{\link{gamboost}} for boosted additive models and
  \code{\link{blackboost}} for boosted trees. See \code{\link{cvrisk}} for
  cross-validated stopping iteration. Furthermore see
  \code{\link{boost_control}}, \code{\link{Family}} and  \code{\link[mboost]{methods}}}
\examples{

    ### a simple two-dimensional example: cars data
    cars.gb <- glmboost(dist ~ speed, data = cars,
                        control = boost_control(mstop = 5000),
                        center = FALSE)
    cars.gb

    ### coefficients should coincide
    coef(cars.gb) + c(cars.gb$offset, 0)
    coef(lm(dist ~ speed, data = cars))

    ### plot fit
    layout(matrix(1:2, ncol = 2))
    plot(dist ~ speed, data = cars)
    lines(cars$speed, predict(cars.gb), col = "red")

    ### now we center the design matrix for
    ### much quicker "convergence"
    cars.gb_centered <- glmboost(dist ~ speed, data = cars,
                                 control = boost_control(mstop = 2000),
                                 center = TRUE)
    par(mfrow=c(1,2))
    plot(cars.gb, main="without centering")
    plot(cars.gb_centered, main="with centering")

    ### alternative loss function: absolute loss
    cars.gbl <- glmboost(dist ~ speed, data = cars,
                         control = boost_control(mstop = 5000),
                         family = Laplace())
    cars.gbl

    coef(cars.gbl) + c(cars.gbl$offset, 0)
    lines(cars$speed, predict(cars.gbl), col = "green")

    ### Huber loss with adaptive choice of delta
    cars.gbh <- glmboost(dist ~ speed, data = cars,
                         control = boost_control(mstop = 5000),
                         family = Huber())

    lines(cars$speed, predict(cars.gbh), col = "blue")
    legend("topleft", col = c("red", "green", "blue"), lty = 1,
           legend = c("Gaussian", "Laplace", "Huber"), bty = "n")

    ### plot coefficient path of glmboost
    par(mai = par("mai") * c(1, 1, 1, 2.5))
    plot(cars.gb)

    ### sparse high-dimensional example
    library("Matrix")
    n <- 100
    p <- 10000
    ptrue <- 10
    X <- Matrix(0, nrow = n, ncol = p)
    X[sample(1:(n * p), floor(n * p / 20))] <- runif(floor(n * p / 20))
    beta <- numeric(p)
    beta[sample(1:p, ptrue)] <- 10
    y <- drop(X \%*\% beta + rnorm(n, sd = 0.1))
    mod <- glmboost(y = y, x = X, center = TRUE) ### mstop needs tuning
    coef(mod, which = which(beta > 0))

}
\keyword{models}
\keyword{regression}
