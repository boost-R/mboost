\name{blackboost}
\alias{blackboost}
\alias{blackboost.formula}
\alias{blackboost.matrix}
\alias{blackboost_fit}
\title{ Gradient Boosting with Regression Trees }
\description{
  Gradient boosting for optimizing arbitrary loss functions where regression
  trees are utilized as base learners.
}
\usage{
\method{blackboost}{formula}(formula, data = list(), weights = NULL, ...)
\method{blackboost}{matrix}(x, y, weights = NULL, ...)
blackboost_fit(object, tree_controls =
    ctree_control(teststat = "max",
                  testtype = "Teststatistic",
                  mincriterion = 0,
                  maxdepth = 2),
    fitmem = ctree_memory(object, TRUE), family = GaussReg(),
    control = boost_control(), weights = NULL)
}
\arguments{
  \item{formula}{ a symbolic description of the model to be fit. }
  \item{data}{ a data frame containing the variables in the model. }
  \item{weights}{ an optional vector of weights to be used in the fitting
                  process. }
  \item{x}{ design matrix. }
  \item{y}{ vector of responses. }
  \item{object}{ an object of class \code{boost_data}, see \code{\link{boost_dpp}}. }
  \item{tree_controls}{ an object of class \code{"\linkS4class{TreeControl}"}, which can be
                  obtained using \code{\link[party]{ctree_control}}.
                  Defines hyper-parameters for the trees which are used as base learners.
                  It is wise
                  to make sure to understand the consequences of altering any of its
                  arguments.}
  \item{fitmem}{ an object of class \code{TreeFitMemory}.}
  \item{family}{ an object of class \code{"\linkS4class{boost_family}"},
                implementing the negative gradient corresponding
                to the loss function to be optimized. By default,
                squared error loss
                for continuous responses is used.}
  \item{control}{ an object of class \code{\link{boost_control}}
                  which defines the hyper-parameters of the
                  boosting algorithm.}
  \item{\dots}{ additional arguments passed to callies. }
}
\details{

  This function implements the `classical'
  gradient boosting utilizing regression trees as base learners.
  Essentially, the same algorithm is implemented in package
  \code{\link[gbm]{gbm}}.  The
  main difference is that arbitrary loss functions to be optimized
  can be specified via the \code{family} argument to \code{blackboost} whereas
  \code{\link[gbm]{gbm}} uses hard-coded loss functions.
  Moreover, the base learners (conditional
  inference trees, see \code{\link[party]{ctree}}) are a little bit more flexible.

  The regression fit is a black box prediction machine and thus
  hardly interpretable.

  Usually, the formula based interface \code{blackboost} should be used.
  When necessary (for example for cross-validation), function
  \code{blackboost_fit} operating on objects of class \code{boost_data}
  is faster alternative.

}
\value{
  An object of class \code{blackboost} with \code{\link{print}}
  and \code{\link{predict}} methods being available.
}
\references{

  Jerome H. Friedman (2001),
  Greedy function approximation: A gradient boosting machine.
  \emph{The Annals of Statistics}, \bold{29}, 1189--1232.

  Greg Ridgeway (1999), The state of boosting.
  \emph{Computing Science and Statistics}, \bold{31},
  172--181.

  Peter Buhlmann and Torsten Hothorn (2007),
  Boosting algorithms: regularization, prediction and model fitting.
  \emph{Statistical Science}, \bold{22}(4), 477--505.

  Torsten Hothorn, Kurt Hornik and Achim Zeileis (2006). Unbiased Recursive
  Partitioning: A Conditional Inference Framework. \emph{Journal of
  Computational and Graphical Statistics}, \bold{15}(3), 651--674.

}
\seealso{\code{\link{glmboost}} for boosted linear models and
  \code{\link{gamboost}} for boosted smooth models. See
  \code{\link{cvrisk}} for cross-validated stoping iteration.
  Furthermore see \code{\link{boost_control}},  \code{\link{Family}} and
  \code{\link[mboost]{methods}}}
\examples{

    ### a simple two-dimensional example: cars data
    cars.gb <- blackboost(dist ~ speed, data = cars,
                          control = boost_control(mstop = 50))
    cars.gb

    ### plot fit
    plot(dist ~ speed, data = cars)
    lines(cars$speed, predict(cars.gb), col = "red")

}
\keyword{models}
\keyword{regression}
