\name{blackboost}
\alias{blackboost}
\title{ Gradient Boosting with Regression Trees }
\description{
  Gradient boosting for optimizing arbitrary loss functions where regression
  trees are utilized as base-learners.
}
\usage{
blackboost(formula, data = list(),
           tree_controls = ctree_control(
               teststat = "max",
               testtype = "Teststatistic",
               mincriterion = 0,
               maxdepth = 2),
           ...)
}
\arguments{
  \item{formula}{ a symbolic description of the model to be fit. }
  \item{data}{ a data frame containing the variables in the model. }
  \item{tree_controls}{ an object of class \code{"TreeControl"}, which can be
                  obtained using \code{\link[party]{ctree_control}}.
                  Defines hyper-parameters for the trees which are used as base-learners.
                  It is wise
                  to make sure to understand the consequences of altering any of its
                  arguments.}
  \item{\dots}{  additional arguments passed to \code{\link{mboost_fit}},
    including \code{weights}, \code{offset}, \code{family} and
    \code{control}. For default values see \code{\link{mboost_fit}}. }
}
\details{

  This function implements the `classical'
  gradient boosting utilizing regression trees as base-learners.
  Essentially, the same algorithm is implemented in package
  \code{\link[gbm]{gbm}}. The
  main difference is that arbitrary loss functions to be optimized
  can be specified via the \code{family} argument to \code{blackboost} whereas
  \code{\link[gbm]{gbm}} uses hard-coded loss functions.
  Moreover, the base-learners (conditional
  inference trees, see \code{\link[party]{ctree}}) are a little bit more flexible.

  The regression fit is a black box prediction machine and thus
  hardly interpretable.

}
\value{
  An object of class \code{mboost} with \code{\link{print}}
  and \code{\link{predict}} methods being available.
}
\references{

  Peter Buehlmann and Torsten Hothorn (2007),
  Boosting algorithms: regularization, prediction and model fitting.
  \emph{Statistical Science}, \bold{22}(4), 477--505.

  Torsten Hothorn, Kurt Hornik and Achim Zeileis (2006). Unbiased recursive
  partitioning: A conditional inference framework. \emph{Journal of
  Computational and Graphical Statistics}, \bold{15}(3), 651--674.

  Yoav Freund and Robert E. Schapire (1996),
  Experiments with a new boosting algorithm.
  In \emph{Machine Learning: Proc. Thirteenth International Conference},
  148--156.

  Jerome H. Friedman (2001),
  Greedy function approximation: A gradient boosting machine.
  \emph{The Annals of Statistics}, \bold{29}, 1189--1232.

  Greg Ridgeway (1999), The state of boosting.
  \emph{Computing Science and Statistics}, \bold{31},
  172--181.

}
\seealso{\code{\link{mboost}} for the generic boosting function and
  \code{\link{glmboost}} for boosted linear models and
  \code{\link{gamboost}} for boosted additive models. See
  \code{\link{cvrisk}} for cross-validated stopping iteration.
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
