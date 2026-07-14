\name{mboost_fit}
\alias{mboost_fit}
\title{ Model-based Gradient Boosting }
\description{
  Work-horse for gradient boosting for optimizing arbitrary loss functions, 
  where component-wise models are utilized as base-learners. Usually, this function 
  is not called directly by the user.
}
\usage{
mboost_fit(blg, response, weights = rep(1, NROW(response)), offset = NULL,
           family = Gaussian(), control = boost_control(), oobweights =
           as.numeric(weights == 0))
}
\arguments{
  \item{blg}{ a list of objects of elements of class \code{blg}, as returned by 
    all base-learners.}
  \item{response}{ the response variable.}
  \item{weights}{ (optional) a numeric vector of weights to be used in 
    the fitting process.}
  \item{offset}{ a numeric vector to be used as offset (optional).}
  \item{family}{a \code{\link{Family}} object.}
  \item{control}{ a list of parameters controlling the algorithm. For
    more details see \code{\link{boost_control}}. }
  \item{oobweights}{ an additional vector of out-of-bag weights, which is
    used for the out-of-bag risk (i.e., if \code{boost_control(risk =
      "oobag")}). This argument is also used internally by
    \code{cvrisk}. }
}
\details{
  The function implements component-wise functional gradient boosting in
  a generic way. This function is the main work horse and used as back-end by
  all boosting algorithms in a unified way. Usually, this function is not
  called directly.  Note that the more convenient modelling interfaces 
  \code{\link{gamboost}}, \code{\link{glmboost}} and \code{\link{blackboost}} 
  all call \code{mboost_fit}.
  
  Basically, the algorithm is initialized with a function
  for computing the negative gradient of the loss function (via its
  \code{family} argument) and one or more base-learners (given as
  \code{blg}). Usually \code{blg} and \code{response} are computed in
  the functions \code{\link{gamboost}}, \code{\link{glmboost}},
  \code{\link{blackboost}} or \code{\link{mboost}}. See there for details 
  on the specification of base-learners.

  The algorithm minimized the in-sample empirical risk defined as
  the weighted sum (by \code{weights}) of the loss function (corresponding
  to the negative gradient) evaluated at the data.

  The structure of the model is determined by the structure
  of the base-learners. If more than one base-learner is given,
  the model is additive in these components.

  Base-learners can be specified via a formula interface
  (function \code{mboost}) or as a list of objects of class \code{bl},
  see, e.g., \code{\link{bols}}.

  \code{oobweights} is a vector used internally by \code{cvrisk}. When carrying
  out cross-validation to determine the optimal stopping iteration of a boosting
  model, the default value of \code{oobweights} (out-of-bag weights) assures
  that the cross-validated risk is computed using the same observation weights
  as those used for fitting the boosting model. It is strongly recommended to
  leave this argument unspecified.
}
\value{
  An object of class \code{mboost} with \code{\link{print}},
  \code{\link{AIC}}, \code{\link{plot}} and \code{\link{predict}}
  methods being available.
}
\seealso{\code{\link{glmboost}} for boosted linear models and
  \code{\link{blackboost}} for boosted trees. See e.g. \code{\link{bbs}}
  for possible base-learners. See \code{\link{cvrisk}} for
  cross-validated stopping iteration. Furthermore see
  \code{\link{boost_control}}, \code{\link{Family}} and
  \code{\link[mboost]{methods}}. }
\examples{
  data("bodyfat", package = "TH.data")

  ### formula interface: additive Gaussian model with
  ### a non-linear step-function in `age', a linear function in `waistcirc'
  ### and a smooth non-linear smooth function in `hipcirc'
  mod <- mboost(DEXfat ~ btree(age) + bols(waistcirc) + bbs(hipcirc),
                data = bodyfat)
  layout(matrix(1:6, nc = 3, byrow = TRUE))
  plot(mod, main = "formula")

  ### the same
  with(bodyfat,
       mod <- mboost_fit(list(btree(age), bols(waistcirc), bbs(hipcirc)),
                         response = DEXfat))
  plot(mod, main = "base-learner")
}
\keyword{models}
\keyword{nonlinear}
