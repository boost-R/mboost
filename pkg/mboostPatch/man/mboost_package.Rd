\name{mboost-package}
\alias{mboost-package}
\docType{package}
\title{
  mboost: Model-Based Boosting
}
\description{
  Functional gradient descent algorithm
  (boosting) for optimizing general risk functions utilizing
  component-wise (penalised) least squares estimates or regression 
  trees as base-learners for fitting generalized linear, additive 
  and interaction models to potentially high-dimensional data.
}
\details{
\tabular{ll}{
Package: \tab mboost\cr
Type: \tab Package\cr
Version: \tab 2.1-3\cr
Date: \tab 2012-XX-XX\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
LazyData: \tab yes\cr
}

  This package is intended for modern regression modelling and stands
  in-between classical generalized linear and additive models, as for example
  implemented by \code{\link{lm}}, \code{\link{glm}}, or \code{\link[mgcv]{gam}},
  and machine-learning approaches for complex interactions models,
  most prominently represented by \code{\link[gbm]{gbm}} and 
  \code{\link[randomForest]{randomForest}}.

  All functionality in this package is based on the generic
  implementation of the optimization algorithm  (function
  \code{\link{mboost_fit}}) that allows for fitting linear, additive,
  and interaction models (and mixtures of those) in low and
  high dimensions. The response may be numeric, binary, ordered,  
  censored or count data.

  Both theory and applications are discussed by Buehlmann and Hothorn (2007).
  UseRs without a basic knowledge of boosting methods are asked
  to read this introduction before analyzing data using this package. 
  The examples presented in this paper are available as package vignette
  \code{mboost_illustrations}.

  Note that the model fitting procedures in this package DO NOT automatically
  determine an appropriate model complexity. This task is the responsibility
  of the data analyst.
}
\section{NEWS in 2.0-series}{

  Version 2.0 comes with new features, is faster and more accurate
  in some aspects. In addition, some changes to the user interface
  were necessary: Subsetting \code{mboost} objects changes the object.
  At each time, a model is associated with a number of boosting iterations
  which can be changed (increased or decreased) using the subset operator.

  The \code{center} argument in \code{\link{bols}} was renamed
  to \code{intercept}. Argument \code{z} renamed to \code{by}.

  The base-learners \code{bns} and \code{bss} are deprecated
  and replaced by \code{bbs} (which results in qualitatively the
  same models but is computationally much more attractive).

  New features include new families (for example for ordinal regression)
  and the \code{which} argument to the \code{coef} and \code{predict}
  methods for selecting interesting base-learners. Predict
  methods are much faster now.

  The memory consumption could be reduced considerably, 
  thanks to sparse matrix technology in package \code{Matrix}.
  Resampling procedures run automatically in parallel 
  if package \code{multicore} is available.

  The most important advancement is a generic implementation
  of the optimizer in function \code{\link{mboost_fit}}.

}
\author{
  Torsten Hothorn \email{Torsten.Hothorn@R-project.org},\cr
  Peter Buehlmann, Thomas Kneib, Matthias Schmid and
  Benjamin Hofner
}

\references{
  Peter Buehlmann and Torsten Hothorn (2007),
  Boosting algorithms: regularization, prediction and model fitting.
  \emph{Statistical Science}, \bold{22}(4), 477--505.

  Torsten Hothorn, Peter Buehlmann, Thomas Kneib, Mattthias Schmid and
  Benjamin Hofner (2010), Model-based Boosting 2.0. \emph{Journal of
  Machine Learning Research}, \bold{11}, 2109--2113.

  Benjamin Hofner, Andreas Mayr, Nikolay Robinzonov and Matthias Schmid
  (2012). Model-based Boosting in R: A Hands-on Tutorial Using the R
  Package mboost. \emph{Department of Statistics, Technical Report No. 120}.\cr
  \url{http://epub.ub.uni-muenchen.de/12754/}

  Available as vignette via: vignette(package = "mboost", "mboost_tutorial")

}

\seealso{
  The main fitting functions include:\cr
  \code{\link{gamboost}} for boosted (generalized) additive models,
  \code{\link{glmboost}} for boosted linear models and
  \code{\link{blackboost}} for boosted trees.\cr
  See there for more details and further links.
}

\examples{

  
  data("bodyfat")
  set.seed(290875)

  ### model conditional expectation of DEXfat given
  model <- mboost(DEXfat ~ 
      bols(age) +                 ### a linear function of age
      btree(hipcirc, waistcirc) + ### a non-linear interaction of
                                  ### hip and waist circumference
      bbs(kneebreadth),           ### a smooth function of kneebreadth
      data = bodyfat, control = boost_control(mstop = 100))

  ### bootstrap for assessing `optimal' number of boosting iterations
  cvm <- cvrisk(model, papply = lapply)

  ### restrict model to mstop(cvm)
  model[mstop(cvm), return = FALSE]
  mstop(model)

  ### plot age and kneebreadth
  layout(matrix(1:2, nc = 2))
  plot(model, which = c("age", "kneebreadth"))

  ### plot interaction of hip and waist circumference
  attach(bodyfat) 
  nd <- expand.grid(hipcirc = h <- seq(from = min(hipcirc),
                                  to = max(hipcirc),
                                  length = 100),
                    waistcirc = w <- seq(from = min(waistcirc),
                                  to = max(waistcirc),
                                  length = 100))
  plot(model, which = 2, newdata = nd)
  detach(bodyfat)

  ### customized plot
  layout(1)
  pr <- predict(model, which = "hip", newdata = nd)
  persp(x = h, y = w, z = matrix(pr, nrow = 100, ncol = 100))

}

\keyword{ package }
\keyword{ smooth }
\keyword{ nonparametric }
\keyword{ models }
