\name{baselearners}
\alias{baselearner}
\alias{baselearners}
\alias{base-learner}
\alias{base-learners}

\title{ Base-learners for Gradient Boosting }
\description{
  Base-learners for fitting base-models in the generic implementation of
  component-wise gradient boosting in function \code{mboost}.
}
\usage{
}
\arguments{
}
\details{
  Three global options affect the base-learners:

  \describe{

    \item{\code{options("mboost_useMatrix")}}{ defaulting to \code{TRUE}
      indicates that the base-learner may use sparse matrix techniques
      for its computations. This reduces the memory consumption but
      might (for smaller sample sizes) require more computing time.}

    \item{\code{options("mboost_indexmin")}}{is an integer that
      specifies the minimum sample size needed to optimize model fitting
      by automatically taking ties into account (default = 10000).}

    \item{\code{options("mboost_dftraceS")}}{\code{FALSE} by default,
      indicating how the degrees of freedom should be computed. Per
      default \deqn{\mathrm{df}(\lambda) = \mathrm{trace}(2S -
      S^{\top}S),}{df(\lambda) = trace(2S - S'S),} with smoother matrix
      \eqn{S = X(X^{\top}X + \lambda K)^{-1} X}{S = X(X'X + \lambda
      K)^(-1) X} is used (see Hofner et al., 2011a). If \code{TRUE}, the
      trace of the smoother matrix \eqn{\mathrm{df}(\lambda) =
      \mathrm{trace}(S)} is used as degrees of freedom.

      Note that these formulae specify the relation of \code{df} and
      \code{lambda} as the smoother matrix \eqn{S} depends only on
      \eqn{\lambda} (and the (fixed) design matrix \eqn{X}, the (fixed)
      penalty matrix \eqn{K}).}

  }
}

\value{
 An object of class \code{blg} (base-learner generator) with a
 \code{dpp} function.

 The call of \code{dpp} returns an object of class
 \code{bl} (base-learner) with a \code{fit} function. The call to
 \code{fit} finally returns an object of class \code{bm} (base-model).
}
\references{

  Iain D. Currie, Maria Durban, and Paul H. C. Eilers (2006),
  Generalized linear array models with applications to
  multidimensional smoothing. \emph{Journal of the Royal
  Statistical Society, Series B--Statistical Methodology},
  \bold{68}(2), 259--280.

  Paul H. C. Eilers (2005), Unimodal smoothing. \emph{Journal of
  Chemometrics}, \bold{19}, 317--328.

  Paul H. C. Eilers and Brian D. Marx (1996), Flexible smoothing with B-splines
  and penalties. \emph{Statistical Science}, \bold{11}(2), 89-121.

  Ludwig Fahrmeir, Thomas Kneib and Stefan Lang (2004), Penalized structured
  additive regression for space-time data: a Bayesian perspective.
  \emph{Statistica Sinica}, \bold{14}, 731-761.

  Jan Gertheiss and Gerhard Tutz (2009), Penalized regression with ordinal
  predictors, \emph{International Statistical Review}, \bold{77}(3), 345--365.

  D. Goldfarb and A. Idnani (1982),  Dual and Primal-Dual Methods
  for Solving Strictly Convex Quadratic Programs.  In J. P. Hennart
  (ed.), Numerical Analysis, Springer-Verlag, Berlin, pp. 226-239.

  D. Goldfarb and A. Idnani (1983),  A numerically stable dual
  method for solving strictly convex quadratic programs.
  \emph{Mathematical Programming}, \bold{27}, 1--33.

  Benjamin Hofner, Torsten Hothorn, Thomas Kneib, and Matthias Schmid (2011a),
  A framework for unbiased model selection based on boosting.
  \emph{Journal of Computational and Graphical Statistics}, \bold{20}, 956--971.

  Benjamin Hofner, Joerg Mueller, and Torsten Hothorn (2011b),
  Monotonicity-Constrained Species Distribution Models,
  \emph{Ecology}, \bold{92}, 1895--1901.

  Benjamin Hofner, Thomas Kneib and Torsten Hothorn (2014), A Unified
  Framework of Constrained Regression. \emph{Statistics & Computing}.
  Online first. DOI:10.1007/s11222-014-9520-y.\cr
  Preliminary version: \url{http://arxiv.org/abs/1403.7118}

  Thomas Kneib, Torsten Hothorn and Gerhard Tutz (2009), Variable
  selection and model choice in geoadditive regression models,
  \emph{Biometrics}, \bold{65}(2), 626--634.

  Torsten Hothorn, Kurt Hornik, Achim Zeileis (2006), Unbiased recursive
  partitioning: A conditional inference framework. \emph{Journal of
    Computational and Graphical Statistics}, \bold{15}, 651--674.

  Torsten Hothorn, Peter Buehlmann, Thomas Kneib, Matthias Schmid and
  Benjamin Hofner (2010), Model-based Boosting 2.0, \emph{Journal of
    Machine Learning Research}, \bold{11}, 2109--2113.

  G. M. Beliakov (2000), Shape Preserving Approximation using Least Squares
  Splines, \emph{Approximation Theory and its Applications},
  \bold{16}(4), 80--98.
}

\seealso{\code{\link{mboost}}}

\examples{

  set.seed(290875)

  n <- 100
  x1 <- rnorm(n)
  x2 <- rnorm(n) + 0.25 * x1
  x3 <- as.factor(sample(0:1, 100, replace = TRUE))
  x4 <- gl(4, 25)
  y <- 3 * sin(x1) + x2^2 + rnorm(n)
  weights <- drop(rmultinom(1, n, rep.int(1, n) / n))

  ### set up base-learners
  spline1 <- bbs(x1, knots = 20, df = 4)
  extract(spline1, "design")[1:10, 1:10]
  extract(spline1, "penalty")
  knots.x2 <- quantile(x2, c(0.25, 0.5, 0.75))
  spline2 <- bbs(x2, knots = knots.x2, df = 5)
  ols3 <- bols(x3)
  extract(ols3)
  ols4 <- bols(x4)

  ### compute base-models
  drop(ols3$dpp(weights)$fit(y)$model) ## same as:
  coef(lm(y ~ x3, weights = weights))

  drop(ols4$dpp(weights)$fit(y)$model) ## same as:
  coef(lm(y ~ x4, weights = weights))

  ### fit model, component-wise
  mod1 <- mboost_fit(list(spline1, spline2, ols3, ols4), y, weights)

  ### more convenient formula interface
  mod2 <- mboost(y ~ bbs(x1, knots = 20, df = 4) +
                     bbs(x2, knots = knots.x2, df = 5) +
                     bols(x3) + bols(x4), weights = weights)
  all.equal(coef(mod1), coef(mod2))


  ### grouped linear effects
  # center x1 and x2 first
  x1 <- scale(x1, center = TRUE, scale = FALSE)
  x2 <- scale(x2, center = TRUE, scale = FALSE)
  model <- gamboost(y ~ bols(x1, x2, intercept = FALSE) +
                        bols(x1, intercept = FALSE) +
                        bols(x2, intercept = FALSE),
                        control = boost_control(mstop = 50))
  coef(model, which = 1)   # one base-learner for x1 and x2
  coef(model, which = 2:3) # two separate base-learners for x1 and x2
                           # zero because they were (not yet) selected.

  ### example for bspatial
  x1 <- runif(250,-pi,pi)
  x2 <- runif(250,-pi,pi)

  y <- sin(x1) * sin(x2) + rnorm(250, sd = 0.4)

  spline3 <- bspatial(x1, x2, knots = 12)
  Xmat <- extract(spline3, "design")
  ## 12 inner knots + 4 boundary knots = 16 knots per direction
  ## THUS: 16 * 16 = 256 columns
  dim(Xmat)
  extract(spline3, "penalty")[1:10, 1:10]

  ## specify number of knots separately
  form1 <- y ~ bspatial(x1, x2, knots = list(x1 = 12, x2 = 14))

  ## decompose spatial effect into parametric part and
  ## deviation with one df
  form2 <- y ~ bols(x1) + bols(x2) + bols(x1, by = x2, intercept = FALSE) +
               bspatial(x1, x2, knots = 12, center = TRUE, df = 1)

\donttest{############################################################
## Do not run and check these examples automatically as
## they take some time

  mod1 <- gamboost(form1)
  plot(mod1)

  mod2 <- gamboost(form2)
  ## automated plot function:
  plot(mod2)
  ## plot sum of linear and smooth effects:
  library(lattice)
  df <- expand.grid(x1 = unique(x1), x2 = unique(x2))
  df$pred <- predict(mod2, newdata = df)
  levelplot(pred ~ x1 * x2, data = df)

## End(Not run and test)
}

  ## specify radial basis function base-learner for spatial effect
  ## and use data-adaptive effective range (theta = NULL, see 'args')
  form3 <- y ~ brad(x1, x2)
  ## Now use different settings, e.g. 50 knots and theta fixed to 0.4
  ## (not really a good setting)
  form4 <- y ~ brad(x1, x2, knots = 50, args = list(theta = 0.4))

\donttest{############################################################
## Do not run and check these examples automatically as
## they take some time
  mod3 <- gamboost(form3)
  plot(mod3)
  dim(extract(mod3, what = "design", which = "brad")[[1]])
  knots <- attr(extract(mod3, what = "design", which = "brad")[[1]], "knots")

  mod4 <- gamboost(form4)
  dim(extract(mod4, what = "design", which = "brad")[[1]])
  plot(mod4)

## End(Not run and test)
}

  ### random intercept
  id <- factor(rep(1:10, each = 5))
  raneff <- brandom(id)
  extract(raneff, "design")
  extract(raneff, "penalty")

  ## random intercept with non-observed category
  set.seed(1907)
  y <- rnorm(50, mean = rep(rnorm(10), each = 5), sd = 0.1)
  plot(y ~ id)
  # category 10 not observed
  obs <- c(rep(1, 45), rep(0, 5))
  model <- gamboost(y ~ brandom(id), weights = obs)
  coef(model)
  fitted(model)[46:50] # just the grand mean as usual for
                       # random effects models


  ### random slope
  z <- runif(50)
  raneff <- brandom(id, by = z)
  extract(raneff, "design")
  extract(raneff, "penalty")

  ### specify simple interaction model (with main effect)
  n <- 210
  x <- rnorm(n)
  X <- model.matrix(~ x)
  z <- gl(3, n/3)
  Z <- model.matrix(~z)
  beta <- list(c(0,1), c(-3,4), c(2, -4))
  y <- rnorm(length(x), mean = (X * Z[,1]) \%*\% beta[[1]] +
                               (X * Z[,2]) \%*\% beta[[2]] +
                               (X * Z[,3]) \%*\% beta[[3]])
  plot(y ~ x, col = z)
  ## specify main effect and interaction
  mod_glm <- gamboost(y ~ bols(x) + bols(x, by = z),
                  control = boost_control(mstop = 100))
  nd <- data.frame(x, z)
  nd <- nd[order(x),]
  nd$pred_glm <- predict(mod_glm, newdata = nd)
  for (i in seq(along = levels(z)))
      with(nd[nd$z == i,], lines(x, pred_glm, col = z))
  mod_gam <- gamboost(y ~ bbs(x) + bbs(x, by = z, df = 8),
                      control = boost_control(mstop = 100))
  nd$pred_gam <- predict(mod_gam, newdata = nd)
  for (i in seq(along = levels(z)))
      with(nd[nd$z == i,], lines(x, pred_gam, col = z, lty = "dashed"))
  ### convenience function for plotting
  par(mfrow = c(1,3))
  plot(mod_gam)


  ### remove intercept from base-learner
  ### and add explicit intercept to the model
  tmpdata <- data.frame(x = 1:100, y = rnorm(1:100), int = rep(1, 100))
  mod <- gamboost(y ~ bols(int, intercept = FALSE) +
                      bols(x, intercept = FALSE),
                  data = tmpdata,
                  control = boost_control(mstop = 1000))
  cf <- unlist(coef(mod))
  ## add offset
  cf[1] <- cf[1] + mod$offset
  signif(cf, 3)
  signif(coef(lm(y ~ x, data = tmpdata)), 3)

  ### much quicker and better with (mean-) centering
  tmpdata$x_center <- tmpdata$x - mean(tmpdata$x)
  mod_center <- gamboost(y ~ bols(int, intercept = FALSE) +
                             bols(x_center, intercept = FALSE),
                         data = tmpdata,
                         control = boost_control(mstop = 100))
  cf_center <- unlist(coef(mod_center, which=1:2))
  ## due to the shift in x direction we need to subtract
  ## beta_1 * mean(x) to get the correct intercept
  cf_center[1] <- cf_center[1] + mod_center$offset -
                  cf_center[2] * mean(tmpdata$x)
  signif(cf_center, 3)
  signif(coef(lm(y ~ x, data = tmpdata)), 3)

\donttest{############################################################
## Do not run and check these examples automatically as
## they take some time

  ### large data set with ties
  nunique <- 100
  xindex <- sample(1:nunique, 1000000, replace = TRUE)
  x <- runif(nunique)
  y <- rnorm(length(xindex))
  w <- rep.int(1, length(xindex))

  ### brute force computations
  op <- options()
  options(mboost_indexmin = Inf, mboost_useMatrix = FALSE)
  ## data pre-processing
  b1 <- bbs(x[xindex])$dpp(w)
  ## model fitting
  c1 <- b1$fit(y)$model
  options(op)

  ### automatic search for ties, faster
  b2 <- bbs(x[xindex])$dpp(w)
  c2 <- b2$fit(y)$model

  ### manual specification of ties, even faster
  b3 <- bbs(x, index = xindex)$dpp(w)
  c3 <- b3$fit(y)$model

  all.equal(c1, c2)
  all.equal(c1, c3)

## End(Not run and test)
}

  ### cyclic P-splines
  set.seed(781)
  x <- runif(200, 0,(2*pi))
  y <- rnorm(200, mean=sin(x), sd=0.2)
  newX <- seq(0,2*pi, length=100)
  ### model without cyclic constraints
  mod <- gamboost(y ~ bbs(x, knots = 20))
  ### model with cyclic constraints
  mod_cyclic <- gamboost(y ~ bbs(x, cyclic=TRUE, knots = 20,
                                 boundary.knots=c(0, 2*pi)))
  par(mfrow = c(1,2))
  plot(x,y, main="bbs (non-cyclic)", cex=0.5)
  lines(newX, sin(newX), lty="dotted")
  lines(newX + 2 * pi, sin(newX), lty="dashed")
  lines(newX, predict(mod, data.frame(x = newX)),
        col="red", lwd = 1.5)
  lines(newX + 2 * pi, predict(mod, data.frame(x = newX)),
        col="blue", lwd=1.5)
  plot(x,y, main="bbs (cyclic)", cex=0.5)
  lines(newX, sin(newX), lty="dotted")
  lines(newX + 2 * pi, sin(newX), lty="dashed")
  lines(newX, predict(mod_cyclic, data.frame(x = newX)),
        col="red", lwd = 1.5)
  lines(newX + 2 * pi, predict(mod_cyclic, data.frame(x = newX)),
        col="blue", lwd = 1.5)

  ### use buser() to mimic p-spline base-learner:
  set.seed(1907)
  x <- rnorm(100)
  y <- rnorm(100, mean = x^2, sd = 0.1)
  mod1 <- gamboost(y ~ bbs(x))
  ## now extract design and penalty matrix
  X <- extract(bbs(x), "design")
  K <- extract(bbs(x), "penalty")
  ## use X and K in buser()
  mod2 <- gamboost(y ~ buser(X, K))
  max(abs(predict(mod1) - predict(mod2)))  # same results

  ### use buser() to mimic penalized ordinal base-learner:
  z <- as.ordered(sample(1:3, 100, replace=TRUE))
  y <- rnorm(100, mean = as.numeric(z), sd = 0.1)
  X <- extract(bols(z))
  K <- extract(bols(z), "penalty")
  index <- extract(bols(z), "index")
  mod1 <- gamboost(y ~  buser(X, K, df = 1, index = index))
  mod2 <- gamboost(y ~  bols(z, df = 1))
  max(abs(predict(mod1) - predict(mod2)))  # same results

  ### kronecker product for matrix-valued responses
  data("volcano", package = "datasets")
  layout(matrix(1:2, ncol = 2))

  ## estimate mean of image treating image as matrix
  image(volcano, main = "data")
  x1 <- 1:nrow(volcano)
  x2 <- 1:ncol(volcano)

  vol <- as.vector(volcano)
  mod <- mboost(vol ~ bbs(x1, df = 3, knots = 10)\%O\%
                      bbs(x2, df = 3, knots = 10),
                      control = boost_control(nu = 0.25))
  mod[250]

  volf <- matrix(fitted(mod), nrow = nrow(volcano))
  image(volf, main = "fitted")

\donttest{############################################################
## Do not run and check these examples automatically as
## they take some time

  ## the old-fashioned way, a waste of space and time
  x <- expand.grid(x1, x2)
  modx <- mboost(vol ~ bbs(Var2, df = 3, knots = 10) \%X\%
                       bbs(Var1, df = 3, knots = 10), data = x,
                       control = boost_control(nu = 0.25))
  modx[250]

  max(abs(fitted(mod) - fitted(modx)))

## End(Not run and test)
}

  ### setting contrasts via contrasts.arg
  x <- as.factor(sample(1:4, 100, replace = TRUE))

  ## compute base-learners with different reference categories
  BL1 <- bols(x, contrasts.arg = contr.treatment(4, base = 1)) # default
  BL2 <- bols(x, contrasts.arg = contr.treatment(4, base = 2))
  ## compute 'sum to zero contrasts' using character string
  BL3 <- bols(x, contrasts.arg = "contr.sum")

  ## extract model matrices to check if it works
  extract(BL1)
  extract(BL2)
  extract(BL3)

  ### setting contrasts using named lists in contrasts.arg
  x2 <- as.factor(sample(1:4, 100, replace = TRUE))

  BL4 <- bols(x, x2,
              contrasts.arg = list(x = contr.treatment(4, base = 2),
                                   x2 = "contr.helmert"))
  extract(BL4)

  ### using special contrast: "contr.dummy":
  BL5 <- bols(x, contrasts.arg = "contr.dummy")
  extract(BL5)
}
\keyword{models}
