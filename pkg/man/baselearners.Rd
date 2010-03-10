\name{baselearners}
\alias{bols}
\alias{bbs}
\alias{bspatial}
\alias{brad}
\alias{brandom}
\alias{btree}
\alias{bmono}
\alias{bmrf}
\alias{bns}
\alias{bss}
\alias{\%+\%}
\alias{\%X\%}
\title{ Base-learners for Gradient Boosting }
\description{
  Base-learners for fitting base-models in the generic implementation of
  component-wise gradient boosting in function \code{mboost}.
}
\usage{
bols(..., by = NULL, index = NULL, intercept = TRUE, df = NULL,
     lambda = 0, contrasts.arg = "contr.treatment")
bbs(..., by = NULL, index = NULL, knots = 20, degree = 3,
    differences = 2, df = 4, lambda = NULL, center = FALSE)
bspatial(...)
brad(..., by = NULL, index = NULL, knots = 100, df = 4, lambda = NULL,
     covFun = stationary.cov,
     args = list(Covariance="Matern", smoothness = 1.5, theta=NULL))
brandom(..., df = 4)
btree(..., tree_controls = ctree_control(stump = TRUE,
                                         mincriterion = 0,
                                         savesplitstats = FALSE))
bmono(..., constraint = c("increasing", "decreasing", "convex", "concave"),
      by = NULL, index = NULL, knots = 20, degree = 3, differences = 2,
      df = 4, lambda = NULL, lambda2 = 1e6, niter=10)
bmrf(..., by = NULL, index = NULL, bnd = NULL, df = 4, lambda = NULL,
    center = FALSE)
bl1 \%+\% bl2
bl1 \%X\% bl2
}
\arguments{
  \item{...}{ one or more predictor variables or one data frame of
              predictor variables.}
  \item{by}{ an optional variable defining varying coefficients,
             either a factor or numeric variable.
	     If \code{by} is a factor, the coding is determined by
	     the global \code{options("contrasts")} or as specified
	     "locally" for the factor (see \code{\link{contrasts}}). Per
	     default treatment coding is used. Note that the main
	     effect needs to be specified in a separate base-learner.}
  \item{index}{ a vector of integers for expanding the variables in \code{...}.
                For example, \code{bols(x, index = index)} is equal to \code{bols(x[index])},
                where \code{index} is an integer of length greater or equal to \code{length(x)}.}
  \item{df}{ trace of the hat matrix for the base-learner defining the base-learner
             complexity. Low values of \code{df} correspond to a large amount of smoothing and
             thus to "weaker" base-learners. Certain restrictions have to be kept for the
             specification of \code{df} since most of the base-learners rely on penalization
             approaches with a non-trivial null space. For example, for P-splines fitted with
             \code{bbs}, \code{df} has to be larger than the order of differences employed in
             the construction of the penalty term. However, when option \code{center=TRUE},
             the effect is centered around its unpenalized part and therefore any positive number
             is admissible for \code{df}.}
  \item{lambda}{ smoothing penalty, computed from \code{df} when \code{df} is specified.}
  \item{knots}{ either the number of (equidistant) interior knots to be used for
                the regression spline fit or a vector including the positions of the interior
                knots. For multiple predictor variables, \code{knots}
		may be a named list where the names in the list are the
		variable names.}
  \item{degree}{ degree of the regression spline.}
  \item{differences}{ 1, 2, or 3. If \code{differences} =
                     \emph{k}, \emph{k}-th-order differences are used as a penalty.}
  \item{intercept}{ if \code{intercept = TRUE} an intercept is added to the design matrix
                    of a linear base-learner. If \code{intercept =
		      FALSE}, continuous covariates should be (mean-) centered.}
  \item{center}{ if \code{center = TRUE} the corresponding effect is
                 re-parameterized such that the unpenalized part of the fit is subtracted and
                 only the deviation effect is fitted. The unpenalized, parametric part has then
         to be included in separate base-learners using \code{bols} (see the examples below).}
  \item{covFun}{ the covariance function (i.e. radial basis)
                 needed to compute the basis functions. Per
                 default the \code{\link[fields]{stationary.cov}} function
                 (package \code{fields}) is used. }
  \item{args}{ a named list of arguments to be passed to
               \code{cov.function}. Thus strongly dependend on the
           specified \code{cov.function}. }
  \item{contrasts.arg}{ a character suitable for input to the \code{\link{contrasts}} replacement function,
                        see \code{\link{model.matrix}}.}
  \item{tree_controls}{ an object of class \code{"\linkS4class{TreeControl}"}, which can be
                  obtained using \code{\link[party]{ctree_control}}.
                  Defines hyper-parameters for the trees which are used as base-learners,
                  stumps are fitted by default.}
  \item{constraint}{ type of constraint to be used. The constraint can
    be either monotonic "increasing" (default), "decreasing" or "convex"
    or "concave". }
  \item{lambda2}{ penalty parameter for the (monotonicity) constraint. }
  \item{niter}{ maximum number of iterations used to compute constraint
    estimates. Increase this number if a warning is displayed. }
  \item{bnd}{
    Object of class \code{bnd}, in which the boundaries of a map are defined and from which neighbourhood
    relations can be construed. See \code{\link[BayesX]{read.bnd}}. If a boundary object is not available,
    the neighbourhood matrix can also be given directly. }
  \item{bl1}{a linear base-learner or a list of linear base-learners.}
  \item{bl2}{a linear base-learner or a list of linear base-learners.}
}
\details{
  \code{bols} refers to linear base-learners (potentially estimated with a ridge penalty), while
  \code{bbs} provide penalized regression splines. \code{bspatial}
  fits bivariate surfaces and \code{brandom} defines random effects base-learners.
  In combination with option \code{by}, these base-learners can be turned into varying
  coefficient terms. The linear base-learners are fitted using Ridge Regression
  where the penalty parameter \code{lambda} is either computed from \code{df}
  (default for \code{bbs}, \code{bspatial}, and \code{brandom}) or specified directly
  (\code{lambda = 0} means no penalization as default for \code{bols}).

  In \code{bols(x)}, \code{x} may be a numeric vector or factor.
  Alternatively, \code{x} can be a data frame containing numeric or
  factor variables. In this case, or when multiple predictor variables
  are specified, e.g., using \code{bols(x1, x2)}, the model is
  equivalent to \code{lm(y ~ ., data = x)} or \code{lm(y ~ x1 + x2)},
  respectively. By default, an intercept term is added to the
  corresponding design matrix (which can be omitted using
  \code{intercept = FALSE}). It is \emph{strongly} adviced to (mean-)
  center continuous covariates, if no intercept is used in \code{bols} (see Hofner
  et al., 2009). When \code{df} is given, a ridge estimator with
  \code{df} degrees of freedom (trace of hat matrix) is used as
  base-learner. Note that all variables are treated as a group, i.e.,
  they enter the model together if the corresponding base-learner is
  selected.

  With \code{bbs}, the P-spline approach of Eilers and Marx (1996) is
  used. P-splines use a squared \emph{k}-th-order difference penalty
  which can be interpreted as an approximation of the integrated squared
  \emph{k}-th derivative of the spline.

  \code{bspatial} implements bivariate tensor product P-splines for the
  estimation of either spatial effects or interaction surfaces. Note
  that \code{bspatial(x, y)} is equivalent to \code{bbs(x, y)}. For
  possible arguments and defaults see there.
  The penalty term is constructed based on bivariate extensions of the
  univariate penalties in \code{x} and \code{y} directions, see Kneib,
  Hothorn and Tutz (2009) for details. Note that the dimensions of the
  penalty matrix increase (quickly) with the number of knots with strong
  impact on computational time. Thus, both should not be chosen to
  large. Different knots for \code{x} and \code{y} can be specified
  by a named list.

  \code{brandom(x)} specifies a random effects base-learner based on a
  factor variable \code{x} that defines the grouping structure of the
  data set. For each level of \code{x}, a separate random intercept is
  fitted, where the random effects variance is governed by the
  specification of the degrees of freedom \code{df}.

  For all linear base-learners the amount of smoothing is determined by the
  trace of the hat matrix, as indicated by \code{df}. If \code{df} is
  specified in \code{bols} a ridge penalty with the according degrees of
  freedom is used. For ordinal variables, a ridge penalty for the
  differences of the adjacent categories (Gertheiss and Tutz 2009) is applied.

  If \code{by} is specified as an additional argument, a
  varying coefficients term is estimated, where \code{by} is the
  interaction variable and the effect modifier is given by either
  \code{x} or \code{x} and \code{y} (specified via \code{...}).
  If \code{bbs} is used, this corresponds to the
  classical situation of varying coefficients, where the effect of
  \code{by} varies over the co-domain of \code{x}. In case of \code{bspatial} as
  base-learner, the effect of \code{by} varies with respect to both
  \code{x} and \code{y}, i.e. an interaction surface between \code{x} and
  \code{y} is specified as effect modifier. For \code{brandom} specification of \code{by}
  leads to the estimation of random slopes for covariate \code{by} with grouping structure
  defined by factor \code{x} instead of a simple random intercept.

  For \code{bbs} and \code{bspatial}, option \code{center} requests that the
  fitted effect is centered around its parametric, unpenalized part. For
  example, with second order difference penalty, a linear effect of \code{x}
  remains unpenalized by \code{bbs} and therefore the degrees of freedom for the base-learner
  have to be larger than two. To avoid this restriction, option \code{center = TRUE}
  subtracts the unpenalized linear effect from the fit, allowing to specify any
  positive number as \code{df}. Note that in this case the linear effect
  \code{x} should generally be specified as an additional base-learner
  \code{bols(x)}. For \code{bspatial} and, for example, second order
  differences, a linear effect of \code{x} (\code{bols(x)}), a linear effect of
  \code{y} (\code{bols(y)}), and their interaction (\code{bols(x*y)}) are
  subtracted from the effect and have to be added separately to the model
  equation. More details on centering can be found in Kneib, Hothorn and Tutz
  (2009) and Fahrmeir, Kneib and Lang (2004).

  \code{brad(x)} specifies penalized radial basis functions as used in
  Kriging. If \code{knots} is used to specify the number of knots, the
  function \code{\link[fields]{cover.design}} is used to specify the
  location of the knots such that they minimize a geometric
  space-filling criterion. Furthermore, knots can be specified directly
  via a matrix. The \code{cov.function} allows to specify the
  radial basis functions. Per default, the flexible Matern correlation
  function is used. This is specified using \code{cov.function =
  stationary.cov} with \code{Covariance = "Matern"} specified via
  \code{args}. If an effective range \code{theta} is applicable for the
  correlation function (e.g., the Matern family) the user can specify
  this value. Per default (if \code{theta = NULL}) the effective range is
  chosen as \eqn{\theta = max(||x_i - x_j||)/c}{theta =
  max(||x_i - x_j||)/c} such that the correlation function
  \deqn{\rho(c; \theta = 1) = \varepsilon,}{rho(c; theta = 1) = epsilon,}
  where \eqn{\varepsilon = 0.001}{epsilon = 0.001}.

  \code{bmrf} builds a base of a markov random field consisting of several regions
  with a neighbourhood structure. The input variable is the observed region.
  The penalty matrix is either construed from a boundary object or must be
  given directly via the option \code{bnd}. In that case the \code{dimnames} of
  the matrix have to be the region names, on the diagonal the number of neighbours
  have to be given for each region, and for each neighbourhood relation the value
  in the matrix has to be -1, else 0. With a boundary object at hand, the fitted
  or predicted values can be directly plotted into the map using \code{\link[BayesX]{drawmap}}.

  For a categorical covariate with non-observed categories
  \code{bols(x)} and \code{brandom(x)} both assign a zero effect
  these categories. However, the non-observed categories must be
  listed in \code{levels(x)}. Thus, predictions are possible
  for new observations if they correspond to this category.

  By default, all linear base-learners include an intercept term (which can
  be removed using \code{intercept = FALSE} for \code{bols} or
  \code{center = TRUE} for \code{bbs}). In this case, an explicit global
  intercept term should be added to \code{gamboost} via \code{bols} (see
  example below). With \code{bols(x, intercept = FALSE)} with categorical
  covariate \code{x} a separate effect for each group (mean effect) is
  estimated (see examples for resulting design matrizes).

  Three global options affect the base-learners: \code{option("mboost_useMatrix")}
  defaulting to \code{TRUE} indicates that the base-learner may use
  sparse matrix techniques for its computations. This reduces the memory
  consumption but might (for smaller sample sizes) require more computing
  time. \code{option("mboost_indexmin")} is an integer for the sample
  size required to optimize model fitting by taking ties into account.
  \code{option("mboost_dftraceS")}, which is also \code{TRUE} by default,
  indicates that the trace of the smoother matrix is used as degrees
  of freedom. If \code{FALSE}, an alternative is used (see
  Hofner et al., 2009).

  Smooth estimates with constraints can be computed using the base-learner
  \code{bmono()} which specifies P-spline base-learners with an
  additional asymmetric penalty enforcing monotonicity or
  convexity/concavity (see Eilers, 2005).

  Two or more linear base-learners can be joined using \code{\%+\%}. A tensor product
  of two or more linear base-learners is returned by \code{\%X\%}.
  These two features are experimental and for expert use only.

  \code{btree} fits a stump to one or more variables. Note that
  \code{\link{blackboost}} is more efficient for boosting stumps.
}
\value{
 An object of class \code{bl} (base-learner) with a \code{dpp} function. The call of \code{dpp} returns
 an object of class \code{bm} (base-model).
}
\references{
  Paul H. C. Eilers (2005), Unimodal smoothing. \emph{Journal of
  Chemometrics}, \bold{19}, 317--328.

  Paul H. C. Eilers and Brian D. Marx (1996), Flexible smoothing with B-splines
  and penalties. \emph{Statistical Science}, \bold{11}(2), 89-121.

  Ludwig Fahrmeir, Thomas Kneib and Stefan Lang (2004), Penalized structured
  additive regression for space-time data: a Bayesian perspective.
  \emph{Statistica Sinica}, \bold{14}, 731-761.

  Jan Gertheiss and Gerhard Tutz (2009), Penalized regression with ordinal
  predictors, \emph{International Statistical Review}, \bold{77}(3), 345--365.

  Benjamin Hofner, Torsten Hothorn, Thomas Kneib, and Matthias Schmid (2009),
  A framework for unbiased model selection based on boosting.
  Technical Report Nr. 72, Institut fuer Statistik, LMU Muenchen.
  \url{http://epub.ub.uni-muenchen.de/11243/}

  Thomas Kneib, Torsten Hothorn and Gerhard Tutz (2009), Variable selection and
  model choice in geoadditive regression models, \emph{Biometrics}, \bold{65}(2),
  626--634.

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
  attributes(spline1)

  knots.x2 <- quantile(x2, c(0.25, 0.5, 0.75))
  spline2 <- bbs(x2, knots = knots.x2, df = 5)
  attributes(spline2)

  attributes(ols3 <- bols(x3))
  attributes(ols4 <- bols(x4))

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
                     bols(x3) + bols(x4))
  all.equal(coef(mod1), coef(mod2))


  ### grouped linear effects
  model <- gamboost(y ~ bols(x1, x2, intercept = FALSE) +
                        bols(x1, intercept = FALSE) +
                        bols(x2, intercept = FALSE),
                        control = boost_control(mstop = 400))
  coef(model, which = 1)   # one base-learner for x1 and x2
  coef(model, which = 2:3) # two separate base-learners for x1 and x2

  ### example for bspatial
  x1 <- runif(250,-pi,pi)
  x2 <- runif(250,-pi,pi)

  y <- sin(x1) * sin(x2) + rnorm(250, sd = 0.4)

  spline3 <- bspatial(x1, x2, knots = 12)
  attributes(spline3)

  ## specify number of knots separately
  form2 <- y ~ bspatial(x1, x2, knots = list(x1 = 12, x2 = 12))

  ## decompose spatial effect into parametric part and
  ## deviation with one df
  form2 <- y ~ bols(x1) + bols(x2) + bols(x1*x2) +
               bspatial(x1, x2, knots = 12, center = TRUE, df = 1)


  ### random intercept
  id <- factor(rep(1:10, each = 5))
  raneff <- brandom(id)
  attributes(raneff)

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
  attributes(raneff)

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
                  control = boost_control(mstop = 1000))
  nd <- data.frame(x, z)
  nd <- nd[order(x),]
  nd$pred_glm <- predict(mod_glm, newdata = nd)
  for (i in seq(along = levels(z)))
      with(nd[nd$z == i,], lines(x, pred_glm, col = z))
  mod_gam <- gamboost(y ~ bbs(x) + bbs(x, by = z),
                      control = boost_control(mstop = 1000))
  nd$pred_gam <- predict(mod_gam, newdata = nd)
  for (i in seq(along = levels(z)))
      with(nd[nd$z == i,], lines(x, pred_gam, col = z, lty = "dashed"))
  ### convinience function for plotting
  par(mfrow = c(1,3))
  plot(mod_gam)


  ### remove intercept from base-learner
  ### and add explicit intercept to the model
  tmpdata <- data.frame(x = 1:100, y = rnorm(1:100), int = rep(1, 100))
  mod <- gamboost(y ~ bols(int, intercept = FALSE) +
                      bols(x, intercept = FALSE),
                  data = tmpdata,
                  control = boost_control(mstop = 2500))
  cf <- unlist(coef(mod))
  cf[1] <- cf[1] + mod$offset
  cf
  coef(lm(y ~ x, data = tmpdata))

  ### quicker and better with (mean-) centering
  tmpdata$x_center <- tmpdata$x - mean(tmpdata$x)
  mod_center <- gamboost(y ~ bols(int, intercept = FALSE) +
                             bols(x_center, intercept = FALSE),
                         data = tmpdata,
                         control = boost_control(mstop = 500))
  cf_center <- unlist(coef(mod_center, which=1:2))
  ## due to the shift in x direction we need to subtract
  ## beta_1 * mean(x) to get the correct intercept
  cf_center[1] <- cf_center[1] + mod_center$offset -
                  cf_center[2] * mean(tmpdata$x) 
  cf_center
  coef(lm(y ~ x, data = tmpdata))

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

}
\keyword{models}
