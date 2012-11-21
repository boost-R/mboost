\name{baselearners}
\alias{bols}
\alias{bbs}
\alias{bspatial}
\alias{brad}
\alias{brandom}
\alias{btree}
\alias{bmono}
\alias{bmrf}
\alias{buser}
\alias{bns}
\alias{bss}
\alias{\%+\%}
\alias{\%X\%}
\alias{\%O\%}
\title{ Base-learners for Gradient Boosting }
\description{
  Base-learners for fitting base-models in the generic implementation of
  component-wise gradient boosting in function \code{mboost}.
}
\usage{
bols(..., by = NULL, index = NULL, intercept = TRUE, df = NULL,
     lambda = 0, contrasts.arg = "contr.treatment")
bbs(..., by = NULL, index = NULL, knots = 20, boundary.knots = NULL,
    degree = 3, differences = 2, df = 4, lambda = NULL, center = FALSE,
    cyclic = FALSE)
bspatial(..., df = 6)
brad(..., by = NULL, index = NULL, knots = 100, df = 4, lambda = NULL,
     covFun = stationary.cov,
     args = list(Covariance="Matern", smoothness = 1.5, theta=NULL))
brandom(..., contrasts.arg = "contr.dummy", df = 4)
btree(..., tree_controls = ctree_control(stump = TRUE,
                                         mincriterion = 0,
                                         savesplitstats = FALSE))
bmono(...,
      constraint = c("increasing", "decreasing", "convex", "concave", "none"),
      by = NULL, index = NULL, knots = 20, boundary.knots = NULL,
      degree = 3, differences = 2, df = 4, lambda = NULL,
      lambda2 = 1e6, niter=10, intercept = TRUE,
      contrasts.arg = "contr.treatment",
      boundary.constraints = FALSE,
      cons.arg = list(lambda = 1e+06, n = NULL, diff_order = NULL))
bmrf(..., by = NULL, index = NULL, bnd = NULL, df = 4, lambda = NULL,
    center = FALSE)
buser(X, K = NULL, by = NULL, index = NULL, df = 4, lambda = NULL)
bl1 \%+\% bl2
bl1 \%X\% bl2
bl1 \%O\% bl2
}
\arguments{
  \item{...}{ one or more predictor variables or one matrix or data
              frame of predictor variables. For smooth base-learners,
	      the number of predictor variables and the number of
	      columns in the data frame / matrix must be less than or
	      equal to 2. If a data frame or matrix (with at least 2
	      columns) is given to \code{bols} or \code{brandom}, it
	      is directly used as the design matrix. Especially, no
	      intercept term is added reagrdless of argument
	      \code{intercept}. If the argument has only one column, it
	      is simplified to a vector and an intercept is added or not
	      according to the argmuent \code{intercept}.}
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
  \item{knots}{ either the number of knots or a vector of the positions
    of the interior knots (for more details see below). For multiple
    predictor variables, \code{knots} may be a named list where the
    names in the list are the variable names.}
  \item{boundary.knots}{ boundary points at which to anchor the B-spline basis
    (default the range of the data). A vector (of length 2)
    for the lower and the upper boundary knot can be specified.This is
    only advised for \code{bbs(..., cyclic = TRUE)}, where the boundary
    knots specify the points at which the cyclic function should be joined. In
    analogy to \code{knots} a names list can be specified. }
  \item{degree}{ degree of the regression spline.}
  \item{differences}{ a non-negative integer, typically 1, 2 or 3. If \code{differences} =
                     \emph{k}, \emph{k}-th-order differences are used as
		     a penalty (\emph{0}-th order differences specify a ridge penalty).}
  \item{intercept}{ if \code{intercept = TRUE} an intercept is added to the design matrix
                    of a linear base-learner. If \code{intercept =
		      FALSE}, continuous covariates should be (mean-) centered.}
  \item{center}{ if \code{center = TRUE} the corresponding effect is
                 re-parameterized such that the unpenalized part of the fit is subtracted and
                 only the deviation effect is fitted. The unpenalized, parametric part has then
         to be included in separate base-learners using \code{bols} (see the examples below).}
  \item{cyclic}{ if \code{cyclic = TRUE} the fitted values coincide at the boundaries
                 (useful for cyclic covariates such as day time etc.).}
  \item{covFun}{ the covariance function (i.e. radial basis)
                 needed to compute the basis functions. Per
                 default the \code{\link[fields]{stationary.cov}} function
                 (package \code{fields}) is used. }
  \item{args}{ a named list of arguments to be passed to
               \code{cov.function}. Thus strongly dependent on the
               specified \code{cov.function}. }
  \item{contrasts.arg}{ a named list of characters suitable for input to
    the \code{\link{contrasts}} replacement function, or the contrast
    matrix itself, see \code{\link{model.matrix}}, or a single character
    string (or contrast matrix) which is then used as contrasts for all
    factors in this base-learner (with the exception of factors in
    \code{by}). See also example below for setting contrasts. Note that
    a special \code{contrasts.arg} exists in package \code{mboost},
    namely "contr.dummy". This contrast is used per default in
    \code{brandom} and can also be used in \code{bols}. It leads to a
    dummy coding as returned by \code{model.matrix(~ x - 1)} were the
    intercept is implicitly included but each factor level gets a
    seperate effect estimate (see example below).}
  \item{tree_controls}{ an object of class \code{"TreeControl"}, which can be
                  obtained using \code{\link[party]{ctree_control}}.
                  Defines hyper-parameters for the trees which are used as base-learners,
                  stumps are fitted by default.}
  \item{constraint}{ type of constraint to be used. The constraint can
    be either monotonic "increasing" (default), "decreasing" or "convex"
    or "concave". Additionally, "none" can be used to specify
    unconstrained P-splines. This is especially of interest in
    conjunction with \code{boundary.constraints = TRUE}.}
  \item{lambda2}{ penalty parameter for the (monotonicity) constraint. }
  \item{niter}{ maximum number of iterations used to compute constraint
    estimates. Increase this number if a warning is displayed. }
  \item{boundary.constraints}{ a logical indicating whether additional
    constraints on the boundaries of the spline should be applied
    (default: FALSE). This is still experimental.}
  \item{cons.arg}{ a named list with additional arguments for boundary
    constraints. The element \code{lambda} specifies the
    penalty parameter that is used for the additional boundary
    constraint. The element \code{n} specifies the number of knots to be
    subject to the constraint and  can be either a scalar (use same
    number of constrained knots on each side) or a vector. Per default
    10\% of the knots on each side are used. The element
    \code{diff_order} can be used to specify the order of the boundary
    penalty: 1 (constant; default for monotonically constrained effects)
    or 2 (linear; default for all other effects).}
  \item{bnd}{
    Object of class \code{bnd}, in which the boundaries of a map are defined and from which neighborhood
    relations can be construed. See \code{\link[BayesX]{read.bnd}}. If a boundary object is not available,
    the neighborhood matrix can also be given directly. }
  \item{X}{ design matrix as it should be used in the penalized least
    squares estimation. Effect modifiers do not need to be included here
    (\code{by} can be used for convenience). }
  \item{K}{ penalty matrix as it should be used in the penalized least
    squares estimation. If \code{NULL} (default), unpenalized estimation
    is used. }
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
  \code{intercept = FALSE}). It is \emph{strongly} advised to (mean-)
  center continuous covariates, if no intercept is used in \code{bols} (see Hofner
  et al., 2011). When \code{df} is given, a ridge estimator with
  \code{df} degrees of freedom (trace of hat matrix) is used as
  base-learner. Note that all variables are treated as a group, i.e.,
  they enter the model together if the corresponding base-learner is
  selected. For ordinal variables, a ridge penalty for the differences
  of the adjacent categories (Gertheiss and Tutz 2009, Hofner et al.
  2011) is applied.

  With \code{bbs}, the P-spline approach of Eilers and Marx (1996) is
  used. P-splines use a squared \emph{k}-th-order difference penalty
  which can be interpreted as an approximation of the integrated squared
  \emph{k}-th derivative of the spline. In \code{bbs} the argument
  \code{knots} specifies either the number of (equidistant) \emph{interior}
  knots to be used for the regression spline fit or a vector including
  the positions of the \emph{interior} knots. Additionally,
  \code{boundary.knots} can be specified. However, this is only advised
  if one uses cyclic constraints, where the \code{boundary.knots}
  specify the points where the function is joined (e.g.,
  \code{boundary.knots = c(0, 2 * pi)} for angles as in a sine function
  or \code{boundary.knots = c(0, 24)} for hours during the day).

  \code{bspatial} implements bivariate tensor product P-splines for the
  estimation of either spatial effects or interaction surfaces. Note
  that \code{bspatial(x, y)} is equivalent to \code{bbs(x, y, df = 6)}.
  For possible arguments and defaults see there. The penalty term is
  constructed based on bivariate extensions of the univariate penalties
  in \code{x} and \code{y} directions, see Kneib, Hothorn and Tutz
  (2009) for details. Note that the dimensions of the penalty matrix
  increase (quickly) with the number of knots with strong impact on
  computational time. Thus, both should not be chosen to large.
  Different knots for \code{x} and \code{y} can be specified by a named
  list.

  \code{brandom(x)} specifies a random effects base-learner based on a
  factor variable \code{x} that defines the grouping structure of the
  data set. For each level of \code{x}, a separate random intercept is
  fitted, where the random effects variance is governed by the
  specification of the degrees of freedom \code{df}. Note that
  \code{brandom(...)} is essentially a wrapper to \code{bols(..., df =
  4)}, i.e., a wrapper that utilizes ridge-penalized categorical
  effects. For possible arguments and defaults see \code{bols}.

  For all linear base-learners the amount of smoothing is determined by
  the trace of the hat matrix, as indicated by \code{df}.

  If \code{by} is specified as an additional argument, a varying
  coefficients term is estimated, where \code{by} is the interaction
  variable and the effect modifier is given by either \code{x} or
  \code{x} and \code{y} (specified via \code{...}). If \code{bbs} is
  used, this corresponds to the classical situation of varying
  coefficients, where the effect of \code{by} varies over the co-domain
  of \code{x}. In case of \code{bspatial} as base-learner, the effect of
  \code{by} varies with respect to both \code{x} and \code{y}, i.e. an
  interaction surface between \code{x} and \code{y} is specified as
  effect modifier. For \code{brandom} specification of \code{by} leads
  to the estimation of random slopes for covariate \code{by} with
  grouping structure defined by factor \code{x} instead of a simple
  random intercept. In \code{bbs}, \code{bspatial} and \code{brandom}
  the computation of the smoothing parameter \code{lambda} for given
  \code{df}, or vice versa, might become (numerically) instable if the
  values of the interaction variable \code{by} become too large. In this
  case, we recommend to rescale the interaction covariate e.g. by
  dividing by \code{max(abs(by))}. If \code{bbs} or \code{bspatial} is
  specified with an factor variable \code{by} with more than two
  factors, the degrees of freedom are shared for the complete
  base-learner (i.e., spread over all factor levels). Note that the null
  space (see next paragraph) increases, as a separate null space for
  each factor level is present. Thus, the minimum degrees of freedom
  increase with increasing number of levels of \code{by} (if
  \code{center = FALSE}).

  For \code{bbs} and \code{bspatial}, option \code{center} requests that
  the fitted effect is centered around its parametric, unpenalized part
  (the so called null space). For example, with second order difference
  penalty, a linear effect of \code{x} remains unpenalized by \code{bbs}
  and therefore the degrees of freedom for the base-learner have to be
  larger than two. To avoid this restriction, option \code{center =
  TRUE} subtracts the unpenalized linear effect from the fit, allowing
  to specify any positive number as \code{df}. Note that in this case
  the linear effect \code{x} should generally be specified as an
  additional base-learner \code{bols(x)}. For \code{bspatial} and, for
  example, second order differences, a linear effect of \code{x}
  (\code{bols(x)}), a linear effect of \code{y} (\code{bols(y)}), and
  their interaction (\code{bols(x*y)}) are subtracted from the effect
  and have to be added separately to the model equation. More details on
  centering can be found in Kneib, Hothorn and Tutz (2009) and Fahrmeir,
  Kneib and Lang (2004).

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

  \code{bmrf} builds a base of a Markov random field consisting of several regions
  with a neighborhood structure. The input variable is the observed region.
  The penalty matrix is either construed from a boundary object or must be
  given directly via the option \code{bnd}. In that case the \code{dimnames} of
  the matrix have to be the region names, on the diagonal the number of neighbors
  have to be given for each region, and for each neighborhood relation the value
  in the matrix has to be -1, else 0. With a boundary object at hand, the fitted
  or predicted values can be directly plotted into the map using \code{\link[BayesX]{drawmap}}.

  \code{buser(X, K)} specifies a base-learner with user-specified design
  matrix \code{X} and penalty matrix \code{K}, where \code{X} and
  \code{K} are used to minimize a (penalized) least squares
  criterion with quadratic penalty. This can be used to easily specify
  base-learners that are not implemented (yet). See examples
  below for details how \code{buser} can be used to mimic existing
  base-learners. Note that for predictions you need to set up the
  design matrix for the new data manually.

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
  estimated (see examples for resulting design matrices).

  Three global options affect the base-learners: \code{option("mboost_useMatrix")}
  defaulting to \code{TRUE} indicates that the base-learner may use
  sparse matrix techniques for its computations. This reduces the memory
  consumption but might (for smaller sample sizes) require more computing
  time. \code{option("mboost_indexmin")} is an integer for the sample
  size required to optimize model fitting by taking ties into account.
  \code{option("mboost_dftraceS")}, which is also \code{TRUE} by default,
  indicates that the trace of the smoother matrix is used as degrees
  of freedom. If \code{FALSE}, an alternative is used (see
  Hofner et al., 2011).

  Smooth estimates with constraints can be computed using the base-learner
  \code{bmono()} which specifies P-spline base-learners with an
  additional asymmetric penalty enforcing monotonicity or
  convexity/concavity (see and Eilers, 2005). For more details in the boosting
  context and monotonic effects of ordinal factors see Hofner, Mueller and Hothorn (2011).

  Two or more linear base-learners can be joined using \code{\%+\%}. A
  tensor product of two or more linear base-learners is returned by
  \code{\%X\%}. When the design matrix can be written as the
  Kronecker product of two matrices \code{X = kronecker(X2, X1)},
  then \code{bl1 \%O\% bl2} with design matrices X1 and X2, respectively,
  can be used to efficiently compute Ridge-estimates
  following Currie, Durban, Eilers (2006). In cases the overall degrees of freedom of the
  combined base-learner increase (additive or multiplicative,
  respectively). These three features are experimental and for expert use
  only.

  \code{btree} fits a stump to one or more variables. Note that
  \code{\link{blackboost}} is more efficient for boosting stumps. For
  further references see Hothorn, Hornik, Zeileis (2006) and Hothorn et
  al. (2010).

  Note that the base-learners \code{bns} and \code{bss} are deprecated
  (and no longer available). Please use \code{bbs} instead, which
  results in qualitatively the same models but is computationally much
  more attractive.

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

  Benjamin Hofner, Torsten Hothorn, Thomas Kneib, and Matthias Schmid (2011),
  A framework for unbiased model selection based on boosting.
  \emph{Journal of Computational and Graphical Statistics}, \bold{20}, 956--971.

  Benjamin Hofner, Joerg Mueller, and Torsten Hothorn (2011),
  Monotonicity-Constrained Species Distribution Models,
  \emph{Ecology}, 92:1895-1901.

  Thomas Kneib, Torsten Hothorn and Gerhard Tutz (2009), Variable
  selection and model choice in geoadditive regression models,
  \emph{Biometrics}, \bold{65}(2), 626--634.

  Torsten Hothorn, Kurt Hornik, Achim Zeileis (2006), Unbiased recursive
  partitioning: A conditional inference framework. \emph{Journal of
    Computational and Graphical Statistics}, \bold{15}, 651--674.

  Torsten Hothorn, Peter Buehlmann, Thomas Kneib, Matthias Schmid and
  Benjamin Hofner (2010), Model-based Boosting 2.0, \emph{Journal of
    Machine Learning Research}, \bold{11}, 2109--2113.


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
                     bols(x3) + bols(x4), weights = weights)
  all.equal(coef(mod1), coef(mod2))


  ### grouped linear effects
  # center x1 and x2 first
  x1 <- scale(x1, center = TRUE, scale = FALSE)
  x2 <- scale(x2, center = TRUE, scale = FALSE)
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
  form2 <- y ~ bols(x1) + bols(x2) + bols(x1, by = x2, intercept = FALSE) +
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
  ### convenience function for plotting
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

  ## the old-fashioned way, a waste of space and time
  x <- expand.grid(x1, x2)
  modx <- mboost(vol ~ bbs(Var2, df = 3, knots = 10)\%X\%
                       bbs(Var1, df = 3, knots = 10), data = x,
                       control = boost_control(nu = 0.25))
  modx[250]

  max(abs(fitted(mod) - fitted(modx)))


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
