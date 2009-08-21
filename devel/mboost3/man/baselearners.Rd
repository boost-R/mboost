\name{baselearners}
\alias{bols}
\alias{bbs}
\alias{bspatial}
\alias{brandom}
\alias{btree}
\title{ Base-learners for Gradient Boosting }
\description{
  Base learners to be utilized in the formula specification of \code{mboost}.
}
\usage{
bols(..., z = NULL, index = NULL, intercept = TRUE, df = NULL, lambda = NULL,
     contrasts.arg = "contr.treatment")
bbs(..., z = NULL, index = NULL, knots = 20, degree = 3,
    differences = 2, df = 4, lambda = NULL, center = FALSE)
bspatial(...)
brandom(..., df = 4)
btree(..., tree_controls = ctree_control(stump = TRUE, mincriterion = 0))
}
\arguments{
  \item{...}{ one or more predictor variables or one data frame of
              predictor variables.}
  \item{z}{ an optional variable defining varying coefficients. 
            For \code{bols} it might be categorical or numeric,
            for \code{bbs} only numeric and binary variables are
            allowed at the moment.}
  \item{index}{ a vector of integers for expanding the variables in \code{...}.
                For example, \code{bols(x, index)} is equal to \code{bols(x[index]}.}
  \item{df}{ trace of the hat matrix for the base learner defining the base learner
             complexity. Low values of \code{df} correspond to a large amount of smoothing and
             thus to "weaker" base learners. Certain restrictions have to be kept for the
             specification of \code{df} since most of the base learners rely on penalisation
             approaches with a non-trivial null space. For example, for P-splines fitted with
             \code{bbs}, \code{df} has to be larger than the order of differences employed in
             the construction of the penalty term. However, when option \code{center=TRUE},
             the effect is centered around its unpenalized part and therefore any positive number
             is admissible for \code{df}.}
  \item{lambda}{ smoothing penality, computed from \code{df} when missing.}
  \item{knots}{ either the number of (equidistant) interior knots to be used for
                the regression spline fit or a vector including the positions of the interior
                knots. For multiple predictor variables, \code{knots} may be a named list.}
  \item{degree}{ degree of the regression spline.}
  \item{differences}{ natural number between 1 and 3. If \code{differences} =
                     \emph{k}, \emph{k}-th-order differences are used as a penalty.}
  \item{intercept}{ if \code{intercept=TRUE} an intercept is added to the design matrix
                    of a linear base-learner.}
  \item{center}{ if \code{center=TRUE} the corresponding effect is
                 re-parameterized such that the unpenalized part of the fit is substracted and
                 only the deviation effect is fitted. The unpenalized, parametric part has then
		 to be included in separate base learners using \code{bols} (see the examples below).}
  \item{contrasts.arg}{ a character suitable for input to the \code{\link{contrasts}} replacement function.}
  \item{tree_controls}{ an object of class \code{\link[party]{TreeControl}}, which can be
                  obtained using \code{\link[party]{ctree_control}}.
                  Defines hyper-parameters for the trees which are used as base learners,
                  stumps are fitted by default.}
}
\details{

  \code{bols} refers to linear base learners (ordinary least squares fit), while
  \code{bbs} refer to penalized regression splines. \code{bspatial}
  fits bivariate surfaces and \code{brandom} defines random effects base learners.
  In combination with option \code{z}, all base learners can be turned into varying
  coefficient terms.

  Linear base learners can be set up using \code{bols}. The function
  can deal with both numeric and factor variables \code{x}.
  By default, an intercept term is added to the corresponding design matrix
  (which can be omitted using \code{intercept = FALSE}). When \code{df} is
  given, a ridge estimator with \code{df} degrees of freedom (trace of hat matrix)
  is used as base learner.

  With \code{bbs}, the P-spline approach of Eilers and Marx (1996) is used.
  P-splines use a squared \emph{k}-th-order difference
  penalty which can be interpreted as an approximation of the integrated squared
  \emph{k}-th derivative of the spline. 

  \code{bspatial} implements bivariate tensor product P-splines for the estimation
  of either spatial effects or interaction surfaces. Note that \code{bspatial(x, y)}
  is equivalent to \code{bbs(x, y)}.
  The penalty term is constructed based on bivariate extensions
  of the univariate penalties in \code{x} and \code{y} directions, see Kneib, Hothorn
  and Tutz (2009) for details. Note that the dimensions of the penalty matrix
  increase (quickly) with the number of knots.
  with strong impact on computational time. Thus, both should not be choosen to large.

  \code{brandom} specifies a random effects base learner based on a factor variable
  \code{x} that defines the grouping structure of the data set. For each level of
  \code{x}, a separate random intercept is fitted, where the random effects variance is
  governed by the specification of the degrees of freedom \code{df}.

  For all base learners the amount of smoothing is determined by the
  trace of the hat matrix, as indicated by \code{df}. If \code{df} is
  specified in \code{bols} a ridge penalty with the according degrees of
  freedom is used.

  If \code{z} is specified as an additional argument, a
  varying coefficients term is estimated, where \code{z} is the
  interaction variable and the effect modifier is given by either
  \code{x} or \code{x} and \code{y}. If \code{bbs} used, this corresponds to the
  classical situation of varying coefficients, where the effect of
  \code{z} varies over the co-domain of \code{x}. In case of \code{bspatial} as
  base learner, the effect of \code{z} varies with respect to both
  \code{x} and \code{y}, i.e. an interaction surface between \code{x} and
  \code{y} is specified as effect modifier. For \code{brandom} specification of \code{z}
  leads to the estimation of random slopes for covariate \code{z} with grouping structure
  defined by factor \code{x} instead of a simple random intercept.

  For \code{bbs} and \code{bspatial}, option \code{center} requests that the
  fitted effect is centered around its parametric, unpenalized part. For
  example, with second order difference penalty, a linear effect of \code{x}
  remains unpenalized by \code{bbs} and therefore the degrees of freedom for the base learner
  have to be larger than 2. To avoid this restriction, option \code{center=TRUE}
  substracts the unpenalized linear effect from the fit, allowing to specify any
  positive number as \code{df}. Note that in this case the linear effect
  \code{x} should generally be specified as an additional base learner
  \code{bols(x)}. For \code{bspatial} and, for example, second order
  differences, a linear effect of \code{x} (\code{bols(x)}), a linear effect of
  \code{y} (\code{bols(y)}), and their interaction (\code{bols(x*y)}) are
  substracted from the effect and have to be added seperately to the model
  equation. More details on centering can be found in Kneib, Hothorn and Tutz
  (2009) and Fahrmeir, Kneib and Lang (2004).

  By default, all base learners include an intercept term (which can 
  be removed using \code{intercept = FALSE} for \code{bols} or \code{center = TRUE}
  for \code{bbs}).
  In this case, an explicit global intercept term should be added to
  \code{gamboost} via \code{bols} (see example below).

  \code{btree} fits a stump to one or two variables. Note that
  \code{\link{blackboost}} is more efficient for boosting stumps.

}
\value{
 An object of class \code{bl} with a \code{dpp} function. The call of \code{dpp} returns
 an object of class \code{bm}.
}
\references{

  Paul H. C. Eilers and Brian D. Marx (1996), Flexible smoothing with B-splines
  and penalties. \emph{Statistical Science}, \bold{11}(2), 89-121.

  Ludwig Fahrmeir, Thomas Kneib and Stefan Lang (2004), Penalized structured
  additive regression for space-time data: a Bayesian perspective.
  \emph{Statistica Sinica}, \bold{14}, 731-761.

  Thomas Kneib, Torsten Hothorn and Gerhard Tutz (2009), Variable selection and
  model choice in geoadditive regression models, \emph{Biometrics}, \bold{65}(2),
  626--634.

}
\seealso{\code{\link{gamboost}}}
\examples{

  x1 <- rnorm(100)
  x2 <- rnorm(100) + 0.25*x1
  x3 <- as.factor(sample(0:1, 100, replace = TRUE))
  x4 <- gl(4, 25)
  y <- 3*sin(x1) + x2^2 + rnorm(100)

  knots.x2 <- quantile(x2, c(0.25,0.5,0.75))

  spline1 <- bbs(x1, knots = 20, df = 4)
attributes(spline1)
spline2 <- bns(x2,knots=knots.x2,df=5)
attributes(spline2)
olsfit <- bols(x3)
attributes(olsfit)

form1 <- y ~ bbs(x1,knots=20,df=4) + bns(x2,knots=knots.x2,df=5)

# example for factors
attributes(bols(x4))

# example for bspatial

x1 <- runif(250,-pi,pi)
x2 <- runif(250,-pi,pi)

y <- sin(x1)*sin(x2) + rnorm(250, sd = 0.4)

spline3 <- bspatial(x1, x2, xknots=12, yknots=12)
attributes(spline3)

form2 <- y ~ bspatial(x1, x2, xknots=12, yknots=12)

# decompose spatial effect into parametric part and deviation with 1 df

form2 <- y ~ bols(x1) + bols(x2) + bols(x1*x2) +
             bspatial(x1, x2, xknots=12, yknots=12, center = TRUE, df=1)

# random intercept

id <- factor(rep(1:10, each=5))
raneff <- brandom(id)
attributes(raneff)

# random slope

z <- runif(50)
raneff <- brandom(id, z=z)
attributes(raneff)

# remove intercept from base learner
# and add explicit intercept to the model

tmpdata <- data.frame(x = 1:100, y = rnorm(1:100), int = rep(1, 100))
mod <- gamboost(y ~ bols(int, center = TRUE) + bols(x, center = TRUE),
                data = tmpdata, control = boost_control(mstop = 2500))
cf <- unlist(coef(mod))
cf[1] <- cf[1] + mod$offset
cf
coef(lm(y ~ x, data = tmpdata))

}
\keyword{models}