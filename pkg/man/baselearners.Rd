\name{baselearners}
\alias{bols}
\alias{bbs}
\alias{bns}
\alias{bss}
\alias{bspatial}
\title{ Base-learners for Gradient Boosting with Smooth Components }
\description{
  Base-learners to be utilized in the formula specification of \code{gamboost()}.
}
\usage{

bols(x, z = NULL, xname = NULL, zname = NULL)
bbs(x, z = NULL, df = 4, knots = NULL, degree = 3, differences = 2, 
    center = FALSE, xname = NULL, zname = NULL)
bns(x, z = NULL, df = 4, knots = NULL, differences = 2, 
    xname = NULL, zname = NULL)
bss(x, df = 4, xname = NULL)
bspatial(x, y, z = NULL, df = 5, xknots = NULL, yknots = NULL,
         degree = 3, differences = 2, center = FALSE, xname = NULL,
         yname = NULL, zname = NULL)

}
\arguments{
  \item{x}{ a vector containing data, either numeric or a factor. }
  \item{y}{ a vector containing data, either numeric or a factor. }
  \item{z}{ an optional vector containing data, either numeric or a factor.}
  \item{xname}{ an optional string indicating the name of the variable whose
  data values are included in the vector \code{x}.}
  \item{yname}{ an optional string indicating the name of the variable whose
  data values are included in the vector \code{y}.}
  \item{zname}{ an optional string indicating the name of the variable whose
  data values are included in the vector \code{z}.}
  \item{df}{ trace of hat matrix used for the (penalized) spline smooth. Low
  values of \code{df} correspond to a large amount of smoothing and thus to
  "weaker" baselearners. Per default, \code{df} has to be greater than 2 but in
  combination with option \code{center} any positive value is admissible.}
  \item{knots}{ either the number of (equidistant) interior knots to be used for
  the regression spline fit or a vector including the positions of the interior
  knots. If \code{knots=NULL}, the interior knots are chosen to be equidistant,
  where the number of interior knots is defined in the same way as in
  \code{smooth.spline}, cf. Hastie and Tibshirani (1990).}
  \item{xknots}{ knots in \code{x}-direction when fitting a bivariate surface
  with \code{bspatial}. See \code{knots} for details.}
  \item{yknots}{ knots in \code{y}-direction when fitting a bivariate surface
  with \code{bspatial}. See \code{knots} for details.}
  \item{degree}{ degree of the regression spline.}
  \item{differences}{ natural number between 1 and 3. If \code{differences} =
  \emph{k}, \emph{k}-th-order differences are used as a penalty. Note that
  \code{df} has to be at least as large as \code{differences}.}
  \item{center}{ If \code{center=TRUE}, the nonparametric effects is
  reparameterised such that the unpenalised part of the fit is substracted and
  only the deviation effect is fitted. The unpenalised, parametric part has then
  to be included in separate base learners using \code{bols}.}
}
\details{

  \code{bols} refers to linear base-learners (ordinary least squares fit), while
  \code{bbs}, \code{bns}, and \code{bss} refer to regression and smoothing
  splines, respectively. With \code{bbs}, the P-spline approach of Eilers and
  Marx (1996) is used. \code{bns} uses the same penalty and interior knots as
  \code{bbs}, but operates with a constrained natural spline basis instead of an
  unconstrained B-spline basis. P-splines use a \emph{k}-th-order difference
  penalty which can be interpreted as an approximation of the integrated squared
  \emph{k}-th derivative of the spline. This approximation is only appropriate
  if the knots are equidistant, so is not recommended to use non-equidistant
  knots for \code{bbs} and \code{bns}. \code{bss} refers to a smoothing spline
  based on the \code{smooth.spline} function. The amount of smoothing is
  determined by the trace of the hat matrix, as indicated by \code{df}. If
  \code{x} is a factor, an ordinary least squares fit is computed. If \code{y}
  is specified as additional argument, a bivariate tensor product penalised
  spline is fitted to estimate an interaction surface between \code{x} and
  \code{y}. If \code{z} is specified, a varying coefficients term is estimated,
  where \code{z} is the interaction variable and the effect modifier is given by
  either \code{x} or \code{x} and \code{y}. If only \code{x} is specified, this
  corresponds to the classical situation of varying coefficients, where the
  effect of \code{z} varies over the domain of \code{x}. If both \code{x} and
  \code{y} are specified, the effect of \code{z} varies with respect to both
  \code{x} and \code{y}, i.e. an interaction surface between \code{x} and
  \code{y} is specified as effect modifier.
  
  For \code{bbs} and \code{bspatial}, option \code{center} requests that the
  fitted effect is centered around its parametric, unpenalised part. For
  example, with second order difference penalty, a linear effect of \code{x}
  remains unpenalised and therefore the degrees of freedom for the base-learner
  have to be larger than 2. To avoid this restriction, option \code{center=TRUE}
  substracts the unpenalised linear effect from the fit, allowing to specify any
  positive number as \code{df}. Note that in this case the linear effect
  \code{x} should generally be specified as an additional base-learner
  \code{bols(x)}. For \code{bspatial} and, for example, second order
  differences, a linear effect of \code{x} (\code{bols(x)}), a linear effect of
  \code{y} (\code{bols(y)}), and their interaction (\code{bols(x*y)}) are
  substracted from the effect and have to be added seperately to the model
  equation. More details on centering can be found in Kneib, Hothorn and Tutz
  (2007) and Fahrmeir, Kneib and Lang (2004).  
}
\value{
 Either a matrix (in case of an ordinary least squares fit) or an object of
 class \code{basis} (in case of a regression or smoothing spline fit) with a
 \code{dpp} function as an additional attribute. The call of \code{dpp} returns
 an object of class \code{basisdpp}.                 
}
\references{

  Paul H. C. Eilers and Brian D. Marx (1996), Flexible smoothing with B-splines
  and penalties. \emph{Statistical Science}, {\bf11}(2), 89-121.

  Ludwig Fahrmeir, Thomas Kneib and Stefan Lang (2004), Penalized structured
  additive regression for space-time data: a Bayesian perspective.
  \emph{Statistica Sinica}, {\bf14}, 731-761.

  Trevor Hastie and Robert Tibshirani (1990), \emph{Generalized Additive Models}. Chapman
  and Hall.

  Thomas Kneib, Torsten Hothorn, and G. Tutz (2007), Variable selection in
  geoadditive regression based on boosting. Unpublished manuscript.

}
\examples{
x1 <- rnorm(100)
x2 <- rnorm(100) + 0.25*x1
x3 <- as.factor(sample(0:1, 100, replace = TRUE))
y <- 3*sin(x1) + x2^2 + rnorm(100)

knots.x2 <- quantile(x2, c(0.25,0.5,0.75))

spline1 <- bbs(x1,knots=20,df=4)
attributes(spline1)
spline2 <- bns(x2,knots=knots.x2,df=5)
attributes(spline2)
olsfit <- bols(x3)
attributes(olsfit)

form1 <- y ~ bbs(x1,knots=20,df=4) + bns(x2,knots=knots.x2,df=5)

# example for bspatial

x1 <- runif(250,-pi,pi)
x2 <- runif(250,-pi,pi)

y <- sin(x1)*sin(x2) + rnorm(250, sd = 0.4)

form2 <- y ~ bspatial(x1, x2, xknots=12, yknots=12)

bs1 <- gamboost(form2)

# decompose spatial effect into parametric part and deviation with 1 df

form2 <- y ~ bols(x1) + bols(x2) + bols(x1*x2) + 
             bspatial(x1, x2, xknots=12, yknots=12, center = TRUE, df=1)

bs1 <- gamboost(form2)

}
\keyword{models}