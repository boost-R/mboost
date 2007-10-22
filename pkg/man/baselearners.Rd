\name{baselearners}
\alias{bols}
\alias{bbs}
\alias{bns}
\alias{bss}
\alias{bspatial}
\alias{brandom}
\title{ Base learners for Gradient Boosting with Smooth Components }
\description{
  Base learners to be utilized in the formula specification of \code{gamboost()}.
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
brandom(x, z = NULL, df = 4, xname = NULL, zname = NULL)

}
\arguments{
  \item{x}{ a vector containing data, either numeric or a factor. }
  \item{y}{ a vector containing numeric data }
  \item{z}{ an optional vector containing numeric data.}
  \item{xname}{ an optional string indicating the name of the variable whose
                data values are given by the vector \code{x}.}
  \item{yname}{ an optional string indicating the name of the variable whose
                data values are given by the vector \code{y}.}
  \item{zname}{ an optional string indicating the name of the variable whose
                data values are given by the vector \code{z}.}
  \item{df}{ trace of the hat matrix for the base learner defining the base learner 
             complexity. Low values of \code{df} correspond to a large amount of smoothing and 
             thus to "weaker" base learners. Certain restrictions have to be kept for the 
             specification of \code{df} since most of the base learners rely on penalisation 
             approaches with a non-trivial null space. For example, for p-splines fitted with 
             \code{bbs}, \code{df} has to be larger than the order of differences employed in
             the construction of the penalty term. However, when option \code{center=TRUE}, 
             the effect is centered around its unpenalised part and therefore any positive number
             is admissible for \code{df}.}
  \item{knots}{ either the number of (equidistant) interior knots to be used for
                the regression spline fit or a vector including the positions of the interior
                knots. If \code{knots=NULL}, the interior knots are chosen to be equidistant,
                where the number of interior knots is defined in the same way as in
                \code{smooth.spline}.}
  \item{xknots}{ knots in \code{x}-direction when fitting a bivariate surface
                 with \code{bspatial}. See \code{knots} for details.}
  \item{yknots}{ knots in \code{y}-direction when fitting a bivariate surface
                 with \code{bspatial}. See \code{knots} for details.}
  \item{degree}{ degree of the regression spline.}
  \item{differences}{ natural number between 1 and 3. If \code{differences} =
                     \emph{k}, \emph{k}-th-order differences are used as a penalty.}
  \item{center}{ If \code{center=TRUE}, the corresponding effect is
                reparameterised such that the unpenalized part of the fit is substracted and
                only the deviation effect is fitted. The unpenalized, parametric part has then
               to be included in separate base learners using \code{bols} (see the examples below).}
}
\details{

  \code{bols} refers to linear base learners (ordinary least squares fit), while
  \code{bbs}, \code{bns}, and \code{bss} refer to penalised regression splines, 
  penalised natural splines, and smoothing splines, respectively. \code{bspatial} 
  fits bivariate surfaces and \code{brandom} defines random effects base learners.
  In combination with option \code{z}, all base learners can be turned into varying 
  coefficient terms.


  With \code{bbs}, the P-spline approach of Eilers and Marx (1996) is used. 
  \code{bns} uses the same penalty and interior knots as \code{bbs}, 
  but operates with a constrained natural spline basis instead of an
  unconstrained B-spline basis. P-splines use a squared \emph{k}-th-order difference
  penalty which can be interpreted as an approximation of the integrated squared
  \emph{k}-th derivative of the spline. This approximation is only valid
  if the knots are equidistant, so is not recommended to use non-equidistant
  knots for \code{bbs} and \code{bns}. \code{bss} refers to a smoothing spline
  based on the \code{smooth.spline} function. 

  \code{bspatial} implements bivariate tensor product P-splines for the estimation
  of either spatial effects (if \code{x} and \code{y} correspond to coordinates)
  or interaction surfaces. The penalty term is constructed based on bivariate extensions
  of the univariate penalties in \code{x} and \code{y} directions, see Kneib, Hothorn 
  and Tutz (2007) for details.

  \code{brandom} specifies a random effects base learner based on a factor variable 
  \code{x} that defines the grouping structure of the data set. For each level of
  \code{x}, a separate random intercept is fitted, where the random effects variance is 
  governed by the specification of the degrees of freedom \code{df}.

  For all base learners except \code{bols}, the amount of smoothing is
  determined by the trace of the hat matrix, as indicated by \code{df}. If \code{z} 
  is specified as an additional argument, a varying coefficients term is estimated,
  where \code{z} is the interaction variable and the effect modifier is given by
  either \code{x} or \code{x} and \code{y}. If only \code{x} is specified and one of the
  nonparametric base learners \code{bbs}, \code{bns} or \code{bss} is used, this
  corresponds to the classical situation of varying coefficients, where the
  effect of \code{z} varies over the domain of \code{x}. In case of \code{bspatial} as 
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
  and penalties. \emph{Statistical Science}, \bold{11}(2), 89-121.

  Ludwig Fahrmeir, Thomas Kneib and Stefan Lang (2004), Penalized structured
  additive regression for space-time data: a Bayesian perspective.
  \emph{Statistica Sinica}, \bold{14}, 731-761.

  Thomas Kneib, Torsten Hothorn and Gerhard Tutz (2007), Variable selection and
  Model Choice in Geoadditive Regression Models. Unpublished manuscript.

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

}
\keyword{models}