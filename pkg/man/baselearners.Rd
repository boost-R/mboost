\name{baselearners}
\alias{bbs}
\alias{bns}
\alias{bss}
\title{ Base-learners for Gradient Boosting with Smooth Components }
\description{
  Base-learners to be utilized in the formula specification of \code{gamboost()}.
}
\usage{

bols(x, z = NULL, xname = NULL, zname = NULL)

bbs(x, z = NULL, df = 4, knots = NULL,
degree = 3, differences = 2, center = FALSE, xname = NULL, zname = NULL)

bns(x, z = NULL, df = 4, knots = NULL, differences = 2,
xname = NULL, zname = NULL)

bss(x, df = 4, xname = NULL)

bspatial(x, y, z = NULL, df = 5, xknots = NULL, yknots = NULL,
degree = 3, differences = 2, center = FALSE, xname = NULL, yname = NULL, zname = NULL)

}
\arguments{
  \item{x,y}{ variable names contained in the \code{data} frame specified in \code{gamboost()}. }
  \item{z}{}
  \item{xname}{ }
  \item{zname}{ }
  \item{df}{ trace of hat matrix used for the (penalized) spline smooth. Low values of \code{df}
  correspond to a large amount of smoothing and thus to "weaker" baselearners. Per default, \code{df}
  has to be greater than 2.}
  \item{knots,xknots,yknots}{ Either the number of (equidistant) interior knots to be used for the spline fit
  or a vector indicating the positions of the interior knots. If \code{knots=NULL}, the interior knots are
  defined in the same way as in \code{smooth.spline}, cf. Hastie and Tibshirani (1990).}
  \item{degree}{ degree of the regression spline to be used as the base-learner.}
  \item{differences}{ natural number greater than or equal to 2. If \code{differences}=k, penalized
  k-th-order differences are used for smoothing.}
  \item{center}{ }
}
\details{

  \code{bols} refers to a linear base-learner (ordinary least squares fit), while \code{bbs},
  \code{bns}, and \code{bss} refer to regression and smoothing splines, respectively. With \code{bbs}, the P-spline
  approach of Eilers and Marx (1996) is used. \code{bns} uses the same penalty as
  \code{bbs}, but operates with a constrained natural spline basis instead of an unconstrained B-spline basis.
  \code{bss} refers to a smoothing spline based on the \code{smooth.spline} function. The
  amount of smoothing is determined by the trace of the hat matrix, as indicated by \code{df}.
  \code{bspatial} ... If \code{x} and/or \code{y} are factors, an ordinary least squares fit is
  computed.
}
\value{
 All the functions above produce an object of class "formula" which contains a symbolic model formula.
}
\references{

Eilers, P.H.C. and Marx, B.D. (1996) Flexible smoothing with B-splines and penalties. Statistical Science,
11(2): 89-121.

Hastie, T. J. and Tibshirani, R. J. (1990) Generalized Additive Models. Chapman and Hall.
}
\examples{
x1 <- rnorm(100)
x2 <- rnorm(100) + 0.25*x1
y <- 3*sin(x1) + x2^2 + rnorm(100)

knots.x2 <- quantile(x2,c(0.25,0.5,0.75))

form1 <- y ~ bbs(x1,knots=20,df=4) + bns(x2,knots=knots.x2,df=5)
}
\keyword{models}