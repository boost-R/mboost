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

bbs(x, z = NULL, df = 4, knots = NULL,
degree = 3, differences = 2, center = FALSE, xname = NULL, zname = NULL)

bns(x, z = NULL, df = 4, knots = NULL, differences = 2,
xname = NULL, zname = NULL)

bss(x, df = 4, xname = NULL)

bspatial(x, y, z = NULL, df = 5, xknots = NULL, yknots = NULL,
degree = 3, differences = 2, center = FALSE, xname = NULL, yname = NULL, zname = NULL)

}
\arguments{
  \item{x}{ a vector containing data, either numeric or a factor. }
  \item{y}{ a vector containing data, either numeric or a factor. }
  \item{z}{ an optional vector containing data, either numeric or a factor.}
  \item{xname}{ an optional string indicating the name of the variable whose data values are included in
  the vector \code{x}.}
  \item{yname}{ an optional string indicating the name of the variable whose data values are included in
  the vector \code{y}.}
  \item{zname}{ an optional string indicating the name of the variable whose data values are included in
  the vector \code{z}.}
  \item{df}{ trace of hat matrix used for the (penalized) spline smooth. Low values of \code{df}
  correspond to a large amount of smoothing and thus to "weaker" baselearners. Per default, \code{df}
  has to be greater than 2.}
  \item{knots}{ either the number of (equidistant) interior knots to be used for the regression spline fit
  or a vector including the positions of the interior knots. If \code{knots=NULL}, the interior knots
  are chosen to be equidistant, where the number of interior knots is defined in the same way as in
  \code{smooth.spline}, cf. Hastie and Tibshirani (1990).}
  \item{xknots}{}
  \item{yknots}{}
  \item{degree}{ degree of the regression spline.}
  \item{differences}{ natural number between 1 and 3. If \code{differences}=\emph{k},
  \emph{k}-th-order differences are used as a penalty.}
  \item{center}{ }
}
\details{

  \code{bols} refers to linear base-learners (ordinary least squares fit), while \code{bbs},
  \code{bns}, and \code{bss} refer to regression and smoothing splines, respectively. With \code{bbs}, the P-spline
  approach of Eilers and Marx (1996) is used. \code{bns} uses the same penalty and interior knots as
  \code{bbs}, but operates with a constrained natural spline basis instead of an unconstrained B-spline basis.
  \code{bss} refers to a smoothing spline based on the \code{smooth.spline} function. The
  amount of smoothing is determined by the trace of the hat matrix, as indicated by \code{df}.
  If \code{x} and/or \code{y} are factors, an ordinary least squares fit is
  computed. If \code{z} is specified, an interaction term between \code{x} and \code{z} is fitted to the data.
  \code{bspatial} ... 
}
\value{
 Either a matrix (in case of an ordinary least squares fit) or an object of class \code{basis} (in case of a
 regression or smoothing spline fit) with a \code{dpp} function as an additional attribute. The call of
 \code{dpp} returns an object of class \code{basisdpp}.
}
\references{

Eilers, P.H.C. and Marx, B.D. (1996) Flexible smoothing with B-splines and penalties. Statistical Science,
11(2): 89-121.

Hastie, T. J. and Tibshirani, R. J. (1990) Generalized Additive Models. Chapman and Hall.
}
\examples{
x1 <- rnorm(100)
x2 <- rnorm(100) + 0.25*x1
x3 <- as.factor(sample(0:1,100,replace=T))
y <- 3*sin(x1) + x2^2 + rnorm(100)

knots.x2 <- quantile(x2,c(0.25,0.5,0.75))

spline1 <- bbs(x1,knots=20,df=4)
attributes(spline1)
spline2 <- bns(x2,knots=knots.x2,df=5)
attributes(spline2)
olsfit <- bols(x3)
attributes(olsfit)

form1 <- y ~ bbs(x1,knots=20,df=4) + bns(x2,knots=knots.x2,df=5)
}
\keyword{models}