\name{helpers}
\alias{df2lambda}
\alias{bl_lin}
\alias{hyper_bbs}

\title{ Helper functions }
\description{
  These helper functions are used internally in \pkg{mboost} and are
  exported only to allow the package \pkg{FDboost} to use them.
}
\usage{
## compute Ridge shrinkage parameter lambda from df or the other way round
df2lambda(X, df = 4, lambda = NULL, dmat = NULL, weights, XtX = NULL)

## hyper parameters for P-splines baselearner (including tensor product P-splines)
hyper_bbs(mf, vary, knots = 20, boundary.knots = NULL, degree = 3,
          differences = 2, df = 4, lambda = NULL, center = FALSE,
          cyclic = FALSE, constraint = "none", deriv = 0L)

## workhorse for fitting (ridge-penalized) baselearners
bl_lin(blg, Xfun, args)
}
\arguments{
  \item{X}{ the design matrix. }
  \item{df}{ degrees of freedom. See \code{\link{bbs}}. }
  \item{lambda}{ smoothing parameter. See \code{\link{bbs}}. }
  \item{dmat}{ penalty matrix. }
  \item{weights}{ regression weights. }
  \item{XtX}{ (weighted) crossproduct of the design matrix. }

  \item{mf}{ model frame. }
  \item{vary}{ names of variables that specify varying coefficients. See
    argument \code{by} of \code{\link{bbs}}. }
  \item{knots, boundary.knots}{ knots. See \code{\link{bbs}}. }
  \item{degree}{ degree of the regression spline. See \code{\link{bbs}}. }
  \item{differences}{ differences used in the penalty. See \code{\link{bbs}}. }
  \item{center}{ use reparameterization? See \code{\link{bbs}}. }
  \item{cyclic}{ use cyclic effects? See \code{\link{bbs}}. }
  \item{constraint}{ type of constraint. See \code{\link{bbs}}. }
  \item{deriv}{ See \code{\link{bbs}}. }

  \item{blg}{ object of class \code{"blg"} that contains the model
    frame, etc.}
  \item{Xfun}{ function to set up the model matrix given the arguments
    in \code{args}.}
  \item{args}{ arguments. E.g. the result of a call to
    \code{hyper_bbs}.}
}
\details{

  Do not call these functions directly. They are only exported to make
  the package \pkg{FDboost} happy.

}
\seealso{\code{\link{mboost}}}
\keyword{misc}
