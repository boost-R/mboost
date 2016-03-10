\name{mboost_intern}
\alias{mboost_intern}
\title{
  Call internal functions.
}
\description{
  Call one of the internal mboost functions.
}
\usage{
mboost_intern(..., fun = c("df2lambda", "hyper_bbs", "hyper_ols",
                  "bl_lin", "bl_lin_matrix",
                  "Complete.cases", "get_index", "isMATRIX",
                  "cbs", "bsplines", "model.frame.blg", 
                  "check_newdata", "do_trace", "rescale_weights", "nnls2D"))
}
\arguments{
  \item{\dots}{ Arguments to \code{fun}. }
  \item{fun}{ The name on an internal mboost function. }
}
\details{
  This function must not be called under any circumstances.
}
\keyword{internal}
