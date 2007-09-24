\name{boost_dpp}
\alias{boost_dpp}
\title{ Data Preprocessing for Gradient Boosting }
\description{
  Formula parsing and model matrix computations for gradient boosting functions.
}
\usage{
boost_dpp(formula, data, weights = NULL,  na.action = na.omit, ...)
}
\arguments{
  \item{formula}{ a symbolic description of the model to be fit. }
  \item{data}{ a data frame containing the variables in the model. }
  \item{weights}{ an optional vector of weights to be used in the fitting
                  process. }
  \item{na.action}{a function which indicates what should happen when the data
          contain \code{NA}s.}
  \item{\dots}{ additional arguments passed to callies. }
}
\details{
  Some very basic formula parsing and model matrix computations.
}
\value{
  An object of class \code{boost_data}.
}
\keyword{datagen}
