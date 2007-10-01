\name{boost_family-class}
\docType{class}
\alias{boost_family-class}
\alias{show,boost_family-method}
\title{Class "boost\_family": Gradient Boosting Family}
\description{Objects of class "boost\_family" define negative gradients of 
             loss functions to be optimized.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{\link{Family}(...)}
}
\section{Slots}{
	 \describe{
    \item{\code{ngradient}:}{ a function with arguments \code{y} and \code{f} 
                              implementing the \emph{negative} gradient of 
                              the \code{loss} function. }
    \item{\code{loss}:}{ a loss function with arguments \code{y} and \code{f} (to be minimized!). }
    \item{\code{risk}:}{ a risk function with arguments \code{y}, \code{f} and \code{w},
               the weighted mean of the loss function by default. }
    \item{\code{offset}:}{ a function with argument \code{y} and \code{w} (weights)
                 for computing a \emph{scalar} offset. }  
    \item{\code{weights}:}{ a logical indicating if weights are allowed. }
    \item{\code{fW}:}{ transformation of the fit for the diagonal weights matrix for an
             approximation of the boosting hat matrix for loss functions other than
             squared error.}
      \item{check_y}{ a function for checking the class / mode of a response variable.}
    \item{\code{name}:}{ a character giving the name of the loss function 
                         for pretty printing. }
    \item{\code{charloss}:}{ a character, the deparsed loss function.}
  }
}
\examples{

    Laplace()

}
\keyword{classes}
