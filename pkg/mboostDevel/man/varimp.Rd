\name{varimp_mboost}
\alias{varimp_mboost}

\title{ Variable Importance }
\description{
In-bag risk reduction per base learner as variable importance in boosting.
}
\usage{
varimp_mboost(object)
}
\arguments{
  \item{object}{ an object of class \code{mboost}.}
  \item{percent}{ logical indicating whether risk reduction should be returned in precent of the total reduction during boosting. }
}
\details{
For each base learner this function cumulates the in-bag risk reduction achieved over all boosting iterations the particular base learner was selected in. This quantifies their specific contribution to risk reduction which can be used to compare the importance of different base learners in the model. 
}
\value{ An object of class \code{varimp_mboost} with an available \code{plot} method. Risk reductions per base learner are ordered according to the \code{mboost} object. }
\references{
}
\author{
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
cars.gb <- gamboost(dist ~ bols(speed) + bbs(speed, center=TRUE), data = cars, control = boost_control(mstop = 50))
plot(varimp_mboost(cars.gb))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
