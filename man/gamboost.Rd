\name{mboost}
\alias{mboost}
\alias{gamboost}
\title{ Gradient Boosting for Additive Models }
\description{
  Gradient boosting for optimizing arbitrary loss functions, where component-wise
  arbitrary base-learners, e.g., smoothing procedures,  are utilized as additive
  base-learners.
}
\usage{
mboost(formula, data = list(), na.action = na.omit, weights = NULL, 
       offset = NULL, family = Gaussian(), control = boost_control(),
       oobweights = NULL, baselearner = c("bbs", "bols", "btree", "bss", "bns"), 
       ...)

gamboost(formula, data = list(), na.action = na.omit, weights = NULL, 
         offset = NULL, family = Gaussian(), control = boost_control(),
         oobweights = NULL, baselearner = c("bbs", "bols", "btree", "bss", "bns"),
         dfbase = 4, ...)
}
\arguments{
  \item{formula}{ a symbolic description of the model to be fit. }
  
  \item{data}{ a data frame containing the variables in the model. }
  \item{na.action}{a function which indicates what should happen when
    the data contain \code{NA}s.}
    \item{weights}{ (optional) a numeric vector of weights to be used in 
    the fitting process.}
  \item{offset}{ a numeric vector to be used as offset (optional).}
  \item{family}{a \code{\link{Family}} object.}
  \item{control}{ a list of parameters controlling the algorithm. For
    more details see \code{\link{boost_control}}. }
  \item{oobweights}{ an additional vector of out-of-bag weights, which is
    used for the out-of-bag risk (i.e., if \code{boost_control(risk =
      "oobag")}). This argument is also used internally by
    \code{cvrisk}. }
  \item{baselearner}{ a character specifying the component-wise base
    learner to be used: \code{\link{bbs}} means P-splines with a
    B-spline basis (see Schmid and Hothorn 2008), \code{\link{bols}}
    linear models and \code{\link{btree}} boosts stumps. 
    \code{bss} and \code{bns} are deprecated.
    Component-wise smoothing splines have been considered in Buehlmann
    and Yu (2003) and Schmid and Hothorn (2008) investigate P-splines
    with a B-spline basis. Kneib, Hothorn and Tutz (2009) also utilize
    P-splines with a B-spline basis, supplement them with their
    bivariate tensor product version to estimate interaction surfaces
    and spatial effects and also consider random effects base
    learners.}
  \item{dfbase}{ a single integer giving the degrees of freedom for P-spline 
    base-learners (\code{\link{bbs}}) globally. }
  \item{\dots}{ additional arguments passed to \code{\link{mboost_fit}}; currently none.}
}
\details{

  A (generalized) additive model is fitted using a boosting algorithm based on
  component-wise base-learners. 
  
  The base-learners can either be specified via the \code{formula} object or via 
  the \code{baselearner} argument. The latter argument is the default base-learner 
  which is used for all variables in the formula, whithout explicit base-learner 
  specification (i.e., if the base-learners are explicitly specified in \code{formula}, 
  the \code{baselearner} argument will be ignored for this variable). 
  
  Of note, \code{"bss"} and \code{"bns"} are deprecated and only in the list for 
  backward compatibility.
  
  Note that more base-learners (i.e., in addition to the ones provided
  via \code{baselearner}) can be specified in \code{formula}. See 
  \code{\link{baselearners}} for details.
  
  The only difference when calling \code{mboost} and \code{gamboost} is that the
  latter function allows one to specify default degrees of freedom for smooth 
  effects specified via \code{baselearner = "bbs"}. In all other cases, 
  degrees of freedom need to be set manually via a specific definition of the 
  corresponding base-learner.
}
\value{
  An object of class \code{mboost} with \code{\link{print}},
  \code{\link{AIC}}, \code{\link{plot}} and \code{\link{predict}}
  methods being available.
}
\references{

  Peter Buehlmann and Bin Yu (2003),
  Boosting with the L2 loss: regression and classification.
  \emph{Journal of the American Statistical Association}, \bold{98},
  324--339.

  Peter Buehlmann and Torsten Hothorn (2007),
  Boosting algorithms: regularization, prediction and model fitting.
  \emph{Statistical Science}, \bold{22}(4), 477--505.

  Thomas Kneib, Torsten Hothorn and Gerhard Tutz (2009), Variable selection and
  model choice in geoadditive regression models, \emph{Biometrics}, \bold{65}(2),
  626--634.

  Matthias Schmid and Torsten Hothorn (2008),
  Boosting additive models using component-wise P-splines as
  base-learners. \emph{Computational Statistics & Data Analysis},
  \bold{53}(2), 298--311.

  Torsten Hothorn, Peter Buehlmann, Thomas Kneib, Mattthias Schmid
  and Benjamin Hofner (2010),
  Model-based Boosting 2.0.
  \emph{Journal of Machine Learning Research}, \bold{11}, 2109 -- 2113.

  Benjamin Hofner, Andreas Mayr, Nikolay Robinzonov and Matthias Schmid
  (2014). Model-based Boosting in R: A Hands-on Tutorial Using the R
  Package mboost. \emph{Computational Statistics}, \bold{29}, 3--35.\cr
  \doi{10.1007/s00180-012-0382-5}

  Available as vignette via: \code{vignette(package = "mboost", "mboost_tutorial")}
}
\seealso{
  See \code{\link{mboost_fit}} for the generic boosting function, 
  \code{\link{glmboost}} for boosted linear models, and
  \code{\link{blackboost}} for boosted trees. 
  
  See \code{\link{baselearners}} for possible base-learners. 
  
  See \code{\link{cvrisk}} for cross-validated stopping iteration. 
  
  Furthermore see \code{\link{boost_control}}, \code{\link{Family}} and
  \code{\link[mboost]{methods}}.
}
\examples{

    ### a simple two-dimensional example: cars data
    cars.gb <- gamboost(dist ~ speed, data = cars, dfbase = 4,
                        control = boost_control(mstop = 50))
    cars.gb
    AIC(cars.gb, method = "corrected")

    ### plot fit for mstop = 1, ..., 50
    plot(dist ~ speed, data = cars)
    tmp <- sapply(1:mstop(AIC(cars.gb)), function(i)
        lines(cars$speed, predict(cars.gb[i]), col = "red"))
    lines(cars$speed, predict(smooth.spline(cars$speed, cars$dist),
                              cars$speed)$y, col = "green")

    ### artificial example: sinus transformation
    x <- sort(runif(100)) * 10
    y <- sin(x) + rnorm(length(x), sd = 0.25)
    plot(x, y)
    ### linear model
    lines(x, fitted(lm(y ~ sin(x) - 1)), col = "red")
    ### GAM
    lines(x, fitted(gamboost(y ~ x,
                    control = boost_control(mstop = 500))),
          col = "green")

}
\keyword{models}
\keyword{nonlinear}
