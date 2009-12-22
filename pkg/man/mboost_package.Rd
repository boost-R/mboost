\name{mboost-package}
\alias{mboost-package}
\docType{package}
\title{
  mboost: Model-Based Boosting
}
\description{
  Functional gradient descent algorithms (boosting) for optimizing
  general loss functions utilizing componentwise least squares, either
  of parametric linear form or smoothing splines, or regression trees as
  base learners for fitting generalized linear, additive and interaction
  models to potentially high-dimensional data.
}
\details{
\tabular{ll}{
Package: \tab mboost\cr
Type: \tab Package\cr
Version: \tab 2.0-0\cr
Date: \tab 2010-01-XX\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
LazyData: \tab yes\cr
}
}
\author{
  Torsten Hothorn \email{Torsten.Hothorn@R-project.org},\cr
  Peter Buhlmann, Thomas Kneib, Matthias Schmid and
  Benjamin Hofner
}

\references{
  Peter Buhlmann and Torsten Hothorn (2007),
  Boosting algorithms: regularization, prediction and model fitting.
  \emph{Statistical Science}, \bold{22}(4), 477--505.
}

\seealso{
  The main fitting functions include:\cr
  \code{\link{gamboost}} for boosted (generalized) additive models,
  \code{\link{glmboost}} for boosted linear models and
  \code{\link{blackboost}} for boosted trees.\cr
  See there for more details and further links.
}

\examples{
data(bodyfat)
model <- gamboost(DEXfat ~ .,
         data = bodyfat, baselearner = "bbs",
         control=boost_control(trace=TRUE))

## extract coefficients of selected base-learners
coef(model)

## extract *all* coefficients
coef(model, which=1:9)

## plot
par(mfrow=c(3,3))
plot(model)

par(mfrow=c(3,3))
plot(model, which=1:9) ##???

## continue until mstop = 200
model[200]
## Note that the original object "model" is altered

## show numbers of selected base-learners
selected(model)

## use cvrisk for 5-fold crossvalidation
weightMat <- cv(model.weights(model), type = "kfold", B = 5)
cv5fold <- cvrisk(model, weightMat)
model[mstop(cv5fold)]
## Note that the original object "model" is altered again
par(mfrow=c(2,3))
plot(model)



## now fit model with intrinsic model selection
## hipcirc + anthro3a + kneebreadth are fitted as
## bols base-learner with intercept = FALSE

model <- gamboost(DEXfat ~ rep(1,nrow(bodyfat)) +
                           hipcirc + anthro3a + kneebreadth +
                           bbs(hipcirc, center=TRUE, df=1) +
                           bbs(anthro3a, center=TRUE, df=1) +
                           bbs(kneebreadth, center=TRUE, df=1),
         data = bodyfat,
         baselearner = function(x) bols(x, intercept=FALSE))

## plot
coef(model)

## anthro3a is selected as linear effect and smooth effect
## to combine this in one plot we can use fitted or predict
par(mfrow=c(2,2))
plot(model, which="anthro3a")
fitted <- rowSums(predict(model, which="anthro3a"))
plot(bodyfat$anthro3a, fitted, main="fitted(anthro3a)")
## or you can use predict to plot on a grid
grid <- seq(min(bodyfat$anthro3a), max(bodyfat$anthro3a), length= 100)
DF <- data.frame(anthro3a=grid)
predict <- rowSums(predict(model, newdata = DF, which="anthro3a"))
plot(grid, predict, main="predict(anthro3a)")
}

\keyword{ package }
\keyword{ smooth }
\keyword{ models }