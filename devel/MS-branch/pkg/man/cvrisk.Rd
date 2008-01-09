\name{cvrisk}
\alias{cvrisk}
\title{ Cross-Validation }
\description{
  Cross-validated estimation of the empirical risk 
  for hyper-parameter selection.
}
\usage{
cvrisk(object, folds, 
       grid = floor(seq(from = floor(mstop(object)/10), 
                       to = mstop(object), length = 10)))
}
\arguments{
  \item{object}{ an object of class \code{gb}. }
  \item{folds}{ a weight matrix with number of rows equal to the number
                of observations. The number of columns corresponds to
                the number of cross-validation runs.}
  \item{grid}{ a vector of stopping parameters the empirical risk
                is to be evaluated for. }
}
\details{

  The number of boosting iterations is a hyper-parameter of the
  boosting algorithms implemented in this package. Honest, 
  i.e., cross-validated, estimates of the empirical risk
  for different stopping parameters \code{mstop} are computed by
  this function which can be utilized to choose an appropriate
  number of boosting iterations to be applied.

  Different forms of cross-validation can be applied, for example
  10-fold cross-validation or bootstrapping. The weights (zero weights
  correspond to test cases) are defined via the \code{folds} matrix.

}
\value{
  An object of class \code{cvrisk}, basically a matrix 
  containing estimates of the empirical risk for a varying number
  of bootstrap iterations. \code{plot} and \code{print} methods
  are available as well as a \code{mstop} method.
}
\references{

  Torsten Hothorn, Friedrich Leisch, Achim Zeileis and Kurt Hornik (2006),
  The design and analysis of benchmark experiments.
  \emph{Journal of Computational and Graphical Statistics}, \bold{14}(3), 
  675--699.
}
\examples{

  data("bodyfat", package = "mboost")
  
  ### fit linear model to data
  model <- glmboost(DEXfat ~ ., data = bodyfat, 
                    control = boost_control(center = TRUE))

  ### AIC-based selection of number of boosting iterations
  maic <- AIC(model)
  maic

  ### inspect coefficient path and AIC-based stopping criterion
  par(mai = par("mai") * c(1, 1, 1, 1.8))
  plot(model)
  abline(v = mstop(maic), col = "lightgray")

  ### 10-fold cross-validation
  n <- nrow(bodyfat)
  k <- 10
  ntest <- floor(n / k)
  cv10f <- matrix(c(rep(c(rep(0, ntest), rep(1, n)), k - 1), 
                    rep(0, n * k - (k - 1) * (n + ntest))), nrow = n)
  cvm <- cvrisk(model, folds = cv10f)
  print(cvm)
  mstop(cvm)
  plot(cvm)

  ### 25 bootstrap iterations
  set.seed(290875)
  bs25 <- rmultinom(25, n, rep(1, n)/n)
  cvm <- cvrisk(model, folds = bs25)
  print(cvm)
  mstop(cvm)

  layout(matrix(1:2, ncol = 2))
  plot(cvm)

  ### trees
  blackbox <- blackboost(DEXfat ~ ., data = bodyfat)
  cvtree <- cvrisk(blackbox, folds = bs25)
  plot(cvtree)

}
\keyword{models}
\keyword{regression}
