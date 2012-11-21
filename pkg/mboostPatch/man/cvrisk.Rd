\name{cvrisk}
\alias{cvrisk}
\alias{cv}
\title{ Cross-Validation }
\description{
  Cross-validated estimation of the empirical risk for hyper-parameter selection.
}
\usage{
cvrisk(object, folds = cv(model.weights(object)),
       grid = 1:mstop(object),
       papply = mclapply,
       fun = NULL, ...)
cv(weights, type = c("bootstrap", "kfold", "subsampling"),
   B = ifelse(type == "kfold", 10, 25), prob = 0.5, strata = NULL)
}
\arguments{
  \item{object}{ an object of class \code{mboost}.}
  \item{folds}{ a weight matrix with number of rows equal to the number
                of observations. The number of columns corresponds to
                the number of cross-validation runs. Can be computed
                using function \code{cv} and defaults to 25 bootstrap samples.}
  \item{grid}{ a vector of stopping parameters the empirical risk
                is to be evaluated for. }
  \item{papply}{ (parallel) apply function, defaults to  \code{\link[parallel]{mclapply}}.
    Alternatively, \code{parLapply} can be used. In the
    latter case, usually more setup is needed (see example for some details).}
  \item{fun}{ if \code{fun} is NULL, the out-of-sample risk is returned. \code{fun},
              as a function of \code{object}, may extract any other characteristic
              of the cross-validated models. These are returned as is.}
  \item{weights}{ a numeric vector of weights for the model to be cross-validated.}
  \item{type}{ character argument for specifying the cross-validation
               method. Currently (stratified) bootstrap, k-fold cross-validation
               and subsampling are implemented.}
  \item{B}{ number of folds, per default 25 for \code{bootstrap} and
            \code{subsampling} and 10 for \code{kfold}.}
  \item{prob}{ percentage of observations to be included in the learning samples
               for subsampling.}
  \item{strata}{ a factor of the same length as \code{weights} for stratification.}
  \item{...}{additional arguments passed to \code{\link[parallel]{mclapply}}
             eventually.}
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

  \code{cvrisk} runs in parallel on OSes where forking is possible
  (i.e., not on Windows) and multiple cores/processors are available.
  The scheduling
  can be changed by the corresponding arguments of
  \code{\link[parallel]{mclapply}} (via the dot arguments).

  The function \code{cv} can be used to build an appropriate
  weight matrix to be used with \code{cvrisk}. If \code{strata} is defined
  sampling is performed in each stratum separately thus preserving
  the distribution of the \code{strata} variable in each fold.
}
\value{
  An object of class \code{cvrisk} (when \code{fun} wasn't specified), basically a matrix
  containing estimates of the empirical risk for a varying number
  of bootstrap iterations. \code{plot} and \code{print} methods
  are available as well as a \code{mstop} method.
}
\references{

  Torsten Hothorn, Friedrich Leisch, Achim Zeileis and Kurt Hornik (2006),
  The design and analysis of benchmark experiments.
  \emph{Journal of Computational and Graphical Statistics}, \bold{14}(3),
  675--699.

  Andreas Mayr, Benjamin Hofner, and Matthias Schmid (2012). The
  importance of knowing when to stop - a sequential stopping rule for
  component-wise gradient boosting. \emph{Methods of Information in
  Medicine}, DOI: \url{http://dx.doi.org/10.3414/ME11-02-0030}
}
\seealso{\code{\link{AIC.mboost}} for
  \code{AIC} based selection of the stopping iteration. Use \code{mstop}
  to extract the optimal stopping iteration from \code{cvrisk}
  object.}
\examples{

  data("bodyfat", package = "mboost")

  ### fit linear model to data
  model <- glmboost(DEXfat ~ ., data = bodyfat, center = TRUE)

  ### AIC-based selection of number of boosting iterations
  maic <- AIC(model)
  maic

  ### inspect coefficient path and AIC-based stopping criterion
  par(mai = par("mai") * c(1, 1, 1, 1.8))
  plot(model)
  abline(v = mstop(maic), col = "lightgray")

  ### 10-fold cross-validation
  cv10f <- cv(model.weights(model), type = "kfold")
  cvm <- cvrisk(model, folds = cv10f, papply = lapply)
  print(cvm)
  mstop(cvm)
  plot(cvm)

  ### 25 bootstrap iterations (manually)
  set.seed(290875)
  n <- nrow(bodyfat)
  bs25 <- rmultinom(25, n, rep(1, n)/n)
  cvm <- cvrisk(model, folds = bs25, papply = lapply)
  print(cvm)
  mstop(cvm)
  plot(cvm)

  ### same by default
  set.seed(290875)
  cvrisk(model, papply = lapply)

  ### 25 bootstrap iterations (using cv)
  set.seed(290875)
  bs25_2 <- cv(model.weights(model), type="bootstrap")
  all(bs25 == bs25_2)

  ### trees
  blackbox <- blackboost(DEXfat ~ ., data = bodyfat)
  cvtree <- cvrisk(blackbox, papply = lapply)
  plot(cvtree)


  ### cvrisk in parallel modes:

  \dontrun{## parallel:::mclapply only runs properly on unix systems
    cvrisk(model)
  }

  \dontrun{## infrastructure needs to be set up in advance
    cl <- makeCluster(25) # e.g. to run cvrisk on 25 nodes via PVM
    myApply <- function(X, FUN, ...) {
      myFun <- function(...) {
          library("mboost") # load mboost on nodes
          FUN(...)
      }
      ## further set up steps as required
      parLapply(cl = cl, X, myFun, ...)
    }
    cvrisk(model, papply = myApply)
    stopCluster(cl)
  }

}
\keyword{models}
\keyword{regression}
