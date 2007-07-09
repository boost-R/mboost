\name{FP}
\alias{FP}
\title{ Fractional Polynomials }
\description{
  Fractional polynomials transformation for continuous covariates.
}
\usage{
FP(x, p = c(-2, -1, -0.5, 0.5, 1, 2, 3))
}
\arguments{
  \item{x}{ a numeric vector. }
  \item{p}{ all powers of \code{x} to be included. }
}
\details{

    A fractional polynomial refers to a model 
    \eqn{\sum_{j = 1}^k (\beta_j x^{p_j} + \gamma_j x^{p_j} \log(x)) + 
    \beta_{k + 1} \log(x)  + \gamma_{k + 1} \log(x)^2}
    where the degree of the fractional polynomial is the number of non-zero regression coefficients
    \eqn{\beta} and \eqn{\gamma}. See \code{\link[mfp]{mfp}} for the
    reference implementation.
  
    Currently, no scaling of \code{x} is implemented. However, one may wish to standardize
    the inputs prior to fitting the model.

}
\value{
  A matrix including all powers \code{p} of \code{x}, all powers \code{p} of \code{log(x)} and
  \code{log(x)}.
}
\references{ 

    Willi Sauerbrei and Patrick Royston (1999), Building multivariable prognostic and 
    diagnostic models: transformation of the predictors by using fractional polynomials. 
    \emph{Journal of the Royal Statistical Society A}, \bold{162}, 71--94.

}
\examples{

    data("bodyfat", package = "mboost")
    tbodyfat <- bodyfat
 
    ### map covariates into [1, 2]
    indep <- names(tbodyfat)[-2]
    tbodyfat[indep] <- lapply(bodyfat[indep], function(x) {
        x <- x - min(x)
        x / max(x) + 1
    })
 
    ### generate formula
    fpfm <- as.formula(paste("DEXfat ~ ", paste("FP(", indep, ")", 
                             collapse = "+")))
    fpfm

    ### fit linear model
    bf_fp <- glmboost(fpfm, data = tbodyfat, 
                      control = boost_control(mstop = 3000))

    ### when to stop
    mstop(aic <- AIC(bf_fp))
    plot(aic)

    ### coefficients
    cf <- coef(bf_fp[mstop(aic)])
    length(cf)
    cf[abs(cf) > 0]

}
\keyword{datagen}
