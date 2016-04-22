
setClass("boost_family", representation = representation(
    ngradient  = "function",
    loss       = "function",
    risk       = "function",
    offset     = "function",
    fW         = "function",
    check_y    = "function",
    weights    = "logical",
    name       = "character",
    charloss   = "character",
    sigmaTF    = "logical"
))

setMethod("show", "boost_family", function(object) {
    cat("\n\t", object@name, "\n\n")
    cat("Loss function:", object@charloss, "\n")
})

Family <- function(ngradient, loss = NULL, risk = NULL, 
                   offset = function(y, w) 0, 
                   fW = function(f) rep(1, length(f)), check_y = function(y) TRUE,
                   weights = TRUE, name = "user-specified", sigmaTF=FALSE) {

    if (is.null(loss))
        loss <- function(y, f) NA
    if (is.null(risk))
        risk <- function(sigma=1, y, f, w = 1) sum(w * loss(y, f))

    RET <- new("boost_family", ngradient = ngradient, loss = loss, 
               risk = risk, offset = offset, fW = fW, check_y = check_y, 
               weights = weights, 
               name = name, charloss = paste(deparse(body(loss)), "\n"),
               sigmaTF = sigmaTF)
    RET
}

### Gaussian (Regression)
GaussReg <- function()
    Family(ngradient = function(sigma=1, y, f, w = 1) y - f,
           loss = function(y, f) (y - f)^2,
           offset = weighted.mean,
           check_y = function(y) {
               if (!is.numeric(y) || !is.null(dim(y)))
                   stop("response is not a numeric vector but ",
                        sQuote("family = GaussReg()"))
               TRUE
           },
           name = "Squared Error (Regression)",
           sigmaTF=FALSE)

### Gaussian (-1 / 1 Binary Classification)
GaussClass <- function()
    Family(ngradient = function(sigma=1, y, f, w = 1) - 2 * y + 2 * y * f,
           loss = function(y, f) 1 - 2 * y * f + (y * f)^2,
           check_y = function(y) {
               if (!is.factor(y))
                   stop("response is not a factor but ",
                           sQuote("family = GaussClass()"))  
               if (nlevels(y) != 2)
                   stop("response is not a factor at two levels but ",
                           sQuote("family = GaussClass()"))
               TRUE
           },

           name = "Squared Error (Classification)",
           sigmaTF=FALSE )

### Laplace
Laplace <- function()
    Family(ngradient = function(sigma=1, y, f, w = 1) sign(y - f),
           loss = function(y, f) abs(y - f),
           offset = function(y, w) median(y),
           check_y = function(y) {
               if (!is.numeric(y) || !is.null(dim(y)))
                   stop("response is not a numeric vector but ",
                        sQuote("family = Laplace()"))
               TRUE
           },
           name = "Absolute Error",
           sigmaTF=FALSE )

### Binomial
Binomial <- function()
    Family(ngradient = function(sigma=1, y, f, w = 1) {
               exp2yf <- exp(-2 * y * f)
               -(-2 * y * exp2yf) / (log(2) * (1 + exp2yf))
           },
           loss = function(y, f) {
               p <- exp(f) / (exp(f) + exp(-f))
               y <- (y + 1) / 2
               -y * log(p) - (1 - y) * log(1 - p)
           },
           offset = function(y, w) {
               p <- weighted.mean(y > 0, w)
               1/2 * log(p / (1 - p))
           },
           fW = function(f) {
               p <- exp(f) / (exp(f) + exp(-f))
               4 * p * (1 - p)
           },
           check_y = function(y) {
               if (!is.factor(y))
                   stop("response is not a factor but ",
                           sQuote("family = Binomial()"))  
               if (nlevels(y) != 2)
                   stop("response is not a factor at two levels but ",
                           sQuote("family = Binomial()"))
               TRUE
           },
           name = "Negative Binomial Likelihood",
           sigmaTF=FALSE )

### Poisson
Poisson <- function()
    Family(ngradient = function(sigma=1, y, f, w = 1) y - exp(f),
           loss = function(y, f) -dpois(y, exp(f), log = TRUE),
           offset = function(y, w) log(weighted.mean(y, w)),
           fW = function(f) exp(f),
           check_y = function(y) {
               if (any(y < 0) || any((y - round(y)) > 0))
                   stop("response is not an integer variable but ",
                        sQuote("family = Poisson()"))
               TRUE
           },
           name = "Poisson Likelihood",
           sigmaTF=FALSE )

### L1Huber
Huber <- function(d = NULL) {
    mc <- match.call()
    if (length(mc) == 2)
        dtxt <- deparse(mc[[2]])
    else
        dtxt <- NULL
    fit <- 0
    Family(ngradient = function(sigma=1, y, f, w = 1) {
               if (is.null(d)) d <- median(abs(y - fit))
               fit <<- f
               ifelse(abs(y - f) < d, y - f, d * sign(y - f))
           },
           loss = function(y, f) {
               if (is.null(d)) d <- median(abs(y - fit))
               ifelse((a <- abs(y - f)) < d, a^2/2, d*(a - d/2))
           },
           offset = function(y, w) median(y),
           check_y = function(y) {
               if (!is.numeric(y) || !is.null(dim(y)))
                   stop("response is not a numeric vector but ",
                        sQuote("family = Huber()"))
               TRUE
           },
           name = paste("Huber Error", 
               ifelse(is.null(d), "(with adaptive d)", 
                                  paste("(with d = ", dtxt, ")", sep = ""))),
           sigmaTF=FALSE )
}

### Adaboost
AdaExp <- function()
    Family(ngradient = function(sigma=1, y, f, w = 1) y * exp(-y * f),
           loss = function(y, f) exp(-y * f),
           offset = function(y, w) {
               p <- weighted.mean(y > 0, w)
               1/2 * log(p / (1 - p))
           },
           check_y = function(y) {
               if (!is.factor(y))
                   stop("response is not a factor but ",
                           sQuote("family = AdaExp()"))  
               if (nlevels(y) != 2)
                   stop("response is not a factor at two levels but ",
                           sQuote("family = AdaExp()"))
               TRUE
           },
           name = "Adaboost Exponential Error",
           sigmaTF=FALSE )

### Cox proportional hazards model (partial likelihood)
CoxPH <- function()
    Family(ngradient = function(sigma=1, y, f, w) {
               time <- y[,1]
               storage.mode(time) <- "double"
               event <- y[,2]
               storage.mode(event) <- "integer"
               if (length(w) == 1) w <- rep(w, length(time))
               storage.mode(w) <- "double"
               if (length(f) == 1)
                   f <- rep(f, length(time))
               storage.mode(f) <- "double"
               .Call("ngradientCoxPLik", time, event, f, w)
           },
           loss = plloss <- function(y, f, w) {
               time <- y[,1]
               event <- y[,2]
               n <- length(time)
               if (length(f) == 1) f <- rep(f, n)
               if (length(w) == 1) w <- rep(w, n)
               indx <- rep(1:n, w)
               time <- time[indx]
               event <- event[indx]
               ef <- exp(f)[indx]
               f <- f[indx]
               n <- length(time)
               risk <- rep(0, n)
               for (i in 1:n)
                   risk[i] <- sum((time >= time[i])*ef)
               event * (f - log(risk))
           },
           risk = function(sigma=1, y, f, w = 1) -sum(plloss(y, f, w)),
           check_y = function(y) {
               if (!inherits(y, "Surv"))
                   stop("response is not an object of class ", sQuote("Surv"),
                        " but ", sQuote("family = CoxPH()"))
               TRUE
           },
           weights = TRUE, 
           name = "Partial Likelihood",
           sigmaTF=FALSE )
           

## Weibull log likelihood for boosting aft models
Weib <- function()
    Family(ngradient = function(sigma=1, y, f, w = 1) {
               lnt <- log(y[,1])
               event <- y[,2]
               eta <- (lnt-f) / sigma
               (event*(exp(eta)-1) + (1-event)*exp(eta)) / sigma               
           },                      
           
           loss = plloss <- function(sigma=1, y, f){
               fw <- function(w){
               exp(w-exp(w))
               }
               Sw <- function(w){
               exp(-exp(w))
               }
               S.stime <- y[,1]
               S.event <- y[,2]
               lnt <- log(S.stime)
               eta <- (lnt-f)/sigma
               - S.event * log( fw( eta ) / sigma ) -
               (1-S.event) * log( Sw( eta ) )
               },
                      
           offset = function(y, w=NULL) {           
               return(1)
           },
           
           risk = function(sigma=1, y, f, w = 1) sum(w * plloss(sigma, y, f)),
           
           check_y = function(y) {
               if (!inherits(y, "Surv"))
                   stop("response is not an object of class ", sQuote("Surv"))
               TRUE
           },
               
           name = "Negative Weibull Likelihood",
           sigmaTF=TRUE )           


## log logistic log likelihood for boosting aft models
Loglog <- function()
    Family(ngradient = function(sigma=1, y, f, w = 1) {
               lnt <- log(y[,1])
               event <- y[,2]
               eta <- (lnt-f) / sigma
               nom <- (exp(-eta)+1)
               (event*(2/nom-1) + (1-event)/nom) / sigma               
           },                      
           
           loss = plloss <- function(sigma=1, y, f){
               fw <- function(w){
               exp(w) / (1+exp(w))^2
               }
               Sw <- function(w){
               1 / (1+exp(w))
               }
               S.stime <- y[,1]
               S.event <- y[,2]
               lnt <- log(S.stime)
               eta <- (lnt-f)/sigma
               - S.event * log( fw( eta ) / sigma ) -
               (1-S.event) * log( Sw( eta ) )
               },
                      
           offset = function(y, w=NULL) {           
               return(1)
           },
           
           risk = function(sigma=1, y, f, w = 1) sum(w * plloss(sigma, y, f)),
           
           check_y = function(y) {
               if (!inherits(y, "Surv"))
                   stop("response is not an object of class ", sQuote("Surv"))
               TRUE
           },
               
           name = "Negative Log Logistic Likelihood",
           sigmaTF=TRUE )
           
           
## log normal log likelihood for boosting aft models
LogNormal <- function()
    Family(ngradient = function(sigma=1, y, f, w = 1) {
               lnt <- log(y[,1])
               event <- y[,2]
               eta <- (lnt-f) / sigma
               (event*eta + (1-event) * dnorm(eta) / (1-pnorm(eta)) ) / sigma               
           },                      
           
           loss = plloss <- function(sigma=1, y, f){
               fw <- function(w){
               dnorm(w)
               }
               Sw <- function(w){
               1 - pnorm(w)
               }
               S.stime <- y[,1]
               S.event <- y[,2]
               lnt <- log(S.stime)
               eta <- (lnt-f)/sigma
               - S.event * log( fw( eta ) / sigma ) -
               (1-S.event) * log( Sw( eta ) )
               },
                      
           offset = function(y, w=NULL) {           
               return(3)
           },
           
           risk = function(sigma=1, y, f, w = 1) sum(w * plloss(sigma, y, f)),
           
           check_y = function(y) {
               if (!inherits(y, "Surv"))
                   stop("response is not an object of class ", sQuote("Surv"))
               TRUE
           },
               
           name = "Negative Log Normal Likelihood",
           sigmaTF=TRUE )
           

Hurdle <- function()
    Family(ngradient = function(sigma=1, y, f, w = 1)
               y/(1+exp(f)) - 1/(sigma*(1+exp(-f))) -
               exp(f)/(sigma*((1+exp(f))^{1/sigma}-1)*(1+exp(f))),
           loss = plloss <- function(sigma=1, y, f, w=1)
               - y*log(exp(f)/(1+exp(f))) + log(1+exp(f)) / sigma -
               log(gamma(y+1/sigma)) + log(gamma(1/sigma)) +
               log(1 - (1+exp(f))^{-1/sigma}),
           risk = function(sigma=1, y, f, w = 1) sum(w * plloss(y=y,
               f=f, sigma=sigma)),
           check_y = function(y) TRUE,
           offset = function(y, w=NULL) {
               return(3)
           },
           name = "Hurdle Likelihood",
           sigmaTF=T) 
           
NBinomial <- function()
    Family(ngradient = function(sigma=1, y, f, w = 1)
               y - (y+sigma)/(exp(f)+sigma)*exp(f),
           loss = plloss <- function(sigma=1, y, f, w=1)
             - (log(gamma(y+sigma)) - log(gamma(sigma)) - log(gamma(y+1)) +
                sigma*log(sigma) - sigma*log(exp(f)+sigma) + y*f -
                y*log(exp(f)+sigma)),
           risk = function(sigma=1, y, f, w = 1) sum(w * plloss(y=y,
               f=f, sigma=sigma)),
           check_y = function(y) TRUE,
           offset = function(y, w=NULL) {
               return(3)
           },
           name = "Negative Negative Binomial Likelihood",
           sigmaTF=T)
           
#ZeroInfNull <- function(off, cp, thetacp)
#    Family(ngradient = function(sigma=1, y, f, w = 1)
#               - exp(f) / (1+exp(f)) + (y==0) *
#               exp(f) / ( exp(f) + ( (exp(cp)+thetacp)/thetacp )^{-thetacp}),
#           loss = plloss <- function(sigma=1, y, f, w=1)
#               log(1+exp(f)) - (y==0) * log( exp(f) + ( (exp(cp)+thetacp)/
#               thetacp)^{-thetacp} ) +
#               (y>0)* ( thetacp* log( (exp(cp)+thetacp)/thetacp) + y*
#               log( 1+exp(-cp)*thetacp) ) + (y>0) * ( log(gamma(thetacp)) -
#               log(gamma(thetacp+y)) ),
#           risk = function(sigma=1, y, f, w = 1) sum(w * plloss(y=y,
#               f=f, sigma=sigma)),
#           check_y = function(y) TRUE,
#           offset = function(y, w=NULL) {
#               off
#           },
#           name = "Zero-inflated NegBin Likelihood, null part",
#           sigmaTF=F) 
           
#ZeroInfCount <- function(off, np, thetanp)
#    Family(ngradient = function(sigma=1, y, f, w = 1)
#               - (y==0) * exp(f) * ( (exp(f)+thetanp)/thetanp )^{-thetanp-1} /
#               ( exp(np) + ( (exp(f)+thetanp)/thetanp )^{-thetanp} ) -
#               (y>0) * thetanp * ( exp(f)/(thetanp+exp(f)) - y* exp(-f)/
#               (1+thetanp*exp(-f)) ),
#           loss = plloss <- function(sigma=1, y, f, w=1)
#               log(1+exp(np)) - (y==0) * log( exp(np) + ( (exp(f)+thetanp)/
#               thetanp)^{-thetanp} ) +
#               (y>0)* ( thetanp*log( (exp(f)+thetanp)/thetanp) + y*
#               log( 1+exp(-f)*thetanp) ) + (y>0) * ( log(gamma(thetanp)) -
#               log(gamma(thetanp+y)) ),
#           risk = function(sigma=1, y, f, w = 1) sum(w * plloss(y=y,
#               f=f, sigma=sigma)),
#           check_y = function(y) TRUE,
#           offset = function(y, w=NULL) {
#               off
#           },
#           name = "Zero-inflated NegBin Likelihood, count part",
#           sigmaTF=F) 
           
#ZeroInfTheta <- function(off, ff, gg)
#    Family(ngradient = function(sigma=1, y, f, w = 1)
#               (y==0) * 1 / ( ((exp(ff)+f)/f)^f * exp(gg) + 1 ) *
#               ( exp(ff) / (exp(ff)+f) - log( (exp(ff)+f)/f ) ) -
#               (y>0) * ( log((exp(ff)+f)/f) - exp(ff)/(exp(ff)+f) +
#               y * exp(-ff) / (1+exp(-ff)*f) ) -
#               (y>0) * (digamma(f)-digamma(y+f)),
#           loss = plloss <- function(sigma=1, y, f, w=1)
#               log(1+exp(gg)) - (y==0) * log( exp(gg) + ( (exp(ff)+f)/
#               f)^{-f} ) +
#               (y>0)* ( f*log( (exp(ff)+f)/f) + y*
#               log( 1+exp(-ff)*f) ) + (y>0) * ( log(gamma(f)) +
#               - log(gamma(f+y)) ),
#           risk = function(sigma=1, y, f, w = 1) sum(w * plloss(y=y,
#               f=f, sigma=sigma)),
#           check_y = function(y) TRUE,
#           offset = function(y, w=NULL) {
#               off
#           },
#           name = "Zero-inflated NegBin Likelihood, theta part",
#           sigmaTF=F) 
           

#ZeroInfNullN <- function(off, cp, thetacp)
#    Family(ngradient = function(sigma=1, y, f, w = 1){
#               zaehler <- - exp(f) / (1+exp(f)) + (y==0) *
#               exp(f) / ( exp(f) + ( (exp(cp)+thetacp)/thetacp )^{-thetacp})
#               TB <- (exp(cp)+thetacp) / thetacp
#               nenner <- exp(f)/(1+exp(f))^2 - (y==0) * exp(f) * TB^{-thetacp} /
#               ( exp(f) + TB^{-thetacp} )^2
#               return(zaehler / nenner)
#               },
#           loss = plloss <- function(sigma=1, y, f, w=1)
#               log(1+exp(f)) - (y==0) * log( exp(f) + ( (exp(cp)+thetacp)/
#               thetacp)^{-thetacp} ) +
#               (y>0)* ( thetacp* log( (exp(cp)+thetacp)/thetacp) + y*
#               log( 1+exp(-cp)*thetacp) ) + (y>0) * ( log(gamma(thetacp)) -
#               log(gamma(thetacp+y)) ),
#           risk = function(sigma=1, y, f, w = 1) sum(w * plloss(y=y,
#               f=f, sigma=sigma)),
#           check_y = function(y) TRUE,
#           offset = function(y, w=NULL) {
#               off
#           },
#           name = "Zero-inflated NegBin Likelihood, Newton null part",
#           sigmaTF=F) 
           
#ZeroInfCountN <- function(off, np, thetanp)
#    Family(ngradient = function(sigma=1, y, f, w = 1){
#               zaehler <- - (y==0) * exp(f) * ( (exp(f)+thetanp)/thetanp )^{-thetanp-1} /
#               ( exp(np) + ( (exp(f)+thetanp)/thetanp )^{-thetanp} ) -
#               (y>0) * thetanp * ( exp(f)/(thetanp+exp(f)) - y* exp(-f)/
#               (1+thetanp*exp(-f)) )
#               TB <- (exp(f)+thetanp)/thetanp
#               nenner1 <- (y==0) * exp(f) * TB^{-thetanp-1} /
#               (exp(np) + TB^{-thetanp})
#               nenner2 <- ( exp(np) + TB^{-thetanp} + exp(f) * TB^{-thetanp-1} ) /
#               ( exp(np) + TB^{-thetanp} )
#               nenner3 <- exp(f) * (-thetanp-1) / (exp(f) + thetanp)
#               nenner4 <- (y>0) * thetanp * ( thetanp*exp(f) / (exp(f)+thetanp)^2 +
#               y*exp(-f) / (1 + thetanp* exp(-f))^2 )
#               return(zaehler / (nenner1 * (nenner2+nenner3) + nenner4 ))
#               },
#           loss = plloss <- function(sigma=1, y, f, w=1)
#               log(1+exp(np)) - (y==0) * log( exp(np) + ( (exp(f)+thetanp)/
#               thetanp)^{-thetanp} ) +
#               (y>0)* ( thetanp*log( (exp(f)+thetanp)/thetanp) + y*
#               log( 1+exp(-f)*thetanp) ) + (y>0) * ( log(gamma(thetanp)) -
#               log(gamma(thetanp+y)) ),
#           risk = function(sigma=1, y, f, w = 1) sum(w * plloss(y=y,
#               f=f, sigma=sigma)),
#           check_y = function(y) TRUE,
#           offset = function(y, w=NULL) {
#               off
#           },
#           name = "Zero-inflated NegBin Likelihood, Newton count part",
#           sigmaTF=F)

#ZeroInfThetaN <- function(off, ff, gg)
#    Family(ngradient = function(sigma=1, y, f, w = 1){
#               zaehler <- (y==0) * 1 / ( ((exp(ff)+f)/f)^f * exp(gg) + 1 ) *
#               ( exp(ff) / (exp(ff)+f) - log( (exp(ff)+f)/f ) ) -
#               (y>0) * ( log((exp(ff)+f)/f) - exp(ff)/(exp(ff)+f) +
#               y * exp(-ff) / (1+exp(-ff)*f) ) -
#               (y>0) * (digamma(f)-digamma(y+f))
#               TB <- (exp(ff)+f)/f
#               nenner1 <- (y==0) / ( TB^f * exp(gg) + 1 )
#               nenner2 <- TB^f * exp(gg) / ( TB^f * exp(gg) + 1)
#               nenner3 <- exp(ff) / (exp(ff)+f) - log( (exp(ff)+f)/f )
#               nenner4 <- log( (exp(ff)+f)/f ) - exp(ff) / (f+exp(ff))
#               nenner5 <- - exp(ff) / ( f*(exp(ff)+f) ) + exp(ff) /
#               ( exp(ff)+f )^2
#               M1 <- nenner1 * (nenner2*nenner3*nenner4 + nenner5)
#               M2 <- (y>0) * ( -exp(ff) / ( f*(exp(ff)+f) ) + exp(ff) /
#               ( exp(ff)+f )^2 + y * f*exp(-2*ff) / ( 1 + exp(-ff)*f )^2 )
#               M3 <- (y>0) * ( trigamma(f) - trigamma(y+f) )
#               return(zaehler / (M1+M2+M3))
#               },
#           loss = plloss <- function(sigma=1, y, f, w=1)
#               log(1+exp(gg)) - (y==0) * log( exp(gg) + ( (exp(ff)+f)/
#               f)^{-f} ) +
#               (y>0)* ( f*log( (exp(ff)+f)/f) + y*
#               log( 1+exp(-ff)*f) ) + (y>0) * ( log(gamma(f)) +
#               - log(gamma(f+y)) ),
#           risk = function(sigma=1, y, f, w = 1) sum(w * plloss(y=y,
#               f=f, sigma=sigma)),
#           check_y = function(y) TRUE,
#           offset = function(y, w=NULL) {
#               off
#           },
#           name = "Zero-inflated NegBin Likelihood, Newton theta part",
#           sigmaTF=F) 
           

#ZeroInfNullN <- function(off, cp, thetacp)
#    Family(ngradient = function(sigma=1, y, f, w = 1)
#               ( - exp(f) / (1+exp(f)) + (y==0) *
#               exp(f) / ( exp(f) + ( (exp(cp)+thetacp)/thetacp )^{-thetacp}) ) /
#               ( exp(f)/(1+exp(f))^2 - (y==0) * exp(f) * ((exp(cp)+
#               thetacp) / thetacp)^{-thetacp} / ( exp(f) + ((exp(cp)+thetacp)
#               / thetacp)^{-thetacp} )^2 ),
#           loss = plloss <- function(sigma=1, y, f, w=1)
#               log(1+exp(f)) - (y==0) * log( pmax( exp(f) + ( (exp(cp)+thetacp)/
#               thetacp)^{-thetacp} , 0.01 ) ) +
#               (y>0)* ( thetacp* log( pmax((exp(cp)+thetacp)/thetacp , 0.01) ) + y*
#               log( pmax( 1+exp(-cp)*thetacp), 0.01) ) + (y>0) * ( log(gamma( pmin(thetacp, 20) )) -
#               log(gamma( pmin( thetacp+y, 20) )) ),
#           risk = function(sigma=1, y, f, w = 1) sum(w * plloss(y=y,
#               f=f, sigma=sigma)),
#           check_y = function(y) TRUE,
#           offset = function(y, w=NULL) {
#               off
#           },
#           name = "Zero-inflated NegBin Likelihood, Newton null part",
#           sigmaTF=F) 
           
#ZeroInfCountN <- function(off, np, thetanp)
#    Family(ngradient = function(sigma=1, y, f, w = 1)
#               ( - (y==0) * exp(f) * ( (exp(f)+thetanp)/thetanp )^{-thetanp-1} /
#               ( exp(np) + ( (exp(f)+thetanp)/thetanp )^{-thetanp} ) -
#               (y>0) * thetanp * ( exp(f)/(thetanp+exp(f)) - y* exp(-f)/
#               (1+thetanp*exp(-f)) ) ) /
#               ( (y==0) * exp(f) * ((exp(f)+thetanp)/
#               thetanp)^{-thetanp-1} / (exp(np) + ((exp(f)+thetanp)/
#               thetanp)^{-thetanp}) *
#               ( ( exp(np) + ((exp(f)+thetanp)/thetanp)^{-thetanp} +
#               exp(f) * ((exp(f)+thetanp)/thetanp)^{-thetanp-1} ) /
#               ( exp(np) + ((exp(f)+thetanp)/thetanp)^{-thetanp} ) +
#               exp(f) * (-thetanp-1) / (exp(f) + thetanp) ) +
#               (y>0) * thetanp * ( thetanp*exp(f) / (exp(f)+thetanp)^2 +
#               y*exp(-f) / (1 + thetanp* exp(-f))^2 ) ),
#           loss = plloss <- function(sigma=1, y, f, w=1)
#               log(1+exp(np)) - (y==0) * log( pmax( exp(np) + ( (exp(f)+thetanp)/
#               thetanp)^{-thetanp} , 0.01) ) +
#               (y>0)* ( thetanp*log( pmax( (exp(f)+thetanp)/thetanp, 0.01) ) + y*
#               log( pmax( 1+exp(-f)*thetanp, 0.01) ) ) + (y>0) * ( log(gamma( pmin(thetanp,20) )) -
#               log(gamma( pmin( thetanp+y, 20) )) ),
#           risk = function(sigma=1, y, f, w = 1) sum(w * plloss(y=y,
#               f=f, sigma=sigma)),
#           check_y = function(y) TRUE,
#           offset = function(y, w=NULL) {
#               off
#           },
#           name = "Zero-inflated NegBin Likelihood, Newton count part",
#           sigmaTF=F)

#ZeroInfThetaN <- function(off, ff, gg)
#    Family(ngradient = function(sigma=1, y, f, w = 1)
#               ( (y==0) * 1 / ( ((exp(ff)+f)/f)^f * exp(gg) + 1 ) *
#               ( exp(ff) / (exp(ff)+f) - log( pmax((exp(ff)+f)/f,0.01) ) ) -
#               (y>0) * ( log( pmax((exp(ff)+f)/f,0.01) ) - exp(ff)/(exp(ff)+f) +
#               y * exp(-ff) / (1+exp(-ff)*f) ) -
#               (y>0) * (digamma(pmax(f,0.1))-digamma(pmax(y+f,0.1))) ) /
#               (
#               (y==0) / ( ((exp(ff)+f)/f)^f * exp(gg) + 1 ) * (
#               (((exp(ff)+f)/f)^f * exp(gg) / ( ((exp(ff)+f)/f)^f * exp(gg) + 1)) *
#               (exp(ff) / (exp(ff)+f) - log( pmax((exp(ff)+f)/f,0.01) )) *
#               (log( pmax((exp(ff)+f)/f,0.01) ) - exp(ff) / (f+exp(ff))) +
#               - exp(ff) / ( f*(exp(ff)+f) ) + exp(ff) /
#               ( exp(ff)+f )^2
#               ) +
#               ((y>0) * ( -exp(ff) / ( f*(exp(ff)+f) ) + exp(ff) /
#               ( exp(ff)+f )^2 + y * f*exp(-2*ff) / ( 1 + exp(-ff)*f )^2 )) +
#               ((y>0) * ( trigamma(pmax(f,0.01)) - trigamma( pmax(y+f,0.01)) ))),
#           loss = plloss <- function(sigma=1, y, f, w=1)
#               log(1+exp(gg)) - (y==0) * log( pmax( exp(gg) + ( (exp(ff)+f)/
#               f)^{-f} , 0.01) ) +
#               (y>0)* ( f*log( pmax( (exp(ff)+f)/f, 0.01) ) + y*
#               log( pmax( 1+exp(-ff)*f, 0.01) ) ) + (y>0) * ( log( gamma( pmin(f, 20) )) +
#               - log(gamma( pmin(f+y, 20) )) ),
#           risk = function(sigma=1, y, f, w = 1) sum(w * plloss(y=y,
#               f=f, sigma=sigma)),
#           check_y = function(y) TRUE,
#           offset = function(y, w=NULL) {
#               off
#           },
#           name = "Zero-inflated NegBin Likelihood, Newton theta part",
#           sigmaTF=F)
           
           
ZeroInfNull <- function(off, cp, thetacp)
    Family(ngradient = function(sigma=1, y, f, w = 1)
               - exp(f) / (1+exp(f)) + (y<0.5) *
               exp(f) / ( exp(f) + ( thetacp/(exp(cp)+thetacp) )^{thetacp}),
           loss = plloss <- function(sigma=1, y, f, w=1)
               log(1+exp(f)) - (y<0.5) * log( pmax( exp(f) + ( thetacp/(exp(cp)+thetacp) )^{thetacp}, 0.000001) ) +
               (y>0.5)* ( thetacp* log( pmax( (exp(cp)+thetacp)/thetacp, 0.000001) ) + y*
               log( pmax( 1+exp(-cp)*thetacp), 0.000001) ) + (y>0.5) * ( log( gamma( pmax(thetacp, 0.000001) )) -
               log( gamma( pmax( thetacp+y, 0.000001) )) ),
           risk = function(sigma=1, y, f, w = 1) sum(w * plloss(y=y,
               f=f, sigma=sigma)),
           check_y = function(y) TRUE,
           offset = function(y, w=NULL) {
               off
           },
           name = "Zero-inflated NegBin Likelihood, null part",
           sigmaTF=F)
           
ZeroInfCount <- function(off, np, thetanp)
    Family(ngradient = function(sigma=1, y, f, w = 1)
               - (y<0.5) * exp(f) * ( thetanp/(exp(f)+thetanp) )^{thetanp+1} /
               ( exp(np) + ( thetanp/(exp(f)+thetanp) )^{thetanp} ) -
               (y>0.5) * thetanp * ( exp(f)/(thetanp+exp(f)) - y* exp(-f)/
               (1+thetanp*exp(-f)) ),
           loss = plloss <- function(sigma=1, y, f, w=1)
               log(1+exp(np)) - (y<0.5) * log( pmax( exp(np) + ( thetanp/(exp(f)+thetanp))^{thetanp}, 0.000001) ) +
               (y>0.5)* ( thetanp*log( pmax((exp(f)+thetanp)/thetanp, 0.000001) ) + y*
               log( pmax( 1+exp(-f)*thetanp, 0.000001)) ) + (y>0.5) * ( log(gamma(pmax(thetanp, 0.000001))) -
               log(gamma(pmax(thetanp+y, 0.000001))) ),
           risk = function(sigma=1, y, f, w = 1) sum(w * plloss(y=y,
               f=f, sigma=sigma)),
           check_y = function(y) TRUE,
           offset = function(y, w=NULL) {
               off
           },
           name = "Zero-inflated NegBin Likelihood, count part",
           sigmaTF=F) 
           
ZeroInfTheta <- function(off, ff, gg)
    Family(loss = plloss <- function(sigma=1, y, f, w=1)
               log(1+exp(gg)) - (y<0.5) * log( pmax( exp(gg) + ( f/(exp(ff)+f))^{f} , 0.000001) ) +
               (y>0.5)* ( f*log( pmax( (exp(ff)+f)/f, 0.000001) ) + y*
               log( pmax( 1+exp(-ff)*f, 0.000001) ) ) + (y>0.5) * ( log( gamma( pmax(f, 0.000001))) +
               - log( pmax( gamma(f+y), 0.000001) ) ),
           ngradient = function(sigma=1, y, f, w = 1)
               (y<0.5) * 1 / ( ((exp(ff)+f)/f)^f * exp(gg) + 1 ) *
               ( exp(ff) / (exp(ff)+f) - log( pmax( (exp(ff)+f)/f, 0.000001) ) ) -
               (y>0.5) * ( log( pmax((exp(ff)+f)/f, 0.000001) ) - exp(ff)/(exp(ff)+f) +
               y * exp(-ff) / (1+exp(-ff)*f) ) -
               (y>0.5) * (digamma( pmax(f, 0.000001) )-digamma( pmax(y+f, 0.000001))),
           risk = function(sigma=1, y, f, w = 1) sum(w * plloss(y=y,
               f=f, sigma=sigma)),
           check_y = function(y) TRUE,
           offset = function(y, w=NULL) {
               off
           },
           name = "Zero-inflated NegBin Likelihood, theta part",
           sigmaTF=F)
           
           
ZeroInfNullP <- function(off, cp)
    Family(ngradient = function(sigma=1, y, f, w = 1)
              (y==0) * exp(f) / (exp(f)+exp(-exp(cp))) - exp(f)/(1+exp(f)),
           loss = plloss <- function(sigma=1, y, f, w=1)
               - (y==0) * ( log(exp(f)+exp(-exp(cp))) - log(1+exp(f)) ) -
               (y>0) * ( -log(1+exp(f)) + y*cp - log(factorial(y)) - exp(cp) ),
           risk = function(sigma=1, y, f, w = 1) sum(w * plloss(y=y,
               f=f, sigma=sigma)),
           check_y = function(y) TRUE,
           offset = function(y, w=NULL) {
               off
           },
           name = "Zero-inflated Poisson Likelihood, null part",
           sigmaTF=F)
           
ZeroInfCountP <- function(off, np)
    Family(ngradient = function(sigma=1, y, f, w = 1)
               - (y==0) * exp(f-exp(f)) / (exp(np)+exp(-exp(f))) +
               (y>0) * (y-exp(f)),
           loss = plloss <- function(sigma=1, y, f, w=1)
               - (y==0) * ( log(exp(np)+exp(-exp(f))) - log(1+exp(np)) ) -
               (y>0) * ( -log(1+exp(np)) + y*f - log(factorial(y)) - exp(f) ),
           risk = function(sigma=1, y, f, w = 1) sum(w * plloss(y=y,
               f=f, sigma=sigma)),
           check_y = function(y) TRUE,
           offset = function(y, w=NULL) {
               off
           },
           name = "Zero-inflated Poisson Likelihood, count part",
           sigmaTF=F) 
           

