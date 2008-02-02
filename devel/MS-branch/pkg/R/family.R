
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
