
setClass("boost_family", representation = representation(
    ngradient  = "function",
    loss       = "function",
    risk       = "function",
    offset     = "function",
    fW         = "function",
    check_y    = "function",
    weights    = "logical",
    name       = "character",
    charloss   = "character"
))

setMethod("show", "boost_family", function(object) {
    cat("\n\t", object@name, "\n\n")
    cat("Loss function:", object@charloss, "\n")
})

Family <- function(ngradient, loss = NULL, risk = NULL, 
                   offset = function(y, w) 0, 
                   fW = function(f) rep(1, length(f)), check_y = function(y) TRUE,
                   weights = TRUE, name = "user-specified") {

    if (is.null(loss))
        loss <- function(y, f) NA
    if (is.null(risk))
        risk <- function(y, f, w = 1) sum(w * loss(y, f))
    RET <- new("boost_family", ngradient = ngradient, loss = loss, 
               risk = risk, offset = offset, fW = fW, check_y = check_y, 
               weights = weights, 
               name = name, charloss = paste(deparse(body(loss)), "\n"))
    RET
}

### Gaussian (Regression)
GaussReg <- function()
    Family(ngradient = function(y, f, w = 1) y - f,
           loss = function(y, f) (y - f)^2,
           offset = weighted.mean,
           check_y = function(y) {
               if (!is.numeric(y) || !is.null(dim(y)))
                   stop("response is not a numeric vector but ",
                        sQuote("family = GaussReg()"))
               TRUE
           },
           name = "Squared Error (Regression)")

### Gaussian (-1 / 1 Binary Classification)
GaussClass <- function()
    Family(ngradient = function(y, f, w = 1) 2 * y - 2 * y * f,
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

           name = "Squared Error (Classification)")

### Laplace
Laplace <- function()
    Family(ngradient = function(y, f, w = 1) sign(y - f),
           loss = function(y, f) abs(y - f),
           offset = function(y, w) median(y),
           check_y = function(y) {
               if (!is.numeric(y) || !is.null(dim(y)))
                   stop("response is not a numeric vector but ",
                        sQuote("family = Laplace()"))
               TRUE
           },
           name = "Absolute Error")

### Binomial
Binomial <- function()
    Family(ngradient = function(y, f, w = 1) {
               exp2yf <- exp(-2 * y * f)
               -(-2 * y * exp2yf) / (log(2) * (1 + exp2yf))
           },
           loss = function(y, f) {
               ### FIXME: this is unstable
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
           name = "Negative Binomial Likelihood")

### Poisson
Poisson <- function()
    Family(ngradient = function(y, f, w = 1) y - exp(f),
           loss = function(y, f) -dpois(y, exp(f), log = TRUE),
           offset = function(y, w) log(weighted.mean(y, w)),
           fW = function(f) exp(f),
           check_y = function(y) {
               if (any(y < 0) || any((y - round(y)) > 0))
                   stop("response is not an integer variable but ",
                        sQuote("family = Poisson()"))
               TRUE
           },
           name = "Poisson Likelihood")

### L1Huber
Huber <- function(d = NULL) {
    mc <- match.call()
    if (length(mc) == 2)
        dtxt <- deparse(mc[[2]])
    else
        dtxt <- NULL
    fit <- 0
    Family(ngradient = function(y, f, w = 1) {
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
                                  paste("(with d = ", dtxt, ")", sep = ""))))
}

### Adaboost
AdaExp <- function()
    Family(ngradient = function(y, f, w = 1) y * exp(-y * f),
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
           name = "Adaboost Exponential Error")

### Cox proportional hazards model (partial likelihood)
CoxPH <- function()
    Family(ngradient = function(y, f, w) {
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
           risk = function(y, f, w = 1) -sum(plloss(y, f, w)),
           check_y = function(y) {
               if (!inherits(y, "Surv"))
                   stop("response is not an object of class ", sQuote("Surv"),
                        " but ", sQuote("family = CoxPH()"))
               TRUE
           },
           weights = TRUE, 
           name = "Partial Likelihood")

QuantReg <- function(tau = 0.5, qoffset = 0.5) {
    stopifnot(tau > 0 && tau < 1)
    stopifnot(qoffset > 0 && qoffset < 1)
    Family(
        ngradient = function(y, f, w = 1) 
            tau*((y - f) >= 0) - (1 - tau)*((y - f)<0) ,
        loss = function(y, f) tau*(y-f)*((y-f)>=0) - (1-tau)*(y-f)*((y-f)<0) ,
        offset = function(y, w = rep(1, length(y))) 
            quantile(y[rep(1:length(y), w)], qoffset),
        check_y = function(y) {
            if (!is.numeric(y) || !is.null(dim(y)))
                stop("response is not a numeric vector but ", 
                     sQuote("family = QuantReg()"))
            return(TRUE)
        },
        name = "Quantile Regression")
}

ExpectileReg <- function (tau = 0.5) {
    stopifnot(tau > 0 && tau < 1)
    Family(
        ngradient = function(y, f, w = 1) 
            2 * tau * (y - f) * ((y - f) > 0) - 2 * (1 - tau) * 
            (f - y) * ((y - f) < 0) + 0 * ((y - f) == 0), 
        loss = function(y, f) tau * (y - f)^2 * 
            ((y - f) >= 0) + (1 - tau) * (y - f)^2 * ((y - f) < 0), 
        offset = function(y, w = rep(1, length(y))) 
            mean(y[w == 1]), 
        check_y = function(y) {
            if (!is.numeric(y) || !is.null(dim(y))) 
                stop("response is not a numeric vector but ", 
                  sQuote("family = ExpectileReg()"))
            TRUE
        }, 
        name = "Expectile Regression")
}
