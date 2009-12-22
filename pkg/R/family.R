
setClass("boost_family", representation = representation(
    ngradient  = "function",
    risk       = "function",
    offset     = "function",
    check_y    = "function",
    weights    = "function",
    nuisance   = "function",
    response   = "function",
    name       = "character",
    charloss   = "character"
))

setClass("boost_family_glm", contains = "boost_family",
    representation = representation(
        fW      = "function"
))

setMethod("show", "boost_family", function(object) {
    cat("\n\t", object@name, "\n\n")
    cat("Loss function:", object@charloss, "\n")
})

Family <- function(ngradient, loss = NULL, risk = NULL,
                   offset = function(y, w) 
                       optimize(risk, interval = range(y), y = y, w = w)$minimum, 
                   check_y = function(y) y,
                   weights = c("any", "none", "zeroone", "case"), 
                   nuisance = function() return(NULL),
                   name = "user-specified", fW = NULL, response = function(f) NA)
{

    if (is.null(loss)) {
        charloss <- ""
        stopifnot(!is.null(risk))
    }
    if (is.null(risk)) {
        stopifnot(!is.null(loss))
        charloss <- paste(deparse(body(loss)), "\n")
        risk <- function(y, f, w = 1) sum(w * loss(y, f), na.rm = TRUE)
    }
    weights <- match.arg(weights)
    check_w <- function(w) {
        switch(weights, 
            "any" = TRUE,
            "none" = isTRUE(all.equal(unique(w), 1)),
            "zeroone" = isTRUE(all.equal(unique(w + abs(w - 1)), 1)),
            "case" = isTRUE(all.equal(unique(w - floor(w)), 0)))
    }
    RET <- new("boost_family", ngradient = ngradient, nuisance = nuisance,
               risk = risk, offset = offset, check_y = check_y, response = response,
               weights = check_w, name = name, charloss = charloss)
    if (!is.null(fW))
        RET <- new("boost_family_glm", ngradient = ngradient, nuisance = nuisance,
                   risk = risk, offset = offset, fW = fW, check_y = check_y, 
                   weights = check_w, name = name, response = response,
                   charloss= charloss)
    RET
}

### Gaussian (Regression)
Gaussian <- function()
    Family(ngradient = function(y, f, w = 1) y - f,
           loss = function(y, f) (y - f)^2,
           offset = weighted.mean,
           check_y = function(y) {
               if (!is.numeric(y) || !is.null(dim(y)))
                   stop("response is not a numeric vector but ",
                        sQuote("family = GaussReg()"))
               y
           },
           name = "Squared Error (Regression)",
           fW = function(f) return(rep(1, length = length(f))),
           response = function(f) f)
GaussReg <- Gaussian

### Gaussian (-1 / 1 Binary Classification)
GaussClass <- function()
    stop("Family", sQuote("GaussClass"), "has been removed")

### Laplace
Laplace <- function()
    Family(ngradient = function(y, f, w = 1) sign(y - f),
           loss = function(y, f) abs(y - f),
           check_y = function(y) {
               if (!is.numeric(y) || !is.null(dim(y)))
                   stop("response is not a numeric vector but ",
                        sQuote("family = Laplace()"))
               y
           },
           name = "Absolute Error")

### Binomial
# lfinv <- binomial()$linkinv
Binomial <- function()
    Family(ngradient = function(y, f, w = 1) {
               exp2yf <- exp(-2 * y * f)
               -(-2 * y * exp2yf) / (log(2) * (1 + exp2yf))
           },
           loss = function(y, f) {
               f <- pmin(abs(f), 36) * sign(f)
               p <- exp(f) / (exp(f) + exp(-f))
               y <- (y + 1) / 2
               -y * log(p) - (1 - y) * log(1 - p)
           },
           offset = function(y, w) {
               p <- weighted.mean(y > 0, w)
               1/2 * log(p / (1 - p))
           },
           fW = function(f) {
               f <- pmin(abs(f), 36) * sign(f)
               p <- exp(f) / (exp(f) + exp(-f))
               4 * p * (1 - p)
           },
           response = function(f) {
               f <- pmin(abs(f), 36) * sign(f)
               p <- exp(f) / (exp(f) + exp(-f))
               return(p)
           },
           check_y = function(y) {
               if (!is.factor(y))
                   stop("response is not a factor but ",
                           sQuote("family = Binomial()"))  
               if (nlevels(y) != 2)
                   stop("response is not a factor at two levels but ",
                           sQuote("family = Binomial()"))
               c(-1, 1)[as.integer(y)]
           },
           name = "Negative Binomial Likelihood")

### Poisson
Poisson <- function()
    Family(ngradient = function(y, f, w = 1) y - exp(f),
           loss = function(y, f) -dpois(y, exp(f), log = TRUE),
           offset = function(y, w) log(weighted.mean(y, w)),
           fW = function(f) exp(f),
           response = function(f) exp(f),
           check_y = function(y) {
               if (any(y < 0) || any((y - round(y)) > 0))
                   stop("response is not an integer variable but ",
                        sQuote("family = Poisson()"))
               y
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
           check_y = function(y) {
               if (!is.numeric(y) || !is.null(dim(y)))
                   stop("response is not a numeric vector but ",
                        sQuote("family = Huber()"))
               y
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
               c(-1, 1)[as.integer(y)]
           },
           name = "Adaboost Exponential Error")

### Cox proportional hazards model (partial likelihood)
CoxPH <- function() {
    plloss <- function(y, f, w) {
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
    }
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
               .Call("ngradientCoxPLik", time, event, f, w, package = "mboost")
           },
           risk = risk <- function(y, f, w = 1) -sum(plloss(y, f, w), na.rm = TRUE),
           offset = function(y, w) 
               optimize(risk, interval = c(0, max(y[,1], na.rm = TRUE)), 
                        y = y, w = w)$minimum, 
           check_y = function(y) {
               if (!inherits(y, "Surv"))
                   stop("response is not an object of class ", sQuote("Surv"),
                        " but ", sQuote("family = CoxPH()"))
               y
           },
           name = "Cox Partial Likelihood")
}

QuantReg <- function(tau = 0.5, qoffset = 0.5) {
    stopifnot(tau > 0 && tau < 1)
    stopifnot(qoffset > 0 && qoffset < 1)
    Family(
        ngradient = function(y, f, w = 1) 
            tau*((y - f) >= 0) - (1 - tau)*((y - f)<0) ,
        loss = function(y, f) tau*(y-f)*((y-f)>=0) - (1-tau)*(y-f)*((y-f)<0) ,
        offset = function(y, w = rep(1, length(y))) 
            quantile(y[rep(1:length(y), w)], qoffset),
        weights = "case",
        check_y = function(y) {
            if (!is.numeric(y) || !is.null(dim(y)))
                stop("response is not a numeric vector but ", 
                     sQuote("family = QuantReg()"))
            y
        },
        name = "Quantile Regression")
}

NBinomial <- function(nuirange = c(0, 100)) {
    sigma <- 1

    plloss <- function(sigma, y, f)
        - (lgamma(y + sigma) - lgamma(sigma) - lgamma(y + 1) +
           sigma * log(sigma) - sigma*log(exp(f) + sigma) + y * f -
           y * log(exp(f) + sigma))

    riskS <- function(sigma, y, fit, w = 1) 
        sum(w * plloss(y = y, f = fit, sigma = sigma))
    risk <- function(y, f, w = 1) 
       sum(w * plloss(y = y, f = f, sigma = sigma))

    ngradient <- function(y, f, w = 1) {
        sigma <<- optimize(riskS, interval = nuirange, 
                           y = y, fit = f, w = w)$minimum
        y - (y + sigma)/(exp(f) + sigma) * exp(f)
    }
	
    Family(ngradient = ngradient, risk = risk,
           check_y = function(y) {
               stopifnot(all.equal(unique(y - floor(y)), 0))
               y
           },
           nuisance = function() return(sigma),
           name = "Negative Negative-Binomial Likelihood")
}

PropOdds <- function(nuirange = c(-0.5, -1), offrange = c(-5, 5)) {

    sigma <- 0
    delta <- 0

    d2s <- function(delta)
        delta[1] + cumsum(c(0, exp(delta[-1])))

    plloss <- function(sigma, y, f, w = 1) {
        if (length(f) == 1) f <- rep(f, length(y))
        tmp <- lapply(1:(length(sigma) + 1), function(i) {
            if (i == 1) return(1 + exp(f - sigma[i]))
            if (i == (length(sigma) + 1)) 
                return(1 - 1/(1 + exp(f - sigma[i - 1])))
            return(1 / (1 + exp(f - sigma[i])) -  
                   1 / (1 + exp(f - sigma[i - 1])))
        })
        loss <- log(tmp[[1]]) * (y == levels(y)[1])
        for (i in 2:nlevels(y))
            loss <- loss - log(tmp[[i]]) * (y == levels(y)[i])
        return(loss)
    }

    riskS <- function(delta, y, fit, w = 1) 
        sum(w * plloss(y = y, f = fit, sigma = d2s(delta)))
    risk <- function(y, f, w = 1) 
        sum(w * plloss(y = y, f = f, sigma = sigma))

    ngradient <- function(y, f, w = 1) {
        delta <<- optim(par = delta, fn = riskS, y = y, 
                        fit = f, w = w, method = "BFGS")$par
        sigma <<- d2s(delta)
        if (length(f) == 1) f <- rep(f, length(y))
        ng <- sapply(1:(length(sigma) + 1), function(i) {
            if (i > 1 & i < (length(sigma) + 1)) {
                ret <- (1 - exp(2 * f - sigma[i - 1] - sigma[i])) / 
                   (1 + exp(f - sigma[i - 1]) + 
                        exp(f - sigma[i]) + 
                        exp(2 * f - sigma[i - 1] - sigma[i]))
            } else {
                if (i == 1) {
                    ret <- -1/(1 + exp(sigma[i] - f))
                } else {
                    ret <- 1 / (1 + exp(f - sigma[i - 1]))
                }
            }
            return(ret * (y == levels(y)[i]))
        })
        rowSums(ng)
    }

    offset <- function(y, w = 1) {
        delta <<- seq(from = nuirange[1], to = nuirange[2], 
                      length = nlevels(y) - 1)
        sigma <<- d2s(delta)
        optimize(risk, interval = offrange, y = y, w = w)$minimum
    }

    Family(ngradient = ngradient, 
           risk = risk, offset = offset,
           check_y = function(y) {
               stopifnot(is.ordered(y))
               y
           },
           nuisance = function() return(sigma),
           response = function(f) {
               ret <- sapply(1:(length(sigma) + 1), function(i) {
                   if (i == 1) return(1 / (1 + exp(f - sigma[i])))
                   if (i == (length(sigma) + 1))
                       return(1 - 1/(1 + exp(f - sigma[i - 1])))
                   return(1 / (1 + exp(f - sigma[i])) -  
                       1 / (1 + exp(f - sigma[i - 1])))
                   })
               ret
               })
}
