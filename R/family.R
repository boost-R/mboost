
setClass("boost_family", representation = representation(
    ngradient  = "function",
    risk       = "function",
    offset     = "function",
    check_y    = "function",
    weights    = "function",
    nuisance   = "function",
    response   = "function",
    rclass     = "function",
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
                   nuisance = function() return(NA),
                   name = "user-specified", fW = NULL,
                   response = function(f) NA,
                   rclass = function(f) NA)
{

    if (is.null(loss)) {
        charloss <- ""
        stopifnot(!is.null(risk))
    } else {
        charloss <- paste(deparse(body(loss)), "\n")
    }
    if (is.null(risk)) {
        stopifnot(!is.null(loss))
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
               risk = risk, offset = offset, check_y = check_y,
               response = response, rclass = rclass, weights = check_w,
               name = name, charloss = charloss)
    if (!is.null(fW))
        RET <- new("boost_family_glm", ngradient = ngradient, nuisance = nuisance,
                   risk = risk, offset = offset, fW = fW, check_y = check_y,
                   weights = check_w, name = name, response = response,
                   rclass = rclass, charloss= charloss)
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
                        sQuote("family = Gaussian()"))
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
           name = "Absolute Error",
           response = function(f) f)

link2dist <- function(link,  ...) {
    p <- get(paste("p", link, sep = ""))
    d <- get(paste("d", link, sep = ""))
    q <- get(paste("q", link, sep = ""))
    ret <- list(p = function(x) p(x, ...),
                d = function(x) d(x, ...),
                q = function(x) q(x, ...))
    attr(ret, "link") <- link
    ret
}

Binomial <- function(type = c("adaboost", "glm"),
                     link = c("logit", "probit", "cloglog", "cauchit", "log"), ...) {
    
    type <- match.arg(type)
    if (type == "adaboost") {
        return(Binomial_adaboost(link = link, ...))
    } else {
        return(Binomial_glm(link = link))
    }
}

### Binomial
# lfinv <- binomial()$linkinv
Binomial_adaboost <- function(link = c("logit", "probit", "cloglog", "cauchit", "log"), ...) {
    
    tmp <- try(match.arg(link), silent = TRUE)
    if (inherits(tmp, "try-error") ) {
        ## use old interface if link is not one of the above (i.e. if link is a distribution)
        link <- link2dist(link, ...)    
    } else {
        ## use new interface otherwise
        link <- make.link(tmp)   
    }
    
    biny <- function(y) {
        if (!is.factor(y))
            stop("response is not a factor but ",
                 sQuote("family = Binomial()"))
        if (nlevels(y) != 2)
            stop("response is not a factor at two levels but ",
                 sQuote("family = Binomial()"))
        return(c(-1, 1)[as.integer(y)])
    }
    
    if (inherits(link, "link-glm") && isTRUE(all.equal(link$name, "logit"))) {
        return(  
          # This is the classic logit family
          # ngradient refers to loss
          # log_2(1 + exp(-2yf))
          # with f = log(p/(1-p))/2 and y in (-1,1)
          # leading to coefs half the size as usually 
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
            rclass = function(f) (f > 0) + 1 ,
            check_y = biny,
            name = "Negative Binomial Likelihood (logit link)")
        )
    }
    
    if (inherits(link, "link-glm")) {
      # now for glm type links that are not logit
      # loss is now 
      # -y * log(p) - (1 - y) * log(1 - p)
      # coefficients hace the usual size! 
      return(Family(ngradient = function(y, f, w = 1) {
        y <- (y + 1) / 2
        p <- link$linkinv(f)
        link$mu.eta(f) * (y / p - (1 - y) / (1 - p))
      },
      loss = function(y, f) {
        p <- link$linkinv(f)
        y <- (y + 1) / 2
        -y * log(p) - (1 - y) * log(1 - p)
      },
      offset = function(y, w) {
        p <- weighted.mean(y > 0, w)
        link$linkfun(p)
      },
      response = function(f) {
        p <- link$linkinv(f)
        return(p)
      },
      rclass = function(f) (f > 0) + 1 ,
      check_y = biny,
      name = paste("Negative Binomial Likelihood --",
                   link$name, "link")))
      
      
    }
    
    trf <- function(f) {
      thresh <- -link$q(.Machine$double.eps)
      pmin(pmax(f, -thresh), thresh)
    }
    
    # now with link2dist (distributions)
    return(Family(ngradient = function(y, f, w = 1) {
        y <- (y + 1) / 2
        p <- link$p(trf(f))
        d <- link$d(trf(f))
        d * (y / p - (1 - y) / (1 - p))
    },
    loss = function(y, f) {
        p <- link$p(trf(f))
        y <- (y + 1) / 2
        -y * log(p) - (1 - y) * log(1 - p)
    },
    offset = function(y, w) {
        p <- weighted.mean(y > 0, w)
        link$q(p)
    },
    response = function(f) {
        p <- link$p(trf(f))
        return(p)
    },
    rclass = function(f) (f > 0) + 1 ,
    check_y = biny,
    name = paste("Negative Binomial Likelihood --",
                 attr(link, "link"), "Link")))
}

### Additional Binomial family 
### Works not only with two-level factors but also 
### with a two-column matrix of number of successes and
### number of failures (see github issue #34 by fabian-s). 
### In case of factors or binary vectors it uses the standard
### 0 and 1 coding; the coefficients hence have the same level
### as the ones resulting from glm() with family = "binomial".
Binomial_glm <- function(link = c("logit", "probit", "cloglog", "cauchit", "log")) {
    
    link <- match.arg(link)
    link <- make.link(link)
    
    y_check <- function(y) {
        if ((is.matrix(y) && NCOL(y)!=2)){
            stop("response should be either a two-column matrix (no. successes ",
                 "and no. failures) or a two level factor or a vector of 0 and 1's ",
                 "for this family")
        }
        if(is.factor(y)){
            if (nlevels(y) != 2) 
                stop("response should be either a two-column matrix ",
                     "(no. successes  and  no. failures) or a two level ", 
                     "factor or a vector of 0 and 1's for this family")
            y <- c(0, 1)[as.integer(y)]
        }
        if(!is.matrix(y)){
            if(!all(y %in% c(0,1))) 
                stop("response should be either a two-column matrix ",
                     "(no. successes  and  no. failures) or a two level ",
                     "factor or a vector of 0 and 1's for this family")
            y <- cbind(y, 1-y)
        }
        return(y)  
    }
    
    loss <- function(y, f, w = 1) {
        ntrials <- rowSums(y)
        y <- y[,1]
        p <- link$linkinv(f)
        -dbinom(x = y, size = ntrials, prob = p, log = log)
    }
    risk <- function(y, f, w = 1) {
        sum(w * loss(y = y, f = f))
    }
    
    ngradient <- function(y, f, w = 1) {
        ntrials <- rowSums(y)
        y <- y[,1]
        p <- link$linkinv(f)
        ngr <-  link$mu.eta(f)*(y - ntrials * p)/(p * (1 - p))
    }
    
    offset <-function(y, w = 1) {
        ntrials <- rowSums(y)
        y <- y[,1]
        p <- mean((y + 0.5)/(ntrials + 1))
        return(link$linkfun(p))
    }
    
    Family(ngradient = ngradient, risk = risk, loss = loss, check_y = y_check, 
           response = function(f) link$linkinv(f), offset = offset,
           rclass = function(f) (f > 0) + 1 ,
           name = "Binomial Distribution (similar to glm)")
}

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
                                  paste("(with d = ", dtxt, ")", sep = ""))),
           response = function(f) f)
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
           rclass = function(f) (f > 0) + 1 ,
           name = "Adaboost Exponential Error")

### Cox proportional hazards model (partial likelihood)
CoxPH <- function() {
    plloss <- function(y, f, w) {
        time <- y[,1]
        event <- y[,2]
        n <- length(time)
        if (length(f) == 1) 
            f <- rep(f, n)
        if (length(w) == 1) 
            w <- rep(w, n)
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
               w[is.na(f)] <- 0.0
               f[is.na(f)] <- 0.0
               .Call("ngradientCoxPLik", time, event, f, w, PACKAGE = "mboost")
           },
           risk = risk <- function(y, f, w = 1) -sum(plloss(y, f, w), na.rm = TRUE),
           offset = function(y, w = 1) 0, ## perhaps use something different
                       ## Note: offset cannot be computed from Cox Partial LH as
                       ## PLH doesn't depend on constant
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

    Family(ngradient = ngradient, 
           risk = risk,
           offset = function(y, w) log(weighted.mean(y, w)),
           check_y = function(y) {
               stopifnot(all.equal(unique(y - floor(y)), 0))
               y
           },
           nuisance = function() return(sigma),
           name = "Negative Negative-Binomial Likelihood",
           response = function(f) exp(f))
}

PropOdds <- function(nuirange = c(-0.5, -1), offrange = c(-5, 5)) {

    sigma <- 0
    delta <- 0

    d2s <- function(delta)
        delta[1] + cumsum(c(0, exp(delta[-1])))

    plloss <- function(sigma, y, f, w = 1) {
        if (length(f) == 1) 
            f <- rep(f, length(y))
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
        if (length(f) == 1) 
            f <- rep(f, length(y))
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
        optimize(risk, interval = offrange, y = y, w = w)$minimum
    }

    response <- function(f) {
        ret <- sapply(1:(length(sigma) + 1), function(i) {
            if (i == 1) return(1 / (1 + exp(f - sigma[i])))
            if (i == (length(sigma) + 1))
            return(1 - 1/(1 + exp(f - sigma[i - 1])))
            return(1 / (1 + exp(f - sigma[i])) -
               1 / (1 + exp(f - sigma[i - 1])))
            })
            ret
    }

    check_y <- function(y) {
        if (!is.ordered(y))
            stop("response must be an ordered factor")
        ## initialize thresholds:
        delta <<- seq(from = nuirange[1], to = nuirange[2],
                      length = nlevels(y) - 1)
        sigma <<- d2s(delta)
        y
    }

    Family(ngradient = ngradient,
           risk = risk, offset = offset,
           check_y = check_y,
           nuisance = function() return(sigma),
           response = response,
           rclass = function(f) apply(response(f), 1, which.max))
}

Weibull <- function(nuirange = c(0, 100)) {
    sigma <- 1

    plloss <- function(sigma, y, f){
        fw <- function(pred){
            exp(pred - exp(pred))
            }
        Sw <- function(pred){
            exp(-exp(pred))
            }
        time <- y[,1]
        ### log(0) won't fly
        time[time == 0] <- sort(unique(time))[2] / 10
        lnt <- log(time)
        Sevent <- y[,2]
        eta <- (lnt - f) / sigma
        - Sevent * log(fw(eta) / sigma) - (1 - Sevent) * log(Sw(eta))
        }

    riskS <- function(sigma, y, fit, w = 1)
        sum(w * plloss(y = y, f = fit, sigma = sigma))
    risk <- function(y, f, w = 1)
       sum(w * plloss(y = y, f = f, sigma = sigma))

    ngradient <- function(y, f, w = 1) {
        sigma <<- optimize(riskS, interval = nuirange,
                           y = y, fit = f, w = w)$minimum

        time <- y[,1]
        ### log(0) won't fly
        time[time == 0] <- sort(unique(time))[2] / 10
        lnt <- log(time)
        event <- y[,2]
        eta <- (lnt - f) / sigma
        (event * (exp(eta) - 1) + (1 - event) * exp(eta)) / sigma
        }

    Family(ngradient = ngradient, risk = risk,
           offset = function(y, w)
               optimize(risk, interval = c(0, max(log(y[,1]), na.rm = TRUE)),
                        y = y, w = w)$minimum,
           check_y = function(y) {
               if (!inherits(y, "Surv"))
                   stop("response is not an object of class ", sQuote("Surv"),
                        " but ", sQuote("family = Weibull()"))
               y
           },
           nuisance = function() return(sigma),
           name = "Negative Weibull Likelihood",
           response = function(f) exp(f))
}

Loglog <- function(nuirange = c(0, 100)) {
    sigma <- 1

    plloss <- function(sigma, y, f){
        fw <- function(pred){
            exp(pred) / (1 + exp(pred))^2
            }
        Sw <- function(pred){
            1 / (1 + exp(pred))
            }

        time <- y[,1]
        ### log(0) won't fly
        time[time == 0] <- sort(time)[2] / 10
        lnt <- log(time)
        Sevent <- y[,2]
        eta <- (lnt - f) / sigma
        - Sevent * log(fw(eta) / sigma) - (1 - Sevent) * log(Sw(eta))
        }

    riskS <- function(sigma, y, fit, w = 1)
        sum(w * plloss(y = y, f = fit, sigma = sigma))
    risk <- function(y, f, w = 1)
       sum(w * plloss(y = y, f = f, sigma = sigma))

    ngradient <- function(y, f, w = 1) {
        sigma <<- optimize(riskS, interval = nuirange,
                           y = y, fit = f, w = w)$minimum
        time <- y[,1]
        ### log(0) won't fly
        time[time == 0] <- sort(unique(time))[2] / 10
        lnt <- log(time)
        event <- y[,2]
        eta <- (lnt - f) / sigma
        nom <- (exp(-eta) + 1)
        (event * (2/nom - 1) + (1 - event) / nom) / sigma
        }

    Family(ngradient = ngradient, risk = risk,
           offset = function(y, w)
               optimize(risk, interval = c(0, max(log(y[,1]), na.rm = TRUE)),
                        y = y, w = w)$minimum,
           check_y = function(y) {
               if (!inherits(y, "Surv"))
                   stop("response is not an object of class ", sQuote("Surv"),
                        " but ", sQuote("family = Loglog()"))
               y
           },
           nuisance = function() return(sigma),
           name = "Negative Log Logistic Likelihood",
           response = function(f) exp(f))
}

Lognormal <- function(nuirange = c(0, 100)) {
    sigma <- 1

    plloss <- function(sigma, y, f){
        fw <- function(pred){
            dnorm(pred)
            }
        Sw <- function(pred){
            1 - pnorm(pred)
            }
        time <- y[,1]
        ### log(0) won't fly
        time[time == 0] <- sort(time)[2] / 10
        lnt <- log(time)
        Sevent <- y[,2]
        eta <- (lnt - f) / sigma
        - Sevent * log(fw(eta) / sigma) - (1 - Sevent) * log(Sw(eta))
        }

    riskS <- function(sigma, y, fit, w = 1)
        sum(w * plloss(y = y, f = fit, sigma = sigma))
    risk <- function(y, f, w = 1)
       sum(w * plloss(y = y, f = f, sigma = sigma))

    ngradient <- function(y, f, w = 1) {
        sigma <<- optimize(riskS, interval = nuirange,
                           y = y, fit = f, w = w)$minimum
        time <- y[,1]
        ### log(0) won't fly
        time[time == 0] <- sort(unique(time))[2] / 10
        lnt <- log(time)
        event <- y[,2]
        eta <- (lnt - f) / sigma
        (event * eta + (1 - event) * dnorm(eta) / (1 - pnorm(eta))) / sigma
        }

    Family(ngradient = ngradient, risk = risk,
           offset = function(y, w)
               optimize(risk, interval = c(0, max(log(y[,1]), na.rm = TRUE)),
                        y = y, w = w)$minimum,
           check_y = function(y) {
               if (!inherits(y, "Surv"))
                   stop("response is not an object of class ", sQuote("Surv"),
                        " but ", sQuote("family = Lognormal()"))
               y
           },
           nuisance = function() return(sigma),
           name = "Negative Lognormal Likelihood",
           response = function(f) exp(f))
}

ExpectReg <- function (tau = 0.5) {
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
                  sQuote("family = ExpectReg()"))
            y
        },
        name = "Expectile Regression")
}

## (1 - AUC)-Loss
AUC <- function() {

	approxGrad <- function(x) {
		ifelse(abs(x) >= 1, 0, ifelse(x >= 0, (1 - x), (x + 1)))
	}
	approxLoss <- function(x) {
		ret <- ifelse(x >= 0, 1 - (1 - x)^2/2, (x + 1)^2/2)
		ifelse(x < -1, 0, ifelse(x > 1, 1, ret))
	}

	ind1 <- 0
	ind0 <- 0
	n0 <- 0
	n1 <- 0
	M0 <- matrix(0,0,0)


	Family(ngradient = function(y, f, w = 1) {
				# if there are weights!=1, observations have to be
				# repeated/omitted accordingly
				if (any(w != 1)) {
					ind1 <- which(rep(y, w) == 1)
					ind0 <- which(rep(y, w) == -1)
					n1 <- length(ind1)
					n0 <- length(ind0)
				}
				
	            if (length(f) == 1)
	                f <- rep(f, n1 + n0)
	            # skip this in the first iteration
	            if (length(unique(f)) != 1) 
	                # scale scores s.t. a gradient of zero makes sense for
	                # differences in f that are bigger than +/-1
	                f <- f/sd(f)
				
				M0 <<- (matrix(f[ind1], nrow = n0, ncol = n1, byrow = TRUE) -
							f[ind0])
				M1 <- approxGrad(M0)
				ng <- vector(n1 + n0, mode = "numeric")
				ng[ind1] <- colSums(M1)
				ng[ind0] <- rowSums(-M1)
				return(ng)
			}, loss = function(y, f, w = 1) {
				# if there are weights!=1, observations have to be
				# repeated/omitted accordingly and M0 must be recomputed
				if (any(w != 1)) {
					ind1 <- which(rep(y, w) == 1)
					ind0 <- which(rep(y, w) == -1)
					n1 <- length(ind1)
					n0 <- length(ind0)
					if (length(f) == 1)
					    f <- rep(f, n1 + n0)
					# skip this in the first iteration
					if (length(unique(f)) != 1) 
					    f <- f/sd(f)
					M0 <- (matrix(f[ind1], nrow = n0, ncol = n1, byrow = TRUE) -
								f[ind0])
				}
				# 1 - approxAUC
				1 - (sum(approxLoss(M0)))/(n1 * n0)
			}, risk = function(y, f, w = 1) {
				# if there are weights!=1, observations have to be
				# repeated/omitted accordingly and M0 must be recomputed
				if (any(w != 1)) {
					ind1 <- which(rep(y, w) == 1)
					ind0 <- which(rep(y, w) == -1)
					n1 <- length(ind1)
					n0 <- length(ind0)
					if (any(c(n0, n1) == 0)) {
						warning("response is constant - AUC is 1.")
						return(0)
					}
					M0 <- (matrix(f[ind1], nrow = n0, ncol = n1, byrow = TRUE) -
								f[ind0])
				}
				# 1 - AUC
				return(1 - (sum(as.numeric(M0 > 0)))/(n1 * n0))
			}, weights = "case",
			offset = function(y, w) {
				0
			}, check_y = function(y) {
				if (!is.factor(y))
					stop("response is not a factor but ", sQuote("family = AUC()"))
				if (nlevels(y) != 2)
					stop("response is not a factor at two levels but ",
							sQuote("family = AUC()"))
				if (length(unique(y)) != 2)
					stop("only one class is present in response.")
				# precompute to speed up ngradient (and warn about conflicting assignments)
				ind1 <<- which(y == levels(y)[2])
				ind0 <<- which(y == levels(y)[1])
				n1 <<- length(ind1)
				n0 <<- length(ind0)
				c(-1, 1)[as.integer(y)]
			},
      rclass = function(f) (f > 0) + 1,
      name = paste("(1 - AUC)-Loss") )
}



GammaReg <- function(nuirange = c(0, 100)) {
    sigma <- 1
    plloss <- function(sigma, y, f){
        lgamma(sigma) + sigma * y * exp(-f) - sigma * log(y) -
              sigma * log(sigma) + sigma * f
    }
    riskS <- function(sigma, y, fit, w = 1)
        sum(w * plloss(y = y, f = fit, sigma = sigma))
    risk <- function(y, f, w = 1)
       sum(w * plloss(y = y, f = f, sigma = sigma))
    ngradient <- function(y, f, w = 1) {
        sigma <<- optimize(riskS, interval = nuirange,
                           y = y, fit = f, w = w)$minimum
        sigma * y * exp(-f) - sigma
    }
    Family(ngradient = ngradient,
           risk = risk,
           offset = function(y, w){
               optimize(risk, interval = c(0, max(y^2, na.rm = TRUE)),
                        y = y, w = w)$minimum
           },
           check_y = function(y){
               if (!is.numeric(y) || !is.null(dim(y)))
                   stop("response is not a numeric vector but ",
                        sQuote("family = GammaReg()"))
               if (any(y < 0))
                   stop("response is not positive but ",
                        sQuote("family = GammaReg()"))
               y
           },
           nuisance = function() return(sigma),
           name = "Negative Gamma Likelihood",
           response = function(f) exp(f))
}

### Hinge Loss
HingeLoss <- function(alphaCost = 0.5) {
    if (alphaCost <= 0 | alphaCost >= 1 | !is.numeric(alphaCost) |
        length(alphaCost) != 1)
        stop("cost parameter must be a number within the range (0, 1)")
    Family(ngradient = function(y, f, w = 1){
               ngr <- ifelse(y > 0, 1 - alphaCost, alphaCost) * y
               ngr[(1 - y * f) < 0] <- 0
               return(ngr)
               },
           loss = function(y, f) ifelse(y > 0, 1 - alphaCost, alphaCost) *
               pmax(0, 1 - y * f),
           offset = function(y, w) {
               p <- weighted.mean(y > 0, w)
               1/2 * log(p / (1 - p))
           },
           check_y = function(y) {
               if (!is.factor(y))
                   stop("response is not a factor but ",
                           sQuote("family = HingeLoss()"))
               if (nlevels(y) != 2)
                   stop("response is not a factor at two levels but ",
                           sQuote("family = HingeLoss()"))
               c(-1, 1)[as.integer(y)]
           },
           rclass = function(f) (f > (2 * alphaCost - 1)) + 1,
           name = "Hinge Loss")
}

# Gehan family for 'mboost' package
#   Brent Johnson, with comments from Benjamin Hofner
#   Reference : Johnson BA and Long Q (2011) Survival ensembles by the sum of pairwise differences with application to
#		lung cancer microarray studies, Annals of Applied Statistics, 5:1081-1101.
Gehan <- function() {
    Family(ngradient = function(y, f, w) {
               Z <- as.vector(log(y[,1]))
               Del <- as.vector(y[,2])
               n1 <- length(Z)
               if(length(f)==1) f <- rep(f,n1)
               if(length(w)==1) w <- rep(w,n1)
               idx <- rep(1:n1,w)
               deltaSubset <- Del[idx]
               resSubset <- Z[idx] - f[idx]  		# e_j
               n2 <- length(deltaSubset)		# n2=sum(w)
               tmpg.fun <- function(I, res, del) {
                   ediff <- res[I] - res
                   arg2 <- sum(del * (ediff >= 0) )
                   arg1 <- ifelse(del[I] > 0.5, sum(ediff <= 0), 0)
                   arg1 - arg2
               }

               g <- rep(0,n1)		## here, the gradient is n1-by-1 with 0s for obs out of the Subset
               g[idx] <- unlist(lapply(1:n2, tmpg.fun, resSubset, deltaSubset))
               return(-1 * g / n2)      ## here, divide by n2=sum(w)
           },
           loss = gehan.loss <- function(y, f, w) {
               Z <- as.vector(log(y[,1]))
               Del <- as.vector(y[,2])
               n1 <- length(Z)
               if(length(f)==1) rep(f, n1)
               if(length(w)==1) rep(w, n1)
               idx <- rep(1:n1,w)
               deltaSubset <- Del[idx]
               resSubset <- Z[idx] - f[idx]
               n2 <- length(deltaSubset)
               tmpl.fun <- function(I, res, del) {
                   if(del[I] < 0.5) out1 <- 0
                   else {
                       ediff <- res[I] - res
                       out1 <- -1 * sum(ediff[ediff <= 0])
                   }
                   out1
               }
               out2 <- rep(0,n1)
               out2[idx] <- unlist(lapply(1:n2, tmpl.fun, resSubset, deltaSubset))/n2
               out2
           },
           risk = risk <- function(y, f, w = 1){
               sum(gehan.loss(y, f, w), na.rm=TRUE)
           },
           weights = "case",
           offset = function(y,w){
               optimize(risk,
                        interval = c(0, max(log(y[,1]), na.rm=TRUE)), y = y, w = w)$minimum
           },
           check_y = function(y) {
               if (!inherits(y,"Surv"))
                   stop("response is not an object of class ",
                        sQuote("Surv"), " but ",
                        sQuote("family = Gehan()"))
               y
           },
           name = "Gehan loss")
}


### Hurdle model, negative binomial non-zero component with log link
Hurdle <- function(nuirange = c(0, 100)){
    sigma <- 1
    plloss <- function(sigma, y, f)
                  - (y * (f - log(1 / sigma + exp(f))) -
                  (log(1 / sigma + exp(f)) + log(sigma)) / sigma +
                  lgamma(y + 1 / sigma) - lgamma(y + 1) - lgamma(1 / sigma) -
                  log(1 - (1 + exp(f) * sigma)^{-1 / sigma}))
    riskS <- function(sigma, y, fit, w = 1)
                 sum(w * plloss(y = y, f = fit, sigma = sigma))
    risk <- function(y, f, w = 1) sum(w * plloss(y = y, f = f, sigma = sigma))
    ngradient = function(y, f, w = 1){
                    sigma <<- optimize(riskS, interval = nuirange, y = y,
                                       fit = f, w = w)$minimum
                    y / sigma / (exp(f) + 1 / sigma) - 1 / sigma /
                    (exp(-f) / sigma + 1) - exp(f) * (1 + exp(f) *
                    sigma)^{-1 / sigma - 1} / (1 - (1 + exp(f) *
                    sigma)^{-1 / sigma})
                    }
    Family(ngradient = ngradient, risk = risk, 
           offset = function(y, w) log(weighted.mean(y, w)),
           check_y = function(y) {
               stopifnot(all.equal(unique(y - floor(y)), 0))
               y
           }, 
           nuisance = function() return(sigma),
           name = "Hurdle model, negative binomial non-zero part",
           response = function(f) exp(f))
}

### multinomial logit model
### NOTE: this family can't be applied out-of-the box
### the rhs for formula needs to read
### ~ bl1 %O% bl0 + bl2 %O% bl1 + ...
### where bl1 is some linear baselearner
### and bl0 is a baselearner with design and penalty term
### equal to the unit matrix for number of levels - 1
### See ?Multinom
Multinomial <- function() {
    lev <- NULL
    biny <- function(y) {
        if (!is.factor(y))
            stop("response is not a factor but ",
                  sQuote("family = Multinomial()"))
        lev <<- levels(y)
        as.vector(model.matrix(~ y - 1)[,-length(lev)])
    }
    return(Family(ngradient = function(y, f, w = 1) {
               if (length(f) != length(y))
                   stop("predictor doesn't correspond to multinomial logit model; see ?Multinomial")
               f <- pmin(abs(f), 36) * sign(f)
               p <- matrix(exp(f), ncol = length(lev) - 1)
               p <- as.vector(p / (1 + rowSums(p)))
               y - p
           },
           loss = function(y, f) {
               f <- pmin(abs(f), 36) * sign(f)
               p <- matrix(exp(f), ncol = length(lev) - 1)
               p <- as.vector(p / (1 + rowSums(p)))
               -y * log(p)
           },
           offset = function(y, w) {
               return(rep(0, length(y)))
           },
           response = function(f) {
               f <- pmin(abs(f), 36) * sign(f)
               p <- matrix(exp(f), ncol = length(lev) - 1)
               p <- cbind(p, 1) / (1 + rowSums(p))
               colnames(p) <- lev
               return(p)
           },
           rclass = function(f) {
               f <- pmin(abs(f), 36) * sign(f)
               p <- matrix(exp(f), ncol = length(lev) - 1)
               p <- cbind(p, 1) / (1 + rowSums(p))
               apply(p, 1, which.max)
           },
           check_y = biny,
           name = "Negative Multinomial Likelihood"))
}

### Cindex family to fit the concordance index for survival data
### The family is based on an estimator proposed by Uno et al. which
### incorporates "inverse probability of censoring weights" which can be 
### either computed inside the family (based on all observations with 
### non-zero weights) or computed before fitting and passed as an argument 
### ("ipcw"), which could be better in cases of cross-validation or resampling.


Cindex <- function (sigma = 0.1, ipcw = 1) {
  
  approxGrad <- function(x) {    ## sigmoid function for gradient
    exp(x/sigma) / (sigma * (1 + exp(x/sigma))^2) 
  }
  
  approxLoss <- function(x) { ## sigmoid function for loss
    1 / (1 + exp(x / sigma))
  }
  
  compute_weights <- function(y, w = 1, ipcw = rep(1,nrow(y))){  ## compute weights
    
    if(all(ipcw == 1))
    {  
      ipcw_wow <- IPCweights(y[w != 0,])
      ipcw <- numeric(nrow(y))
      ipcw[w!=0] <- ipcw_wow
    }
    
    survtime <- y[,1]
    n <- nrow(y)
    wweights <- matrix( (ipcw)^2, nrow = n, ncol = n)
    weightsj <- matrix(survtime, nrow = n, ncol = n) 
    weightsk <- matrix(survtime, nrow = n, ncol = n, byrow = TRUE) 
    weightsI <- ifelse(weightsj == weightsk, .5, (weightsj < weightsk) + 0) - diag(.5, n,n)
    wweights <- wweights * weightsI
    Wmat <- w %o% w
    wweights <- wweights * Wmat
    wweights <- wweights / sum(wweights)
    rm(weightsI); rm(weightsk); rm(weightsj)
    return(wweights)
  }
  
  ngradient = function(y, f, w = 1) {   ## negative gradient
    if (!all(w %in% c(0,1)))
      stop(sQuote("weights"), " must be either 0 or 1 for family ",
           sQuote("UnoC"))
    
    survtime <- y[,1]
    event <- y[,2]
    n <- length(survtime)
    
    if (length(w) == 1) 
        w <- rep(1, n)
    if (length(f) == 1)
        f <- rep(f, n)
    
    etaj <- matrix(f, nrow = n, ncol = n, byrow = TRUE)
    etak <- matrix(f, nrow = n, ncol = n)
    etaMat <- etak - etaj
    rm(etaj); rm(etak);
    
    weights_out <- compute_weights(y, w, ipcw = ipcw)
    M1 <- approxGrad(etaMat) * weights_out
    ng <- colSums(M1) - rowSums(M1)
    
    return(ng)
  }
  
  risk = function(y, f, w = 1) { ## empirical risk
    survtime <- y[,1]
    event <- y[,2]
    n <- length(survtime)
    
    if (length(f) == 1)
      f <- rep(f, n)
    
    etaj <- matrix(f, nrow = n, ncol = n, byrow = TRUE)
    etak <- matrix(f, nrow = n, ncol = n)
    etaMat <- (etak - etaj)
    rm(etaj); rm(etak);
    weights_out <- compute_weights(y, w, ipcw = ipcw )
    M1 <- approxLoss(etaMat) * weights_out
    return(- sum(M1))
  }
  
  Family(   ## build the family object
    ngradient = ngradient,
    risk = risk,
    weights = "zeroone",
    offset = function(y, w = 1) {0},
    check_y = function(y) {
      if (!inherits(y, "Surv"))
        stop("response is not an object of class ", sQuote("Surv"),
             " but ", sQuote("family = UnoC()"))
      y},
    rclass = function(f){},
    name = paste("Concordance Probability by Uno"))
}


# RCG model for bounded response variables, log link
RCG <- function (nuirange = c(0, 1), offrange = c(-5, 5)) {
    sigma <- c(1, 0.2)
    plloss <- function(sigma, y, f, w = 1) {
        alpha <- sigma[1]
        rho <- sigma[2]
        if (length(f) == 1)
            f <- rep(f, length(y))
        my <- 1 + (exp(f) - 1) * y
        dmy <- my ^ 2 - 4 * rho * exp(f) * y * (1 - y)
        fy <- lgamma(2 * alpha) - 2 * lgamma(alpha) + alpha * f + 
            alpha * log(1 - rho) + log(my) + (alpha - 1) * log(y * (1 - y)) - 
            (alpha + 0.5) * log(dmy)
        return(-fy)
    }
    riskS <- function(sigma, y, fit, w = 1)
        sum(w * plloss(y = y, f = fit, sigma = sigma))
    risk <- function(y, f, w = 1)
        sum(w * plloss(y = y, f = f, sigma = sigma))
    ngradient <- function(y, f, w = 1) {
        sigma <<- optim(par = sigma, fn = riskS, y = y, fit = f,
            lower = c(0.001, 0.001), upper = c(100, 0.99),
            w = w, method = "L-BFGS-B")$par
        if (length(f) == 1)
            f <- rep(f, length(y))
        alpha <- sigma[1]
        rho <- sigma[2]
        my <- 1 + (exp(f) - 1) * y
        dmy <- my ^ 2 - 4 * rho * exp(f) * y * (1 - y)
        return(alpha + (y * exp(f)) / my - (2 * alpha + 1) * exp(f) *
                   y * (my - 2 * rho * (1 - y)) / dmy)
    }
    offset <- function(y, w = 1)
        0
    dRCG <- function(f, sigma) {
        alpha <- sigma[1]
        rho <- sigma[2]
        densw <- function(w, alpha, rho, f)
            w * gamma(2 * alpha) / (gamma(alpha)) ^ 2 * (1 - rho) ^ alpha *
            (exp(f))^alpha * ((exp(f) - 1) * w + 1) * (w * (1 - w))^(alpha - 1) /
            ((((exp(f) - 1) * w + 1) ^ 2 - 4 * rho * exp(f) * w * (1 - w)) ^ (alpha + 0.5))
        res <- rep(0, length(f))
        for (k in 1:length(f))
            res[k] <- integrate(densw, lower = 0, upper = 1, alpha, rho, f[k])$value
        return(res)
    }
    response <- function(f, sigma)
        dRCG(f, sigma)
    Fisher <- function(y, alpha, rho, X, coefs) {
        X <- data.frame(X)
        f.i <- matrix(0, ncol = length(coefs), nrow = length(coefs))
        for (i in 1:length(y)) {
            wobs <- as.numeric(y[i])
            x <- as.matrix(X[i,])
            xt <- t(x)
            xtx <- xt %*% x
            exb <- as.numeric(exp(x %*% coefs))
            D1 <- (exb - 1) * wobs + 1
            D2 <- D1 ^ 2 - 4 * rho * exb * wobs * (1 - wobs)
            f.i <- f.i + (D1 * wobs * exb - (wobs * exb) ^ 2) * xtx / (D1 ^ 2) -
                (alpha + 0.5) / (D2 ^ 2) * 
                (2 * wobs * exb * xtx * (2 * wobs * exb - wobs + 1 - 2 * rho * (1 - wobs)) 
                 * D2 - xtx * (2 * wobs * exb * (D1 - 2 * rho * (1 - wobs)))^2)
        }
        Fisher.information <- -f.i
        return(Fisher.information)
    }
    
    Family(
        ngradient = ngradient,
        risk = risk,
        offset = offset,
        # Fisher = Fisher,
        check_y = function(y) {
            if(!all(y < 1 & y > 0))
                stop("y must be < 1 and > 0")
            y
        },
        nuisance = function()
            return(sigma),
        response = function(f)
            response(f, sigma = sigma),
        rclass = function(f)
            f,
        name = "Ratio of Correlated Gammas (RCG)"
    )
}
