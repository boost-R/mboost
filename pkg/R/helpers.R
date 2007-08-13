
### two possible interpretations of weights:
### 1) case counts: observation i is w_i times in the sample
### 2) relative weights: observation i is given weight w_i
rescale_weights <- function(w) {
    if (max(abs(w - floor(w))) < sqrt(.Machine$double.eps))
        return(w)
    return(w / sum(w) * sum(w > 0))
}

### data preprocessing
boost_dpp <- function(formula, data, weights = NULL, na.action = na.omit, ...) {

    if (is.null(weights)) {
        env <- ModelEnvFormula(formula, data, na.action = na.action, ...)
        weights <- rep.int(1, NROW(env@get("response")))
    } else {
        env <- ModelEnvFormula(formula, data, na.action = na.action, 
                               other = list(weights = ~ weights), ...)
        weights <- env@get("weights")[[1]]
    }
    y <- env@get("response")
    if (length(y) != 1)
        stop("cannot deal with multivariate response variables")
    y <- y[[1]]
    x <- env@get("designMatrix")

    if (is.factor(y)) {
        if (nlevels(y) != 2)
            stop("not a binary classification problem")
        yfit <- as.numeric(y) - 1
        yfit[yfit == 0] <- -1
    } else {
        yfit <- y
    }

    RET <- gb_xyw(x, y, weights)
    RET$formula <- formula
    RET$menv <- env
    class(RET) <- "boost_data"
    RET
}

gb_xyw <- function(x, y, w) {

    if (is.null(w))
        w <- rep.int(1, NROW(x))

    if (is.factor(y)) {
        if (nlevels(y) != 2)
            stop("not a binary classification problem")
        yfit <- as.numeric(y) - 1
        yfit[yfit == 0] <- -1
    } else {
        yfit <- y
    }

    list(x = x, y = y, yfit = yfit, w = w)
}

### check measurement scale of response for some losses
check_y_family <- function(y, family) {

    if (isTRUE(all.equal(attributes(family), 
                  attributes(Binomial())))) {
        if (!is.factor(y))
            warning("response is not a factor but ", 
                    sQuote("family = Binomial()"))
        if (nlevels(y) != 2)
            warning("response is not a factor at two levels but ", 
                    sQuote("family = Binomial()"))
    }
    if (isTRUE(all.equal(attributes(family), 
                  attributes(CoxPH())))) {
        if (!inherits(y, "Surv"))
            stop("response is not an object of class ", sQuote("Surv"), 
                 " but ", sQuote("family = CoxPH()"))
    }
    if (isTRUE(all.equal(attributes(family), 
                  attributes(Poisson())))) {
        if (any(y < 0) || any((y - round(y)) > 0))
            stop("response is not an integer variable but ", 
                 sQuote("family = Poisson()"))
    }
}

### check for negative gradient corresponding to L2 loss
checkL2 <- function(object)
    isTRUE(all.equal(attributes(object$family)[-c(4:8)], 
                     attributes(GaussReg())[-c(4:8)]))

gm <- function(object) UseMethod("gm")

gm.glmboost <- function(object) {

    mstop <- nrow(object$ensemble)
    x <- object$data$x
    if (object$control$center) x <- object$data$center(x)
    RET <- matrix(0, nrow = NROW(x), ncol = mstop)

    jsel <- object$ensemble[,"xselect"]
    cf <- object$ensemble[,"coef"] * object$control$nu 

    for (m in 1:mstop)
        RET[,m] <- cf[m] * x[,jsel[m]]

    RET[,1] <- RET[,1] + object$offset

    return(RET)
}

gm.gamboost <- function(object) {

    mstop <- nrow(object$ensemble)
    x <- object$data$input
    RET <- matrix(0, nrow = NROW(x[[1]]), ncol = mstop)
    nu <- object$control$nu

    for (m in 1:mstop)
        RET[,m] <- nu * fitted(object$ensembless[[m]])

    RET[,1] <- RET[,1] + object$offset

    return(RET) 
}

### partial fits
gamplot <- function(object, newdata = NULL) {

     if (is.null(newdata)) {
         x <- object$data$input
         pr <- function(obj) fitted(obj)
     } else {
         x <- newdata
         pr <- function(obj) predict(obj, newdata = x)
     }
     lp <- matrix(0, ncol = length(x), nrow = NROW(x[[1]]))
     ens <- object$ensemble
     ensss <- object$ensembless
     nu <- object$control$nu
     mstop <- nrow(ens)
     for (m in 1:mstop) {
         xselect <- ens[m,"xselect"]
         lp[,xselect] <- lp[,xselect] + nu * pr(ensss[[m]])
     }
     colnames(lp) <- colnames(x)
     lp
}

bhatmat <- function(n, H, xselect, fitm, fW) {
    
    B <- matrix(0, nrow = n, ncol = n)
    I <- diag(n)
    tr <- numeric(length(xselect))

    for (m in 1:length(xselect)) {
        B <- B + (H[[xselect[m]]] * fW(fitm[,m])) %*% (I - B)
        tr[m] <- sum(diag(B))
    }
    list(hatmatrix = B, trace = tr)
}

hatglm <- function(model) {

     x <- model$data$x
     if (model$control$center) x <- model$data$center(x)
     nr <- colSums(x^2)
     nu <- model$control$nu
     x <- t(t(x) / sqrt(nr)) * sqrt(nu)
     xselect <- model$ensemble[,"xselect"]
     fW <- model$family@fW

     g <- mboost:::gm.glmboost(model)
     fitm <- t(apply(g, 1, function(a) cumsum(a)))

     B <- matrix(0, nrow = nrow(x), ncol = nrow(x))
     tr <- numeric(length(xselect))

     for (m in 1:length(xselect)) {
        xs <- x[,xselect[m]]
        xw <- xs * fW(fitm[,m])
        B <- B + xs %*% (xw - crossprod(xw, B))
        tr[m] <- .Call("sumdiag", B)
    }
    op <- list(hatmatrix = B, trace = tr)
    RET <- diag(op[[1]])       
    attr(RET, "hatmatrix") <- op[[1]]
    attr(RET, "trace") <- op[[2]]
    RET
}


### fractional polynomials transformation
### all powers `p' of `x', all powers `p' of `x' times `log(x)' and `log(x)'
### see Sauerbrei & Royston (1999), JRSS A (162), 71--94
FP <- function(x, p = c(-2, -1, -0.5, 0.5, 1, 2, 3)) {
    xname <- deparse(substitute(x))
    ### <FIXME> do we need this? 
    ### map x into [1, 2]
    ##if (scale) {
    ##    x <- (x - mean(x))
    ##    x <- x / (max(abs(x)) * 2)
    ##    x <- x - min(x) + 1
    ##}
    ### </FIXME>
    if (any(x < sqrt(.Machine$double.eps)))
        stop("negative values in ", sQuote(xname), "\n")
    Xp <- sapply(p, function(p) x^p)
    Xp <- cbind(Xp, log(x))
    Xp <- cbind(Xp, Xp * log(x))   
    colnames(Xp) <- c(paste(xname, "^", p, sep = ""),
                      paste("log(", xname, ")", sep = ""),
                      paste("log(", xname, ")", xname, "^", p, sep = ""),
                      paste("log(", xname, ")^2", sep = ""))
    Xp
}

### inverse probability of censoring weights
### see van der Laan & Robins (2003)
IPCweights <- function(x, maxweight = 5) {

    if (!extends(class(x), "Surv"))
        stop(sQuote("x"), " is not a Surv object")

    event <- x[,"status"]
    x[,"status"] <- 1 - event
    km <- survfit(x)   
    Ghat <- getsurv(km, time = x[,"time"])
    Ghat[event == 0] <- 1
    w <- event / Ghat
    w[w > maxweight] <- maxweight
    w
}

### extract survival probabilities
### taken from ipred:::getsurv
### DO NOT TOUCH HERE
getsurv <- function(obj, times)
{
    # get the survival probability for times from KM curve j'

    if (!inherits(obj, "survfit")) stop("obj is not of class survfit")
    # <FIXME: methods may have problems with that>
    class(obj) <- NULL
    # </FIXME>
    lt <- length(times)
    nsurv <- times

    # if the times are the same, return the km-curve

    if(length(times) == length(obj$time)) {
        if (all(times == obj$time)) return(obj$surv)
    }

    # otherwise get the km-value for every element of times separatly

    inside <- times %in% obj$time
    for (i in (1:lt)) {
        if (inside[i])
            nsurv[i] <- obj$surv[obj$time == times[i]]
        else  {
            less <- obj$time[obj$time < times[i]]
            if (length(less) == 0)
                nsurv[i] <- 1
            else
                nsurv[i] <- obj$surv[obj$time == max(less)]
        }
    }
    nsurv
}


### modified and speeded up version of `stats::smooth.spline'
### everything is for internal use in `mboost' only!

sknotl <- function(x, nk = NULL) {
    ## if (!all.knots)
    ## return reasonable sized knot sequence for INcreasing x[]:
    n.kn <- function(n) {
        ## Number of inner knots
        if(n < 50) n
        else trunc({
            a1 <- log( 50, 2)
	    a2 <- log(100, 2)
            a3 <- log(140, 2)
	    a4 <- log(200, 2)
	    if	(n < 200) 2^(a1+(a2-a1)*(n-50)/150)
	    else if (n < 800) 2^(a2+(a3-a2)*(n-200)/600)
	    else if (n < 3200)2^(a3+(a4-a3)*(n-800)/2400)
	    else  200 + (n-3200)^0.2
        })
    }
    n <- length(x)
    if(is.null(nk)) nk <- n.kn(n)
    else if(!is.numeric(nk)) stop("'nknots' must be numeric <= n")
    else if(nk > n)
        stop("cannot use more inner knots than unique 'x' values")
    c(rep(x[1], 3), x[seq(1,n, len= nk)], rep(x[n], 3))
}

smoothbase <- function(x, ux, y, w, df) {

    ### handle dummy codings and linear fits, essentially
    if (length(ux) < 4 || df == 1) {
        if (all(x %in% c(0, 1))) {
            X <- matrix(x, ncol = 1)
        } else {
            X <- cbind(1, x)
        }
        object <- lm.wfit(x = X, 
                          y = matrix(y, ncol = 1), w = w)
        object$yfit <- object$fitted.values
        class(object) <- "lmfit"
        return(object)
    }

    ### original `smooth.spline' parameters
    spar <- NULL
    cv <- FALSE
    all.knots <- FALSE
    df.offset <- 0
    penalty <- 1
    nknots <- NULL
    control.spar = list()

    contr.sp <- list(low = -1.5,## low = 0.      was default till R 1.3.x
                     high = 1.5,
                     tol = 1e-4,## tol = 0.001   was default till R 1.3.x
                     eps = 2e-8,## eps = 0.00244 was default till R 1.3.x
                     maxit = 500, trace = getOption("verbose"))
    contr.sp[(names(control.spar))] <- control.spar
    if(!all(sapply(contr.sp[1:4],is.double)) ||
       contr.sp$tol < 0 || contr.sp$eps <= 0 || contr.sp$maxit <= 0)
        stop("invalid 'control.spar'")

    n <- length(x)

    ## Replace y[] for same x[] (to 6 digits precision) by their mean :
    ## BUT ONLY WHEN NEEDED !!!
    if (length(ux) == length(x)) {
        tmp <- cbind(w, w*y, w*y^2)[order(x),]
    } else {
        ox <- match(x, ux)
        tmp <- cbind(w, w*y, w*y^2)
        out <- matrix(0, nrow = length(ux), ncol = 3)
        tmp <- .Call("wybar", ox, sort(unique(ox)), tmp, out, 
                     PACKAGE = "mboost")
    }
    wbar <- tmp[, 1]
    ybar <- tmp[, 2]/ifelse(wbar > 0, wbar, 1)
    yssw <- sum(tmp[, 3] - wbar*ybar^2) # will be added to RSS for GCV
    nx <- length(ux)
    if(nx <= 3) stop("need at least four unique 'x' values")
    if(cv && nx < n)
        warning("crossvalidation with non-unique 'x' values seems doubtful")
    r.ux <- ux[nx] - ux[1]
    xbar <- (ux - ux[1])/r.ux           # scaled to [0,1]
    if(all.knots) {
	knot <- c(rep(xbar[1], 3), xbar, rep(xbar[nx], 3))
	nk <- nx + 2
    } else {
	knot <- sknotl(xbar, nknots)
	nk <- length(knot) - 4
    }

    ## ispar != 1 : compute spar (later)
    ispar <-
        if(is.null(spar) || missing(spar)) { ## || spar == 0
            if(contr.sp$trace) -1 else 0
        } else 1
    spar <- if(ispar == 1) as.double(spar) else double(1)
    ## was <- if(missing(spar)) 0 else if(spar < 1.01e-15) 0 else  1
    ## icrit {../src/sslvrg.f}:
    ##		(0 = no crit,  1 = GCV ,  2 = ord.CV , 3 = df-matching)
    icrit <- if(cv) 2 else  1
    dofoff <- df.offset
    if(!missing(df)) {
	if(df > 1 && df <= nx) {
	    icrit <- 3
	    dofoff <- df
	} else warning("you must supply 1 < df <= n,  n = #{unique x} = ", nx)
    }
    iparms <- as.integer(c(icrit,ispar, contr.sp$maxit))
    names(iparms) <- c("icrit", "ispar", "iter")

    fit <- .Fortran("qsbart",		# code in ../src/qsbart.f
		    as.double(penalty),
		    as.double(dofoff),
		    x = as.double(xbar),
		    y = as.double(ybar),
		    w = as.double(wbar),
		    ssw = as.double(yssw),
		    as.integer(nx),
		    as.double(knot),
		    as.integer(nk),
		    coef = double(nk),
		    ty = double(nx),
		    lev = double(nx),
		    crit = double(1),
		    iparms = iparms,
		    spar = spar,
		    parms = unlist(contr.sp[1:4]),
		    isetup = as.integer(0),
		    scrtch = double((17 + 0) * nk + 1),
		    ld4  = as.integer(4),
		    ldnk = as.integer(1),
		    ier = integer(1),
		    DUP = FALSE, PACKAGE = "stats"
		    )[c("coef", "ty")]

    fit.object <- list(knot = knot, nk = nk, min = ux[1], range = r.ux,
		       coef = fit$coef)
    class(fit.object) <- "smooth.spline.fit"
    if (length(ux) == length(x)) {
        fit.object$yfit <- fit$ty[rank(x)]
    } else {
        fit.object$yfit <- predict(fit.object, x = x)$y
    }
    fit.object
}


predict.lmfit <- function(object, newdata) {
    x <- newdata
    if (length(object$coef) == 2)
        return(as.vector(cbind(1, x) %*% object$coef))
    return(as.vector(x * object$coef))
}
