
mboost_fit <- function(blg, response, weights = NULL, offset = NULL, 
                       family = GaussReg(), control = boost_control()) {

    ### hyper parameters
    mstop <- 0
    risk <- control$risk
    stopifnot(!control$constraint)
    nu <- control$nu
    trace <- control$trace
    tracestep <- options("width")$width / 2

    ### extract negative gradient and risk functions
    ngradient <- family@ngradient
    riskfct <- family@risk

    ### handle missing responses (via zero weights)
    yna <- is.na(response)
    y <- response
    if (any(yna)) {
        weights[] <- 0
        y[is.na(y)] <- y[1]
    }
    check_y_family(response, family)

    ### recode to -1, +1
    if (is.factor(response)) {
        y <- as.numeric(y) - 1
        y[y == 0] <- -1
    }

    ### unweighted problem
    if (is.null(weights)) weights <- rep.int(1, length(y))
    WONE <- (max(abs(weights - 1)) < .Machine$double.eps)
    if (!family@weights && !WONE)
        stop(sQuote("family"), " is not able to deal with weights")

    ### rescale weights (because of the AIC criterion)
    ### <FIXME> is this correct with zero weights??? </FIXME>
    weights <- rescale_weights(weights)
    oobweights <- as.numeric(weights == 0)

    ### set up the fitting functions
    bl <- lapply(blg, dpp, weights = weights)
    blfit <- lapply(bl, function(x) x$fit)
    fit1 <- blfit[[1]]

    xselect <- NA
    ens <- vector(mode = "list", length = control$mstop)

    ### vector of empirical risks for all boosting iterations
    ### (either in-bag or out-of-bag)
    mrisk <- NA
    tsums <- numeric(length(bl))
    ss <- vector(mode = "list", length = length(bl))

    ### initialized the boosting algorithm
    fit <- offset
    if (is.null(offset))
        fit <- offset <- family@offset(y, weights)
    u <- ustart <- ngradient(y, fit, weights)

    ### set up a function for boosting
    boost <- function(niter) {
        for (m in (mstop + 1):(mstop + niter)) {

            ### is there is more than one baselearner
            if (length(bl) > 1) {
                tsums <- rep(-1, length(bl))
                ### fit least squares to residuals _componentwise_
                for (i in 1:length(bl)) {
                    ss[[i]] <- sstmp <- blfit[[i]](y = u) ### try(fit(bl[[i]], y = u))
                    ### if (inherits(ss[[i]], "try-error")) next
                    tsums[i] <- mean.default(weights * ((sstmp$fitted()) - u)^2, 
                                             na.rm = TRUE)
                }
                if (all(tsums < 0))
                    stop("could not fit base learner in boosting iteration ", m)
                xselect[m] <<- which.min(tsums)
                basess <- ss[[xselect[m]]]	
            } else {
                ### a little faster for one baselearner only
                basess <- fit1(y = u)
                xselect[m] <<- 1
            }

            ### update step
            fit <<- fit + nu * basess$fitted()

            ### negative gradient vector, the new `residuals'
            u <<- ngradient(y, fit, weights)

            ### evaluate risk, either for the learning sample (inbag)
            ### or the test sample (oobag)
            if (risk == "inbag") mrisk[m] <<- riskfct(y, fit, weights)
            if (risk == "oobag") mrisk[m] <<- riskfct(y, fit, oobweights)

            ### save the model
            ens[[m]] <<- basess

            ### print status information
            if (trace)
                mboost:::do_trace(m, risk = mrisk, step = tracestep, width = niter)
        }
        mstop <<- mstop + niter 
        return(TRUE)
    }
    ### actually go for initial mstop iterations!
    tmp <- boost(control$mstop)

    ### prepare a (very) rich objects
    RET <- list(baselearner = blg,      ### the baselearners (without weights)
                basemodel = bl,         ### the basemodels (with weights)
                offset = offset,        ### offset
                ustart = ustart,        ### first negative gradients
                control = control,      ### control parameters
                family = family,        ### family object
                response = response,    ### the response variable
                "(weights)" = weights       ### weights used for fitting
    )

    ### update to new weights; just a fresh start
    RET$update <- function(weights = NULL, risk = "oobag") {
        control$mstop <- mstop
        control$risk <- risk
        mboost_fit(blg = blg, response = response, weights = weights, 
                   offset = offset, family = family, control = control)
    }

    ### number of iterations performed so far
    RET$mstop <- function() mstop

    ### which basemodels have been selected so far?
    RET$xselect <- function(cw = FALSE) {
        ### <FIXME> is this the right place?
        if (inherits(bl[[1]], "bl_cwlin") && (length(bl) == 1 && cw))
            return(as.integer(sapply(ens[1:mstop], function(m) m$model[2])))
        ### </FIXME>
        return(xselect[1:mstop])
    }

    ### current fitted values
    RET$fitted <- function() fit

    ### current negative gradient
    RET$resid <- function() u

    ### current risk fct.
    RET$risk <- function() mrisk[1:mstop]

    ### negative risk (at current iteration)
    RET$logLik <- function() -mrisk[mstop]

    ### figure out which baselearners are requested
    thiswhich <- function(which = NULL, usedonly = FALSE) {
        if (is.null(which)) which <- 1:length(bl)
        if (is.character(which)) {
            i <- sapply(which, function(w) {
                wi <- grep(w, names(bl))
                ifelse(length(wi) > 0, wi, NA)
            })
            if (any(is.na(i)))
                warning(paste(which[is.na(i)], collapse = ","), " not found")
            which <- i
        } 
        ### return only those selected so far
        if (usedonly) which <- which[which %in% RET$xselect()]
        return(which)
    }

    ### prepare for computing predictions in the following ways
    ### - for all baselearners (which = NULL) or selected onces
    ### - for each baselearners separately (components = TRUE)
    ### - aggregated ("sum"), the complete path over all
    ###   boosting iterations done so far ("cumsum") or
    ###   not aggregated at all ("none")
    RET$predict <- function(newdata = NULL, which = NULL, components = FALSE,
                            aggregate = c("sum", "cumsum", "none")) {

        indx <- ((1:length(xselect)) <= mstop)
        which <- thiswhich(which, usedonly = FALSE)
        if (length(which) == 0) return(NULL)

        aggregate <- match.arg(aggregate)

        n <- ifelse(!is.null(newdata), nrow(newdata), length(y))
        pfun <- function(w, agg) {
            ix <- xselect == w & indx
            n <- 
            if (!any(ix)) return(rep.int(0, n))
            nu * bl[[w]]$predict(ens[ix],         
                newdata = newdata, aggregate = agg)
        }

        pr <- switch(aggregate, "sum" = {
            pr <- sapply(which, pfun, agg = "sum")
            colnames(pr) <- names(bl)[which]
            if (components) return(pr)
            offset + rowSums(pr)
        }, "cumsum" = {
            if (components) {
                pr <- lapply(which, pfun, agg = "cumsum")
                names(pr) <- names(bl)[which]
                return(pr)
            } else {
                pr <- lapply(which, pfun, agg = "none")
                ret <- matrix(0, nrow = n, ncol = sum(indx))
                tmp <- integer(length(bl)) + 1
                for (i in 1:sum(indx)) { 
                    ret[,i] <- pr[[xselect[i]]][,tmp[xselect[i]]] 
                    if (i > 1) ret[,i] <- ret[,i] + ret[,i-1]
                    tmp[xselect[i]] <- tmp[xselect[i]] + 1
                }
                return(ret)
            }
         }, "none" = {
                pr <- lapply(which, pfun, agg = "none")
                ret <- matrix(0, nrow = n, ncol = sum(indx))
                tmp <- integer(length(bl)) + 1
                for (i in 1:sum(indx)) { 
                    ret[,i] <- pr[[xselect[i]]][,tmp[xselect[i]]]
                    tmp[xselect[i]] <- tmp[xselect[i]] + 1   
                }
                return(ret)
         })
        return(pr)
    }

    ### extract a list of the model frames of the single baselearners
    RET$model.frame <- function(which = NULL) {
        which <- thiswhich(which, usedonly = FALSE)
        # if (length(which) == 1) return(model.frame(blg[[which]]))
        tmp <- lapply(blg[which], model.frame)
        ret <- vector(mode = "list", length = length(bl))
        names(ret) <- names(bl)
        ret[which] <- tmp
        ret
    }

    ### update to a new number of boosting iterations mstop
    ### i <= mstop means less iterations than current
    ### i >  mstop needs additional computations
    ### updates take place in THIS ENVIRONMENT,
    ### some models are CHANGED!
    RET$subset <- function(i) {
        if (i <= mstop || length(xselect) > i) {
            mstop <<- i
            fit <<- RET$predict()
            u <<- ngradient(y, fit, weights)
        } else {
            tmp <- boost(i - mstop)
        }
    } 

    ### if baselearners have a notion of coefficients,
    ### extract these either aggregated ("sum"),
    ### their coefficient path ("cumsum") or not
    ### aggregated at all ("none")
    RET$coef <- function(which = NULL, aggregate = c("sum", "cumsum", "none")) {
        
        indx <- ((1:length(xselect)) <= mstop)
        which <- thiswhich(which, usedonly = FALSE)
        if (length(which) == 0) return(NULL)

        aggregate <- match.arg(aggregate)
        cfun <- function(w) {
            ix <- (xselect == w & indx)
            if (!any(ix)) return(NULL)
            cf <- sapply(ens[ix], coef)
            ret <- switch(aggregate, 
                "sum" = rowSums(cf) * nu,
                "cumsum" = {
                    M <- triu(crossprod(Matrix(1, nc = ncol(cf))))
                    as.matrix(nu * (cf %*% M))
                },
                "none" = nu * cf
            )
        }
        if (length(which) == 1) {
            ret <- cfun(which)
            attr(ret, "offset") <- offset
            names(ret) <- bl[[which]]$Xnames
            return(ret)
        }
        tmp <- lapply(which, cfun)
        ret <- vector(mode = "list", length = length(bl))
        names(ret) <- names(bl)
        ret[which] <- tmp
        attr(ret, "offset") <- offset
        return(ret)
    }

    ### function for computing hat matrices of individual predictors
    RET$hatvalues <- function(which = NULL) {
        which <- thiswhich(which, usedonly = TRUE)
        tmp <- lapply(bl[which], function(b) hatvalues(b) * nu)
        ret <- vector(mode = "list", length = length(bl))
        names(ret) <- names(bl)
        ret[which] <- tmp
        ret
    }

    class(RET) <- "mboost"
    return(RET)
}

### compute predictions
### <FIXME>: add arguments (for documentation)
###          add link argument (family needs to be touched)
### </FIXME>
predict.mboost <- function(object, type = c("lp", "response"), ...) {
    pr <- object$predict(...)
    type <- match.arg(type)
    if (is.factor(y <- object$response) && type == "response")
        return(factor(levels(y)[(lp > 0) + 1], levels = levels(y)))
    return(pr)
}

### extract coefficients
coef.mboost <- function(object, ...)
    object$coef(...)

### compute boosting hat matrix and its trace
hatvalues.mboost <- function(model, ...) {
    H <- model$hatvalues(...)
    n <- length(model$response)
    if (checkL2(model)) {
        op <- .Call("R_trace_gamboost", as.integer(n), H,
                    as.integer(model$xselect()), PACKAGE = "mboost")
    } else {
        fitm <- predict(model, aggregate = "cumsum")
        op <- bhatmat(n, H, model$xselect(), fitm, model$family@fW)
    }
    RET <- diag(op[[1]])  
    attr(RET, "hatmatrix") <- op[[1]]
    attr(RET, "trace") <- op[[2]]
    RET
}

AIC.mboost <- function(object, method = c("corrected", "classical", "gMDL"),
                       df = c("trace", "actset"), ..., k = 2) {

    df <- match.arg(df)
    if (df == "trace") {
        hatval <- hatvalues(object)
        RET <- AICboost3(object, method = method,
                         df = attr(hatval, "trace"), k = k)
    } 
    if (df == "actset") {
        ### compute active set: number of non-zero coefficients
        ### for each boosting iteration
        xs <- object$xselect(cw = TRUE)
        xu <- sort(sapply(unique(xs), function(i) which(xs == i)[1]))
        xu <- c(xu, mstop(object) + 1)
        df <- rep(1:(length(xu) - 1), diff(xu))
        ### <FIXME>: offset = 0 may mean hat(offset) = 0 or
        ### no offset computed at all!
        if (object$offset != 0) df <- df + 1
        ### </FIXME>
        RET <- AICboost3(object, method = method, 
                         df = df, k = k)
    }
    return(RET)
}

AICboost3 <- function(object, method = c("corrected", "classical", "gMDL"), df, k = 2) {

    if (object$control$risk != "inbag")
        return(NA)
    method <- match.arg(method)

    if (checkL2(object) && method == "classical")
        stop("classical AIC method not implemented for Gaussian family")
    if (!checkL2(object) && method == "corrected")
        stop("corrected AIC method not implemented for non-Gaussian family")

    sumw <- sum(object$weights)
    if (method == "corrected") 
        AIC <- log(object$risk() / sumw) +
               (1 + df/sumw) / (1 - (df + 2)/sumw)

    ### loss-function is to be MINIMIZED, take -2 * logLik == 2 * risk
    if (method == "classical")
        AIC <- 2 * object$risk() + k * df
    if (method == "gMDL"){
        s <- object$risk()/(sumw - df)
        AIC <- log(s) + df/sumw * log((sum(object$response^2) - object$risk())
                      /(df * s))
        }
    mstop <- which.min(AIC)
    RET <- AIC[mstop]

    attr(RET, "mstop") <- which.min(AIC)
    attr(RET, "df") <- df  
    attr(RET, "AIC") <- AIC
    attr(RET, "corrected") <- method == "corrected"

    class(RET) <- "gbAIC"
    return(RET)
}


### compute fitted values
fitted.mboost <- function(object, ...) {
    args <- list(...)
    if (length(args) == 0)
        return(object$fitted())
    object$predict(...)
}

### residuals (the current negative gradient)
resid.mboost <- function(object, ...) 
    object$resid()

logLik.mboost <- function(object, ...)
    object$logLik()

### restrict or enhance models to less/more
### boosting iterations.
### ATTENTION: x gets CHANGED!
"[.mboost" <- function(x, i, ...) {
    stopifnot(length(i) == 1 && i > 0)
    x$subset(i)
    return(NULL)
}

mstop.mboost <- function(object, ...) object$mstop()

update.mboost <- function(object, weights, ...)
    object$update(weights)

model.frame.mboost <- function(formula, ...)
    formula$model.frame(...)

response.mboost <- function(object, ...)
    object$response

### main user interface function
### formula may contain unevaluated baselearners as well
### as additional variables
###     y ~ bols3(x1) + x2 + btree(x3)
### is evaluated as
###     y ~ bols3(x1) + baselearner(x2) + btree(x3)
### see mboost_fit for the dots
mboost <- function(formula, data = list(), baselearner = bbs3, ...) {

    ### OK, we need at least variable names to go ahead
    if (as.name(formula[[3]]) == ".") {
        formula <- as.formula(paste(as.character(formula[[2]]),
            "~", paste(names(data)[names(data) != all.vars(formula[[2]])], 
                       collapse = "+"), collapse = ""))
    }
    ### instead of evaluating a model.frame, we evaluate
    ### the expressions on the lhs of formula directly
    "+" <- function(a,b) {
        ### got baselearner, fine!
        if (inherits(a, "blg")) a <- list(a)
        ### a single variable; compute baselearner
        if (!is.list(a)) a <- list(baselearner(a))
        ### got baselearner, fine!
        if (inherits(b, "blg")) b <- list(b)
        ### a single variable, compute baselearner
        if (!is.list(b)) b <- list(baselearner(b))
        ### join both baselearners in a list
        c(a, b)
    }
    ### set up all baselearners
    bl <- eval(as.expression(formula[[3]]), envir = data)
    ### if there is just one, assign it to a list anyway
    if (inherits(bl, "blg")) bl <- list(bl)
    ### just a check
    stopifnot(all(sapply(bl, inherits, what = "blg")))
    ### we need identifiers for the baselearners, 
    ### split the formula at `+'
    nm <- strsplit(as.character(as.expression(formula[[3]])), "\\+")[[1]]
    nm <- gsub(" ", "", nm)
    names(bl) <- nm
    ### baselearners constructed indirectly via `baselearner'
    ### don't know the variable names yet
    ### check for parts of the lhs of formula
    ### that didn't specify a baselearner
    funs <- grep("\\(", nm)
    missing <- 1:length(nm)
    if (length(funs) > 0)
        missing <- missing[-funs]
    if (length(missing) > 0) {
        ### assign variable names
        for (m in missing) bl[[m]]$set_names(nm[m])
        ### assign names containing `baselearner'
        names(bl)[missing] <- paste(deparse(substitute(baselearner)), 
                                    "(", nm[missing], ")", sep = "")
    }
    ### get the response
    response <- eval(as.expression(formula[[2]]), envir = data)
    mboost_fit(bl, response = response, ...)
}

### nothing to do there
Gamboost <- mboost

### just one single tree-based baselearner
Blackboost <- function(formula, data = list(), ...) {

    if (formula[[3]] == as.name(".")) {
        xvars <- names(data)[names(data) != all.vars(formula)]
    } else {
        xvars <- all.vars(formula[[3]])
    }
    formula <- as.formula(paste(formula[[2]], "~ btree3(", 
        paste(xvars, collapse = ","), ")", collapse = ""))
    mboost(formula = formula, data = data, ...)
}

### fit a linear model componentwise
Glmboost <- function(formula, data = list(), weights = NULL, 
                     na.action = na.pass, contrasts.arg = NULL, 
                     center = FALSE, control = boost_control(), ...) {

    ### get the model frame first
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "weights"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$na.action <- na.action
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    ### center argument moved to this function
    if (control$center) {
        center <- TRUE
        warning("boost_control center deprecated")
    }
    ### set up the model.matrix and center (if requested)
    X <- model.matrix(attr(mf, "terms"), data = mf, 
                      contrasts.args = contrasts.args)
    cm <- rep(0, ncol(X))
    if (center) {
        cm <- colMeans(X, na.rm = TRUE)
        ### center numeric variables only
        center <- attr(X, "assign") %in% which(sapply(mf, is.numeric)[-1])
        cm[!center] <- 0
        X <- scale(X, center = cm, scale = FALSE)
    }
    ### this function will be used for predictions later
    newX <- function(newdata) {
        X <- model.matrix(delete.response(attr(mf, "terms")), data = newdata,
                          contrasts.args = contrasts.args)
        scale(X, center = cm, scale = FALSE)
    }

    ### component-wise linear models baselearner
    bl <- list(bolscw(X))
    response <- model.response(mf)
    weights <- model.weights(mf)
    ret <- mboost_fit(bl, response = response, weights = weights, 
                      control = control, ...)
    ret$newX <- newX
    ### need specialized method (hatvalues etc. anyway)
    class(ret) <- c("Glmboost", "mboost")
    return(ret)
}

Glmboost.matrix <- function(x, y, center = FALSE, 
                            control = boost_control(), ...) {

    X <- x
    if (control$center) {
        center <- TRUE
        warning("boost_control center deprecated")
    }
    if (center) {
        cm <- colMeans(X, na.rm = TRUE)
        ### center numeric variables only
        center <- attr(X, "assign") %in% which(sapply(mf, is.numeric)[-1])
        cm[!center] <- 0
        X <- scale(X, center = cm, scale = FALSE)
    }
    newX <- function(newdata) {
        if (isMATRIX(newdata)) {
            if (all(colnames(X) == colnames(newdata)))
                return(newdata)
        }
        return(NULL)
    }
    bl <- list(bolscw(X))
    ret <- mboost_fit(bl, response = y, control = control, ...)
    ret$newX <- newX
    ### need specialized method (hatvalues etc. anyway)
    class(ret) <- c("Glmboost", "mboost")
    return(ret)
}


    

predict.Glmboost <- function(object, newdata = NULL, ...) {

    if (!is.null(newdata)) {
        newdata <- object$newX(newdata)
    }
    object$predict(newdata = newdata, ...)
}

hatvalues.Glmboost <- function(model, ...) {

    if (!checkL2(model)) return(hatvalues.mboost(model))
    Xf <- t(model$basemodel[[1]]$MPinv()) * model$control$nu
    X <- model$baselearner[[1]]$get_data()
    op <- .Call("R_trace_glmboost", X, Xf,
                as.integer(model$xselect(cw = TRUE)),
                PACKAGE = "mboost")
    RET <- diag(op[[1]])
    attr(RET, "hatmatrix") <- op[[1]]  
    attr(RET, "trace") <- op[[2]] 
    RET
}
