
mboost_fit <- function(blg, response, weights = NULL, offset = NULL, 
                       family = GaussReg(), control = boost_control()) {

    ### hyper parameters
    mstop <- 0
    risk <- control$risk
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
    if (is.null(weights)) weights <- rep.int(1, NROW(y))
    WONE <- (max(abs(weights - 1)) < .Machine$double.eps)
    if (!family@weights && !WONE)
        stop(sQuote("family"), " is not able to deal with weights")

    ### rescale weights (because of the AIC criterion)
    ### <FIXME> is this correct with zero weights??? </FIXME>
    weights <- rescale_weights(weights)
    oobweights <- as.numeric(weights == 0)
    if (control$risk == "oobag") {
        triskfct <- function(y, f) riskfct(y, f, oobweights)
    } else {
        triskfct <- function(y, f) riskfct(y, f, weights)
    }

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
    offsetarg <- offset
    if (is.null(offset))
        fit <- offset <- family@offset(y, weights)
    u <- ustart <- ngradient(y, fit, weights)

    ### set up function for fitting baselearner(s)
    if (length(bl) > 1) {
        basefit <- function(u, m) {
           tsums <- rep(-1, length(bl))
           ### fit least squares to residuals _componentwise_
           for (i in 1:length(bl)) {
               ss[[i]] <<- sstmp <- blfit[[i]](y = u) ### try(fit(bl[[i]], y = u))
               ### if (inherits(ss[[i]], "try-error")) next
               tsums[i] <- mean.default(weights * ((sstmp$fitted()) - u)^2,
                                        na.rm = TRUE)
           }
           if (all(tsums < 0))
               stop("could not fit base learner in boosting iteration ", m)
           xselect[m] <<- which.min(tsums)
           return(ss[[xselect[m]]])
        }
    } else {
        basefit <- function(u, m) {
            xselect[m] <<- 1
            return(fit1(y = u))
        }
    }


    ### set up a function for boosting
    boost <- function(niter) {
        for (m in (mstop + 1):(mstop + niter)) {

            ### fit baselearner(s)
            basess <- basefit(u, m)

            ### update step
            ### <FIXME> handle missing values!
            fit <<- fit + nu * basess$fitted()
            ### <FIXME>

            ### negative gradient vector, the new `residuals'
            u <<- ngradient(y, fit, weights)

            ### evaluate risk, either for the learning sample (inbag)
            ### or the test sample (oobag)
            mrisk[m] <<- triskfct(y, fit)

            ### save the model
            ens[[m]] <<- basess

            ### print status information
            ### print xselect???
            if (trace)
                do_trace(m, risk = mrisk, step = tracestep, width = niter)
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
        ### use user specified offset only (since it depends on weights otherwise)
        if (!is.null(offsetarg)) offsetarg <- offset
        mboost_fit(blg = blg, response = response, weights = weights, 
                   offset = offsetarg, family = family, control = control)
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
                wi <- grep(w, names(bl), fixed = TRUE)
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
    ### - for all selected baselearners (which = NULL) or chosen ones (selected or not)
    ### - aggregated ("sum"), the complete path over all
    ###   boosting iterations done so far ("cumsum") or
    ###   not aggregated at all ("none")
    ### - always returns a matrix
    RET$predict <- function(newdata = NULL, which = NULL, 
                            aggregate = c("sum", "cumsum", "none")) {

        indx <- ((1:length(xselect)) <= mstop)
        which <- thiswhich(which, usedonly = nw <- is.null(which))
        if (length(which) == 0) return(NULL)

        aggregate <- match.arg(aggregate)

        n <- ifelse(!is.null(newdata), nrow(newdata), length(y))
        pfun <- function(w, agg) {
            ix <- xselect == w & indx
            m <- Matrix(0, nrow = n, ncol = sum(indx))
            if (!any(ix)) {
                if (agg == "sum") return(rep.int(0, n))
                return(m)
            }
            ret <- nu * bl[[w]]$predict(ens[ix],         
                   newdata = newdata, aggregate = agg)
            if (agg == "sum") return(ret)
            m[, which(ix)] <- ret
            m
        }

        pr <- switch(aggregate, "sum" = {
            pr <- sapply(which, pfun, agg = "sum")
            if (!nw) return(pr)
            ### only if no selection of baselearners
            ### was made via the `which' argument
            offset + matrix(rowSums(pr), ncol = 1)
        }, "cumsum" = {
            if (!nw) {
                pr <- lapply(which, pfun, agg = "none")
                M <- triu(crossprod(Matrix(1, nc = sum(indx))))
                pr <- lapply(pr, function(x) as(x %*% M, "matrix"))
                names(pr) <- names(bl)[which]
                return(pr)
            } else {
                ret <- Matrix(0, nrow = n, ncol = sum(indx))
                for (i in 1:length(bl)) ret <- ret + pfun(i, agg = "none")
                M <- triu(crossprod(Matrix(1, nc = sum(indx))))
                as(ret %*% M, "matrix") + offset
            }
         }, "none" = {
            if (!nw) {
                pr <- lapply(which, pfun, agg = "none")
                names(pr) <- names(bl)[which]
                return(pr)
            } else {
                ret <- Matrix(0, nrow = n, ncol = sum(indx))
                for (i in 1:length(bl)) ret <- ret + pfun(i, agg = "none")
                as(ret, "matrix")
            }
         })
        return(pr)
    }

    ### extract a list of the model frames of the single baselearners
    RET$model.frame <- function(which = NULL) {
        which <- thiswhich(which, usedonly = is.null(which))
        ret <- lapply(blg[which], model.frame)
        names(ret) <- names(bl)[which]
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
        which <- thiswhich(which, usedonly = is.null(which))
        if (length(which) == 0) return(NULL)

        aggregate <- match.arg(aggregate)
        cfun <- function(w) {
            ix <- (xselect == w & indx)
            if (!any(ix)) return(NULL)
            cf <- sapply(ens[ix], coef)
            if (!is.matrix(cf)) cf <- matrix(cf, nrow = 1)
            ret <- switch(aggregate, 
                "sum" = rowSums(cf) * nu,
                "cumsum" = {
                    M <- triu(crossprod(Matrix(1, nc = ncol(cf))))
                    as.matrix(nu * (cf %*% M))
                },
                "none" = nu * cf
            )
        }
        ret <- lapply(which, cfun)
        names(ret) <- names(bl)[which]
        attr(ret, "offset") <- offset
        return(ret)
    }

    ### function for computing hat matrices of individual predictors
    RET$hatvalues <- function(which = NULL) {
        which <- thiswhich(which, usedonly = is.null(which))
        ### make sure each list element corresponds to one baselearner
        ### non-selected baselearners receive NULL
        ret <- vector(mode = "list", length = length(bl))
        ret[which] <- lapply(bl[which], function(b) hatvalues(b) * nu)
        names(ret) <- names(bl)
        ret
    }

    class(RET) <- "mboost"
    return(RET)
}

### main user interface function
### formula may contain unevaluated baselearners as well
### as additional variables
###     y ~ bols3(x1) + x2 + btree(x3)
### is evaluated as
###     y ~ bols3(x1) + baselearner(x2) + btree(x3)
### see mboost_fit for the dots
mboost <- function(formula, data = list(), baselearner = c("bbs", "bols", "btree", "bss", "bns"), ...) {

    if (is.character(baselearner)) {
        baselearner <- match.arg(baselearner)
        if (baselearner %in% c("bss", "bns")) {
            warning("bss and bns are deprecated, bbs is used instead")
            baselearner <- "bbs"
        }
        baselearner <- get(baselearner, mode = "function", 
                           envir = parent.frame())
    }
    stopifnot(is.function(baselearner))

    ### OK, we need at least variable names to go ahead
    if (length(formula[[3]]) == 1) {
        if (as.name(formula[[3]]) == ".") {
            formula <- as.formula(paste(deparse(formula[[2]]),
                "~", paste(names(data)[names(data) != all.vars(formula[[2]])], 
                           collapse = "+"), collapse = ""))
        }
    }
    ### instead of evaluating a model.frame, we evaluate
    ### the expressions on the rhs of formula directly
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
    if (!is.list(bl)) bl <- list(baselearner(bl))
    ### just a check
    stopifnot(all(sapply(bl, inherits, what = "blg")))
    ### we need identifiers for the baselearners, 
    ### split the formula at `+'
    nm <- strsplit(paste(deparse(formula[[3]]), collapse = ""), "\\+")[[1]]
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
    ret <- mboost_fit(bl, response = response, ...)
    ret$call <- match.call()
    ret
}

### nothing to do there
gamboost <- function(formula, data = list(), baselearner = c("bbs", "bols", "btree", "bss", "bns"), 
                     dfbase = 4, ...) {

    if (is.character(baselearner)) {
        baselearner <- match.arg(baselearner)
        if (baselearner %in% c("bss", "bns")) {
            warning("bss and bns are deprecated, bbs is used instead")
            baselearner <- "bbs"
        }
        baselearner <- get(baselearner, mode = "function",
                           envir = parent.frame())
    }
    stopifnot(is.function(baselearner))
    if (isTRUE(all.equal(baselearner, bbs)))
        baselearner <- function(...) bbs(..., df = dfbase)
    ret <- mboost(formula = formula, data = data, baselearner = baselearner, ...)
    ret$call <- match.call()
    ret
}


### just one single tree-based baselearner
blackboost <- function(formula, data = list(), 
    tree_controls = ctree_control(teststat = "max",
                               testtype = "Teststatistic",
                               mincriterion = 0,
                               maxdepth = 2),
    ...) {

    ### get the model frame first
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "weights"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$na.action <- na.pass
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    response <- model.response(mf)
    mf <- mf[,-1, drop = FALSE]
    bl <- list(btree(mf, tree_controls = tree_controls))
    ret <- mboost_fit(bl, response = response, ...)
    ret$call <- cl
    ret
}

### fit a linear model componentwise
glmboost <- function(x, ...) UseMethod("glmboost", x)

glmboost.formula <- function(formula, data = list(), weights = NULL, 
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
        warning("boost_control(center = TRUE) is deprecated, use glmboost(..., center = TRUE)")
    }
    ### set up the model.matrix and center (if requested)
    X <- model.matrix(attr(mf, "terms"), data = mf, 
                      contrasts.arg = contrasts.arg)
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
        mf <- model.frame(attr(mf, "terms"), data = newdata, na.action = na.pass)
        X <- model.matrix(attr(mf, "terms"), data = mf,
                          contrasts.arg = contrasts.arg)
        scale(X, center = cm, scale = FALSE)
    }

    ### component-wise linear models baselearner
    bl <- list(bolscw(X))
    response <- model.response(mf)
    weights <- model.weights(mf)
    ret <- mboost_fit(bl, response = response, weights = weights, 
                      control = control, ...)
    ret$newX <- newX
    ret$call <- cl
    ### need specialized method (hatvalues etc. anyway)
    class(ret) <- c("glmboost", "mboost")
    return(ret)
}

glmboost.matrix <- function(x, y, center = FALSE, 
                            control = boost_control(), ...) {

    X <- x
    if (control$center) {
        center <- TRUE
        warning("boost_control center deprecated")
    }
    if (center) {
        cm <- colMeans(X, na.rm = TRUE)
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
    ret$call <- match.call()
    ### need specialized method (hatvalues etc. anyway)
    class(ret) <- c("glmboost", "mboost")
    return(ret)
}

glmboost.default <- function(x, ...) {
    if (extends(class(x), "Matrix"))
        return(glmboost.matrix(x = x, ...))
    stop("no method for objects of class ", sQuote(class(x)), 
         " implemented")
}
