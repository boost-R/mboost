
mboost_fit <- function(blg, response, weights = rep(1, NROW(response)),
                       offset = NULL, family = Gaussian(), control =
                       boost_control(), oobweights = as.numeric(weights == 0)) {

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
        weights[yna] <- 0
        y[is.na(y)] <- y[1]
    }
    y <- check_y_family(response, family)

    ### unweighted problem
    if (is.null(weights)) weights <- rep.int(1, NROW(y))
    if (!family@weights(weights))
        stop(sQuote("family"), " is not able to deal with weights")

    ### rescale weights (because of the AIC criterion)
    ### <FIXME> is this correct with zero weights??? </FIXME>
    weights <- rescale_weights(weights)
    if (is.null(oobweights))
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
    nuisance <- vector(mode = "list", length = control$mstop)

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
    cwlin <- FALSE
    if (length(bl) > 1) {
        bnames <- names(bl)
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
        cwlin <- inherits(bl[[1]], "bl_cwlin")
        if (cwlin) {
            bnames <- bl[[1]]$Xnames
            basefit <- function(u, m) {
                mod <- fit1(y = u)
                xselect[m] <<- mod$model["xselect"]
                return(mod)
            }
        } else {
            bnames <- names(bl)
            basefit <- function(u, m) {
                xselect[m] <<- 1L
                return(fit1(y = u))
            }
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
            nuisance[[m]] <<- family@nuisance()

            ### print status information
            ### print xselect???
            if (trace)
                do_trace(m, mstop = mstop, risk = mrisk, 
                         step = tracestep, width = niter)
        }
        mstop <<- mstop + niter
        return(TRUE)
    }
    ### actually go for initial mstop iterations!
    tmp <- boost(control$mstop)

    ### prepare a (very) rich objects
    RET <- list(baselearner = blg,          ### the baselearners (without weights)
                basemodel = bl,             ### the basemodels (with weights)
                offset = offset,            ### offset
                ustart = ustart,            ### first negative gradients
                control = control,          ### control parameters
                family = family,            ### family object
                response = response,        ### the response variable
                rownames = bnames,          ### rownames of learning data
                "(weights)" = weights,      ### weights used for fitting
                nuisance = 
                    function() nuisance     ### list of nuisance parameters
    )

    ### update to new weights; just a fresh start
    RET$update <- function(weights = NULL, oobweights = NULL, risk = "oobag") {
        control$mstop <- mstop
        control$risk <- risk
        ### use user specified offset only (since it depends on weights otherwise)
        if (!is.null(offsetarg)) offsetarg <- offset
        mboost_fit(blg = blg, response = response, weights = weights,
                   offset = offsetarg, family = family, control = control,
                   oobweights = oobweights)
    }

    ### number of iterations performed so far
    RET$mstop <- function() mstop

    ### which basemodels have been selected so far?
    RET$xselect <- function()
        return(xselect[1:mstop])

    ### current fitted values
    RET$fitted <- function() as.vector(fit)

    ### current negative gradient
    RET$resid <- function() u

    ### current risk fct.
    RET$risk <- function() mrisk[1:mstop]

    ### negative risk (at current iteration)
    RET$logLik <- function() -mrisk[mstop]

    ### figure out which baselearners are requested
    thiswhich <- function(which = NULL, usedonly = FALSE) {
        if (is.null(which)) which <- 1:max(RET$xselect())
        if (is.character(which)) {
            i <- sapply(which, function(w) {
                wi <- grep(w, bnames, fixed = TRUE)
                if (length(wi) > 0) return(wi)
                return(NA)
            })
            if (any(is.na(i)))
                warning(paste(which[is.na(i)], collapse = ","), " not found")
            which <- i
        }
        ### return only those selected so far
        if (usedonly) which <- which[which %in% RET$xselect()]
        return(which)
    }
    RET$which <- thiswhich

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
            if (cwlin) w <- 1
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
            if (!is.matrix(pr)) pr <- matrix(pr, nrow = 1)
            return(offset + matrix(rowSums(pr), ncol = 1))
        }, "cumsum" = {
            if (!nw) {
                pr <- lapply(which, pfun, agg = "none")
                pr <- lapply(pr, function(x) .Call("R_mcumsum", as(x, "matrix")))
                names(pr) <- bnames[which]
                return(pr)
            } else {
                ret <- Matrix(0, nrow = n, ncol = sum(indx))
                for (i in 1:max(xselect)) ret <- ret + pfun(i, agg = "none")
                return(.Call("R_mcumsum", as(ret, "matrix")) + offset)
            }
         }, "none" = {
            if (!nw) {
                pr <- lapply(which, pfun, agg = "none")
                for (i in 1:length(pr)) pr[[i]] <- as(pr[[i]], "matrix")
                names(pr) <- bnames[which]
                return(pr)
            } else {
                ret <- Matrix(0, nrow = n, ncol = sum(indx))
                for (i in 1:max(xselect)) ret <- ret + pfun(i, agg = "none")
                return(as(ret, "matrix"))
            }
         })
         return(pr)
    }

    ### extract a list of the model frames of the single baselearners
    RET$model.frame <- function(which = NULL) {
        which <- thiswhich(which, usedonly = is.null(which))
        ret <- lapply(blg[which], model.frame)
        names(ret) <- bnames[which]
        ret
    }

    ### update to a new number of boosting iterations mstop
    ### i <= mstop means less iterations than current
    ### i >  mstop needs additional computations
    ### updates take place in THIS ENVIRONMENT,
    ### some models are CHANGED!
    RET$subset <- function(i) {
        if (i <= mstop || i <= length(xselect)) {
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
            cf <- numeric(mstop)
            if (!any(ix)) {
                if (!cwlin)
                    cf <- matrix(0, nrow = length(bl[[w]]$Xnames), ncol = mstop)
            } else {
                cftmp <- sapply(ens[ix], coef)
                nr <- NROW(cftmp)
                if (!is.matrix(cftmp)) nr <- 1
                cf <- matrix(0, nrow = nr, ncol = mstop)
                cf[, which(ix)] <- cftmp
            }
            if (!is.matrix(cf)) cf <- matrix(cf, nrow = 1)
            ret <- switch(aggregate,
                "sum" = rowSums(cf) * nu,
                "cumsum" = {
                    .Call("R_mcumsum", as(cf, "matrix"))
                },
                "none" = nu * cf
            )
        }
        ret <- lapply(which, cfun)
        names(ret) <- bnames[which]
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
        names(ret) <- bnames
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
mboost <- function(formula, data = list(),
    baselearner = c("bbs", "bols", "btree", "bss", "bns"), ...) {

    if (is.character(baselearner)) {
        baselearner <- match.arg(baselearner)
        bname <- baselearner
        if (baselearner %in% c("bss", "bns")) {
            warning("bss and bns are deprecated, bbs is used instead")
            baselearner <- "bbs"
        }
        baselearner <- get(baselearner, mode = "function",
                           envir = parent.frame())
    } else {
        bname <- deparse(substitute(baselearner))
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
        cl <- match.call()
        ### got baselearner, fine!
        if (inherits(a, "blg")) a <- list(a)
        ### a single variable; compute baselearner
        if (!is.list(a)) {
            a <- list(baselearner(a))
            a[[1]]$set_names(deparse(cl[[2]]))
        }
        ### got baselearner, fine!
        if (inherits(b, "blg")) b <- list(b)
        ### a single variable, compute baselearner
        if (!is.list(b)) {
            b <- list(baselearner(b))
            b[[1]]$set_names(deparse(cl[[3]]))
        }
        ### join both baselearners in a list
        c(a, b)
    }
    ### set up all baselearners
    bl <- eval(as.expression(formula[[3]]), 
               envir = c(as.list(data), list("+" = get("+"))),
               enclos = environment(formula))
    ### rhs was one single baselearner
    if (inherits(bl, "blg")) bl <- list(bl)
    ### rhs was one single variable
    if (!is.list(bl)) { 
        bl <- list(baselearner(bl))
        bl[[1]]$set_names(as.character(formula[[3]]))
    }

    ### just a check
    stopifnot(all(sapply(bl, inherits, what = "blg")))

    ### assign calls as names of base learners
    names(bl) <- sapply(bl, function(x) x$get_call())

    ### get the response
    response <- eval(as.expression(formula[[2]]), envir = data,
                     enclos = environment(formula))
    ret <- mboost_fit(bl, response = response, ...)
    if (is.data.frame(data) && nrow(data) == length(response))
        ret$rownames <- rownames(data)
    else
        ret$rownames <- 1:NROW(response)
    ret$call <- match.call()
    ret
}

### nothing to do there
gamboost <- function(formula, data = list(),
    baselearner = c("bbs", "bols", "btree", "bss", "bns"), dfbase = 4, ...) {

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
        baselearner <- function(...) bbs(as.data.frame(list(...)), df = dfbase)
    ret <- mboost(formula = formula, data = data, baselearner = baselearner, ...)
    ret$call <- match.call()
    class(ret) <- c("gamboost", class(ret))
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
    ret$rownames <- rownames(mf)
    class(ret) <- c("blackboost", class(ret))
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
        mf <- model.frame(delete.response(attr(mf, "terms")),
            data = newdata, na.action = na.pass)
        X <- model.matrix(delete.response(attr(mf, "terms")), data = mf,
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
    ret$hatvalues <- function() {
        H <- vector(mode = "list", length = ncol(X))
        MPinv <- ret$basemodel[[1]]$MPinv()
        for (j in unique(ret$xselect()))
            H[[j]] <- (X[,j] %*% MPinv[j, ,drop = FALSE]) * control$nu
        H
    }
    ret$rownames <- rownames(mf)
    class(ret) <- c("glmboost", "mboost")
    return(ret)
}

glmboost.matrix <- function(x, y, center = FALSE,
                            control = boost_control(), ...) {q

    X <- x
    if (is.null(colnames(X))) colnames(X) <- paste("V", 1:ncol(X), sep = "")
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
        stop(sQuote("newdata"), " is not a matrix with the same variables as ",
              sQuote("x"))
        return(NULL)
    }
    bl <- list(bolscw(X))
    ret <- mboost_fit(bl, response = y, control = control, ...)
    ret$newX <- newX
    ret$call <- match.call()
    ### need specialized method (hatvalues etc. anyway)
    ret$hatvalues <- function() {
        H <- vector(mode = "list", length = ncol(X))
        MPinv <- ret$basemodel[[1]]$MPinv()
        for (j in unique(ret$xselect()))
            H[[j]] <- (X[,j] %*% MPinv[j, ,drop = FALSE]) * control$nu
        H
    }
    ret$rownames <- rownames(X)
    class(ret) <- c("glmboost", "mboost")
    return(ret)
}

glmboost.default <- function(x, ...) {
    if (extends(class(x), "Matrix"))
        return(glmboost.matrix(x = x, ...))
    stop("no method for objects of class ", sQuote(class(x)),
         " implemented")
}
