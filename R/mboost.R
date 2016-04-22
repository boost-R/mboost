
mboost_fit <- function(blg, response, weights = rep(1, NROW(response)),
                       offset = NULL, family = Gaussian(),
                       control = boost_control(),
                       oobweights = as.numeric(weights == 0)) {

    ### hyper parameters
    mstop <- 0
    risk <- control$risk
    nu <- control$nu
    trace <- control$trace
    stopintern <- control$stopintern
    if (is.numeric(stopintern)) {
        stopeps <- stopintern
        stopintern <- TRUE
    } else {
        stopeps <- 0
    }
    tracestep <- options("width")$width / 2

    ### extract negative gradient and risk functions
    ngradient <- family@ngradient
    riskfct <- family@risk

    ### handle missing responses (via zero weights)
    yna <- is.na(response)
    y <- response
    if (any(yna)) {
        weights[yna] <- 0
        y[is.na(y)] <- y[!yna][1] # use first non-missing value
        warning("response contains missing values;\n",
                " weights of corresponding observations are set to zero",
                " and thus these observations are not used for fitting")
    }
    y <- check_y_family(y, family)

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
           if (all(is.na(tsums)) || all(tsums < 0))
               stop("could not fit any base-learner in boosting iteration ", m)
           xselect[m] <<- which.min(tsums)
           return(ss[[xselect[m]]])
        }
    } else {
        cwlin <- inherits(bl[[1]], "bl_cwlin")
        if (cwlin) {
            bnames <- bl[[1]]$Xnames
            basefit <- function(u, m) {
                mod <- fit1(y = u)
                if(all(is.na(coef(mod))))
                    stop("could not fit any base-learner in boosting iteration ", m)
                xselect[m] <<- mod$model["xselect"]
                return(mod)
            }
        } else {
            bnames <- names(bl)
            basefit <- function(u, m) {
                xselect[m] <<- 1L
                mod <- fit1(y = u)
                if((!is.null(coef(mod)) && all(is.na(coef(mod)))) ||
                   all(is.na(mean.default(weights * ((mod$fitted()) - u)^2, na.rm = TRUE))))
                    stop("could not fit base-learner in boosting iteration ", m)
                return(mod)
            }
        }
    }

    ## if names are missing try to get these from the calls
    if (is.null(bnames) && !cwlin)
        names(blg) <- names(bl) <- bnames <- sapply(blg, function(x) x$get_call())


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

            ### happened for family = Poisson() with nu = 0.1
            if (any(!is.finite(u[!is.na(u)])))
                stop("Infinite residuals: please decrease step-size nu in ",
                     sQuote("boost_control"))

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

            ### internal stopping (for oobag risk only)
            if (stopintern & (m > 1)) {
                if ((mrisk[m] - mrisk[m - 1]) > stopeps) break
            }
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
                rownames =                  ### rownames of learning data
                    paste("obs", 1:length(weights), sep = ""),
                "(weights)" = weights,      ### weights used for fitting
                nuisance =
                    function() nuisance     ### list of nuisance parameters
    )

    ### update to new weights; just a fresh start
    RET$update <- function(weights = NULL, oobweights = NULL, risk = "oobag",
                           trace = NULL) {

        control$mstop <- mstop
        if (!is.null(risk))
            control$risk <- risk
        if (!is.null(trace))
            control$trace <- trace
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
        if (is.null(which)) which <- 1:length(bnames)
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

        pfun <- function(w, agg) {
            ix <- xselect == w & indx
            if (!any(ix))
                return(0)
            if (cwlin) w <- 1
            ret <- nu * bl[[w]]$predict(ens[ix],
                   newdata = newdata, aggregate = agg)
            if (agg == "sum") return(ret)
            m <- Matrix(0, nrow = nrow(ret), ncol = sum(indx))
            m[, which(ix)] <- ret
            m
        }

        pr <- switch(aggregate, "sum" = {
            pr <- do.call("cbind", lapply(which, pfun, agg = "sum"))
            if (!nw){
                colnames(pr) <- bnames[which]
                attr(pr, "offset") <- offset
                return(pr)
            } else {
                ## only if no selection of baselearners
                ## was made via the `which' argument
                ret <- matrix(rowSums(pr), ncol = 1)
                if (length(offset) != 1 && !is.null(newdata)) {
                    warning(sQuote("length(offset) > 1"),
                            ": User-specified offset is not a scalar, ",
                            "thus offset not used for prediction when ",
                            sQuote("newdata"), " is specified")
                } else {
                    ret <- ret + offset
                }
                return(ret)
            }
        }, "cumsum" = {
            if (!nw) {
                pr <- lapply(which, pfun, agg = "none")
                pr <- lapply(pr, function(x) {
                    .Call("R_mcumsum", as(x, "matrix"), PACKAGE = "mboost")
                })
                names(pr) <- bnames[which]
                attr(pr, "offset") <- offset
                return(pr)
            } else {
                pr <- 0
                for (i in 1:max(xselect)) pr <- pr + pfun(i, agg = "none")
                pr <- .Call("R_mcumsum", as(pr, "matrix"), PACKAGE = "mboost")
                if (length(offset) != 1 && !is.null(newdata)) {
                    warning(sQuote("length(offset) > 1"),
                            ": User-specified offset is not a scalar, ",
                            "thus offset not used for prediction when ",
                            sQuote("newdata"), " is specified")
                } else {
                    pr <- pr + offset
                }
                return(pr)
            }
         }, "none" = {
            if (!nw) {
                pr <- lapply(which, pfun, agg = "none")
                for (i in 1:length(pr)) pr[[i]] <- as(pr[[i]], "matrix")
                names(pr) <- bnames[which]
                attr(pr, "offset") <- offset
                return(pr)
            } else {
                pr <- 0
                for (i in 1:max(xselect)) pr <- pr + pfun(i, agg = "none")
                pr <- as(pr, "matrix")
                attr(pr, "offset") <- offset
                return(pr)
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
            ## no need to recompute everything if mstop isn't changed
            if (i != mstop) {
                mstop <<- i
                fit <<- RET$predict()
                u <<- ngradient(y, fit, weights)
            }
        } else {
            ## if prior reduction of mstop,
            ## first increase mstop to old value first
            if (mstop != length(xselect)) {
                mstop <<- length(xselect)
                fit <<- RET$predict()
                u <<- ngradient(y, fit, weights)
            }
            ## now fit the rest
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
                if (!cwlin) {
                    nm <- bl[[w]]$Xnames
                    cf <- matrix(0, nrow = length(nm), ncol = mstop)
                }
            } else {
                if (inherits(ens[ix][[1]], "bm_cwlin") && !cwlin) {
                    cftmp <- sapply(ens[ix], coef, all = TRUE)
                } else {
                    cftmp <- sapply(ens[ix], coef)
                }
                nr <- NROW(cftmp)
                if (!is.matrix(cftmp)) nr <- 1
                cf <- matrix(0, nrow = nr, ncol = mstop)
                cf[, which(ix)] <- cftmp
            }
            if (!is.matrix(cf)) cf <- matrix(cf, nrow = 1)
            ## check if base-learner has coefficients
            if(any(sapply(cf, is.null))){
                ret <- NULL
            } else {
                ret <- switch(aggregate,
                              "sum" = rowSums(cf) * nu,
                              "cumsum" = {
                                  .Call("R_mcumsum", as(cf, "matrix") * nu,
                                        PACKAGE = "mboost")
                              },
                              "none" = nu * cf
                              )
            }
            ### set names, but not for bolscw base-learner
            if (!cwlin) {
                nm <- bl[[w]]$Xnames
                if (is.matrix(ret)) {
                    rownames(ret) <- nm
                } else {
                   names(ret) <- nm
               }
            }
            ret
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
mboost <- function(formula, data = list(), na.action = na.omit,
    baselearner = c("bbs", "bols", "btree", "bss", "bns"), ...) {


    ## We need at least variable names to go ahead
    if (length(formula[[3]]) == 1) {
        if (as.name(formula[[3]]) == ".") {
            formula <- as.formula(paste(deparse(formula[[2]]),
                "~", paste(names(data)[names(data) != all.vars(formula[[2]])],
                           collapse = "+"), collapse = ""))
        }
    }

    if (is.data.frame(data)) {
        if (!all(complete.cases(data))) {
            ## drop cases with missing values in any of the specified variables:
            vars <- all.vars(formula)[all.vars(formula) %in% names(data)]
            data <- na.action(data[, vars])
        }
    } else {
        if (any(unlist(lapply(data, is.na))))
            warning(sQuote("data"),
                    " contains missing values. Results might be affected. Consider removing missing values.")
    }

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
        bl[[1]]$set_names(deparse(formula[[3]]))
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
blackboost <- function(formula, data = list(), weights = NULL,
    na.action = na.pass,
    tree_controls = party::ctree_control(teststat = "max",
                               testtype = "Teststatistic",
                               mincriterion = 0,
                               maxdepth = 2, savesplitstats = FALSE),
    ...) {

    ### get the model frame first
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "weights", "na.action"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    response <- model.response(mf)
    weights <- model.weights(mf)
    ## drop outcome
    mf <- mf[,-1, drop = FALSE]
    ## drop weights
    mf$"(weights)" <- NULL
    bl <- list(btree(mf, tree_controls = tree_controls))
    ret <- mboost_fit(bl, response = response, weights = weights, ...)
    ret$call <- cl
    ret$rownames <- rownames(mf)
    class(ret) <- c("blackboost", class(ret))
    ret
}

### fit a linear model componentwise
glmboost <- function(x, ...) UseMethod("glmboost", x)

glmboost.formula <- function(formula, data = list(), weights = NULL,
                             na.action = na.pass, contrasts.arg = NULL,
                             center = TRUE, control = boost_control(), ...) {

    ## We need at least variable names to go ahead
    if (length(formula[[3]]) == 1) {
        if (as.name(formula[[3]]) == ".") {
            formula <- as.formula(paste(deparse(formula[[2]]),
                "~", paste(names(data)[names(data) != all.vars(formula[[2]])],
                           collapse = "+"), collapse = ""))
        }
    }

    ### get the model frame first
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "weights", "na.action"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf$data <- data ## use cc data
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    ### center argument moved to this function
    if (!control$center) {
        center <- FALSE
        warning("boost_control(center = FALSE) is deprecated, use glmboost(..., center = FALSE)")
    }
    ### set up the model.matrix and center (if requested)
    X <- model.matrix(attr(mf, "terms"), data = mf,
                      contrasts.arg = contrasts.arg)
    assign <- attr(X, "assign")
    cm <- numeric(ncol(X))
    if (center) {
        if (!attr(attr(mf, "terms"), "intercept") == 1)
            warning("model with centered covariates does not contain intercept")
        cm <- colMeans(X, na.rm = TRUE)
        cm[assign == 0] <- 0
        X <- scale(X, center = cm, scale = FALSE)
    }
    ### this function will be used for predictions later
    newX <- function(newdata) {
        mf <- model.frame(delete.response(attr(mf, "terms")),
                          data = newdata, na.action = na.action)
        X <- model.matrix(delete.response(attr(mf, "terms")),
                          data = mf, contrasts.arg = contrasts.arg)
        scale(X, center = cm, scale = FALSE)
    }

    ### component-wise linear models baselearner
    bl <- list(bolscw(X))
    response <- model.response(mf)
    weights <- model.weights(mf)
    ret <- mboost_fit(bl, response = response, weights = weights,
                      control = control, ...)
    ret$newX <- newX
    ret$assign <- assign
    ret$center <- cm
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
    ### specialized method for model.frame
    ret$model.frame <- function(which = NULL) {
        if (!is.null(which))
            warning("Argument ", sQuote("which"), " is ignored")
        mf
    }
    ### save standard update function for re-use
    update <- ret$update
    ### needs a specialized update function as well
    ret$update <- function(weights = NULL, oobweights = NULL, risk = "oobag",
                           trace = NULL) {
        ## call standard update function
        res <- update(weights = weights, oobweights = oobweights, risk = risk,
                      trace = trace)
        ## now re-set all special arguments
        res$newX <- newX
        res$assign <- assign
        res$center <- cm
        res$call <- cl
        ### need specialized method (hatvalues etc. anyway)
        res$hatvalues <- function() {
            H <- vector(mode = "list", length = ncol(X))
            MPinv <- res$basemodel[[1]]$MPinv()
            for (j in unique(res$xselect()))
                H[[j]] <- (X[,j] %*% MPinv[j, ,drop = FALSE]) * control$nu
            H
        }
        res$rownames <- rownames(mf)
        ### specialized method for model.frame
        res$model.frame <- function(which = NULL) {
            if (!is.null(which))
                warning("Argument ", sQuote("which"), " is ignored")
            mf
        }
        class(res) <- c("glmboost", "mboost")
        res
    }
    class(ret) <- c("glmboost", "mboost")
    return(ret)
}

glmboost.matrix <- function(x, y, center = TRUE, weights = NULL,
                            na.action = na.pass,
                            control = boost_control(), ...) {

    X <- x
    if (nrow(X) != NROW(y))
        stop("dimensions of ", sQuote("x"), " and ", sQuote("y"),
             " do not match")
    if (is.null(colnames(X)))
        colnames(X) <- paste("V", 1:ncol(X), sep = "")

    ## drop cases with missing values in any of the specified variables:
    if (any(!Complete.cases(cbind(X, y)))) {
        X <- na.action(X)
        if (!is.null(removed <- attr(X, "na.action"))) {
            y <- y[-removed]
        }
    }

    if (!control$center) {
        center <- FALSE
        warning("boost_control(center = FALSE) is deprecated, use glmboost(..., center = FALSE)")
    }
    assign <- numeric(ncol(X))
    cm <- numeric(ncol(X))
    ### guess intercept
    intercept <- which(colSums(abs(scale(X, center = TRUE, scale = FALSE)), na.rm=TRUE)
                                   < .Machine$double.eps)
    if (length(intercept) > 0)
        intercept <- intercept[colSums(abs(X[, intercept, drop = FALSE]), na.rm=TRUE)
                               > .Machine$double.eps]
    INTERCEPT <- length(intercept) == 1
    if (INTERCEPT) {
        assign[-intercept] <- 1:(ncol(X) - 1)
    } else {
        assign <- 1:ncol(X)
    }
    if (center) {
        cm <- colMeans(X, na.rm = TRUE)
        if (!INTERCEPT)
            warning("model with centered covariates does not contain intercept")
        cm[assign == 0] <- 0
        X <- scale(X, center = cm, scale = FALSE)
    }
    newX <- function(newdata) {
        if (isMATRIX(newdata)) {
            if (all(colnames(X) == colnames(newdata)))
                return(scale(newdata, center=cm, scale=FALSE))
        }
        stop(sQuote("newdata"), " is not a matrix with the same variables as ",
             sQuote("x"))
        return(NULL)
    }
    bl <- list(bolscw(X))
    ret <- mboost_fit(bl, response = y, control = control, ...)
    ret$newX <- newX
    ret$assign <- assign
    ret$center <- cm
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
    ### specialized method for model.frame
    ret$model.frame <- function(which = NULL) {
        if (!is.null(which))
            warning("Argument ", sQuote("which"), " is ignored")
        X
    }
    ### save standard update function for re-use
    update <- ret$update
    ### needs a specialized update function as well
    ret$update <- function(weights = NULL, oobweights = NULL, risk = "oobag",
                           trace = NULL) {
        ## call standard update function
        res <- update(weights = weights, oobweights = oobweights, risk = risk,
                      trace = trace)
        ## now re-set all special arguments
        ret$newX <- newX
        res$assign <- assign
        res$center <- cm
        res$call <- match.call()
        ### need specialized method (hatvalues etc. anyway)
        res$hatvalues <- function() {
            H <- vector(mode = "list", length = ncol(X))
            MPinv <- res$basemodel[[1]]$MPinv()
            for (j in unique(res$xselect()))
                H[[j]] <- (X[,j] %*% MPinv[j, ,drop = FALSE]) * control$nu
            H
        }
        res$rownames <- rownames(X)
        ### specialized method for model.frame
        res$model.frame <- function(which = NULL) {
            if (!is.null(which))
                warning("Argument ", sQuote("which"), " is ignored")
            X
        }
        class(res) <- c("glmboost", "mboost")
        res
    }
    class(ret) <- c("glmboost", "mboost")
    return(ret)
}

glmboost.default <- function(x, ...) {
    if (extends(class(x), "Matrix"))
        return(glmboost.matrix(x = x, ...))
    stop("no method for objects of class ", sQuote(class(x)),
         " implemented")
}
