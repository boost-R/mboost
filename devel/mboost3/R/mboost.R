
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

    fit <- offset
    if (is.null(offset))
        fit <- offset <- family@offset(y, weights)
    u <- ustart <- ngradient(y, fit, weights)

    ### start boosting iteration
    boost <- function(niter) {
        for (m in (mstop + 1):(mstop + niter)) {

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
            if (control$saveensss)
                ens[[m]] <<- basess

            ### print status information
            if (trace)
                mboost:::do_trace(m, risk = mrisk, step = tracestep, width = niter)
        }
        mstop <<- mstop + niter 
        return(TRUE)
    }
    tmp <- boost(control$mstop)

    RET <- list(baselearner = blg,
                basemodel = bl,
                offset = offset,        ### offset
                ustart = ustart,        ### first negative gradients
                control = control,      ### control parameters
                family = family,        ### family object
                response = response,    ### the response variable
                weights = weights       ### weights used for fitting
    )

    RET$update <- function(weights = NULL) {
        control$mstop <- mstop
        mboost(blg = blg, response = response, weights = weights, 
               offset = offset, family = family, control = control)
    }

    RET$mstop <- function() mstop

    RET$xselect <- function(cw = FALSE) {
        if (inherits(bl[[1]], "bl_cwlin") && (length(bl) == 1 && cw))
            return(as.integer(sapply(ens[1:mstop], function(m) m$model[2])))
        return(xselect[1:mstop])
    }

    RET$fitted <- function() fit

    RET$resid <- function() u

    RET$risk <- function() mrisk[1:mstop]

    RET$logLik <- function() -mrisk[mstop]

    thiswhich <- function(which = NULL, usedonly = FALSE) {
        if (is.null(which)) which <- 1:length(bl)
        if (is.character(which)) {
            i <- match(which, names(bl))
            if (any(is.na(i)))
                warning(paste(which[is.na(i)], collapse = ","), " not found")
            which <- i
        } 
        if (usedonly) which <- which[which %in% RET$xselect()]
        return(which)
    }

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

    RET$model.frame <- function(which = NULL) {
        which <- thiswhich(which, usedonly = FALSE)
        tmp <- lapply(blg[which], model.frame)
        ret <- vector(mode = "list", length = length(bl))
        names(ret) <- names(bl)
        ret[which] <- tmp
        ret
    }

    RET$subset <- function(i) {
        if (i <= mstop || length(xselect) > i) {
            mstop <<- i
            fit <<- RET$predict()
            u <<- ngradient(y, fit, weights)
        } else {
            tmp <- boost(i - mstop)
        }
    } 

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

predict.mboost <- function(object, type = c("lp", "response"), ...) {
    pr <- object$predict(...)
    type <- match.arg(type)
    if (is.factor(y <- object$response) && type == "response")
        return(factor(levels(y)[(lp > 0) + 1], levels = levels(y)))
    return(pr)
}

coef.mboost <- function(object, ...)
    object$coef(...)

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


fitted.mboost <- function(object, ...) {
    args <- list(...)
    if (length(args) == 0)
        return(object$fitted())
    object$predict(...)
}

resid.mboost <- function(object, ...) 
    object$resid()

logLik.mboost <- function(object, ...)
    object$logLik()

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

mboost <- function(formula, data = list(), ...) {

    "+" <- function(a,b) {
        if (inherits(a, "blg")) a <- list(a)
        if (inherits(b, "blg")) b <- list(b)
        c(a, b)
    }
    bl <- eval(as.expression(formula[[3]]), envir = data)
    if (inherits(bl, "blg")) bl <- list(bl)
    stopifnot(all(sapply(bl, inherits, what = "blg")))
    nm <- strsplit(as.character(as.expression(formula[[3]])), "\\+")[[1]]
    nm <- gsub(" ", "", nm)
    names(bl) <- nm
    response <- eval(as.expression(formula[[2]]), envir = data)
    mboost_fit(bl, response = response, ...)
}

Glmboost <- function(formula, data = list(), weights = NULL, na.action = na.pass, 
                     contrasts.arg = NULL, center = FALSE, control = boost_control(), ...) {

    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "weights"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$na.action <- na.action
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    if (control$center) {
        center <- TRUE
        warning("boost_control center deprecated")
    }
    X <- model.matrix(attr(mf, "terms"), data = mf, 
                      contrasts.args = contrasts.args)
    cm <- rep(0, ncol(X))
    if (center) {
        cm <- colMeans(X, na.rm = TRUE)
        center <- attr(X, "assign") %in% which(sapply(mf, is.numeric)[-1])
        cm[!center] <- 0
        X <- scale(X, center = cm, scale = FALSE)
    }
    newX <- function(newdata) {
        X <- model.matrix(delete.response(attr(mf, "terms")), data = newdata,
                          contrasts.args = contrasts.args)
        scale(X, center = cm, scale = FALSE)
    }

    bl <- list(bolscw(X))
    response <- model.response(mf)
    weights <- model.weights(mf)
    ret <- mboost_fit(bl, response = response, weights = weights, 
                      control = control, ...)
    ret$newX <- newX
    class(ret) <- c("Glmboost", "mboost")
    return(ret)
}

predict.Glmboost <- function(object, newdata = NULL, ...) {

    if (!is.null(newdata)) {
        if (!isMATRIX(newdata)) {
            newdata <- object$newX(newdata)
        }
    }
    object$predict(newdata = newdata, ...)
}

hatvalues.Glmboost <- function(model, ...) {

    if (!checkL2(model)) return(hatvalues.mboost(model))
    Xf <- t(model$basemodel[[1]]$MPinv()) * model$control$nu
    X <- model.frame(model$baselearner[[1]])
    op <- .Call("R_trace_glmboost", X, Xf,
                as.integer(model$xselect(cw = TRUE)),
                PACKAGE = "mboost")
    RET <- diag(op[[1]])
    attr(RET, "hatmatrix") <- op[[1]]  
    attr(RET, "trace") <- op[[2]] 
    RET
}

Gamboost <- function(formula, data = list(), baselearner = bbs3, ...) {

    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    if (is.character(baselearner)) 
        baselearner <- get(baselearner, mode = "function", 
                           envir = parent.frame())
    bl <- vector(mode = "list", length = ncol(mf) - 1)
    names(bl) <- names(mf)[-1]
    for (i in 2:ncol(mf)) {
        tmp <- mf[[i]]
        bl[[i - 1]] <- baselearner(tmp)
        bl[[i - 1]]$set_names(names(mf)[i])
    }
    if (!inherits(bl, "blg") && !is.list(bl)) bl <- list(bbs3(bl))
    stopifnot(all(sapply(bl, inherits, what = "blg")))
    response <- mf[,1]
    mboost_fit(bl, response = response, ...)
}

Blackboost <- function(formula, data = list(), ...) {

    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    bl <- list(btree(mf[,-1, drop = FALSE]))
    names(bl) <- "btree"
    response <- mf[,1]
    mboost_fit(bl, response = response, ...)
}
