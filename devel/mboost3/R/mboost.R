
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
                do_trace(m, risk = mrisk, step = tracestep, width = niter)
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

    RET$xselect <- function() xselect[1:mstop]

    RET$fitted <- function() fit

    RET$risk <- function() mrisk[1:mstop]

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
                ret <- matrix(0, nrow = nrow(pr[[1]]), ncol = length(pr))
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
                ret <- matrix(0, nrow = nrow(pr[[1]]), ncol = length(pr))
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
        which <- thiswhich(which, usedonly = TRUE)
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

    RET$coef <- function(which = NULL, aggregate = TRUE) {
        
        indx <- ((1:length(xselect)) <= mstop)
        which <- thiswhich(which, usedonly = FALSE)
        if (length(which) == 0) return(NULL)

        if (length(which) == 1) {
            ix <- (xselect == which & indx)
            if (!any(ix)) return(NULL)
            cf <- sapply(ens[ix], coef)
            if (aggregate) {
                ret <- (rowSums(cf) * nu)
            } else {
                M <- triu(crossprod(Matrix(1, nc = ncol(cf))))
                ret <- nu * (cf %*% M)
            }
            attr(ret, "offset") <- offset
            return(ret)
        }

        if (!aggregate) 
            stop("aggregate=FALSE only available for one single baselearner")

        tmp <- lapply(which, function(w) {
            ix <- (xselect == w & indx)
            if (!any(ix)) return(NULL)
            tmp <- ens[ix]
            cf <- 0
            for (i in 1:length(tmp))
                cf <- cf + coef(tmp[[i]])
            return(cf * nu)
        })
        ret <- vector(mode = "list", length = length(bl))
        names(ret) <- names(bl)
        ret[which] <- tmp
        attr(ret, "offset") <- offset
        return(ret)
    }

    ### function for computing hat matrices of individual predictors
    RET$hatvalues <- function(which = NULL) {
        which <- thiswhich(which, usedonly = TRUE)
        tmp <- lapply(bl[which], hatvalues)
        ret <- vector(mode = "list", length = length(bl))
        names(ret) <- names(bl)
        ret[which] <- tmp
        ret
    }

    class(RET) <- "mboost"
    return(RET)
}

predict.mboost <- function(object, ...)
    object$predict(...)

coef.mboost <- function(object, ...)
    object$coef(...)

hatvalues.mboost <- function(model, ...)
    model$hatvalues(...)

fitted.mboost <- function(object, ...)
    object$fitted()

"[.mboost" <- function(x, i, ...) {
    x$subset(i)
    return(NULL)
}

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

Glmboost <- function(formula, data = list(), ...) {

    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    X <- model.matrix(attr(mf, "terms"), mf)
    bl <- list(bolscw(X))
    response <- eval(as.expression(formula[[2]]), envir = data)
    mboost_fit(bl, response = response, ...)
}

Gamboost <- function(formula, data = list(), baselearner = bbs3, ...) {

    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
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
