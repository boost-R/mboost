
mboost <- function(blg, y, weights = NULL, offset = NULL, 
                   family = GaussReg(), control = boost_control()) {

    check_y_family(y, family)

    ### hyper parameters
    mstop <- 0
    risk <- control$risk
    constraint <- control$constraint
    nu <- control$nu
    trace <- control$trace
    tracestep <- options("width")$width / 2

    ### extract negative gradient and risk functions
    ngradient <- family@ngradient
    riskfct <- family@risk

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

    xselect <- integer(mstop)
    ens <- vector(mode = "list", length = mstop)

    ### vector of empirical risks for all boosting iterations
    ### (either in-bag or out-of-bag)
    mrisk <- numeric(mstop)
    mrisk[1:mstop] <- NA
    tsums <- numeric(length(bl))
    ss <- vector(mode = "list", length = length(bl))

    fit <- offset
    if (is.null(offset))
        fit <- offset <- family@offset(y, weights)
    u <- ustart <- ngradient(y, fit, weights)

    ### start boosting iteration
    boost <- function(niter) {
        for (m in (mstop + 1):(mstop + niter)) {

            ### fit least squares to residuals _componentwise_
            for (i in 1:length(bl)) {
                tsums[i] <- -1
                ss[[i]] <- try(fit(bl[[i]], y = u))
                if (inherits(ss[[i]], "try-error")) next
                tsums[i] <- mean(weights * (fitted(ss[[i]]) - u)^2, na.rm = TRUE)
            }

            if (all(tsums < 0))
                stop("could not fit base learner in boosting iteration ", m)
            xselect[m] <<- order(tsums)[1]
            basess <- ss[[xselect[m]]]

            ### update step
            fit <<- fit + nu * fitted(basess)

            ### negative gradient vector, the new `residuals'
            u <<- ngradient(y, fit, weights)

            ### evaluate risk, either for the learning sample (inbag)
            ### or the test sample (oobag)
            if (risk == "inbag") mrisk[m] <<- riskfct(y, fit, weights)
            if (risk == "oobag") mrisk[m] <<- riskfct(y, fit, oobweights)

            ### save the model, i.e., the selected coefficient and variance
            if (control$saveensss)
                ens[[m]] <<- basess

            ## free memory
            rm("basess")

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
                response = y,           ### the response variable
                weights = weights       ### weights used for fitting
    )

    RET$update <- function(weights = NULL) {
        control$mstop <- mstop
        mboost(blg = blg, y = y, weights = weights, offset = offset,
               family = family, control = control)
    }

    RET$mstop <- function() mstop

    RET$xselect <- function() xselect[1:mstop]

    RET$fitted <- function() fit

    RET$risk <- function() risk[1:mstop]

    RET$predict <- function(newdata = NULL, which = NULL, 
                            components = FALSE, Sum = TRUE) {

        indx <- (1:length(xselect) <= mstop)
        if (is.null(which))
            which <- sort(unique(xselect[indx]))

        if (length(which) == 1)
            return(bl[[which]]$predict(ens[xselect == which & indx], 
                                       newdata = newdata, Sum = Sum))

        pr <- sapply(which, function(w) 
            bl[[w]]$predict(ens[xselect == w & indx], newdata = newdata))
        colnames(pr) <- names(bl)[which]
        if (components) return(nu * pr)
        offset + nu * rowSums(pr)
    }

    RET$model.frame <- function(which = NULL) {
        if (is.null(which)) which <- sort(unique(RET$xselect()))
        ret <- lapply(blg[which], model.frame)
        names(ret) <- names(bl)[which]
        ret
    }

    RET$subset <- function(i) {
        if (i <= mstop || length(xselect) > i) {
            mstop <<- i
            fit <<- RET$predict()
        } else {
            tmp <- boost(i - mstop)
        }
    } 

    RET$coef <- function(which = NULL, Sum = TRUE) {
        
        indx <- (1:length(xselect) <= mstop)
        if (is.null(which))
            which <- sort(unique(xselect[indx]))

        if (length(which) == 1) {
            cf <- sapply(ens[xselect == which & indx], coef)
            if (Sum) return(rowSums(cf) * nu)
            M <- triu(crossprod(Matrix(1, nc = ncol(cf))))
            return(nu * (cf %*% M))
        }

        ret <- lapply(which, function(w) {
            tmp <- ens[xselect == w & indx]
            cf <- 0
            for (i in 1:length(tmp))
                cf <- cf + coef(ens[[i]])
            return(cf * nu)
        })
        names(ret) <- names(bl)[which]
        return(ret)
    }

    ### function for computing hat matrices of individual predictors
    RET$hatvalues <- function(which = NULL) {
        if (is.null(which)) which <- sort(unique(RET$xselect()))
        ret <- lapply(bl[which], hatvalues)
        names(ret) <- names(bl)[which]
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
