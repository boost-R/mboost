
confint.mboost <- function(object, parm = NULL, level = 0.95,
                           B = 1000, B.mstop = 25, newdata = NULL,
                           which = parm, ...) {

    which <- object$which(which, usedonly = FALSE)
    if (!all(which %in% object$which(NULL, usedonly = FALSE)))
        stop(sQuote("which"), " is wrongly specified")

    ## create new data and/or restructure data
    newdata <- .create_newdata(object, newdata, which)

    outer.folds <- cv(model.weights(object), B = B)
    predictions <- vector("list", B)

    cat("Start computing bootstrap confidence intervals... \n")
    for (i in 1:B) {
        cat("\rB =", i)
        ## update model
        mod <- update(object, weights = outer.folds[, i],
                      risk = "inbag", trace = FALSE)
        if (B.mstop > 0) {
            ## <FIXME> are the weights handled correctly?
            cvr <- cvrisk(mod, folds = cv(model.weights(mod), B = B.mstop))
            mod[mstop(cvr)]
        }
        predictions[[i]] <- .predict_confint(mod, newdata = newdata,
                                             which = which)
    }
    cat("\n")

    ## prepare returned object
    res <- list(level = level, boot_pred = predictions, data = newdata,
                model = object)
    attr(res, "which") <- which
    class(res) <- "mboost.ci"
    return(res)
}

confint.glmboost <- function(object, parm = NULL, level = 0.95,
                             B = 1000, B.mstop = 25,
                             which = parm, ...) {

    outer.folds <- cv(model.weights(object), B = B)
    which <- object$which(which, usedonly = FALSE)
    if (!all(which %in% object$which(NULL, usedonly = FALSE)))
        stop(sQuote("which"), " is wrongly specified")

    coefficients <- matrix(NA, ncol = length(which), nrow = B)
    colnames(coefficients) <- names(coef(object, which = which))

    for (i in 1:B) {
        ## update model
        mod <- update(object, weights = outer.folds[, i],
                      risk = "inbag")
        if (B.mstop > 0) {
            ## <FIXME> are the weights handled correctly?
            cvr <- cvrisk(mod, folds = cv(model.weights(mod), B = B.mstop))
            mod[mstop(cvr)]
        }
        coefficients[i, ] <- unlist(coef(mod, which = which, off2int = TRUE))
    }

    ## prepare returned object
    res <- list(confint = .ci_glmboost(coefficients, level = level, which = which),
                level = level, boot_coefs = coefficients, model = object)
    attr(res, "which") <- which
    class(res) <- "glmboost.ci"
    return(res)
}

.ci_glmboost <- function(coefficients, level, which = NULL) {
    quantiles <- c((1 - level)/2, 1 - (1 - level)/2)

    tmp <- apply(coefficients, 2, FUN = quantile, probs = quantiles)
    CI <- as.data.frame(t(tmp))[which, ]
    return(CI)
}

## pe = add point estimte
print.glmboost.ci <- function(x, which = NULL, level = x$level, pe = FALSE, ...) {

    if (is.null(which)) {
        which <- attr(x, "which")
    } else {
        which <- x$model$which(which, usedonly = FALSE)
        if (!all(which %in% attr(x, "which")))
            stop(sQuote("which"), " is wrongly specified")
    }

    if (!is.null(level) && level != x$level) {
        CI <- .ci_glmboost(x$boot_coefs,  level = level, which = which)
    } else {
        CI <- x$confint[which, ]
    }

    if (pe) {
        tmp <- data.frame(beta = coef(x$model, which))
        CI <- cbind(tmp, CI)
    }
    if (length(which) > 1) {
        cat("\tBootstrap Confidence Intervals\n")
    } else {
        cat("\tBootstrap Confidence Interval\n")
    }
    print(CI, ...)
    return(invisible(x))
}


## ## Aditionally needed: Check for multivariate base-learners (except bols)

## FIXME: check for by variable and bivariate base-learners which both need a
## different data set for prediction
## FIXME: what about factor variables? do we get the correct levels?
.create_newdata <- function(object, newdata = NULL, which, ...) {
    if (is.null(newdata)) {
        data <- newdata <- model.frame(object, which = which)
        nms <- names(object$baselearner)[which]

        for (w in which) {
            ## get data from w-th base-learner
            tmp <- data[[w]][rep(1, 100), , drop = FALSE]

            ## are there varying coefficients (i.e. argument by)
            get_vary <- object$baselearner[[w]]$get_vary
            vary <- ""
            if (!is.null(get_vary)) vary <- get_vary()
            if (vary != "") {
                if (is.factor(tmp[[vary]])) {
                    if (nlevels(tmp[[vary]]) > 2)
                        stop("Atomatic data creation for ", sQuote("by"),
                             " variables with more than two levels is",
                             " currently not supported;",
                             " Specify ", sQuote("newdata"), " instead.")
                    data[[w]][[vary]] <- factor(levels(data[[w]][[vary]])[-1],
                                                levels = levels(data[[w]][[vary]]))
                }
                if (is.numeric(tmp[[vary]]))
                    data[[w]][[vary]] <- 1
            }

            ## now make grid
            grid <- function(x) {
                if (is.numeric(x)) {
                    return(seq(min(x), max(x), length = 100))
                } else {
                    return(rep(unique(x), length.out = 100))
                }
            }

            for (j in 1:ncol(data[[w]]))
                tmp[, colnames(data[[w]])[j]] <- grid(data[[w]][,j])

            ## FIXME: what about btree and bmrf?

            ## check if any base-learner is a bivariate smooth effect, i.e. if
            ## base-learner is multivariate and bbs, bspatial or brad
            which.vary <- colnames(tmp) == vary
            multivar <- grepl("bbs|bspatial|brad", nms[w]) &
                ncol(tmp[!which.vary]) >= 2
            if (multivar) {
                ## make grid
                egrid <- expand.grid(tmp[!which.vary])
                if (vary != "") {
                    x.vary <- tmp[which.vary]
                    rownames(x.vary) <- NULL
                    tmp <- cbind(egrid, x.vary)
                } else {
                    tmp <- egrid
                }
            }
            ## reset rownames
            rownames(tmp) <- NULL
            newdata[[w]] <- tmp
        }
    } else {
        ## restructure new data
        data <- model.frame(object, which = which)
        nms <- lapply(data, colnames)
        tmp <- lapply(nms, function(x) newdata[, x, drop = FALSE])
        newdata <- tmp
    }
    return(newdata)
}

## special prediction function for the construction of confidence intervals:
.predict_confint <- function(object, newdata = NULL, which, ...) {
    nrows <- sapply(newdata, nrow)
    predictions <- sapply(which, function(w)
                          matrix(predict(object, newdata[[w]], which = w),
                                               ncol = 1, nrow = nrows[w]))
    if (is.matrix(predictions))
        predictions <- as.data.frame(predictions)
    names(predictions) <- names(newdata[which])

    ## ###### FIXME FIXME FIXME
    ## WAS ist mit predict == 0?
    ## Wie geht es in .ci_mboost weiter?
    ## ###### FIXME FIXME FIXME
    return(predictions)
}


### plot functions
plot.mboost.ci <- function(x, which, level = x$level,
                           ylim = NULL, type = "l", col = "black",
                           ci.col = rgb(170, 170, 170, alpha = 85,
                                        maxColorValue = 255),
                           raw = FALSE, ...) {

    if (missing(which)) {
        which <- attr(x, "which")
    } else {
        which <- x$model$which(which, usedonly = FALSE)
        if (!all(which %in% attr(x, "which")))
            stop(sQuote("which"), " is wrongly specified")
    }
    if (length(which) > 1)
        stop("Specify a single base-learner using ", sQuote("which"))

    CI <- .ci_mboost(x$boot_pred, level = level, which = which, raw = raw)

    if (is.null(ylim)) {
        ylim <- range(CI)
    }

    plot(x$model, which = which, type = "n", ylim = ylim,
         col = col, ...)
    lines(x, which, level, col = ci.col, raw = raw, ...)
    lines(x$model, which = which, type = "l", col = col, ...)
}

lines.mboost.ci <- function(x, which, level = x$level,
                            col = rgb(170, 170, 170, alpha = 85,
                                      maxColorValue = 255),
                            raw = FALSE, ...) {

    if (missing(which)) {
        which <- attr(x, "which")
    } else {
        which <- x$model$which(which, usedonly = FALSE)
        if (!all(which %in% attr(x, "which")))
            stop(sQuote("which"), " is wrongly specified")
    }
    if (length(which) > 1)
        stop("Specify a single base-learner using ", sQuote("which"))


    CI <- .ci_mboost(x$boot_pred, level = level, which = which, raw = raw)

    x.data <- x$data[[which]]
    if (ncol(x.data) > 1) {
        stop("Cannot plot lines for more than 1 dimension")
    } else {
        x.data <- x.data[, 1]
    }

    if (!raw) {
        polygon(c(x.data, rev(x.data)),
                c(CI[1, ], rev(CI[2,])),
                col = col, border = col)
    } else {
        matlines(x$data[[which]], CI, col = col, lty = "solid", ...)
    }
}

.ci_mboost <- function(predictions, level, which = NULL, raw = FALSE) {

    preds <- sapply(predictions, function(p) p[[which]])
    if (!raw) {
        quantiles <- c((1 - level)/2, 1 - (1 - level)/2)
        preds <- apply(preds, 1, FUN = quantile, probs = quantiles)
    }

    return(preds)
}
