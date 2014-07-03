
confint.mboost <- function(object, B = 1000, newdata = NULL,
                           B.mstop = 25, which = NULL, ...) {

    which <- object$which(which, usedonly = FALSE)

    ## create new data and/or restructure data
    newdata <- .create_newdata(object, newdata, which)

    outer.folds <- cv(model.weights(object), B = B)
    predictions <- vector("list", B)

    for (i in 1:B) {
        ## update model
        mod <- update(object, weights = outer.folds[, i],
                      risk = "inbag")
        if (B.mstop > 0) {
            ## <FIXME> are the weights handled correctly?
            cvr <- cvrisk(mod, folds = cv(model.weights(mod), B = B.mstop))
            mod[mstop(cvr)]
        }
        predictions[[i]] <- .predict_confint(mod, newdata = newdata,
                                             which = which)
    }
    res <- list(boot_pred = predictions, data = newdata, model = object)
    class(res) <- "mboost.ci"
    return(res)
}

confint.glmboost <- function(object, B = 1000, B.mstop = 25,
                             which = NULL, ...) {

    outer.folds <- cv(model.weights(object), B = B)
    which <- object$which(which, usedonly = FALSE)

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
    res <- list(boot_coefs = coefficients, model = object)
    class(res) <- "glmboost.ci"
    return(res)
}

print.glmboost.ci <- function(x, which = NULL, level = 0.95) {
    quantiles <- c((1 - level)/2, 1 - (1 - level)/2)
    which <- x$model$which(which, usedonly = FALSE)
    tmp <- apply(x$boot_coefs[, which], 2, FUN = quantile, probs = quantiles)
    print(t(tmp))
}

## ## check for varing...
## data <- model.frame(x, which = w)[[1]]
## get_vary <- x$baselearner[[w]]$get_vary
## vary <- ""
## if (!is.null(get_vary)) vary <- get_vary()
## if (!is.null(newdata)) data <- newdata[, colnames(data), drop = FALSE]
## if (vary != "") {
##     v <- data[[vary]]
##     if (is.factor(v)) v <- factor(levels(v)[-1], levels = levels(v))
##     if (is.numeric(v)) v <- 1
## }

## ## Aditionally needed: Check for multivariate base-learners (except bols)


## check for by variable and bivariate base-learners which both need a different
## data set for prediction
.create_newdata <- function(object, newdata = NULL, which, ...) {
    if (is.null(newdata)) {
        data <- newdata <- model.frame(object, which = which)
        for (w in which) {
            ## make grid!
            tmp <- data[[w]][rep(1, 100), , drop = FALSE]
            grid <- function(x) {
                if (is.numeric(x)) {
                    return(seq(min(x), max(x), length = 100))
                } else {
                    return(rep(levels(x), length.out = 100))
                }
            }
            for (j in 1:ncol(data[[w]]))
                tmp[, colnames(data[[w]])[j]] <- grid(data[[w]][,j])
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

## .create_newdata.glmboost <- function(object, newdata, ...) {
##     if (is.null(newdata)) {
##         data <- model.frame(object)
##         ## make grid!
##         tmp <- data[rep(1, 100), ]
##         grid <- function(x) {
##             if (is.numeric(x)) {
##                 return(seq(min(x), max(x), length = 100))
##             } else {
##                 return(rep(levels(x), length.out = 100))
##             }
##         }
##         for (j in 1:ncol(data))
##             tmp[, colnames(data)[j]] <- grid(data[,j])
##         newdata <- tmp
##     }
##     return(newdata)
## }

## special prediction function for the construction of confidence intervals:
.predict_confint <- function(object, newdata = NULL, which, ...) {
    predictions <- matrix(NA, ncol = length(which), nrow = nrow(newdata[[1]]))
    for (w in which) {
        predictions[, w] <- predict(object, newdata[[w]], which = w)
    }
    return(predictions)
}

# .predict_confint.glmboost <- function(object, newdata, which, ...) {
#     warning("shouldn't we return confints for coef?")
#     predict(object, newdata = newdata, which = which)
# }

plot.mboost.ci <- function(x, which, level = 0.95,
                           ylim = NULL, type = "l", col = "black",
                           ci.col = "grey",  raw = FALSE, ...) {

    which <- x$model$which(which, usedonly = FALSE)

    if (is.null(ylim)) {
        preds <- sapply(x$boot_pred, function(p) p[, which])
        if (!raw) {
            quantiles <- c((1 - level)/2, 1 - (1 - level)/2)
            tmp <- apply(preds, 1, FUN = quantile, probs = quantiles)
            if (is.null(ylim))
                ylim <- range(tmp)
        } else {
            if (is.null(ylim))
                ylim <- range(preds)
        }
    }

    plot(x$model, which = which, type = "n", ylim = ylim,
         col = col, ...)
    lines(x, which, level, col = ci.col, raw = raw, ...)
    lines(x$model, which = which, type = "l", col = col, ...)
}

lines.mboost.ci <- function(x, which, level = 0.95, col = "grey",
                            raw = FALSE, ...) {
    preds <- sapply(x$boot_pred, function(p) p[, which])
    x.data <- x$data[[which]]
    if (ncol(x.data) > 1) {
        stop("Cannot plot lines for more than 1 dimenstion")
    } else {
        x.data <- x.data[, 1]
    }
    if (!raw) {
        quantiles <- c((1 - level)/2, 1 - (1 - level)/2)
        tmp <- apply(preds, 1, FUN = quantile, probs = quantiles)
        polygon(c(x.data, rev(x.data)),
                c(tmp[1, ], rev(tmp[2,])),
                col = col, border = col)
    } else {
        matlines(x$data[[which]], preds, col = col, lty = "solid", ...)
    }
}
