
##
## cross-validation (bootstrap, k-fold cv etc.) of empirical risk
## for boosting algorithms
##

cvrisk <- function(object, ...)
    UseMethod("cvrisk")

cvrisk.mboost <- function (object, folds = cv(model.weights(object)),
                           grid = 0:mstop(object), papply = mclapply,
                           fun = NULL, mc.preschedule = FALSE,
                           ...) {
    
    papply <- match.fun(papply)
    weights <- model.weights(object)
    if (any(weights == 0))
        warning("zero weights")
    if (is.null(folds)) {
        folds <- rmultinom(25, length(weights), weights/sum(weights))
        attr(folds, "type") <- "25-fold bootstrap"
    } else {
        stopifnot(is.matrix(folds) && nrow(folds) == length(weights))
    }
    fitfct <- object$update
    oobrisk <- matrix(0, nrow = ncol(folds), ncol = length(grid))
    if (!is.null(fun))
        stopifnot(is.function(fun))
    fam_name <- object$family@name
    call <- deparse(object$call)
    if (is.null(fun)) {
        dummyfct <- function(weights, oobweights) {
            mod <- fitfct(weights = weights, oobweights = oobweights)
            mstop(mod) <- max(grid)
            ## return all risk values in grid (+ 1 as 0 is included)
            risk(mod)[grid + 1]
        }
        if (fam_name == "Cox Partial Likelihood" && all(colSums(folds == 0) == 1))
            stop("Leave-one-out cross-validation cannot be used with ", sQuote("family = CoxPH()"))
    } else { ## !is.null(fun)
        dummyfct <- function(weights, oobweights) {
            mod <- fitfct(weights = weights, oobweights = oobweights)
            mod[max(grid)]
            ## make sure dispatch works correctly
            class(mod) <- class(object)
            fun(mod)
        }
    }

    ## use case weights as out-of-bag weights (but set inbag to 0)
    OOBweights <- matrix(rep(weights, ncol(folds)), ncol = ncol(folds))
    OOBweights[folds > 0] <- 0
    if (identical(papply, mclapply)) {
        oobrisk <- papply(1:ncol(folds),
                          function(i) try(dummyfct(weights = folds[, i],
                                                   oobweights = OOBweights[, i]),
                                          silent = TRUE),
                          mc.preschedule = mc.preschedule,
                          ...)
    } else {
        oobrisk <- papply(1:ncol(folds),
                          function(i) try(dummyfct(weights = folds[, i],
                                                   oobweights = OOBweights[, i]),
                                          silent = TRUE),
                          ...)
    }
    ## if any errors occured remove results and issue a warning
    if (any(idx <- sapply(oobrisk, is.character))) {
        warning(sum(idx), " fold(s) encountered an error. ",
                "Results are based on ", ncol(folds) - sum(idx),
                " folds only.\n",
                "Original error message(s):\n",
                sapply(oobrisk[idx], function(x) x))
        oobrisk[idx] <- NA
    }
    if (!is.null(fun))
        return(oobrisk)
    oobrisk <- t(as.data.frame(oobrisk))
    oobrisk <- oobrisk / colSums(OOBweights)
    colnames(oobrisk) <- grid
    rownames(oobrisk) <- 1:nrow(oobrisk)
    attr(oobrisk, "risk") <- fam_name
    attr(oobrisk, "call") <- call
    attr(oobrisk, "mstop") <- grid
    attr(oobrisk, "type") <- ifelse(!is.null(attr(folds, "type")),
                                    attr(folds, "type"), "user-defined")
    class(oobrisk) <- "cvrisk"
    oobrisk
}

print.cvrisk <- function(x, ...) {
    cat("\n\t Cross-validated", attr(x, "risk"), "\n\t",
        attr(x, "call"), "\n\n")
    print(colMeans(x, na.rm = TRUE))
    cat("\n\t Optimal number of boosting iterations:", mstop(x), "\n")
    return(invisible(x))
}

plot.cvrisk <- function(x, xlab = "Number of boosting iterations",
                        ylab = attr(x, "risk"),
                        ylim = range(x), main = attr(x, "type"), ...) {
    
    ## force evaluation of attributes (to not evaluate them in the wrong situation)
    force(ylab)
    force(main)
    
    mstops <- attr(x, "mstop")

    x <- x[, apply(x, 2, function(y) all(!is.na(y))), drop = FALSE]
    cm <- colMeans(x)
    plot(1:ncol(x), cm, ylab = ylab, ylim = ylim,
         type = "n", lwd = 2, xlab = xlab,
         main = main, axes = FALSE, ...)
    out <- apply(x, 1, function(y) lines(1:ncol(x),y, col = "lightgrey"))
    rm(out)
    ms <- which.min(cm)
    lines(c(ms, ms), c(min(c(0, ylim[1] * ifelse(ylim[1] < 0, 2, 0.5))), cm[ms]),
          lty = 2)
    lines(1:ncol(x), cm, type = "l")
    axis(1, at = 1:ncol(x), labels = mstops)
    axis(2)
    box()
}

mstop.cvrisk <- function(object, ...)
    attr(object, "mstop")[which.min(colSums(object, na.rm = TRUE))]

cv <- function(weights, type = c("bootstrap", "kfold", "subsampling"),
               B = ifelse(type == "kfold", 10, 25),
               prob = 0.5, strata = NULL) {

    type <- match.arg(type)
    n <- length(weights)

    if (is.null(strata)) strata <- gl(1, n)
    if (!is.factor(strata)) stop(sQuote("strata"), " must be a factor")
    folds <- matrix(0, nrow = n, ncol = B)

    ### <FIXME> handling of weights needs careful documentation </FIXME>
    for (s in levels(strata)) {
        indx <- which(strata == s)
        folds[indx,] <- switch(type,
                               "bootstrap" = cvboot(length(indx), B = B, weights[indx]),
                               "kfold" = cvkfold(length(indx), k = B) * weights[indx],
                               "subsampling" = cvsub(length(indx), prob = prob, B = B) * weights[indx])
    }
    attr(folds, "type") <- paste(B, "-fold ", type, sep = "")
    return(folds)
}


cvboot <- function(n, B, weights)
    rmultinom(B, n, weights / sum(weights))

cvkfold <- function(n, k) {
    #if (k > n / 2) stop("k > n/2")
    fl <- floor(n/k)
    folds <- c(rep(c(rep(0, fl), rep(1, n)), k - 1),
               rep(0, n * k - (k - 1) * (fl + n)))
    matrix(folds, nrow = n)[sample(1:n),, drop = FALSE]
}

cvsub <- function(n, prob, B) {
    k <- floor(n * prob)
    indx <- rep(c(0, 1), c(n - k, k))
    replicate(B, sample(indx))[sample(1:n),, drop = FALSE]
}
