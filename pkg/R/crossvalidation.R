
##
## cross-validation (bootstrap, k-fold cv etc.) of empirical risk
## for boosting algorithms
##

cvrisk <- function(object, folds = cv(model.weights(object)), grid = 1:mstop(object),
                   parallel = require("multicore"), fun = NULL, ...) {

    weights <- model.weights(object)
    if (any(weights == 0)) warning("zero weights")
    if (is.null(folds)) {
        ### bootstrap by default
        folds <- rmultinom(25, length(weights), weights / sum(weights))
    } else {
        stopifnot(is.matrix(folds) &&
                  nrow(folds) == length(weights))
    }

    fitfct <- object$update
    oobrisk <- matrix(0, nrow = ncol(folds), ncol = length(grid))
    if (!is.null(fun)) stopifnot(is.function(fun))

    fam_name <- object$family@name
    call <- deparse(object$call)

    if (is.null(fun)) {
        dummyfct <- function(weights) {
            mod <- fitfct(weights = weights)
            mod[max(grid)]
            mod$risk()[grid]
        }
    } else {
        dummyfct <- function(weights) {
            mod <- fitfct(weights = weights)
            mod[max(grid)]
            ### make sure dispatch works correctly
            class(mod) <- class(object)
            fun(mod)
        }
    }

    oobrisk <- MYapply(1:ncol(folds), function(i) dummyfct(folds[,i]),
                       parallel = parallel, ...)
    if (!is.null(fun)) return(oobrisk)
    oobrisk <- t(as.data.frame(oobrisk))
    oobrisk <- oobrisk/colSums(folds == 0)
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
    print(colMeans(x))
    cat("\n\t Optimal number of boosting iterations:", mstop(x), "\n")
    return(invisible(x))
}

plot.cvrisk <- function(x, ylab = attr(x, "risk"),
                        xlab = "Number of boosting iterations",
                        ylim = range(x), main = attr(x, "type"), ...) {

    cm <- colMeans(x)
    plot(1:ncol(x), cm, ylab = ylab, ylim = ylim,
         type = "n", lwd = 2, xlab = xlab,
         main = main, axes = FALSE, ...)
    out <- apply(x, 1, function(y) lines(1:ncol(x),y, col = "lightgrey"))
    rm(out)
    ms <- which.min(cm)
    lines(c(ms, ms), c(min(c(0, ylim[1] * ifelse(ylim[1] < 0, 2, 0.5))), cm[ms]), lty = 2)
    lines(1:ncol(x), cm, type = "l")
    axis(1, at = 1:ncol(x), label = attr(x, "mstop"))
    axis(2)
    box()
}

mstop.cvrisk <- function(object, ...)
    attr(object, "mstop")[which.min(colSums(object))]

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
    if (k > n / 2) stop("k > n/2")
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
