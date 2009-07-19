
##
## cross-validation (bootstrap, k-fold cv etc.) of empirical risk
## for boosting algorithms
##

Cvrisk <- function(object, folds = NULL, grid = c(1:mstop(object)), 
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

    myapply <- lapply
    if (parallel) {
        if (!multicore:::isChild())
            myapply <- mclapply
    }

    fam_name <- object$family@name
    call <- deparse(object$call)

    dummyfct <- function(weights){
        model <- fitfct(weights = weights)
        if (!is.null(fun)) return(fun(model))
        ret <- model$risk()[grid]
        ret
    }

    oobrisk <- myapply(1:ncol(folds), function(i) dummyfct(folds[,i]),
                       ...)
    if (!is.null(fun)) return(oobrisk)
    oobrisk <- t(as.data.frame(oobrisk))
    oobrisk <- oobrisk/colSums(folds == 0)
    colnames(oobrisk) <- grid
    rownames(oobrisk) <- 1:nrow(oobrisk)
    attr(oobrisk, "risk") <- fam_name
    attr(oobrisk, "call") <- call
    attr(oobrisk, "mstop") <- grid
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

plot.cvrisk <- function(x, ylab = attr(x, "risk"), ylim = range(x),
                        main = attr(x, "call"), ...) {

    cm <- colMeans(x)
    plot(1:ncol(x), cm, ylab = ylab, ylim = ylim,
         type = "n", lwd = 2,
         xlab = "Number of boosting iterations",
         main = main, ...)
    out <- apply(x, 1, function(y) lines(1:ncol(x),y, col = "lightgrey"))
    rm(out)
    ms <- which.min(cm)
    lines(c(ms, ms), c(min(c(0, ylim[1] * ifelse(ylim[1] < 0, 2, 0.5))), cm[ms]), lty = 2)
    lines(1:ncol(x), cm, type = "l")
}

mstop.cvrisk <- function(object, ...)
    attr(object, "mstop")[which.min(colSums(object))]
