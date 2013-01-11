
stabsel <- function(object, FWER = 0.05, cutoff, q,
                    folds = cv(model.weights(object), type = "subsampling", B = 100),
                    papply = mclapply, ...) {

    p <- length(variable.names(object))
    ibase <- 1:p

    if (!missing(q) && p < q)
        stop("Average number of selected base-learners ", sQuote("q"),
             " must be smaller \n  than the number of base-learners",
             " specified in the model ", sQuote("object"))

    if (!(FWER > 0 && FWER < 0.5))
        stop(sQuote("FWER"), " must be between 0 and 0.5")

    if (! xor(missing(cutoff), missing(q)))
        stop(" Either ", sQuote("cutoff"), " or ", sQuote("q"),
             "must be specified (but not both).")

    if (missing(cutoff)) {
        cutoff <- min(0.9, tmp <- (q^2 / (FWER * p) + 1) / 2)
        upperbound <- q^2 / p / (2 * cutoff - 1)
        if (tmp > 0.9 && upperbound - FWER > FWER/2) {
            warning("Upper bound for FWER >> ", FWER,
                    " for the given value of ", sQuote("q"),
                    " (true upper bound = ", upperbound, ")")
        }
    }
    if (missing(q)){
        stopifnot(cutoff >= 0.5)
        q <- ceiling(sqrt(FWER * (2 * cutoff - 1) * p))
        upperbound <- q^2 / p / (2 * cutoff - 1)
        if (upperbound - FWER > FWER/2)
            warning("Upper bound for FWER >> ", FWER,
                    " for the given value of ", sQuote("cutoff"),
                    " (true upper bound = ", upperbound, ")")
    }

    fun <- function(model) {
        xs <- selected(model)
        qq <- sapply(1:length(xs), function(x) length(unique(xs[1:x])))
        if (qq[length(xs)] < q)
            warning(sQuote("mstop"), " too small to select ", sQuote("q"),
                    " base-learners; Increase ", sQuote("mstop"),
                    " bevor applying ", sQuote("stabsel"))
        xs[qq > q] <- xs[1]
        xs
    }
    ss <- cvrisk(object, fun  = fun,
                 folds = folds,
                 papply = papply)
    ret <- matrix(0, nrow = length(ibase), ncol = m <- mstop(object))
    for (i in 1:length(ss)) {
        tmp <- sapply(ibase, function(x)
            ifelse(x %in% ss[[i]], which(ss[[i]] == x)[1], m + 1))
        ret <- ret + t(sapply(tmp, function(x) c(rep(0, x - 1), rep(1, m - x + 1))))
    }

    phat <- ret / length(ss)
    rownames(phat) <- names(variable.names(object))
    if (extends(class(object), "glmboost"))
        rownames(phat) <- variable.names(object)
    ret <- list(phat = phat, selected = which((mm <- apply(phat, 1, max)) >= cutoff),
                max = mm, cutoff = cutoff, q = q)
    class(ret) <- "stabsel"
    ret
}

print.stabsel <- function(x, ...) {

    cat("\tStability Selection\n")
    if (length(x$selected) > 0) {
        cat("\nSelected base-learners:\n")
        print(x$selected)
    } else {
        cat("\nNo base-learner selected\n")
    }
    cat("\nSelection probabilities:\n")
    print(x$max[x$max > 0])
    cat("\nCutoff: ", x$cutoff, "; ", sep = "")
    cat("q: ", x$q, "\n\n")
    invisible(x)
}

plot.stabsel <- function(x, main = deparse(x$call), col = NULL, ...) {

    h <- x$phat
    h <- h[rowSums(h) > 0, , drop = FALSE]
    if (is.null(col))
        col <- hcl(h = 40, l = 50, c = h[,ncol(h)] / max(h) * 490)
    matplot(t(h), type = "l", lty = 1, xlab = "Number of boosting iterations",
            ylab = "Selection Probability", main = main, col = col, ylim = c(0, 1), ...)
    abline(h = x$cutoff, lty = 1, col = "lightgray")
    axis(4, at = x$phat[rowSums(x$phat) > 0, ncol(x$phat)],
         labels = rownames(x$phat)[rowSums(x$phat) > 0], las = 1)
}


fitsel <- function(object, newdata = NULL, which = NULL, ...) {
    fun <- function(model) {
        tmp <- predict(model, newdata = newdata,
                       which = which, agg = "cumsum")
        ret <- c()
        for (i in 1:length(tmp))
            ret <- rbind(ret, tmp[[i]])
        ret
    }
    ss <- cvrisk(object, fun = fun, ...)
    ret <- matrix(0, nrow = nrow(ss[[1]]), ncol = ncol(ss[[1]]))
    for (i in 1:length(ss))
        ret <- ret + sign(ss[[i]])
    ret <- abs(ret) / length(ss)
    ret
}
