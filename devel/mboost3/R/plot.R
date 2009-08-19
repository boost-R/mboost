
### just a try
plot.mboost <- function(x, which = NULL, newdata = NULL, 
                        type = "b", rug = TRUE, eylim = TRUE, ...) {

    if (is.null(which)) which <- sort(unique(x$xselect()))
    pr <- predict(x, newdata = newdata, which = which)
    mf <- model.frame(x, which = which)
    if (!is.null(newdata)) {
        for (i in 1:ncol(pr))
            mf[[i]] <- newdata[, colnames(mf[[i]]), drop = FALSE]
    }

    if (eylim) ylim <- range(pr)

    for (i in 1:ncol(pr)) {
        dat <- mf[[i]]
        p <- pr[,i]
        if (!eylim) ylim <- range(p)

        if (ncol(dat) == 1) {
            plot(sort(dat[[1]]), p[order(dat[[1]])], type = type, xlab = names(mf[[i]]),
                 ylab = names(mf)[i], ylim = ylim, ...)
            if (rug) rug(dat[[1]])
        }
        if (ncol(dat) == 2) {
            dat$pr <- p
            coordinates(dat) <- as.formula(paste("~", paste(names(mf[[i]]), collapse = "+"), sep = ""))
            print(spplot(dat, "pr", xlab = names(mf[[i]])[1], ylab = names(mf[[i]])[2], ...))
        }
        if (ncol(dat) > 2)
            stop("not yet implemented")
    }
}

plot.glmboost <- function(x, main = deparse(x$call), col = NULL, ...) {

    cp <- coef(x, aggregate = "cumsum")
    cf <- cp[, ncol(cp)]
    if (is.null(col))  
        col <- hcl(h = 40, l = 50, c= abs(cf) / max(abs(cf)) * 490)
    matplot(t(cp), type = "l", lty = 1, xlab = "Number of boosting iterations",
            ylab = "Coefficients", main = main, col = col, ...)
    abline(h = 0, lty = 1, col = "lightgray")
    axis(4, at = cf, labels = variable.names(x), las = 1)

}
