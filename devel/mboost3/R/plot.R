
### just a try
plot.mboost <- function(x, which = NULL, newdata = NULL,
                        type = "b", rug = TRUE, ylim = NULL, 
                        xlab = variable.names(x), ylab = expression(f[partial]), ...) {

    which <- x$which(which, usedonly = TRUE)
    pr <- predict(x, which = which)
    if (is.null(ylim)) ylim <- range(pr)
    stopifnot(length(xlab) %in% c(1, length(variable.names(x))))
    stopifnot(length(ylab) %in% c(1, length(variable.names(x))))

    for (w in which) {

        data <- model.frame(x, which = w)[[1]]
        get_vary <- x$baselearner[[w]]$get_vary
        vary <- ""
        if (!is.null(get_vary)) vary <- get_vary()
        if (!is.null(newdata)) data <- newdata[, colnames(data), drop = FALSE]
        if (vary != "") {
            v <- data[[vary]]
            if (is.factor(v)) v <- factor(rep(levels(v)[2], 1), levels = levels(v))
            if (is.numeric(v)) v <- 1
            data[[vary]] <- v[rep(1, nrow(data))]
        }
        pr <- predict(x, newdata = data, which = w)
        if (!is.null(vary)) data <- data[, colnames(data) != vary, drop = FALSE]

        xl <- ifelse(length(xlab) > max(which), xlab[w], xlab[1])
        yl <- ifelse(length(ylab) > max(which), ylab[w], ylab[1])

        if (ncol(data) == 1) {
            plot(sort(data[[1]]), pr[order(data[[1]])], type = type, 
                 xlab = xl, ylab = yl, ylim = ylim, ...)
            if (rug) rug(data[[1]])
        }
        if (ncol(data) == 2) {
            fm <- as.formula(paste("pr ~ ", paste(colnames(data), collapse = "*"), sep = ""))
            print(levelplot(fm, data = data, ...))
        }
        if (ncol(data) > 2) {
            for (v in colnames(data)) {
                tmp <- data
                pardep <- sapply(data[[v]], function(vv) {
                    tmp[[v]] <- vv
                    mean(predict(x, newdata = tmp, which = w))
                })
                plot(sort(data[[v]]), pardep[order(data[[v]])], type = type,
                     xlab = v, ylab = "Partial Dependency", ylim = ylim, ...)
            }
        }
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
