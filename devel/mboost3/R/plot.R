
### just a try
plot.mboost <- function(x, which = NULL, newdata = NULL, 
                        type = ifelse(is.null(newdata), "b", "l"),  
                        rug = TRUE, eylim = TRUE, ...) {

    pr <- predict(x, newdata = newdata, which = which, 
                  components  = TRUE)
    mf <- model.frame(x, which = which)
    if (!is.null(newdata)) {
        for (i in 1:length(pr)) 
            mf[[i]] <- newdata[, colnames(mf[[i]]), drop = FALSE]
    }

    if (eylim) ylim <- range(pr)

    for (i in 1:ncol(pr)) {
        dat <- mf[[i]]
        p <- pr[,i]
        if (!eylim) ylim <- range(p)

        if (ncol(dat) == 1) {
            plot(sort(dat[[1]]), p[order(dat[[1]])], type = type, xlab = names(mf[[i]]),
                 ylab = colnames(pr)[i], ylim = ylim, ...)
            if (rug) rug(dat[[1]])
        }
        if (ncol(mf[[i]]) == 2) {
            x1 <- mf[[i]][1]
            x2 <- mf[[i]][2]
            plot(x1, x2, xlab = names(mf[[i]])[1], ylab = names(mf[[i]])[2])
        }
        if (ncol(mf[[i]]) > 2) {
            for (j in 1:ncol(mf[[i]]))
                plot(sort(dat[[j]]), p[order(dat[[j]])], type = type, xlab = names(mf[[i]])[j],
                     ylab = colnames(pr)[i], ylim = ylim, ...)
                if (rug) rug(dat[[1]])
        }
    }
}
