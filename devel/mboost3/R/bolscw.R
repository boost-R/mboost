
bolscw <- function(..., z = NULL) {

    mf <- list(...)
    if (length(mf) == 1 && (is.matrix(mf[[1]]) || is.data.frame(mf[[1]]))) {
        mf <- mf[[1]]
    } else {
        mf <- as.data.frame(mf)
        cl <- as.list(match.call(expand.dots = FALSE))[2][[1]]
        colnames(mf) <- sapply(cl, function(x) as.character(x))
    }
    vary <- ""
    if (!is.null(z)) {
        mf <- cbind(mf, z)
        colnames(mf) <- c(colnames(mf), deparse(substitute(z)))
        vary <- colnames(mf)[ncol(mf)]
    }

    ret <- list(model.frame = function() return(mf),
                get_names = function() colnames(mf),
                set_names = function(value) attr(mf, "names") <<- value)
    class(ret) <- "blg"

    newX <- function(newdata) {
        if (is.matrix(newdata)) return(mf)
        X <- model.matrix(~ ., data = newdata[, colnames(mf) != vary, drop = FALSE])
        if (vary != "")
            X <- X * newdata[[vary]]
        return(X)
    }

    X <- newX(mf)

    cfM <- Matrix(0, nrow = ncol(X), ncol = 1)

    ret$dpp <- function(weights) {

        xw <- t(X * weights)
        xtx <- colSums(X^2 * weights)
        sxtx <- sqrt(xtx)
        MPinvS <- (1 / sxtx) * xw

        fit <- function(y) {
            xselect <- which.max(abs(mu <- MPinvS %*% y))  
            coef <- mu[xselect] / sxtx[xselect]
            cfM[xselect,] <- coef
            ret <- list(model = as.matrix(cfM),
                        fitted = function() 
                            coef * X[,xselect])
            class(ret) <- c("bm_lin", "bm")
            return(ret)
        }

        predict <- function(bm, newdata = NULL, aggregate = c("sum", "cumsum", "none")) {

            aggregate <- match.arg(aggregate)
            cf <- switch(aggregate, "sum" = {
                cf <- 0
                for (i in 1:length(bm))
                    cf <- cf + coef(bm[[i]])
                cf
            },
            "cumsum" = {
                cf <- coef(bm[[1]])
                for (i in 2:length(bm))
                    cf <- cbind(cf, cf[[i-1]] + coef(bm[[i]]))
                cf
            },
            "none" = {
                cf <- coef(bm[[1]])
                for (i in 2:length(bm))
                    cf <- cbind(cf, coef(bm[[i]]))
                cf
            })

            if (!is.null(newdata))
                X <- newX(newdata)
            return(X %*% cf)
        }
    
        ret <- list(fit = fit, predict = predict)
        class(ret) <- c("bl_lin", "bl")
        return(ret)
    }
    return(ret)
}

