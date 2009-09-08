
### Linear baselearner, potentially Ridge-penalized
bcobs <- function(..., index = NULL, lambda = 1, knots = 20, 
                  constraint = c("none", "increase", "decrease",
                                 "convex", "concave", "periodic")) {

    cll <- match.call()
    mf <- list(...)
    if (is.data.frame(mf[[1]])) {
        mf <- mf[[1]]
    } else {
        mf <- as.data.frame(mf)
        cl <- as.list(match.call(expand.dots = FALSE))[2][[1]]
        colnames(mf) <- sapply(cl, function(x) as.character(x))
    }

    x <- mf[[1]]
    boundary.knots <- range(x, na.rm = TRUE)
    if (length(knots) == 1) {
        knots <- seq(from = boundary.knots[1],
                     to = boundary.knots[2], length = knots + 2)
    }

    CC <- all(Complete.cases(mf))
    ### option
    DOINDEX <- is.data.frame(mf) && (nrow(mf) > 10000 || is.factor(mf[[1]]))
    if (is.null(index)) {
        ### try to remove duplicated observations or 
        ### observations with missings
        if (!CC || DOINDEX) {
            index <- get_index(mf)
            mf <- mf[index[[1]],,drop = FALSE]
            index <- index[[2]]
        }
    }

    ret <- list(model.frame = function() 
                    if (is.null(index)) return(mf) else return(mf[index,,drop = FALSE]),
                get_call = function() cll,
                get_data = function() mf,
                get_index = function() index,
                get_names = function() colnames(mf),
                get_vary = function() vary,
                set_names = function(value) attr(mf, "names") <<- value)
    class(ret) <- "blg"

    ret$dpp <- function(weights) {

        weights[!Complete.cases(mf)] <- 0
        w <- weights
        if (!is.null(index)) w <- as.vector(tapply(weights, index, sum))

        fit <- function(y) {
            if (!is.null(index)) {
                y <- as.vector(tapply(weights * y, index, sum))
            } else {
                y <- y * weights
            }
            model <- cobs(x = mf[[1]], y = y, w = w, constraint = constraint)
            ret <- list(model = model, 
                        fitted = function() {
                            ret <- fitted(model)
                            if (is.null(index)) return(ret)
                            return(ret[index])
                        })
            class(ret) <- "bm"
            ret
        }

        ### prepare for computing predictions
        predict <- function(bm, newdata = NULL, aggregate = c("sum", "cumsum", "none")) {
            pr <- sapply(bm, function(mod) {
                if (is.null(newdata)) return(mod$fitted())
                cobs:::predict.cobs(mod$model, z = newdata[[colnames(mf)]])[,2]
            })
            aggregate <- match.arg(aggregate)
            pr <- switch(aggregate, "sum" = 
                matrix(rowSums(pr), ncol = 1),
            "cumsum" = {
                M <- triu(crossprod(Matrix(1, nc = ncol(cf))))
                as(pr %*% M, "matrix")
            },
            "none" = pr)
            if (is.null(index)) return(pr[,,drop = FALSE])
            return(pr[index,,drop = FALSE])
        }

        ret <- list(fit = fit,
                    predict = predict)
        class(ret) <- "bl"
        return(ret)

    }
    return(ret)
}

library("mboost3")
library("cobs")
attach(asNamespace("mboost3"))

data("bodyfat", package = "mboost3")

bb <- function(...) bcobs(..., constraint = "increase")
bmod <- mboost(DEXfat ~ . , data = bodyfat, baselearner = bb,
            control = boost_control(mstop = 50))
pdf("bodyfat_constr.pdf")
layout(matrix(1:9, nc = 3))
plot(bmod)
dev.off()

pdf("bodyfat_unconstr.pdf")
mod <- gamboost(DEXfat ~ ., data = bodyfat, control = boost_control(mstop = 50))
layout(matrix(1:9, nc = 3))
plot(mod)
dev.off()

