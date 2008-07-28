
library("splines")

get_vars <- function(data, index, classes = c("numeric", "factor")) {

    x <- data[, index, drop = FALSE]
    stopifnot(all(sapply(x, class) %in% classes))
    x
}

spDes <- function(x, nknots = 20, degree = 3) {

    offset <- diff(range(x, na.rm = TRUE)) / 10000
    diffs <- (diff(range(x, na.rm = TRUE)) + 2 * offset) / (nknots - 1)
    knots <- seq(from = min(x, na.rm = TRUE) - offset - degree * diffs, 
                 to = max(x, na.rm = TRUE) + offset + degree * diffs,
                 by = diffs)

    design <- function(x)
        splineDesign(knots, x, ord = degree + 1, outer.ok = TRUE) 

    design
}

bl_pLS <- function(data, index, by = NULL, weights = rep.int(1, nrow(data)), 
                   setup = spDes, df = 4, ...) {

    x <- get_vars(data, c(index, by), classes = "numeric")

    xu <- unique(x)
    inx <- match(apply(x, 1, function(x) paste(x, collapse = "\r")), 
          apply(xu, 1, function(x) paste(x, collapse = "\r")))
    wu <- as.vector(tapply(weights, inx, sum))

    if (!is.null(by)) {
        byx <- xu[, ncol(xu)]
        xu <- xu[, -ncol(xu)]
    }

    tmp <- spDes(xu, ...)

    Xu <- tmp$design(xu)
    K <- tmp$pen()
    K <- K * mboost:::df2lambda(Xu, df = df, dmat = K, weights = wu)

    Xw <- Xu * wu
    XtX <- crossprod(Xw, Xu)
    Xsolve <- solve(XtX + K, t(Xu))
    
    fit <- function(y) {
        yu <- tapply(y * weights, inx, sum)
        as.vector(Xsolve %*% yu)
    }
    list(fit = fit, design = tmp$design)
}




npp <- function(data, index, setup = spDes, differences = 2, ...) {

    x <- get_vars(data, index, classes = "numeric")
    stopifnot(length(x) == 1)
    x <- x[[1]]

    design <- setup(x, ...)
    Xnp <- design(x)
    ### FIXME: spaltenweise Mittelwert abziehen?
    D <- diff(diag(ncol(Xnp)), differences = differences)
    Xnp <- tcrossprod(Xnp, D) %*% solve(tcrossprod(D))

    Xp <- cbind(1, poly(x, degree = differences - 1))

    return(list(Xp = Xp, Xnp = Xnp))
}

foo <- function(data, index, ...) {

     x <- get_vars(data, index, classes = "numeric")
     xu <- unique(x)
     inx <- match(apply(x, 1, function(x) paste(x, collapse = "\r")),
                  apply(xu, 1, function(x) paste(x, collapse = "\r")))

     dmat <- vector(mode = "list", length = length(x))
     for (i in 1:length(x)) 
         dmat[[i]] <- npp(data, i, ...)

     d <- ncol(dmat[[1]]$Xp) + 1
     ind <- vector(mode = "list", length = ncol(xu))
     for (i in 1:length(x))
         ind[[i]] <- 1:d
     ind <- as.matrix(do.call("expand.grid", ind))[-1,, drop = FALSE]
     tmp <- vector(mode = "list", length = ncol(ind))
     baseX <- vector(mode = "list", length = nrow(ind))
     if (ncol(ind) == 1) {
         for (i in 1:(nrow(ind) - 1))
             baseX[[i]] <- dmat[[1]]$Xp[, ind[i,], drop = FALSE]
         baseX[[nrow(ind)]] <- dmat[[1]]$Xnp
     } else {
         for (i in 1:nrow(ind)) {
             for (j in 1:ncol(ind)) {
                 if (ind[i,j] < d)
                     tmp[[j]] <- dmat[[j]]$Xp[, ind[i,j], drop = FALSE]
                 else
                     tmp[[j]] <- dmat[[j]]$Xnp
             }
             baseX[[i]] <- tensor.prod.model.matrix(tmp)
         }
     }
     baseX
}
