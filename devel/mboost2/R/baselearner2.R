
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

    rm(x)
    design <- function(x)
        splineDesign(knots, x, ord = degree + 1, outer.ok = TRUE) 

    design
}

### mypoly function?

Terms <- function(data, index, by = NULL, differences = 2, ...) {

     ret <- list()
     count <- 0
     if (length(differences) != length(index))
         differences <- rep(differences, length(index))
     for (i in 1:length(index)) {
         for (d in (0:differences[i])) {
             count <- count + 1
             ret[[count]] <- list(varid = index[i], 
                 spdesign = spDes(get_vars(data, index[i]), ...), 
                 which = d, differences = differences[i], by = by)
             class(ret[[count]]) <- "Terms"
         }
     }
     if (length(ret) == (differences[1] + 1)) return(ret)

     ind <- vector(mode = "list", length = length(index))
     for (i in 1:length(index))
         ind[[i]] <- 0:differences[i]
     ind <- as.matrix(do.call("expand.grid", ind))[-1,, drop = FALSE]

     ret <- list(Terms = ret, Tensor = ind)
     class(ret) <- "MTerms"
     return(ret)
}

model.matrix.Terms <- function(object, data, ...) {

    x <- get_vars(data, object$varid, classes = "numeric")
    stopifnot(length(x) == 1)
    nam <- names(x)
    x <- x[[1]]

    xu <- unique(x)
    index <- match(x, xu)

    if (object$which == object$differences) {

         Xnp <- object$spdesign(xu)
         ### FIXME: spaltenweise Mittelwert abziehen?
         D <- diff(diag(ncol(Xnp)), differences = object$differences)
         ret <- tcrossprod(Xnp, D) %*% solve(tcrossprod(D))
         attr(ret, "index") <- index
    } else {

         if (object$which == 0) {
             ret <- matrix(1, nrow = NROW(x), ncol = 1)
         } else {

             Xp <- poly(xu, degree = object$differences - 1)
             ret <- Xp[, object$which, drop = FALSE]
             attributes(ret) <- attributes(Xp)
             class(ret) <- "matrix"
             attr(ret, "index") <- index
         }
    }
    attr(ret, "Terms") <- object

    return(ret)
}

model.matrix.MTerms <- function(object, data, ...) {

    mm <- lapply(object$Terms, model.matrix, data = data, ...)

    tmp <- vector(mode = "list", length = ncol(object$Tensor))
    baseX <- vector(mode = "list", length = nrow(object$Tensor))

    for (i in 1:nrow(object$Tensor)) {
         for (j in 1:ncol(object$Tensor)) {
             tmp[[j]] <- mm[[object$Tensor[i,j] + 1]]
         }
         baseX[[i]] <- Tensor(tmp)
    } 
    c(mm, baseX)
}

### todo: poly failed
### dimensions of Tensor
