
### the classical tree-based baselearner; stumps by default
### (also fits an additive model)
btree <- function(...,
    tree_controls = party::ctree_control(stump = TRUE,
                                  mincriterion = 0,
                                  savesplitstats = FALSE)) {

    if (!requireNamespace("party"))
        stop("cannot load ", sQuote("party"))

    cll <- match.call()
    cll[[1]] <- as.name("btree")

    ctrl <- tree_controls
    mf <- list(...)
    if (length(mf) == 1 && is.data.frame(mf[[1]])) {
        mf <- mf[[1]]
    } else {
        mf <- as.data.frame(mf)
        cl <- as.list(match.call(expand.dots = FALSE))[2][[1]]
        colnames(mf) <- sapply(cl, function(x) as.character(x))
    }

    ret <- list(model.frame = function() return(mf),
                get_call = function(){
                    cll <- deparse(cll, width.cutoff=500L)
                    if (length(cll) > 1)
                        cll <- paste(cll, collapse="")
                    cll
                },
                get_names = function() colnames(mf),
                set_names = function(value) {
                    if(length(value) != length(colnames(mf)))
                        stop(sQuote("value"), " must have same length as ",
                             sQuote("colnames(mf)"))
                    for (i in 1:length(value)){
                        cll[[i+1]] <<- as.name(value[i])
                    }
                    attr(mf, "names") <<- value
                })
    class(ret) <- "blg"


    ret$dpp <- function(weights) {

        ### construct design matrix etc.
        y <- vector(length = nrow(mf), mode = "numeric")
        ### name for working response (different from any x)
        rname <- paste(sample(LETTERS, 25, replace = TRUE), collapse = "")
        fm <- as.formula(paste(rname, " ~ ", paste(colnames(mf), collapse = "+")))
        df <- mf
        df[[rname]] <- y
        object <- party_intern(fm, data = df, fun = "ctreedpp")
        fitmem <- party::ctree_memory(object, TRUE)
        where <- rep.int(0, nrow(mf))
        storage.mode(where) <- "integer"
        storage.mode(weights) <- "double"

        fitfun <- function(y) {

            party_intern(y, object@responses, fun = "R_modify_response")
            tree <- party_intern(object, weights, fitmem, ctrl, where,
                                 fun = "R_TreeGrow")
            party_intern(tree, TRUE, fun = "R_remove_weights")

            fitted <- function() {
                wh <- party_intern(tree, object@inputs, 0.0,
                                   fun = "R_get_nodeID")
                return(unlist(party_intern(tree, wh, fun = "R_getpredictions")))
            }

            ret <- list(model = tree, fitted = fitted)
            class(ret) <- c("bm_tree", "bm")
            ret
        }

        predict <- function(bm, newdata = NULL, aggregate = c("sum", "cumsum", "none")) {
            aggregate <- match.arg(aggregate)

            if (is.null(newdata)) {
                newinp <- object@inputs
            } else {
                newinp <- party_intern(object, newdata, fun = "newinputs")
            }

            pr <- 0
            for (i in 1:length(bm)) {
                wh <- party_intern(bm[[i]]$model, newinp, 0.0,
                                   fun = "R_get_nodeID")
                pri <- unlist(party_intern(bm[[i]]$model, wh,
                                           fun = "R_getpredictions"))
                if (aggregate == "sum") {
                    pr <- pr + pri
                } else {
                    if (i > 1) {
                        pr <- cbind(pr, pri)
                    } else {
                        pr <- pri
                    }
                    if (aggregate == "cumsum")
                        if (i > 1) pr[,i] <- pr[,i] + pr[,i-1]
                }
            }
            return(pr)
        }

        ret <- list(fit = fitfun, predict = predict)
        class(ret) <- c("bl_tree", "bl")
        ret
    }
    return(ret)
}
