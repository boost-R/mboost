mboost_intern <- function(...,
    fun = c("df2lambda", "hyper_bbs", "hyper_ols",
            "bl_lin", "bl_lin_matrix",
            "Complete.cases", "get_index", "isMATRIX",
            "cbs", "bsplines", "model.frame.blg")) {

    fun <- match.arg(fun)
    do.call(fun, list(...))
}
