mboost_intern <- function(...,
    fun = c("df2lambda", "hyper_bbs", "bl_lin", "bl_lin_matrix",
            "Complete.cases", "get_index", "isMATRIX",
            "hyper_bbs", "cbs", "bsplines", "names.blg",
            "model.frame.blg")) {

    fun <- match.arg(fun)
    do.call(fun, list(...))
}
