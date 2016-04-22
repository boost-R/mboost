
boost_control <- function(mstop = 100, nu = 0.1,
                          risk = c("inbag", "oobag", "none"),
                          center = TRUE, trace = FALSE,
                          nupen = NULL) {

   risk <- match.arg(risk)
   RET <- list(mstop = mstop, nu = nu,
               risk = risk, center = center,
               trace = trace, nupen = nupen)
   class(RET) <- c("boost_control")
   RET
}
