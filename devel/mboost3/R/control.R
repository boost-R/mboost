
boost_control <- function(mstop = 100, nu = 0.1, 
                          risk = c("inbag", "oobag", "none"),
                          center = FALSE, trace = FALSE) {

   risk <- match.arg(risk)
   RET <- list(mstop = mstop, nu = nu, 
               risk = risk, center = center,
               trace = trace)
   class(RET) <- c("boost_control")
   RET
}
