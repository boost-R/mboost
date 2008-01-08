
boost_control <- function(mstop = 100, nu = 0.1, constraint = FALSE, 
                          risk = c("inbag", "oobag", "none"), 
                          savedata = TRUE, center = FALSE, trace = FALSE) {

   risk <- match.arg(risk)
   RET <- list(mstop = mstop, nu = nu, constraint = constraint,
               risk = risk, savedata = savedata, center = center,
               trace = trace)
   class(RET) <- c("boost_control")
   RET
}
