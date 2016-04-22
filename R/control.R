
boost_control <- function(mstop = 100, nu = 0.1,
                          risk = c("inbag", "oobag", "none"),
                          stopintern = FALSE,
                          center = TRUE, trace = FALSE) {

   risk <- match.arg(risk)
   stopintern <- stopintern & (risk == "oobag")
   RET <- list(mstop = mstop, nu = nu,
               risk = risk, stopintern = stopintern,
               center = center, trace = trace)
   class(RET) <- c("boost_control")
   RET
}
