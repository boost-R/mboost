
lwd <- 2
cex <- 1
sink("tmpfile")
library("mboost")
sink()
file.remove("tmpfile")
cat("\n\n\t%%%% DON'T EDIT THIS FILE\n\n")
options(prompt = "R> ", width = 60, continue = "     ", digits = 5)
if (!file.exists("figures"))
    dir.create("figures")
cat("\\setkeys{Gin}{width = 0.95\\textwidth}")
set.seed(290875)

perfplot <- function(x, grid, alpha = NULL, border = 1,
                     xlab = "", ylab = "Performance") {
  x <- as.matrix(x)
  nc <- NCOL(x)
  nr <- NROW(x)
  nam <- grid
  if(is.null(alpha)) alpha <- 0.17 - 0.0002 * nr

  plot(rep(1:nc, 2), rep(range(x), nc), type = "n", axes = FALSE,
       ylab = ylab, xlab = xlab, xlim = c(0.6, nc + 0.4))
  axis(1, at = 1:nc, labels = nam)
  axis(2)
  box()
  matlines(t(x), col = rgb(0,0,0, alpha), lty = 1)
  lines(colMeans(x))
}

data("bodyfat", package = "TH.data")
data("wpbc", package = "TH.data")
