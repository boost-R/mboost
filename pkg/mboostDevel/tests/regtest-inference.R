require("mboostDevel")
attach(asNamespace("mboostDevel"))

### (Slightly) modified version of the code accompanying the paper:
###   Shah, R. D. and Samworth, R. J. (2013), Variable selection with error
###   control: Another look at Stability Selection, J. Roy. Statist. Soc., Ser.
###   B, 75, 55-80. DOI: 10.1111/j.1467-9868.2011.01034.x
###
### Original code available from
###   http://www.statslab.cam.ac.uk/~rds37/papers/r_concave_tail.R
### or
###   http://www.statslab.cam.ac.uk/~rjs57/r_concave_tail.R
r.TailProbs <- function(eta, B, r) {
    ## If pi = ceil{ B * 2 * eta} / B + 1/B,..., 1 return the tail probability.
    ## If pi < ceil{ B * 2 * eta} / B return 1

    MAXa <- 100000
    MINa <- 0.00001

    s <- -1/r
    etaB <- eta * B
    k_start <- (ceiling(2 * etaB) + 1)
    output <- rep(1, B)
    if (k_start > B)
        return(output)

    a_vec <- rep(MAXa,B)

    Find.a <- function(prev_a)
        uniroot(Calc.a, lower = MINa, upper = prev_a,
                tol = .Machine$double.eps^0.75)$root

    Calc.a <- function(a) {
        denom <- sum((a + 0:k)^(-s))
        num <- sum((0:k) * (a + 0:k)^(-s))
        num / denom - etaB
    }

    for(k in k_start:B)
        a_vec[k] <- Find.a(a_vec[k-1])

    # NB this function makes use of several gloabl variables
    OptimInt <- function(a) {
        num <- (k + 1 - etaB) * sum((a + 0:(t-1))^(-s))
        denom <- sum((k + 1 - (0:k)) * (a + 0:k)^(-s))
        1 - num / denom
    }

    prev_k <- k_start
    for(t in k_start:B) {
        cur_optim <- rep(0, B)
        cur_optim[B] <- OptimInt(a_vec[B])
        if (prev_k <= (B-1)) {
            for (k in prev_k:(B-1))
                cur_optim[k] <- optimize(f=OptimInt, lower = a_vec[k+1],
                                         upper = a_vec[k], maximum  = TRUE)$objective
        }
        output[t] <- max(cur_optim)
        prev_k <- which.max(cur_optim)
    }
    return(output)
}

pminD <- function(theta, B, r = c(-1/2, -1/4)) {
    pmin(c(rep(1, B), r.TailProbs(theta^2, B, r[1])),
         r.TailProbs(theta, 2*B, r[2]))
}

## test r-concave bound
B <- 50
x <- (1:(2 * B))/(2 * B)
p <- 1000
q <- 50
theta <- q/p
if (FALSE) {
    plot(x, log(pminD(theta, B)), xlab = "pi")
    abline(v = ceiling(2 * theta * 100) / 100)
    Ds <- cbind(c(rep(1, B), r.TailProbs(theta^2, B, -1/2)),
                r.TailProbs(theta, 2*B, -1/4))
    round(log(Ds), 2)
    lines(x, log(Ds[,1]), col = "red", lwd = 2)
    lines(x, log(Ds[,2]), col = "blue", lwd = 2)
}

## r-concave bound of Shah & Samworth (2013)
bound_ss <- (pminD(theta, B) * p)[40:100]
plot(x[40:100], bound_ss, xlab = "pi", ylim = c(0, 50))
## Bound of Meinshausen & Buehlmann (2010)
points(x[40:100], q^2 / (2 * x[40:100] - 1) / p, col = "red")
## now our implementation
bound <- rep(NA, 61)
for (i in 40:100) {
    bound[i - 39] <- minD(q, p, i/100, B) * p
}
points((40:100)/100, bound, col = "green")
stopifnot(all((bound - bound_ss) < sqrt(.Machine$double.eps)))

## test r-concave bound
B <- 50
x <- (1:(2 * B))/(2 * B)
p <- 1000
q <- 490
theta <- q/p

## r-concave bound of Shah & Samworth (2013)
bound_ss <- (pminD(theta, B) * p)[40:100]
plot(x[40:100], bound_ss, xlab = "pi")
## Bound of Meinshausen & Buehlmann (2010)
points(x[40:100], q^2 / (2 * x[40:100] - 1) / p, col = "red")
## now our implementation
bound <- rep(NA, 61)
for (i in 40:100) {
    bound[i - 39] <- minD(q, p, i/100, B) * p
}
points((40:100)/100, bound, col = "green")
stopifnot(all((bound - bound_ss) < sqrt(.Machine$double.eps)))

### computation of q from other values
cutoff <- 0.6
PFER <- 0.2
B <- 50
p <- 200
(q <- optimal_q(p = p, cutoff = cutoff, PFER = PFER, B = B))
# check:
round(minD(q, p, cutoff, B) * p, 3)
round(minD(q + 1, p, cutoff, B) * p, 3)


### computation of cutoff from other values
PFER <- 0.2
B <- 50
p <- 200
q <- 7
(cutoff <- optimal_cutoff(p = p, q = q, PFER = PFER, B = B))
# check:
round(minD(q, p, cutoff, B) * p, 3)
round(minD(q, p, cutoff - 1e-2, B) * p, 3)


### check stabsel interface
data("bodyfat", package = "TH.data")
mod <- glmboost(DEXfat ~ ., data = bodyfat)
(sbody <- stabsel(mod, q = 3, PFER = 0.2))
dim(sbody$phat)
(sbody <- stabsel(mod, q = 3, PFER = 0.2, error.bound = "SS"))
dim(sbody$phat)


## check stabsel_parameters and (theoretical) error control
cutoff <- 0.6
for (i in 1:10) {
    print(stabsel_parameters(cutoff = cutoff, q = i, p = 100, error.bound = "MB"))
}
for (i in 1:10) {
    print(stabsel_parameters(cutoff = cutoff, q = i, p = 100, error.bound = "SS"))
}

## check if missing values are determined correctly (especially at the extreme values)
p <- 100
B <- 50
cutoff <- 0.6
# low PFER
PFER <- 0.001
(res <- stabsel_parameters(p = p, cutoff = cutoff, PFER = PFER, B = B, error.bound = "SS"))
stabsel_parameters(p = p, cutoff = cutoff, q = res$q, B = B, error.bound = "SS")
# high PFER
PFER <- 50
(res <- stabsel_parameters(p = p, cutoff = cutoff, PFER = PFER, B = B, error.bound = "SS"))
stabsel_parameters(p = p, cutoff = cutoff, q = res$q, B = B, error.bound = "SS")
# medium PFER
PFER <- 1
(res <- stabsel_parameters(p = p, cutoff = cutoff, PFER = PFER, B = B, error.bound = "SS"))
stabsel_parameters(p = p, cutoff = cutoff, q = res$q, B = B, error.bound = "SS")
stabsel_parameters(p = p, cutoff = cutoff, q = res$q + 1, B = B, error.bound = "SS")


q <- 10
# high PFER
PFER <- 5
(res <- stabsel_parameters(p = p, q = q, PFER = PFER, B = B, error.bound = "SS"))
stabsel_parameters(p = p, cutoff = res$cutoff, q = q, B = B, error.bound = "SS")
# low PFER
PFER <- 0.001
(res <- stabsel_parameters(p = p, q = q, PFER = PFER, B = B, error.bound = "SS"))
stabsel_parameters(p = p, cutoff = res$cutoff, q = q, B = B, error.bound = "SS")
# medium PFER
PFER <- 1
(res <- stabsel_parameters(p = p, q = q, PFER = PFER, B = B, error.bound = "SS"))
stabsel_parameters(p = p, cutoff = res$cutoff, q = q, B = B, error.bound = "SS")
stabsel_parameters(p = p, cutoff = res$cutoff - 0.01, q = q, B = B, error.bound = "SS")
