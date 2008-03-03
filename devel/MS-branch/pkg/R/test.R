library(mboost)

NBinomial <- function()
    Family(ngradient = function(sigma=1, y, f, w = 1)
               y - (y+sigma)/(exp(f)+sigma)*exp(f),
           loss = plloss <- function(sigma=1, y, f, w=1)
             - (log(gamma(y+sigma)) - log(gamma(sigma)) - log(gamma(y+1)) +
                sigma*log(sigma) - sigma*log(exp(f)+sigma) + y*f -
                y*log(exp(f)+sigma)),
           risk = function(sigma=1, y, f, w = 1) sum(w * plloss(y=y,
               f=f, sigma=sigma)),
           check_y = function(y) TRUE,
           offset = function(y, w=NULL) {
               return(3)
           },
           name = "Negative Negative Binomial Likelihood",
           sigmaTF=T)

ctrl <-  boost_control(mstop=1000,center=T)

m1 <- glm.nb(Days ~ Eth+Lrn+Sex, quine, link = log)
m2 <- glmboost(Days ~ Eth+Lrn+Sex, data=quine, family=NBinomial(),
               control=ctrl)
               
               
               
               
################################################################################

seval <- function(S, stime.enl){

X <- c(1,S$surv)
Y <- c(S$time,10e06)
Q <- function(x) (X[which((x<=Y))])[1]
sapply(stime.enl, Q)

}

seval2 <- function(S, stime.enl){

X <- c(S$surv)
Y <- c(S$time)
Q <- function(x) {
    q <- rev(X[which((x>=Y))])[1]
    if (is.na(q)) return(1) else
        return(q)
        }
sapply(stime.enl, Q)

}


BSLogNormal <- function(t)
    Family(ngradient = function(sigma=1, y, f, w = 1) {
               Sw <- function(w){
               1 - pnorm(w)
               }
               lnt <- log(y[,1])
               eta <- (lnt-f) / sigma
               yt <- y[,1] > t
               2 * dnorm(eta)/sigma * (yt-Sw(eta))
           },

           loss = plloss <- function(sigma=1, y, f){
               Sw <- function(w){
               1 - pnorm(w)
               }
               S.stime <- y[,1]
               lnt <- log(S.stime)
               eta <- (lnt-f)/sigma
               yt <- y[,1] > t
               (yt-Sw(eta))^2
               },

           offset = function(y, w=NULL) {
               return(0)
           },

           risk = function(sigma=1, y, f, w = 1) sum(w * plloss(sigma, y, f)),

           check_y = function(y) {
               if (!inherits(y, "Surv"))
                   stop("response is not an object of class ", sQuote("Surv"))
               TRUE
           },

           name = "Negative Log Normal Likelihood",
           sigmaTF=TRUE )
           
           
sbrier <- function(sigma, y, f, t){
    W1.z <- y[,2]*(y[,1]<=t)
    W1.n <- seval(survfit(Surv(y[,1],(1-y[,2]))),y[,1])
    W1 <- W1.z/W1.n
    W2.z <- y[,1]>t
    W2.n <- seval2(survfit(Surv(y[,1],(1-y[,2]))),t)
    W2 <- W2.z/W2.n
    Sw <- function(w){
        1 - pnorm(w)
        }
    S.stime <- y[,1]
    lnt <- log(S.stime)
    eta <- (lnt-f)/sigma
    yt <- S.stime > t
    Sw(eta)^2 * W1 + (1-Sw(eta))^2 * W2
}

library(mboost)

model1 <- survreg(Surv(futime, fustat) ~ ecog.ps + rx, ovarian, dist='lognormal')
y <- Surv(ovarian$futime, ovarian$fustat)
f <- predict(model1, type="lp")
sigma <- model1$scale
mean(sbrier(sigma, y, f, t=500))
model1$coef


ctrl <- boost_control(center=T, mstop=500)

t <- 500
y <- Surv(ovarian$futime,ovarian$fustat)
W1.z <- y[,2]*(y[,1]<=t)
W1.n <- seval(survfit(Surv(y[,1],(1-y[,2]))),y[,1])
W1 <- W1.z/W1.n
W2.z <- y[,1]>t
W2.n <- seval2(survfit(Surv(y[,1],(1-y[,2]))),t)
W2 <- W2.z/W2.n
weights <- W1+W2


model2 <- glmboost(Surv(futime, fustat) ~ ecog.ps + rx, data=ovarian,
          family=BSLogNormal(500), weights=weights, control=ctrl)