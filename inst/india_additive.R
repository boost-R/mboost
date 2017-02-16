
library("mboost")
library("party")

load("india.Rdata")

##################
# Additive Model #
##################

quantiles <- c(0.05, 0.1, 0.5)

inds <- paste("w", 1:50, sep="")

fmadd <- stunting ~ bols(intercept, intercept=FALSE) + 
    bols(cage) + bbs(cage, center=TRUE, df=1, knots = 20) +
    bols(breastfeeding) + bbs(breastfeeding, center=TRUE, df=1, knots=20) +
    bols(mbmi) + bbs(mbmi, center=TRUE, df=1, knots=20) +
    bols(mage) + bbs(mage, center=TRUE, df=1, knots=20) +
    bols(medu) + bbs(medu, center=TRUE, df=1, knots=20) +
    bols(edupartner) + bbs(edupartner, center=TRUE, df=1, knots=20) +
    bols(csex) +
    bols(ctwin) +
    bols(cbirthorder) +
    bols(munemployed) +
    bols(mreligion) +
    bols(mresidence) +
    bols(deadchildren) +
    bols(wealth) +
    bols(electricity) +
    bols(radio) +
    bols(television) +
    bols(refrigerator) +
    bols(bicycle) +
    bols(motorcycle) +
    bols(car)

for(i in 1:50)
{
    indiafit <- india[india[,inds[i]]<3,]
    indiapred <- india[india[,inds[i]]==3,]
    wfit <- as.numeric(indiafit[,inds[i]]==1)
    
    for(j in 1:length(quantiles))
    {
        it <- 10000
        inc <- 10000
        maxit <- 200000
        bc <- boost_control(mstop = it, nu=0.25, trace = TRUE, risk = "oob")
        
        model <- gamboost(fmadd,
                          data = indiafit,
                          control = bc,
                          weights = wfit,
                          family = QuantReg(tau = quantiles[j])
        )
        
        risk <- model$risk()                    
        while( ((risk[it-inc+1] - risk[it])>=0.05) &
               (it<=maxit) ) 
        {
            it <- it + inc
            model[it]
            risk <- model$risk()                    
        }
        
        save(model, file=paste("additive/model",i,"_",quantiles[j],".Rdata",sep=""))
        rm(model)
        gc()
    }
}

