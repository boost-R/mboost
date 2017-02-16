
library("mboost")
library("party")

load("india.Rdata")

#######
# VCM #
#######

quantiles <- c(0.05, 0.1, 0.5)
inds <- paste("w", 1:50, sep="")

fmvcm <- stunting ~ bols(intercept, intercept=FALSE) + 
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
    bols(car) + 
    bols(cage, by=csex) + bbs(cage, center=TRUE, df=1, knots = 20, by=csex) +
    bols(breastfeeding, by=csex) + bbs(breastfeeding, center=TRUE, df=1, knots=20, by=csex) +
    bols(mbmi, by=csex) + bbs(mbmi, center=TRUE, df=1, knots=20, by=csex) +
    bols(mage, by=csex) + bbs(mage, center=TRUE, df=1, knots=20, by=csex) +
    bols(medu, by=csex) + bbs(medu, center=TRUE, df=1, knots=20, by=csex) +
    bols(edupartner, by=csex) + bbs(edupartner, center=TRUE, df=1, knots=20, by=csex) +
    bols(ctwin, by=csex) +
    bols(cbirthorder, by=csex) +
    bols(munemployed, by=csex) +
    bols(mreligion, by=csex) +
    bols(mresidence, by=csex) +
    bols(deadchildren, by=csex) +
    bols(wealth, by=csex) +
    bols(electricity, by=csex) +
    bols(radio, by=csex) +
    bols(television, by=csex) +
    bols(refrigerator, by=csex) +
    bols(bicycle, by=csex) +
    bols(motorcycle, by=csex) +
    bols(car, by=csex) 

for(i in 1:length(inds))
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
        
        model <- gamboost(fmvcm,
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
        
        save(model, file=paste("VCM/model",i,"_",quantiles[j],".Rdata",sep=""))
        rm(model)
        gc()
    }
}