
library("mboost")
library("party")

load("india.Rdata")

##########
# Stumps #
##########

quantiles <- c(0.05, 0.1, 0.5)
inds <- paste("w", 1:50, sep="")

features <- names(india)[c(11:16, 19:33)]
fmtree <- paste(features, collapse = " + ")
fmtree <- as.formula(paste("stunting", fmtree, sep=" ~ "))

for(i in 1:length(inds))
{
    indiafit <- india[india[,inds[i]]<3,]
    indiapred <- india[india[,inds[i]]==3,]
    wfit <- as.numeric(indiafit[,inds[i]]==1)
    
    for(j in 1:length(quantiles))
    {
        it <- 1000
        inc <- 1000
        maxit <- 50000
        bc <- boost_control(mstop = it, nu=0.25, trace = TRUE, risk = "oob")
        tc <- ctree_control(maxdepth = 1, savesplitstats = FALSE)
        
        model <- blackboost(fmtree, 
                            data = indiafit, 
                            control = bc,
                            tree_controls = tc,
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
        
        pred <- predict(model[which.min(risk)], newdata=indiapred)
        save(pred, file=paste("stumps/model",i,"_",quantiles[j],".Rdata",sep=""))
    }
}
