
library("foreign")
library("car")

indiaraw <- read.dta("IAKR51FL.DTA")
vars <- character(0)
nams <- character(0)

# Response related variables:
# hw2      "Weight in kilograms (1 decimal)"
# hw3      "Height in centimeters (1 decimal)"
# hw70     "Ht/A Standard deviations (according to WHO)" = Stunting
# hw71     "Wt/A Standard deviations (according to WHO)" = Wasting
# hw72     "Wt/Ht Standard deviations (according to WHO)" = Underweight
# hw73     "BMI Standard deviations (according to WHO)"
# hw5      "Ht/A Standard deviations" = old definition of Stunting
# hw8      "Wt/A Standard deviations" = old definition of Wasting
# hw11     "Wt/Ht Standard deviations" = old definition of Underweight

indiaraw$hw2 <- indiaraw$hw2/10
indiaraw$hw3 <- indiaraw$hw3/10
indiaraw$cbmi <- indiaraw$hw2/((indiaraw$hw3/100)^2)

vars <- c(vars, c("hw2", "hw3", "hw70", "hw71", "hw72", "hw73", "hw5", "hw8",
                  "hw11","cbmi"))
nams <- c(nams, c("cweight", "cheight", "stunting", "wasting", "underweight",
                  "cbmiscore", "stuntingold", "wastingold", "underweightold",
                  "cbmi"))

# v008     "Date of interview (CMC)"
# b3       "Date of birth (CMC)"
# m5       "Months of breastfeeding"
# b4       "Sex of child"
# bord     "Birth order number"
# b0       "Child is twin"

indiaraw$cage <- indiaraw$v008-indiaraw$b3
indiaraw$m5[indiaraw$m5==94] <- 0
indiaraw$m5[indiaraw$m5>90] <- NA
indiaraw$ctwin <- recode(indiaraw$b0, "c('1st of multiple','2nd of multiple','3rd of multiple','4th of multiple','5th of multiple')='twin'")

vars <- c(vars, c("cage", "m5", "b4", "ctwin", "bord"))
nams <- c(nams, c("cage", "breastfeeding", "csex", "ctwin", "cbirthorder"))


# v445     "Body mass index for respondent"
# v437     "Respondent's weight (kilos-1d)"
# v438     "Respondent's height (cms-1d)"
# b3       "Date of birth (CMC)"
# v011     "Date of birth (CMC)"
# v133     "Education in single years"
# v715     "Partner's education-single years"
# v717     "Respondent's occupation"
# v130     "Religion"
# v025     "Type of place of residence"
# v206     "Sons who have died"
# v207     "Daughters who have died"

indiaraw$v445 <- indiaraw$v445/100
indiaraw$v437 <- indiaraw$v437/10
indiaraw$v438 <- indiaraw$v438/10
indiaraw$mage <- round((indiaraw$b3-indiaraw$v011)/12, 0)
indiaraw$munemployed <- as.numeric(indiaraw$v717)
indiaraw$munemployed[indiaraw$munemployed>1] <- 2
indiaraw$munemployed <- factor(indiaraw$munemployed, labels=c("unemployed","employed"))
#indiaraw$cunemployed <- recode(indiaraw$v717, "c('prof., tech., manag.','clerical','sales','agric-self employed','agric-employee','household & domestic','services','skilled & unskilled manual','[for nfhs-3 skilled and unskilled manual combined]','unskilled','don\'t know')='working'")
indiaraw$mreligion <- recode(indiaraw$v130, "c('buddhist/neo-buddhist','jain','jewish','parsi/zoroastrian','no religion','donyi polo')='other'")
indiaraw$deadchildren <- indiaraw$v206 + indiaraw$v207

vars <- c(vars, c("v445", "v437", "v438", "mage", "v133", "v715", "munemployed", "mreligion", "v025", "deadchildren"))
nams <- c(nams, c("mbmi", "mweight", "mheight", "mage", "medu", "edupartner", "munemployed", "mreligion", "mresidence", "deadchildren"))


# v190     "Wealth index"
# v119     "Has electricity"
# v120     "Has radio"
# v121     "Has television"
# v122     "Has refrigerator"
# v123     "Has bicycle"
# v124     "Has motorcycle/scooter"
# v125     "Has car"
vars <- c(vars, c("v190", "v119", "v120", "v121", "v122", "v123", "v124", "v125"))
nams <- c(nams, c("wealth", "electricity", "radio", "television", "refrigerator", "bicycle", "motorcycle", "car"))


india <- indiaraw[,vars]
names(india) <- nams
rm(indiaraw)


#######
# plausibility checks

hist(india$cweight)
india$cweight[india$cweight>25] <- NA
hist(india$cheight)
india$cheight[india$cheight>140] <- NA
hist(india$stunting)
india$stunting[india$stunting>2000] <- NA
hist(india$wasting)
india$wasting[india$wasting>2000] <- NA
hist(india$underweight)
india$underweight[india$underweight>2000] <- NA

hist(india$cbmi)
india$cbmi[india$cbmi>26] <- NA
india$cbmi[india$cbmi<10] <- NA
plot(india$cbmi, india$stunting)

hist(india$cage)
plot(india$cage, india$stunting)

table(india$breastfeeding)
plot(india$breastfeeding, india$stunting)
sum(india$breastfeeding > india$cage, na.rm=T)

table(india$csex)
plot(india$csex, india$stunting)

table(india$ctwin)
plot(india$ctwin, india$stunting)

table(india$cbirthorder)
india$cbirthorder[india$cbirthorder>5] <- 5
india$cbirthorder <- factor(india$cbirthorder)
plot(india$cbirthorder, india$stunting)

hist(india$mbmi)
india$mbmi[india$mbmi>40] <- NA
plot(india$mbmi, india$stunting)

hist(india$mweight)
india$mweight[india$mweight>120] <- NA
india$mweight[india$mweight<20] <- NA
plot(india$mweight, india$stunting)

hist(india$mheight)
india$mheight[india$mheight>185] <- NA
india$mheight[india$mheight<123] <- NA
plot(india$mheight, india$stunting)

hist(india$mage)
plot(india$mage, india$stunting)

hist(india$medu)
plot(india$medu, india$stunting)

hist(india$edupartner)
india$edupartner[india$edupartner>80] <- NA
plot(india$edupartner, india$stunting)

table(india$munemployed)
plot(india$munemployed, india$stunting)

table(india$mreligion)
plot(india$mreligion, india$stunting)

table(india$mresidence)
plot(india$mresidence, india$stunting)

table(india$deadchildren)
india$deadchildren[india$deadchildren>3] <- 3
india$deadchildren <- factor(india$deadchildren)
plot(india$deadchildren, india$stunting)

table(india$wealth)
plot(india$wealth, india$stunting)

table(india$electricity)
india$electricity[india$electricity=="not de jure resident"] <- NA
india$electricity<- factor(india$electricity, labels=c("no","yes"))
plot(india$electricity, india$stunting)

table(india$radio)
india$radio[india$radio=="not de jure resident"] <- NA
india$radio<- factor(india$radio, labels=c("no","yes"))
plot(india$radio, india$stunting)

table(india$television)
india$television[india$television=="not de jure resident"] <- NA
india$television<- factor(india$television, labels=c("no","yes"))
plot(india$television, india$stunting)

table(india$refrigerator)
india$refrigerator[india$refrigerator=="not de jure resident"] <- NA
india$refrigerator<- factor(india$refrigerator, labels=c("no","yes"))
plot(india$refrigerator, india$stunting)

table(india$bicycle)
india$bicycle[india$bicycle=="not de jure resident"] <- NA
india$bicycle<- factor(india$bicycle, labels=c("no","yes"))
plot(india$bicycle, india$stunting)

table(india$motorcycle)
india$motorcycle[india$motorcycle=="not de jure resident"] <- NA
india$motorcycle<- factor(india$motorcycle, labels=c("no","yes"))
plot(india$motorcycle, india$stunting)

table(india$car)
india$car[india$car=="not de jure resident"] <- NA
india$car<- factor(india$car, labels=c("no","yes"))
plot(india$car, india$stunting)

india <- india[complete.cases(india),]

india$intercept <- 1

center <- function(x)
  return(x-mean(x))

india$cagec <- center(india$cage)
india$magec <- center(india$mage)
india$mbmic <- center(india$mbmi)
india$meduc <- center(india$medu)
india$breastfeedingc <- center(india$breastfeeding)
india$eduparterc <- center(india$edupartner)

# create random partitions of the data set
set.seed(12345)

help <- rep(1:3, each=nrow(india)/3)
india$w1 <- sample(help)
india$w2 <- sample(help)
india$w3 <- sample(help)
india$w4 <- sample(help)
india$w5 <- sample(help)
india$w6 <- sample(help)
india$w7 <- sample(help)
india$w8 <- sample(help)
india$w9 <- sample(help)
india$w10 <- sample(help)
india$w11 <- sample(help)
india$w12 <- sample(help)
india$w13 <- sample(help)
india$w14 <- sample(help)
india$w15 <- sample(help)
india$w16 <- sample(help)
india$w17 <- sample(help)
india$w18 <- sample(help)
india$w19 <- sample(help)
india$w20 <- sample(help)
india$w21 <- sample(help)
india$w22 <- sample(help)
india$w23 <- sample(help)
india$w24 <- sample(help)
india$w25 <- sample(help)
india$w26 <- sample(help)
india$w27 <- sample(help)
india$w28 <- sample(help)
india$w29 <- sample(help)
india$w30 <- sample(help)
india$w31 <- sample(help)
india$w32 <- sample(help)
india$w33 <- sample(help)
india$w34 <- sample(help)
india$w35 <- sample(help)
india$w36 <- sample(help)
india$w37 <- sample(help)
india$w38 <- sample(help)
india$w39 <- sample(help)
india$w40 <- sample(help)
india$w41 <- sample(help)
india$w42 <- sample(help)
india$w43 <- sample(help)
india$w44 <- sample(help)
india$w45 <- sample(help)
india$w46 <- sample(help)
india$w47 <- sample(help)
india$w48 <- sample(help)
india$w49 <- sample(help)
india$w50 <- sample(help)

write.table(india, "india.raw", col.names=TRUE, row.names=FALSE, sep=" ", quote=FALSE)
save("india", file="india.Rdata")

