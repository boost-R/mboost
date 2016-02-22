require("mboostDevel")

set.seed(2102)

### a simple two-dimensional example
data("cars")

# gamboost with only one baselearner
gam <- gamboost(dist ~ bols(speed), 
  data = cars, control = boost_control(mstop = 50))

# all risk reduction should be assigned to the only baselearner
stopifnot( round((gam$family@risk(y = gam$response, f = gam$offset) - 
  gam$risk()[length(gam$risk())]) / length(gam$response), 
  digits = 4 ) == round(varimp_mboost(gam), digits = 4) )
# in addition check percentage
stopifnot( varimp_mboost(gam, percent = TRUE) == 1 )
# check plots
plot(varimp_mboost(gam))
plot(varimp_mboost(gam, percent = TRUE))


# glmboost with only one variable / without intercept
glm <- glmboost(dist ~ 0 + speed, data = cars, center = FALSE,
  control = boost_control(mstop = 50))
stopifnot( varimp_mboost(glm, percent = TRUE) == 1 )
# check plots
plot(varimp_mboost(glm))
plot(varimp_mboost(glm, percent = TRUE))


### gamboost with multiple baselearners
data("iris")
iris$target <- factor(ifelse(iris$Species == "setosa", 1, 0))
iris$int <- rep(1, nrow(iris)) 

gam = gamboost(target ~ bols(int, intercept = FALSE) +
                        bbs(Sepal.Width, center = TRUE) +
                        bbs(Sepal.Length, center = TRUE) + 
                        bbs(Petal.Width, center = TRUE) + 
                        bbs(Petal.Length, center = TRUE),
  data = iris, control = boost_control(mstop = 100), 
  family = Binomial(link = c("logit")))
# relative risk reduction adds up to 1
stopifnot( sum(varimp_mboost(gam, percent = TRUE)) == 1 )

# check plotting (baselearners always sorted top-down)
plot(varimp_mboost(gam))


### glmboost with multiple variables and intercept
glm <- glmboost(target ~ Sepal.Width + Sepal.Length + Petal.Width + Petal.Length,
  data = iris, control = boost_control(mstop = 50), 
  family = Binomial(link = c("logit")))
# relative risk reduction adds up to 1
stopifnot( sum(varimp_mboost(glm, percent = TRUE)) == 1 )

# check plotting (baselearners always sorted top-down)
plot(varimp_mboost(glm))
plot(varimp_mboost(glm, percent = TRUE))


### some more tests for plotting

# show only limited number of bars, aggregate rest
plot(varimp_mboost(glm, percent = TRUE), nbars = 4)
plot(varimp_mboost(glm, percent = TRUE), nbars = 3)
plot(varimp_mboost(glm, percent = TRUE), nbars = 2)
plot(varimp_mboost(glm, percent = TRUE), nbars = 1) # other <-> 100%

# change number of characters in bar labels
plot(varimp_mboost(gam, percent = TRUE), maxchar = 3)
plot(varimp_mboost(gam, percent = TRUE), maxchar = 100)

# change xlim
plot(varimp_mboost(gam), xlim = c(0,5))


### other larger than single variables, top-down order is overruled
mydata = as.data.frame( matrix(nrow = 100, ncol = 100, data = rnorm(n = 10000)) )
mydata$target = factor(rep(c(1,0), each = 50))
glm = glmboost(target ~ ., data = mydata, control = boost_control(mstop = 100), 
  family = Binomial(link = c("logit")))

plot(varimp_mboost(glm, percent = TRUE), nbars = 5)
