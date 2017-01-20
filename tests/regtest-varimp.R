require("mboost")
require("lattice")

set.seed(2102)


### tests for variable importance functions
### a simple two-dimensional example
data("cars")

# gamboost with only one baselearner
gam <- gamboost(dist ~ bols(speed), 
  data = cars, control = boost_control(mstop = 50))

# risk reduction should completely be assigned to the only baselearner
stopifnot( round((gam$family@risk(y = gam$response, f = gam$offset) - 
  gam$risk()[length(gam$risk())]) / length(gam$response), 
  digits = 4 ) == round(varimp(gam), digits = 4) )
# check plots
plot(varimp(gam))
plot(varimp(gam), percent = FALSE, type = "blearner")

# glmboost with only one variable / without intercept
glm <- glmboost(dist ~ 0 + speed, data = cars, center = FALSE,
  control = boost_control(mstop = 50))
# check plot
plot(varimp(glm))

# glmboost with only one variable / with intercept
glm <- glmboost(dist ~ speed, data = cars, center = FALSE,
  control = boost_control(mstop = 50))
# check plot
plot(varimp(glm), percent = FALSE)


### examples with multiple baselearners
data("iris")
iris$target <- factor(ifelse(iris$Species == "setosa", 1, 0))
# separate intercept
iris$int <- rep(1, nrow(iris)) 

# gamboost with multiple variables
gam <- gamboost(target ~ bols(int, intercept = FALSE) +
                        bols(Sepal.Width, intercept = FALSE) +
                        bbs(Sepal.Width, center = TRUE) +
                        bols(Sepal.Length, intercept = FALSE) +
                        bbs(Sepal.Length, center = TRUE) + 
                        bbs(Petal.Width, center = TRUE) + 
                        bbs(Petal.Length, center = TRUE),
  data = iris, control = boost_control(mstop = 100), 
  family = Binomial(link = c("logit")))
# check distinct baselearners and variables
stopifnot( length(varimp(gam)) == 7 ) # 7 baselearner
stopifnot( length(unique(attr(varimp(gam), "variable_names"))) == 5 ) # 5 vars
# check plotting
plot(varimp(gam), type = "blearner")
plot(varimp(gam), type = "variable")
plot(varimp(gam), type = "variable", nbars = 3)
plot(varimp(gam), type = "variable", nbars = 1)

# glmboost with multiple variables and intercept
glm <- glmboost(target ~ Sepal.Width + Sepal.Length + Petal.Width +Petal.Length,
  data = iris, control = boost_control(mstop = 100), 
  family = Binomial(link = c("logit")))
# check plotting
plot(varimp(glm))
plot(varimp(glm), type = "blearner", percent = FALSE)


### some more tests for plotting
# show only limited number of bars, aggregate rest
plot(varimp(glm), nbars = 4)
plot(varimp(glm), nbars = 3)
plot(varimp(glm), nbars = 2)
plot(varimp(glm), nbars = 1) # other <-> 100%

# change number of characters in bar labels
plot(varimp(gam), maxchar = 3)
plot(varimp(gam), maxchar = 100)

# change xlim
plot(varimp(gam), xlim = c(0,5))
plot(varimp(gam), xlim = c(0,0.2))

# change auto.key
plot(varimp(gam), auto.key = FALSE)

# change order of baselearners
plot(varimp(gam), blorder = "alphabetical")
plot(varimp(gam), blorder = "rev_alphabetical")
plot(varimp(gam), blorder = "formula")

# with type 'blearner'
plot(varimp(gam), type = "blearner")
plot(varimp(gam), type = "blearner", blorder = "alphabetical")
plot(varimp(gam), type = "blearner", blorder = "rev_alphabetical")
plot(varimp(gam), type = "blearner", blorder = "formula")


### again multiple baselearners (with simulated data)
mydata = as.data.frame(matrix(nrow = 100, ncol = 100, data = rnorm(n = 10000)))
mydata$target = factor(rep(c(1,0), each = 50))
glm = glmboost(target ~ ., data = mydata, control = boost_control(mstop = 100), 
  family = Binomial(link = c("logit")))
# check plot
# show ten bars only (default setting)
plot(varimp(glm))
# other larger than single variables, but still displayed at the bottom
plot(varimp(glm), nbars = 5)
