
source("config.R")

x <- readLines("_config.yml")
x[grep("^baseurl", x)] <- paste("baseurl: http://", pkg, ".r-forge.r-project.org", sep = "")

writeLines(x, file.path(dest, "_config.yml"))

wd <- setwd(dest)

system("jekyll build")

setwd(wd)

system(paste("cp -ra", file.path(dest, "_site/*"), publish))

