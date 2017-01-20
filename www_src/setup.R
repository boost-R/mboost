
library("pkg2html")
library("markdown")

source("config.R")

if (!file.exists(dest))
    dir.create(dest)

stopifnot(file.exists(dest)) 

template <- system.file("template", package = "pkg2html")

system(paste("cp -ra", file.path(template, "*"), dest, sep = " "))

wd <- setwd(file.path(dest, "_data"))
R2yaml(pkg)
writeLines(bib2yaml(file.path(wd, bib), 
           c("Buehlmann:2008:StatSci",
"Hothorn+Buehlmann+Kneib+Schmid+Hofner_2010",
"Hothorn:2006:Bioinformatics:16940323",
"Hothorn:2006:Biostatistics:16344280",
"Schmid+Hothorn:2008a",
"Kneib+Hothorn+Tutz:2009",
"Hothorn+Mueller+Schroeder_2011",
"Fenske+Kneib+Hothorn_2011",
"Fenske_Fahrmeir_Hothorn_2013",
"Mayr_Hothorn_Fenske_2012",
"Schmid+Hothorn:2008b",
"Hofner+Hothorn+Kneib_2011",
"Schmid_Hothorn_Krause_2012",
"Hothorn_Kneib_Buehlmann_2014",
"Moest_Schmid_Faschingbauer_2014")),
           con = "cites.yml")

setwd(wd)
setwd(file.path(dest, "_posts"))
NEWS2md(pkg)

setwd(wd)

Rmd <- list.files(pattern = "Rmd$")

for (f in Rmd)
    writeLines(Rmd2html(f), con = file.path(dest, gsub("Rmd$", "html", f)))

x <- readLines(file.path(dest, "_data", "pkg.yml"))
x <- c(x, paste("headpic: /img/", pkg, ".png", sep = ""))
writeLines(x, con = file.path(dest, "_data", "pkg.yml"))

file.copy("_config.yml", dest, overwrite = TRUE)

yml <- list.files(pattern = "yml$")
yml <- yml[-grep("^_", yml)]
sapply(yml, function(f) file.copy(f, file.path(dest, "_data"), overwrite = TRUE))

system(paste("cat site.css >> ", file.path(dest, "css", "main.css")))

system(paste("cp ", file.path(publish, "img", "*"), file.path(dest, "img")))

