
pkg <- "mboost"
dest <- "html"
publish <- "../www"

download.file("http://user.math.uzh.ch/hothorn/TH.bib", dest = "TH.bib")
system("cat mboost.bib >> TH.bib")
bib <- "TH.bib"