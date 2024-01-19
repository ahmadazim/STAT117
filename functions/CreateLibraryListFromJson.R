library("rjson")
dat <- fromJSON(paste(readLines("~/Dropbox/GP_Teaching/STA117/Textbook117/renv.lock"), collapse=""))
cat(capture.output(dput(names(which(sapply(dat$Packages, function(x) x$Source == "Bioconductor") == T)))))
cat(capture.output(dput(names(which(sapply(dat$Packages, function(x) x$Source == "Repository") == T)))))
# creates directly callable vectors. 