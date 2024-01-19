allecho <- c("L7a", "L8a", "L8ca", "L11ca", "L12", "L13a", "L16ca", "L17ca", "L18ca") # All instructor annotated slide names
noecho <- NULL # Remove all instructor slides for student version
.showfigs <- function(x) any(x %in% noecho) # Include all instructor slides
#.showfigs <- function(x) any(x %in% allecho) # Include no instructor slides
# Period is included in front of function so that it doesn't show up in "ls()"

# load these all in for every chapter
library(curatedOvarianData)
library(ROCR)
library(coda)
library(rjags)
library(R2jags)
library(pracma)
library(survival)
library(fdrtool)
library(meta)
library(switchBox)
library(klaR)
library(survcomp)
library(MASS)
library(gam)
library(knitr)

# control sizing for all chapters
library(knitr)
knitr::opts_chunk$set(out.height = "\\textheight",  out.width = "\\textwidth")
knitr::opts_chunk$set(tidy=TRUE,tidy.opts=list(width.cutoff=70))
options(scipen = 999) # Prevents scientific notation.
