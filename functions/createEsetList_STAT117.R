createEsetList <- function(combat=F, rescale=F, affyonly=T, savedata=T, verbose = F){

  ## -----------------------------------------------------------------------------
  ## PARAMETER DESCRIPTIONS:
  ##
  ## combat   - combat = T runs the ComBat function
  ## rescale  - rescale = T rescales the data 
  ## affyonly - affyonly = T only returns Affymetrix datasets
  ## savedata - savedata = T saves esets_STAT117.RData file in working directory
  ##                         or esets_affy_STAT117.RData file if affyonly = T
  ## verbose  - verbose = F supresses all messages
  ## 
  ## OUTPUT:
  ## createEsetList() returns an eset list and/or .RData file.
  ## -----------------------------------------------------------------------------

  ## -----------------------------------------------------------------------------
  ## load libraries needed
  ## -----------------------------------------------------------------------------
  
  library(genefilter)
  library(survival)
  library(logging)
  library("curatedOvarianData", character.only=TRUE)
  # Contains default settings - source in
  source("functions/patientselection_STAT117.config")
  
  if(verbose == T){
    kLogFile <-"createEsetList.log"
    basicConfig()
    addHandler(writeToFile, logger="", file=kLogFile)
    
    loginfo("Inside script createEsetList_STAT117.R - inputArgs =")
    loginfo(paste("Loading", "curatedOvarianData", sessionInfo()$otherPkgs[["curatedOvarianData"]]$Version))
  }
  
  ## -----------------------------------------------------------------------------
  ## needed functions
  ## -----------------------------------------------------------------------------
  filterQuantile <- function(object, q){
    if (!identical(q >=0 && q < 1, TRUE))
      stop("require 0 <= q < 1")
    if (!identical(class(object) == "ExpressionSet", TRUE))
      stop("object must be an ExpressionSet")
    gene.sd <- esApply(object,1,sd, na.rm=TRUE)
    gene.quantile <- quantile(gene.sd, probs=q)
    actual.makescutoff <- sum(gene.sd < gene.quantile) / length(gene.sd)
    ##make sure the correct number of genes are getting filtered:
    if (abs(q - actual.makescutoff) > 0.01){
      stop("Not scaling this object, likely pre-scaled.")
    }else{
      object <- object[gene.sd > gene.quantile, ]
    }
    return(object)
  }
  ##recursive intersect function
  intersectMany <- function(lst){
    ## Find the intersection of multiple vectors stored as elements of a
    ## list, through a tail-recursive function.
    if (length(lst)==2){
      return(intersect(lst[[1]],lst[[2]]))
    }else{
      return(intersectMany(c(list(intersect(lst[[1]],lst[[2]])),lst[-1:-2])))
    }
  }
  
  ##Split out non-specific probe sets
  expandProbesets <- function (eset, sep = "///"){
    x <- lapply(featureNames(eset), function(x) strsplit(x, sep)[[1]])
    eset <- eset[order(sapply(x, length)), ]
    x <- lapply(featureNames(eset), function(x) strsplit(x, sep)[[1]])
    idx <- unlist(sapply(1:length(x), function(i) rep(i, length(x[[i]]))))
    xx <- !duplicated(unlist(x))
    idx <- idx[xx]
    x <- unlist(x)[xx]
    eset <- eset[idx, ]
    featureNames(eset) <- x
    eset
  }
  
  ## -----------------------------------------------------------------------------
  ##load the esets
  ## -----------------------------------------------------------------------------
  #data(list=data(package="curatedOvarianData")[[3]][,3])
  # For speed, saved to data folder
  load("data/allesets.RData")
  strEsets <- ls(pattern="^.*_eset$")
  esets <- list()
  
  ## -----------------------------------------------------------------------------
  ##Explicit removal of samples from specified datasets:
  ## -----------------------------------------------------------------------------
  
  # Reduces strEsets to the 17 Affymetrix Studies if affyonly = T
  if(affyonly){
    affylist <- c("GSE14764_eset","GSE18520_eset","GSE19829.GPL570_eset","GSE19829.GPL8300_eset","GSE20565_eset","GSE2109_eset","GSE26193_eset","GSE26712_eset","GSE30161_eset","GSE44104_eset","GSE6008_eset","GSE6822_eset","GSE9891_eset","PMID15897565_eset","PMID17290060_eset","PMID19318476_eset","TCGA_eset")
    affylist_index <- which(strEsets %in% affylist)
    strEsets <- strEsets[affylist_index]
  }
  delim <- ":"   ##This is the delimiter used to specify dataset:sample,
  ##e.g. TCGA_eset:TCGA.24.1927
  # remove.samples INCORPORATED BY DEFAULT
  if(exists("remove.samples")){
    ##split into those with the delimiter and those without:
    remove.samples.delim <- grep(delim, remove.samples, fixed=TRUE, value=TRUE)
    ##over-write remove.samples for later, keeping this variable in
    ##its historical use:
    remove.samples.orig <- remove.samples
    remove.samples.nodataset <- grep(delim, remove.samples, fixed=TRUE, invert=TRUE, value=TRUE)
    if(length(remove.samples.delim) > 0){
      datasets <- gsub(paste(delim, ".+", sep=""), "", remove.samples.delim)
      samples <- gsub(paste(".+", delim, sep=""), "", remove.samples.delim)
      remove.samples.delim <- lapply(unique(datasets), function(ds){
        samples[datasets %in% ds]
      })
      names(remove.samples.delim) <- unique(datasets)
    }
  }

  if(verbose == T){loginfo("Clean up the esets.")}
  for (strEset in strEsets){
    if(verbose == T){print(strEset)}
    eset <- get(strEset)
    ##Deal with genes which had a single probe mapping to multiple genes
    # DROP BY DEFAULT
    if (exists("probes.not.mapped.uniquely")){
      if(identical(probes.not.mapped.uniquely, "drop")){
        ##Drop rows without unique gene name
        eset <- eset[!grepl("///",featureNames(eset),fixed=TRUE),]
      }else if (identical(probes.not.mapped.uniquely, "split")){
        ##Split out rows without unique gene name
        eset <- expandProbesets(eset)
      }
    }
    ## Run ComBat
    # COMBAT IS NOT RUN BY DEFAULT
    if (exists("combat") && combat) {
      # workaround bug #12
      if (identical(pubMedIds(eset), "21720365")) {
        eset$batch[match("TCGA.23.1023", sampleNames(eset))] <- 12
      }
      if (sum(is.na(eset$batch)) == 0) {
        tmp <- NA # Prevents tmp from being used in next iteration
        if (verbose == F){
          suppressMessages(invisible(capture.output(try(tmp <- sva::ComBat(exprs(eset),mod=model.matrix(~rep(1, ncol(eset))),batch=eset$batch), silent=TRUE))))
          #try(tmp <- sva::ComBat(exprs(eset),mod=model.matrix(~rep(1, ncol(eset))),batch=eset$batch), silent=TRUE)
        } else {
          tmp <- try(sva::ComBat(exprs(eset),mod=model.matrix(~rep(1, ncol(eset))),batch=eset$batch), silent=TRUE)
        }
        if (class(tmp)=="matrix") {
          if(verbose == T){loginfo(paste("Making ComBat correction to", strEset))}
          exprs(eset) <- tmp
        }
      }
    }
    ##filter genes with standard deviation below the required quantile
    # NO FILTER BY DEFAULT
    if(exists("quantile.cutoff") && quantile.cutoff > 0 && quantile.cutoff < 1){
      eset <- filterQuantile(eset, q=quantile.cutoff)
    }
    ##rescale to z-scores
    # FALSE BY DEFAULT, COMMENTED OUT
    if(exists("rescale") && rescale){
      exprs(eset) <- t(scale(t(exprs(eset))))
    }
    ##samples to be removed
    remove <- rep(FALSE, ncol(eset))
    ##remove samples without required metadata
    # days_to_death AND vital_status REQUIRED BY DEFAULT
    if(exists("meta.required") && length(meta.required) > 0){
      for (varname in meta.required){
        if (varname %in% colnames(pData(eset))){
          remove[ is.na(eset[[varname]]) ] <- TRUE
        }
      }
    }
    ##remove samples not matching regexes
    all.rules <- ls(pattern="rule\\.[0-9]+")
    for (one.rule in all.rules){
      this.remove <- !grepl(get(one.rule)[2], eset[[ get(one.rule)[1] ]])
      # FALSE BY DEFAULT
      if(!strict.checking)
        this.remove[ is.na(eset[[ get(one.rule)[1] ]]) ] <- FALSE
      remove[this.remove] <- TRUE
    }
    ##remove samples pre-specified for removal, that have a dataset specified:
    # NOT RUN BY DEFAULT
    if(exists("remove.samples.delim")){
      if (strEset %in% names(remove.samples.delim)){
        remove[sampleNames(eset) %in% remove.samples.delim[[strEset]]] <- TRUE
      }
    }
    ##remove samples pre-specified for removal, that did *not* have a dataset specified:
    # NOT RUN BY DEFAULT
    if(exists("remove.samples.nodataset"))
      remove[sampleNames(eset) %in% remove.samples.nodataset] <- TRUE
    ##do the actual removal
    eset <- eset[, !remove]
    # NOT RUN BY DEFAULT
    if (exists("considered.datasets") && !(strEset %in% considered.datasets))
    {
      if(verbose == T){loginfo(paste("excluding",strEset,"(considered.datasets)"))}
      next
    }
    ##include study if it has enough samples and events:
    # MINIMUM SET AT 15
    if (exists("min.number.of.events") && !is.na(min.number.of.events)
        && exists("min.sample.size") && !is.na(min.sample.size)
        && min.number.of.events > 0
        && sum(eset$vital_status == "deceased") < min.number.of.events
        || ncol(eset) < min.sample.size)
    {
      if(verbose == T){
        loginfo(paste("excluding",strEset,"(min.number.of.events or min.sample.size)"))
      }
      next
    }
    # MINIMUM SET AT 1000
    if (exists("min.number.of.genes") && nrow(eset) < min.number.of.genes) {
      if(verbose == T){loginfo(paste("excluding",strEset,"(min.number.of.genes)"))}
      next
    }
    # FALSE BY DEFAULT
    if (exists("remove.retracted") && remove.retracted && length(grep("retracted", experimentData(eset)@other$warnings$warnings)) > 0){
      if(verbose == T){loginfo(paste("excluding",strEset,"(remove.retracted)"))}
      next
    }
    # TRUE BY DEFAULT
    if (exists("remove.subsets") && remove.subsets && length(grep("subset", experimentData(eset)@other$warnings$warnings)) > 0){
      if(verbose == T){loginfo(paste("excluding",strEset,"(remove.subsets)"))}
      next
    }
    if(verbose == T){loginfo(paste("including",strEset))}
    esets[[strEset]] <- eset
    rm(eset)
  }

  ## END LOOP

  ##optionally take the intersection of genes common to all platforms:
  # FALSE BY DEFAULT
  if(exists("keep.common.only") && keep.common.only){
    features.per.dataset <- lapply(esets, featureNames)
    intersect.genes <- intersectMany(features.per.dataset)
    esets <- lapply(esets, function(eset){
      eset <- eset[intersect.genes, ]
      return(eset)
    })
  }

  ids.with.missing.data <- which(sapply(esets, function(X) sum(!complete.cases(exprs(X))) > 0))
  if(verbose == T){loginfo(paste("Ids with missing data:", paste(names(ids.with.missing.data), collapse=", ")))}

  # FALSE BY DEFAULT
  if (length(ids.with.missing.data) > 0 && exists("impute.missing") && impute.missing) {
    for (i in ids.with.missing.data) {
      require(impute)
      exprs(esets[[i]]) = impute.knn(exprs(esets[[i]]))$data
    }
  }

  # ADDS days_to_death AND deceased BY DEFAULT
  if (exists("add.surv.y") && is.function(add.surv.y)) {
    for (i in 1:length(esets)) {
      esets[[i]]$y = add.surv.y(esets[[i]])
    }
  }
  if (savedata == T && affyonly == F){
    save(esets,file="data/esets_STAT117.RData")
  } 
  if (savedata == T && affyonly == T){
    save(esets,file="data/esets_affy_STAT117.RData")
  }
  
  return(esets)
}

# For meta analysis dataset
# esets_meta <- createEsetList(combat=F,rescale=F,affyonly=T,savedata=F,verbose=F)

# For dataset replacing Tamoxifen (allow students to make choice on rescale)
# esets_bin <- createEsetList(combat=T,rescale=F,affyonly=F,savedata=F,verbose=F)



