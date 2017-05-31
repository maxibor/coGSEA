method_checker <- function(method_vector){
  if(base::length(method_vector) == 0){
      base::print("Please select at least one GSEA method")
  }
  methods <- list()

  methods$camera = "camera" %in% method_vector
  methods$gage = "gage" %in% method_vector
  methods$globaltest = "globaltest" %in% method_vector
  methods$gsva = "gsva" %in% method_vector
  methods$ssgsea = "ssgsea" %in% method_vector
  methods$zscore = "zscore" %in% method_vector
  methods$plage = "plage" %in% method_vector
  methods$ora = "ora" %in% method_vector
  methods$padog = "padog" %in% method_vector
  methods$roast = "roast" %in% method_vector
  methods$safe = "safe" %in% method_vector
  methods$setrank = "setrank" %in% method_vector

  return(methods)
}

prepareTwoGroupsData <- function(voom.results, contrast, genesetIdx,
        min.samples = NULL, verbose = FALSE){
    if (!base::is.matrix(contrast)){
        # find the reference samples
        if (base::is.null(voom.results$targets) ||
                base::is.null(voom.results$targets$group))
            base::stop(paste0("The data frame 'targets' of the object 'voom.results' ",
                         "must have a column named 'group'."))
        ref.group = base::levels(factor(voom.results$targets$group))[1]
        if (verbose)
            cat(paste0("   Reference group is identified as ", ref.group, "\n"))
        cnt.sam.indx = base::which(voom.results$targets$group == ref.group)
    }
    # use gene indexes instead of gene IDs
    groupData = list()
    gsets = list()
    for (j in 1:base::length(genesetIdx)){
        gsets[[j]] = base::as.character(genesetIdx[[j]])
    }
    base::names(gsets) = base::names(genesetIdx)
    groupData[["gsets"]] = gsets
    # Extract a logCPM matrix for each contrast
    data.log = voom.results$E
    base::rownames(data.log) = base::as.character(base::seq(1, base::nrow(data.log)))
    design = voom.results$design
    sam.idx = 1:base::ncol(data.log)
    groupData[["data"]] = list()
    base::set.seed(05081986)
    contr.num = base::ifelse(base::is.matrix(contrast), base::ncol(contrast), base::length(contrast))
    for(i in 1:contr.num){
        if (base::is.matrix(contrast)){
            # find the indexes of the treatment group samples
            d = design[, contrast[,i] > 0]
            if (base::is.null(base::ncol(d))){
                tre.sam.indx = sam.idx[ d == 1]
            }else if (base::ncol(d) > 1){
                tre.sam.indx = base::c()
                for (j in 1:base::ncol(d))
                    tre.sam.indx = base::c(tre.sam.indx, sam.idx[ d[,j]
                                            == 1])
            }
            else
                base::stop("Invalid contrasts selected.")
            # find the indexes of the control group samples
            d = design[, contrast[,i] < 0]
            if (base::is.null(base::ncol(d))){
                cnt.sam.indx = sam.idx[ d == 1]
            }else if (base::ncol(d) > 1){
                cnt.sam.indx = base::c()
                for (j in 1:base::ncol(d))
                    cnt.sam.indx = base::c(cnt.sam.indx, sam.idx[ d[,j]
                                            == 1])
            }
            else
                stop("Invalid contrasts selected.")
        }else{
            tre.sam.indx = sam.idx[ design[, contrast[i]] == 1]
        }
        # Check if a minimum number of samples is required
        if (! base::is.null(min.samples)){
            if (base::length(tre.sam.indx) == 1)
                tre.sam.indx = base::rep(tre.sam.indx, min.samples)
            else if (base::length(tre.sam.indx) < min.samples)
                tre.sam.indx = base::c(tre.sam.indx,
                        base::sample(tre.sam.indx, min.samples - base::length(tre.sam.indx)))

            if (base::length(cnt.sam.indx) == 1)
                cnt.sam.indx = base::rep(cnt.sam.indx, min.samples)
            else if (base::length(cnt.sam.indx) < min.samples)
                cnt.sam.indx = base::c(cnt.sam.indx,
                        base::sample(cnt.sam.indx, min.samples - base::length(cnt.sam.indx)))
        }
        # logCPM matrix has control samples then treatment samples
        data.log.sel = data.log[, base::c(cnt.sam.indx, tre.sam.indx)]
        groupData$data[[i]] = list()
        groupData$data[[i]][["logCPM"]] = data.log.sel
        # group1 is control / reference
        groupData$data[[i]][["group1"]] = base::seq(1,base::length(cnt.sam.indx))
        groupData$data[[i]][["group2"]] = base::seq(base::length(cnt.sam.indx) +
                        1,base::ncol(data.log.sel))
    }
    base::names(groupData$data) = base::colnames(contrast)
    return(groupData)
}



getNumberofSamples <- function(voom.results, contrast){
    if (base::is.null(voom.results$design)){
        return(0)
    }
    if (base::is.matrix(contrast)){
        samples = base::c()
        sam.idx = base::colnames(voom.results$E)
        for(i in 1:base::ncol(contrast)){
            d = voom.results$design[, contrast[,i] != 0]
            if (base::is.null(base::ncol(d))){
                samples = sam.idx[ d == 1]
            }else if (base::ncol(d) > 1){
                for (j in 1:base::ncol(d))
                    samples = base::c(samples, sam.idx[ d[,j] == 1])
            }
            else
                base::stop("Invalid contrasts selected.")
        }
        return(base::length(base::unique(samples)))
    }else{
        return(base::nrow(voom.results$design))
    }
}


genesetToIdx <- function(geneset = "H", specie = "Mus musculus", entrezGenesIds){

    if (specie == "Mus musculus"){
        if (geneset == "H"){
            base::load(base::url("http://bioinf.wehi.edu.au/software/MSigDB/mouse_H_v5p2.rdata"))
            genesetToIdx = limma::ids2indices(Mm.H, entrezGenesIds)
        } else if (geneset == "C2_KEGG"){
            base::load(base::url("http://bioinf.wehi.edu.au/software/MSigDB/mouse_c2_v5p2.rdata"))
            Mm.c2.subset = Mm.c2[base::grep("KEGG",base::attributes(Mm.c2)$names)]
            genesetToIdx = limma::ids2indices(Mm.c2.subset, entrezGenesIds)
        } else if(geneset == "C2_REACTOME"){
            base::load(base::url("http://bioinf.wehi.edu.au/software/MSigDB/mouse_c2_v5p2.rdata"))
            Mm.c2.subset = Mm.c2[base::grep("REACTOME",base::attributes(Mm.c2)$names)]
            genesetToIdx = limma::ids2indices(Mm.c2.subset, entrezGenesIds)
        } else if (base::class(geneset) == "list"){
            genesetToIdx = limma::ids2indices(geneset, entrezGenesIds, remove.empty=TRUE)
        } else {
            cat(paste(geneset , "is not a valid MSigDB geneset. Available genesets : \n - H \n - C2"))
            return(NULL)
        }
    } else if (specie == "Homo sapiens"){
        if (geneset == "H"){
            base::load(base::url("http://bioinf.wehi.edu.au/software/MSigDB/human_H_v5p2.rdata"))
            genesetToIdx = limma::ids2indices(Hs.H, entrezGenesIds)
        } else if (geneset == "C2_KEGG"){
            base::load(base::url("http://bioinf.wehi.edu.au/software/MSigDB/human_c2_v5p2.rdata"))
            Hs.c2.subset = Hs.c2[base::grep("KEGG",base::attributes(Hs.c2)$names)]
            genesetToIdx = limma::ids2indices(Hs.c2.subset, entrezGenesIds)
        } else if (geneset == "C2_REACTOME"){
            base::load(base::url("http://bioinf.wehi.edu.au/software/MSigDB/human_c2_v5p2.rdata"))
            Hs.c2.subset = Hs.c2[base::grep("REACTOME",base::attributes(Hs.c2)$names)]
            genesetToIdx = limma::ids2indices(Hs.c2.subset, entrezGenesIds)
        } else if (base::class(geneset) == "list"){
            genesetToIdx = limma::ids2indices(geneset, entrezGenesIds, remove.empty=TRUE)
        } else {
            cat(paste(geneset , "is not a valid MSigDB geneset. Available genesets : \n - H \n - C2"))
            return(NULL)
        }
    } else {
        cat(paste(specie,"is not a valid specie. Available specie : \n - Mus musculus \n - Homo sapiens"))
        return(NULL)
    }
    return(genesetToIdx)
}

msig_to_setrankDb = function(msigdb, orga, dbname){
  # msigdb : database R object
  # orga : organism name (HUMAN|MOUSE)
  # dbname : Msig database name (C2|H)
  # example : mouse_H = msig_to_setrankDb(Mm.H, "MOUSE", "H")
  i = 1
  df = data.frame(geneID = factor(), termID = factor(), termName = factor(), dbName = factor(), description = factor())
  for (name in(base::names(msigdb))){
    for (gene in msigdb[[name]]){
      df = base::rbind(df, data.frame(geneID = factor(gene), termID = factor(paste(orga,"MSigDB",dbname,i,sep = "_")), termName = factor(name), dbName = factor("MSigDB"), description = factor("")))
    }
    i = i+1
  }
  return(df)
}

dgeToElist = function(dgeobject){
  tmp = methods::new("EList")
  for (i in base::names(dge)){
    tmp[[i]] = dge[[i]]
  }
  if (base::is.null(tmp$E)){
    tmp$E = tmp$counts
  } else if (base::is.null(tmp$counts)){
    tmp$counts = tmp$E
  } else {
    cat("There is a problem with the count matrix")
    return()
  }
  if(tmp$genes$ENTREZ){
    tmp$genes = tmp$genes$ENTREZ
  }
  base::rownames(tmp$E) = tmp$genes
  base::rownames(tmp$counts) = base::rownames(tmp$E)
  if (base::is.null(tmp$common.dispersion)){
    cat("Dispersion estimates needed.\nPlease run estimateDisp() on your dge Object.\nRun '?estimateDisp()' for help")
    return()
  }
  return(DGE)
}

cGSEAcore = function(ElistObject, contrast, geneset = "H", specie = "Mus musculus", entrezGenesIds, camera = TRUE, gage = TRUE, globaltest = TRUE, gsva = TRUE, ssgsea = TRUE, zscore = TRUE, plage = TRUE, ora = TRUE, padog = TRUE, roast = TRUE, safe = TRUE, setrank = FALSE, num.workers = 4, verbose = TRUE){

  genesetIdx = genesetToIdx(geneset = geneset, specie = specie, entrezGenesIds = entrezGenesIds)

  output = list()
  res = list()
  time = list()

  if(camera == TRUE){
    print("RUNNING CAMERA")
    camerastart = base::Sys.time()
    res$camera = runcamera(voom.results = ElistObject, contrast = contrast, genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
    camerastop = base::Sys.time()
    time$camera = base::difftime(camerastop, camerastart, units = c("secs"))
  }

  if (gage == TRUE){
    print("RUNNING GAGE")
    gagestart = base::Sys.time()
    res$gage = rungage(voom.results = ElistObject, contrast = contrast, genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
    gagestop = base::Sys.time()
    time$gage = base::difftime(gagestop,gagestart, units = c("secs"))
  }
  if (globaltest == TRUE){
    print("RUNNING GLOBALTEST")
    globalteststart = base::Sys.time()
    res$globaltest = runglobaltest(voom.results = ElistObject, contrast = contrast, genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
    globalteststop = base::Sys.time()
    time$globaltest = base::difftime(globalteststop,globalteststart, units = c("secs"))
  }
  if (gsva == TRUE){
    print("RUNNING GSVA")
    gsvastart = base::Sys.time()
    res$gsva = rungsva(method = "gsva", voom.results = ElistObject, contrast = contrast,  genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
    gsvastop = base::Sys.time()
    time$gsva = base::difftime(gsvastop,gsvastart, units = c("secs"))
  }
  if (ssgsea == TRUE){
    print("RUNNING SSGSEA")
    ssgseastart = base::Sys.time()
    res$ssgsea = rungsva(method = "ssgsea", voom.results = ElistObject, contrast = contrast,  genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
    ssgseastop = base::Sys.time()
    time$ssgsea = base::difftime(ssgseastop,ssgseastart, units = c("secs"))
  }
  if (zscore == TRUE){
    print("RUNNING ZSCORE")
    zscorestart = base::Sys.time()
    res$zscore = rungsva(method = "zscore", voom.results = ElistObject, contrast = contrast,  genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
    zscorestop = base::Sys.time()
    time$zscore = base::difftime(zscorestop,zscorestart, units = c("secs"))
  }
  if (plage == TRUE){
    print("RUNNING PLAGE")
    plagestart = base::Sys.time()
    res$plage = rungsva(method = "plage", voom.results = ElistObject, contrast = contrast,  genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
    plagestop = base::Sys.time()
    time$plage = base::difftime(plagestop,plagestart, units = c("secs"))
  }
  if (ora == TRUE){
    print("RUNNING ORA")
    orastart = base::Sys.time()
    res$ora = runora(voom.results = ElistObject, contrast = contrast,  genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
    orastop = base::Sys.time()
    time$ora = base::difftime(orastop, orastart, units = c("secs"))
  }
  if (padog == TRUE){
    print("RUNNING PADOG")
    padogstart = base::Sys.time()
    res$padog = runpadog(voom.results = ElistObject, contrast = contrast,  genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
    padogstop = base::Sys.time()
    time$padog = base::difftime(padogstop, padogstart, units = c("secs"))
  }
  if (roast == TRUE){
    print("RUNNING ROAST")
    roaststart = base::Sys.time()
    res$roast = runroast(voom.results = ElistObject, contrast = contrast,  genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
    roaststop = base::Sys.time()
    time$roast = base::difftime(roaststop,roaststart, units = c("secs"))
  }
  if (safe == TRUE){
    print("RUNNING SAFE")
    safestart = base::Sys.time()
    res$safe = runsafe(voom.results = ElistObject, contrast = contrast,  genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
    safestop = base::Sys.time()
    time$safe = base::difftime(safestop, safestart, units = c("secs"))
  }
  if (setrank == TRUE){
    print("RUNNING SETRANK")
    setrankstart = base::Sys.time()
    res$setrank = runsetrank(voom.results = ElistObject, contrast = contrast,  geneset = geneset, specie = specie, num.workers = num.workers, verbose = verbose)
    setrankstop = base::Sys.time()
    time$setrank = base::difftime(setrankstop, setrankstart, units = c("secs"))
  }
  output$time = time
  output$res = res
  output$contrast = contrast
  output$collection = genesetIdx
  output$Elist = ElistObject
  return(output)
}

rank_getter = function(geneSetname, method_result){
  out = list()
  ranks = base::c()
  pvals = base::c()
  for (aset in geneSetname){
    for (i in base::seq(1,base::nrow(method_result))){
      if (base::rownames(method_result)[i] == aset){
        ranks = base::append(ranks, base::as.data.frame(method_result)[["Rank"]][i])
        pvals = base::append(pvals, base::as.data.frame(method_result)[["p.value"]][i])
      }
    }
  }
  out$ranks = ranks
  out$pvals = pvals
  return(out)
}


cGSEAOutputTable = function(cGSEAcoreOutput, contrastLevel = contrastLevel, geneSetLevels){
  dt = list()
  for (method in base::names(cGSEAcoreOutput$res)){
    method_stats = rank_getter(geneSetname = geneSetLevels,method_result = cGSEAcoreOutput[["res"]][[method]][[contrastLevel]])

    dt[[paste(method,"Rank",sep="_")]] = base::as.numeric(method_stats[["ranks"]])
    # ranks[[method]] = rank_getter(geneSetname = geneSetLevels,method_result = cGSEAOutput[["res"]][[method]][[contrastLevel]])
    #p.values
    dt[[paste(method,"p.value",sep="_")]] = base::as.numeric(method_stats[["pvals"]])
  }

  # df = data.frame(matrix(unlist(dt), nrow = length(geneSetLevels), byrow = FALSE),stringsAsFactors=FALSE)
  df = data.frame(matrix(base::unlist(dt), nrow = base::dim(cGSEAcoreOutput$res[[1]][[1]])[1], byrow = FALSE),stringsAsFactors=FALSE)

  base::colnames(df) = base::names(dt)
  base::rownames(df) = geneSetLevels
  return(df)
}

getlogFCFromLMFit <- function(voom.results, contrast,
        logFC.cutoff, fdr.cutoff){
    # to be changed for gene symbols support
    base::stopifnot(base::class(voom.results) == "EList")
    print("limma DE analysis is carried out ... ")
    # fit linear model for each gene using limma package functions
    vfit = limma::lmFit(voom.results, design=voom.results$design) # Fit linear model
# for each gene given a series of arrays
    if (base::is.matrix(contrast)){
        vfit = limma::contrasts.fit(vfit, contrast) # make all pair-wise comparisons
        contr.names = base::colnames(contrast)
        contr.num =  base::ncol(contrast)
        coefs = 1:base::ncol(contrast)
# between the groups
    }else{
        contr.names = base::names(contrast)
        contr.num = base::length(contrast)
        coefs = contrast
    }
    ebayes.results = limma::eBayes(vfit) # compute moderated t-statistics, moderated
#F-statistic, and log-odds of differential expression by empirical
#    Bayes moderation of the standard errors towards a common value
    logFC = matrix(0, base::nrow(ebayes.results), contr.num)
    limma.tops = list()
    for (i in 1:base::length(coefs)){
        top.table = limma::topTable(ebayes.results, coef=coefs[i],
                number=Inf, sort.by="none")
        limma.fc = top.table$logFC
        base::names(limma.fc) = base::rownames(ebayes.results)
        logFC[, i] = limma.fc
        base::rownames(top.table) = base::rownames(ebayes.results)
        limma.tops[[contr.names[i]]] = top.table
        de.genes = top.table[top.table[, "adj.P.Val"] <= fdr.cutoff, ]
        if (base::nrow(de.genes) == 0)
            cat(paste0("WARNING: it seems the contrast ",
                    contr.names[i],
                    " has no DE genes at the selected 'fdr.cutoff'.\n",
                    "The 'fdr.cutoff' was ignored in the calculations.\n"))
        if (base::nrow(de.genes) > 0){
            de.genes = de.genes[base::abs(de.genes[, "logFC"]) >= logFC.cutoff, ]
            if (base::nrow(de.genes) == 0)
                cat(paste0("WARNING: it seems the contrast ",
                    contr.names[i],
                    " has no DE genes at the selected 'logFC.cutoff'.\n",
                    "The 'logFC.cutoff' was ignored in the calculations.\n"))
        }
    }

    base::rownames(logFC) = base::rownames(ebayes.results)
    base::colnames(logFC) = contr.names

    return(list(logFC=logFC, limma.results=ebayes.results, limma.tops=limma.tops))
}

signifCal = function(combiPval, avgLFC){
  # combiPval : combined Pvalue after correction (correction  : BH, combination : fisher)
  #avgLFC : average logFC in geneset from edgerR/limma
  signif = -base::log10(combiPval)*base::abs(avgLFC)
  return(signif)
}

geneSetsLogFC = function(logFcMatrix, genesetCollection, condition){
  geneSetCollectionLogFC = base::c()
  for (geneset in base::names(genesetCollection)){
    setLogFC = base::c()
    for (gene in genesetCollection[[geneset]]){
      setLogFC = base::append(setLogFC, x = logFcMatrix$limma.tops[[condition]]$logFC[base::which(base::rownames(logFcMatrix$limma.tops[[condition]]) == gene)])
    }
    avgLogFC = base::mean(setLogFC)
    geneSetCollectionLogFC = base::append(geneSetCollectionLogFC, avgLogFC)
  }
  base::names(geneSetCollectionLogFC) = base::names(genesetCollection)
  return(geneSetCollectionLogFC)
}

adjustPval =  function(pvalTable, type){
  switch(type,
    holm = base::apply(pvalTable, 2, p.adjust, method = "holm"),
    hochberg = base::apply(pvalTable, 2, p.adjust, method = "hochberg"),
    hommel = base::apply(pvalTable, 2, p.adjust, method = "hommel"),
    bonferroni = base::apply(pvalTable, 2, p.adjust, method = "bonferroni"),
    BH = base::apply(pvalTable, 2, p.adjust, method = "BH"),
    BY = base::apply(pvalTable, 2, p.adjust, method = "BY"),
    fdr = base::apply(pvalTable, 2, p.adjust, method = "fdr"),
    none = base::apply(pvalTable, 2, p.adjust, method = "none"))
}

combinePval =  function(pvalAdjTable, type){

  switch(type,
    sumz = base::apply(pvalAdjTable, 1, sumz),
    votep = base::apply(pvalAdjTable, 1, votep),
    minimump = base::apply(pvalAdjTable, 1, minimump),
    sumlog = base::apply(pvalAdjTable, 1, sumlog),
    sump = base::apply(pvalAdjTable, 1, sump),
    logitp = base::apply(pvalAdjTable, 1, logitp),
    meanp = base::apply(pvalAdjTable, 1, meanp),
    maximump = base::apply(pvalAdjTable, 1, maximump))
}

comparisonSummaryData = function(preparedDataResumPlot2){
  m = base::do.call(cbind, preparedDataResumPlot2)
  base::colnames(m) = base::rep(base::c("x.data","y.data","dir","rank","sig","id","gsSize"),base::length(base::names(preparedDataResumPlot2)))
  m = stats::na.omit(m)
  m.id = base::abbreviate(base::rownames(m), minlength = 6, use.classes = FALSE, dot = FALSE, named = FALSE)
  m = m[,-base::grep("id", colnames(m))]
  base::colnames(m) = base::rep(base::c("x.data","y.data","dir","rank","sig","gsSize"),base::length(base::names(preparedDataResumPlot2)))
  nms = base::colnames(m)
  m2 = base::sapply(base::unique(nms), function(x) rowMeans (m[, nms == x]))
  m2 = base::as.data.frame(m2)
  m2 = base::cbind(m2,m.id)
  m2[,base::c("x.data","y.data","dir","rank","sig","gsSize")] = base::apply(m2[,base::c("x.data","y.data","dir","rank","sig","gsSize")],2,as.numeric)
  base::colnames(m2)[base::ncol(m2)] = "id"
  return(m2)
}

combineRanks = function(resTable){
    rankCol = base::grep("_Rank", base::colnames(resTable))
    meanRank = base::apply(resTable[,rankCol],1,mean)
    return(base::cbind(resTable, base::as.numeric(meanRank)))
}

genes_in_gs = function(geneset_collection, gene_list){

    gs_check = sapply(geneset_collection, function(x) x %in% gene_list)

    for (i in 1:length(names(geneset_collection))){
        geneset_collection[[i]] = geneset_collection[[i]][gs_check[[i]]]
    }

    return(geneset_collection)
}

prepareData = function(cGSEAcoreOutput, alpha = 0.05, directoryPath, pvalAdjMethod = "BH", pvalCombMethod = "sumlog", min.intersection.size = 1, shinyMode = FALSE){
  #Saving genesets and genes to a file

  genecollec = genes_in_gs(geneset_collection = cGSEAcoreOutput$collection, gene_list = base::rownames(cGSEAcoreOutput$Elist$E))
  if (shinyMode == FALSE){
    base::sink(paste(directoryPath,"/geneset_collection_genes.csv", sep = ""))
    writeLines(paste(names(genecollec),unlist(lapply(genecollec, paste, collapse=",")), sep = ","))
    sink()
  }


  logFCTable = getlogFCFromLMFit(voom.results = cGSEAcoreOutput$Elist, contrast = cGSEAcoreOutput$contrast, logFC.cutoff = 0, fdr.cutoff = 1)

  time = base::unlist(cGSEAcoreOutput$time)
  base::names(time) = base::names(cGSEAcoreOutput$time)


  output = list()
  result = list()
  heatmap = list()
  PCA = list()
  abrev = list()
  correlation = list()
  clustering = list()
  snailPlot = list()
  resumPlot1 = list()
  resumPlot2 = list()
  resumPlot2Table = list()
  limma.out = logFCTable



  for (condi in base::colnames(cGSEAcoreOutput$contrast)){
    resTable = cGSEAOutputTable(cGSEAcoreOutput = cGSEAcoreOutput, contrastLevel = condi, geneSetLevels = names(cGSEAcoreOutput$collection))
    resTable = stats::na.omit(resTable)
    # pvalue adjustment and combination
    pvalcol = base::grep("p.value", base::colnames(resTable))
    pvalTable = resTable[,pvalcol]
    #pvalAdjTable = apply(pvalTable, 2, p.adjust, method = "BH") #Benjamini Hochberg for pvalue adjustment
    pvalAdjTable = adjustPval(pvalTable, type = pvalAdjMethod)
    base::colnames(pvalAdjTable) = base::gsub("_p.value","_adj_p.value",base::colnames(pvalAdjTable))
    # pvalComb = apply(pvalAdjTable, 1, sumlog) #fisher method for pvalue combination
    pvalComb = combinePval(pvalAdjTable, type = pvalCombMethod)
    pvalComb = base::unlist(base::lapply(pvalComb, '[[',3))

    #result Table
    resTable2 = base::cbind(resTable, pvalAdjTable)
    resTable3 = base::cbind(resTable2, pvalComb)
    base::colnames(resTable3)[base::ncol(resTable3)] = "combined_p.value"
    result[[condi]] = stats::na.omit(resTable3)
    result[[condi]] = combineRanks(result[[condi]])
    base::colnames(result[[condi]])[base::ncol(result[[condi]])] = "Avg_Rank"
    if (shinyMode == FALSE){
      utils::write.csv(result[[condi]], file = paste(directoryPath,"/result_",condi,".csv", sep = ""))
    }
    #ranks
    rankCol = base::grep("Rank",base::colnames(resTable))
    rankTable = resTable[,rankCol]
    base::colnames(rankTable) = base::gsub("_Rank","", base::colnames(rankTable))

    #data for pVal binary heatmap
    heatmap[[condi]] = base::t(base::apply(pvalAdjTable,1, function(x) base::ifelse(x < alpha, 1, 0)))
    base::colnames(heatmap[[condi]]) = base::gsub("_adj_p.value","",base::colnames(heatmap[[condi]]))
    this.apply = base::apply(heatmap[[condi]], 1, sum)
    max_methods = base::max(this.apply)
    heatmap[[condi]] = heatmap[[condi]][base::c(base::which(this.apply == max_methods),base::which(this.apply == (max_methods-1))),]

    #data for PCA
    PCA[[condi]] = base::t(rankTable)

    #data for clustering
    clustering[[condi]] = stats::hclust(stats::dist(base::t(rankTable), method = "euclidean"), method = "ward.D2")

    #data for correlation plot
    correlation[[condi]] = stats::cor(rankTable)

    #for SnailPlot
    gsFactor = base::as.factor(base::names(cGSEAcoreOutput$collection))
    listInput = base::apply(pvalAdjTable, 2, function(x) base::names(base::which(x < alpha)))
    cond = base::sapply(listInput, function(x) base::length(x) > 0)
    listInput = listInput[cond]
    listInput = base::sapply(listInput, function(x) factor(x, levels = base::levels(gsFactor)))
    base::names(listInput) = base::gsub("_adj_p.value","", base::names(listInput))
    snailPlot[[condi]] = listInput

    abrev_names_vector = base::abbreviate(base::gsub("HALLMARK_","H_",base::names(cGSEAcoreOutput$collection)), minlength = 6, use.classes = FALSE, dot = FALSE, named = FALSE, method = "left.kept")
    abrev_names_mat = base::cbind(cGSEAcoreOutput$collection,abrev_names_vector)
    base::colnames(abrev_names_mat) = base::c("original","abbreviation")
    abrev[[condi]] = base::t(base::as.data.frame(abrev_names_mat[,2]))
    abrev[[condi]] = base::cbind(base::rownames(abrev[[condi]]), abrev[[condi]])
    base::colnames(abrev[[condi]]) = base::c("Gene Set Name","abbreviation")
    if (shinyMode == FALSE){
      utils::write.csv(abrev[[condi]], file = paste(directoryPath,"/abrreviations_",condi,".csv", sep = ""))
    }

    #for resumPlot1 - simple summary plot
    gsLogFC = geneSetsLogFC(logFcMatrix = logFCTable, genesetCollection = cGSEAcoreOutput$collection, condition = condi)
    resumPlot1[[condi]]$x = -base::log10(pvalComb)
    resumPlot1[[condi]]$y = gsLogFC
    resumPlot1[[condi]]$labels = abrev_names_vector

    #for resumPlot2 - advanced summary plot
    signifScore = signifCal(combiPval = pvalComb, avgLFC = gsLogFC)
    avgRank = base::apply(resTable[,rankCol],1,mean)

    resumPlot2[[condi]] = base::cbind(base::as.numeric(-base::log10(pvalComb)), base::abs(base::as.numeric(gsLogFC)), base::as.numeric(gsLogFC/base::abs(gsLogFC)), base::as.numeric(avgRank), base::as.numeric(signifScore), abrev_names_vector, base::unlist(base::lapply(cGSEAcoreOutput$collection, length)))
    base::colnames(resumPlot2[[condi]]) = base::c("x.data","y.data","dir","rank","sig","id", "gsSize")
    base::rownames(resumPlot2[[condi]]) = base::names(cGSEAcoreOutput$collection)
    resumPlot2[[condi]] = base::as.data.frame(resumPlot2[[condi]])
    resumPlot2[[condi]][,base::c("x.data","y.data","dir","rank","sig", "gsSize")] = base::apply(resumPlot2[[condi]][,base::c("x.data","y.data","dir","rank","sig", "gsSize")],2, as.numeric)

    #limma data

    if (shinyMode == FALSE){
        utils::write.csv(logFCTable$limma.tops[[condi]], file = paste(directoryPath,"limma_results_",condi,".csv",sep = ""))
    }

    # resumPlot2Table[[condi]] = cbind(as.numeric(-log10(pvalComb)), abs(as.numeric(gsLogFC)), as.numeric(gsLogFC/abs(gsLogFC)), as.numeric(avgRank), as.numeric(signifScore), abrev_names_vector, unlist(lapply(cGSEAcoreOutput$collection, length)))
    # colnames(resumPlot2Table[[condi]]) = c("x.data","y.data","dir","rank","sig","id", "gsSize")
    # rownames(resumPlot2Table[[condi]]) = names(cGSEAcoreOutput$collection)
    # resumPlot2Table[[condi]] = as.data.frame(resumPlot2[[condi]])
    # resumPlot2Table[[condi]][,c("x.data","y.data","dir","rank","sig", "gsSize")] = apply(resumPlot2[[condi]][,c("x.data","y.data","dir","rank","sig", "gsSize")],2, as.numeric)
  }
  if (base::length(base::colnames(cGSEAcoreOutput$contrast)) > 1){
    output$comparison = comparisonSummaryData(resumPlot2)
  }

  # output$comparisonTable = comparisonSummaryData(resumPlot2Table)

  output$time = time
  output$result = result
  output$heatmap = heatmap
  output$PCA = PCA
  output$contrast = cGSEAcoreOutput$contrast
  output$clustering = clustering
  output$correlation = correlation
  output$snailPlot = snailPlot
  output$resumPlot1 = resumPlot1
  output$resumPlot2 = resumPlot2
  output$minIntersectionSize = min.intersection.size
  output$abbreviation = abrev

  return(output)

}
