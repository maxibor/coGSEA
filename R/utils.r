method_checker <- function(method_vector){
  if(length(method_vector) == 0){
      print("Please select at least one GSEA method")
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
    if (!is.matrix(contrast)){
        # find the reference samples
        if (is.null(voom.results$targets) ||
                is.null(voom.results$targets$group))
            stop(paste0("The data frame 'targets' of the object 'voom.results' ",
                         "must have a column named 'group'."))
        ref.group = levels(factor(voom.results$targets$group))[1]
        if (verbose)
            cat(paste0("   Reference group is identified as ", ref.group, "\n"))
        cnt.sam.indx = which(voom.results$targets$group == ref.group)
    }
    # use gene indexes instead of gene IDs
    groupData = list()
    gsets = list()
    for (j in 1:length(genesetIdx)){
        gsets[[j]] = as.character(genesetIdx[[j]])
    }
    names(gsets) = names(genesetIdx)
    groupData[["gsets"]] = gsets
    # Extract a logCPM matrix for each contrast
    data.log = voom.results$E
    rownames(data.log) = as.character(seq(1, nrow(data.log)))
    design = voom.results$design
    sam.idx = 1:ncol(data.log)
    groupData[["data"]] = list()
    set.seed(05081986)
    contr.num = ifelse(is.matrix(contrast), ncol(contrast), length(contrast))
    for(i in 1:contr.num){
        if (is.matrix(contrast)){
            # find the indexes of the treatment group samples
            d = design[, contrast[,i] > 0]
            if (is.null(ncol(d))){
                tre.sam.indx = sam.idx[ d == 1]
            }else if (ncol(d) > 1){
                tre.sam.indx = c()
                for (j in 1:ncol(d))
                    tre.sam.indx = c(tre.sam.indx, sam.idx[ d[,j]
                                            == 1])
            }
            else
                stop("Invalid contrasts selected.")
            # find the indexes of the control group samples
            d = design[, contrast[,i] < 0]
            if (is.null(ncol(d))){
                cnt.sam.indx = sam.idx[ d == 1]
            }else if (ncol(d) > 1){
                cnt.sam.indx = c()
                for (j in 1:ncol(d))
                    cnt.sam.indx = c(cnt.sam.indx, sam.idx[ d[,j]
                                            == 1])
            }
            else
                stop("Invalid contrasts selected.")
        }else{
            tre.sam.indx = sam.idx[ design[, contrast[i]] == 1]
        }
        # Check if a minimum number of samples is required
        if (! is.null(min.samples)){
            if (length(tre.sam.indx) == 1)
                tre.sam.indx = rep(tre.sam.indx, min.samples)
            else if (length(tre.sam.indx) < min.samples)
                tre.sam.indx = c(tre.sam.indx,
                        sample(tre.sam.indx, min.samples - length(tre.sam.indx)))

            if (length(cnt.sam.indx) == 1)
                cnt.sam.indx = rep(cnt.sam.indx, min.samples)
            else if (length(cnt.sam.indx) < min.samples)
                cnt.sam.indx = c(cnt.sam.indx,
                        sample(cnt.sam.indx, min.samples - length(cnt.sam.indx)))
        }
        # logCPM matrix has control samples then treatment samples
        data.log.sel = data.log[, c(cnt.sam.indx, tre.sam.indx)]
        groupData$data[[i]] = list()
        groupData$data[[i]][["logCPM"]] = data.log.sel
        # group1 is control / reference
        groupData$data[[i]][["group1"]] = seq(1,length(cnt.sam.indx))
        groupData$data[[i]][["group2"]] = seq(length(cnt.sam.indx) +
                        1,ncol(data.log.sel))
    }
    names(groupData$data) = colnames(contrast)
    return(groupData)
}


getNumberofSamples <- function(voom.results, contrast){
    if (is.null(voom.results$design)){
        return(0)
    }
    if (is.matrix(contrast)){
        samples = c()
        sam.idx = colnames(voom.results$E)
        for(i in 1:ncol(contrast)){
            d = voom.results$design[, contrast[,i] != 0]
            if (is.null(ncol(d))){
                samples = sam.idx[ d == 1]
            }else if (ncol(d) > 1){
                for (j in 1:ncol(d))
                    samples = c(samples, sam.idx[ d[,j] == 1])
            }
            else
                stop("Invalid contrasts selected.")
        }
        return(length(unique(samples)))
    }else{
        return(nrow(voom.results$design))
    }
}


genesetToIdx <- function(geneset = "H", specie = "Mus musculus", entrezGenesIds){

	if(!(require(limma))){
		source("http://www.bioconductor.org/biocLite.R")
		biocLite("limma")
		require(limma)
	}

    if (specie == "Mus musculus"){
        if (geneset == "H"){
            load(url("http://bioinf.wehi.edu.au/software/MSigDB/mouse_H_v5p2.rdata"))
            genesetToIdx = ids2indices(Mm.H, entrezGenesIds)
        } else if (geneset == "C2_KEGG"){
            load(url("http://bioinf.wehi.edu.au/software/MSigDB/mouse_c2_v5p2.rdata"))
            Mm.c2.subset = Mm.c2[grep("KEGG",attributes(Mm.c2)$names)]
            genesetToIdx = ids2indices(Mm.c2.subset, entrezGenesIds)
        } else if(geneset == "C2_REACTOME"){
            load(url("http://bioinf.wehi.edu.au/software/MSigDB/mouse_c2_v5p2.rdata"))
            Mm.c2.subset = Mm.c2[grep("REACTOME",attributes(Mm.c2)$names)]
            genesetToIdx = ids2indices(Mm.c2.subset, entrezGenesIds)
        } else if (class(geneset) == "list"){
            genesetToIdx = ids2indices(geneset, entrezGenesIds, remove.empty=TRUE)
        } else {
            cat(paste(geneset , "is not a valid MSigDB geneset. Available genesets : \n - H \n - C2"))
            return(NULL)
        }
    } else if (specie == "Homo sapiens"){
        if (geneset == "H"){
            load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_H_v5p2.rdata"))
            genesetToIdx = ids2indices(Hs.H, entrezGenesIds)
        } else if (geneset == "C2_KEGG"){
            load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c2_v5p2.rdata"))
            Hs.c2.subset = Hs.c2[grep("KEGG",attributes(Hs.c2)$names)]
            genesetToIdx = ids2indices(Hs.c2.subset, entrezGenesIds)
        } else if (geneset == "C2_REACTOME"){
            load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c2_v5p2.rdata"))
            Hs.c2.subset = Hs.c2[grep("REACTOME",attributes(Hs.c2)$names)]
            genesetToIdx = ids2indices(Hs.c2.subset, entrezGenesIds)
        } else if (class(geneset) == "list"){
            genesetToIdx = ids2indices(geneset, entrezGenesIds, remove.empty=TRUE)
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
  for (name in(names(msigdb))){
    for (gene in msigdb[[name]]){
      df = rbind(df, data.frame(geneID = factor(gene), termID = factor(paste(orga,"MSigDB",dbname,i,sep = "_")), termName = factor(name), dbName = factor("MSigDB"), description = factor("")))
    }
    i = i+1
  }
  return(df)
}

dgeToElist = function(dgeobject){
  tmp = new("EList")
  for (i in names(dge)){
    tmp[[i]] = dge[[i]]
  }
  if (is.null(tmp$E)){
    tmp$E = tmp$counts
  } else if (is.null(tmp$counts)){
    tmp$counts = tmp$E
  } else {
    cat("There is a problem with the count matrix")
    return()
  }
  if(tmp$genes$ENTREZ){
    tmp$genes = tmp$genes$ENTREZ
  }
  rownames(tmp$E) = tmp$genes
  rownames(tmp$counts) = rownames(tmp$E)
  if (is.null(tmp$common.dispersion)){
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
    camerastart = Sys.time()
    res$camera = runcamera(voom.results = ElistObject, contrast = contrast, genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
    camerastop = Sys.time()
    time$camera = difftime(camerastop, camerastart, units = c("secs"))
  }

  if (gage == TRUE){
    print("RUNNING GAGE")
    gagestart = Sys.time()
    res$gage = rungage(voom.results = ElistObject, contrast = contrast, genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
    gagestop = Sys.time()
    time$gage = difftime(gagestop,gagestart, units = c("secs"))
  }
  if (globaltest == TRUE){
    print("RUNNING GLOBALTEST")
    globalteststart = Sys.time()
    res$globaltest = runglobaltest(voom.results = ElistObject, contrast = contrast, genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
    globalteststop = Sys.time()
    time$globaltest = difftime(globalteststop,globalteststart, units = c("secs"))
  }
  if (gsva == TRUE){
    print("RUNNING GSVA")
    gsvastart = Sys.time()
    res$gsva = rungsva(method = "gsva", voom.results = ElistObject, contrast = contrast,  genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
    gsvastop = Sys.time()
    time$gsva = difftime(gsvastop,gsvastart, units = c("secs"))
  }
  if (ssgsea == TRUE){
    print("RUNNING SSGSEA")
    ssgseastart = Sys.time()
    res$ssgsea = rungsva(method = "ssgsea", voom.results = ElistObject, contrast = contrast,  genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
    ssgseastop = Sys.time()
    time$ssgsea = difftime(ssgseastop,ssgseastart, units = c("secs"))
  }
  if (zscore == TRUE){
    print("RUNNING ZSCORE")
    zscorestart = Sys.time()
    res$zscore = rungsva(method = "zscore", voom.results = ElistObject, contrast = contrast,  genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
    zscorestop = Sys.time()
    time$zscore = difftime(zscorestop,zscorestart, units = c("secs"))
  }
  if (plage == TRUE){
    print("RUNNING PLAGE")
    plagestart = Sys.time()
    res$plage = rungsva(method = "plage", voom.results = ElistObject, contrast = contrast,  genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
    plagestop = Sys.time()
    time$plage = difftime(plagestop,plagestart, units = c("secs"))
  }
  if (ora == TRUE){
    print("RUNNING ORA")
    orastart = Sys.time()
    res$ora = runora(voom.results = ElistObject, contrast = contrast,  genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
    orastop = Sys.time()
    time$ora = difftime(orastop, orastart, units = c("secs"))
  }
  if (padog == TRUE){
    print("RUNNING PADOG")
    padogstart = Sys.time()
    res$padog = runpadog(voom.results = ElistObject, contrast = contrast,  genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
    padogstop = Sys.time()
    time$padog = difftime(padogstop, padogstart, units = c("secs"))
  }
  if (roast == TRUE){
    print("RUNNING ROAST")
    roaststart = Sys.time()
    res$roast = runroast(voom.results = ElistObject, contrast = contrast,  genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
    roaststop = Sys.time()
    time$roast = difftime(roaststop,roaststart, units = c("secs"))
  }
  if (safe == TRUE){
    print("RUNNING SAFE")
    safestart = Sys.time()
    res$safe = runsafe(voom.results = ElistObject, contrast = contrast,  genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
    safestop = Sys.time()
    time$safe = difftime(safestop, safestart, units = c("secs"))
  }
  if (setrank == TRUE){
    print("RUNNING SETRANK")
    setrankstart = Sys.time()
    res$setrank = runsetrank(voom.results = ElistObject, contrast = contrast,  geneset = geneset, specie = specie, num.workers = num.workers, verbose = verbose)
    setrankstop = Sys.time()
    time$setrank = difftime(setrankstop, setrankstart, units = c("secs"))
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
  ranks = c()
  pvals = c()
  for (aset in geneSetname){
    for (i in seq(1,nrow(method_result))){
      if (rownames(method_result)[i] == aset){
        ranks = append(ranks, as.data.frame(method_result)[["Rank"]][i])
        pvals = append(pvals, as.data.frame(method_result)[["p.value"]][i])
      }
    }
  }
  out$ranks = ranks
  out$pvals = pvals
  return(out)
}


cGSEAOutputTable = function(cGSEAcoreOutput, contrastLevel = contrastLevel, geneSetLevels){
  dt = list()
  for (method in names(cGSEAcoreOutput$res)){
    method_stats = rank_getter(geneSetname = geneSetLevels,method_result = cGSEAcoreOutput[["res"]][[method]][[contrastLevel]])

    dt[[paste(method,"Rank",sep="_")]] = as.numeric(method_stats[["ranks"]])
    # ranks[[method]] = rank_getter(geneSetname = geneSetLevels,method_result = cGSEAOutput[["res"]][[method]][[contrastLevel]])
    #p.values
    dt[[paste(method,"p.value",sep="_")]] = as.numeric(method_stats[["pvals"]])
  }

  # df = data.frame(matrix(unlist(dt), nrow = length(geneSetLevels), byrow = FALSE),stringsAsFactors=FALSE)
  df = data.frame(matrix(unlist(dt), nrow = dim(cGSEAcoreOutput$res[[1]][[1]])[1], byrow = FALSE),stringsAsFactors=FALSE)

  colnames(df) = names(dt)
  rownames(df) = geneSetLevels
  return(df)
}

getlogFCFromLMFit <- function(voom.results, contrast,
        logFC.cutoff, fdr.cutoff){
    # to be changed for gene symbols support
    stopifnot(class(voom.results) == "EList")
    print("limma DE analysis is carried out ... ")
    # fit linear model for each gene using limma package functions
    vfit = lmFit(voom.results, design=voom.results$design) # Fit linear model
# for each gene given a series of arrays
    if (is.matrix(contrast)){
        vfit = contrasts.fit(vfit, contrast) # make all pair-wise comparisons
        contr.names = colnames(contrast)
        contr.num =  ncol(contrast)
        coefs = 1:ncol(contrast)
# between the groups
    }else{
        contr.names = names(contrast)
        contr.num = length(contrast)
        coefs = contrast
    }
    ebayes.results = eBayes(vfit) # compute moderated t-statistics, moderated
#F-statistic, and log-odds of differential expression by empirical
#    Bayes moderation of the standard errors towards a common value
    logFC = matrix(0, nrow(ebayes.results), contr.num)
    limma.tops = list()
    for (i in 1:length(coefs)){
        top.table = topTable(ebayes.results, coef=coefs[i],
                number=Inf, sort.by="none")
        limma.fc = top.table$logFC
        names(limma.fc) = rownames(ebayes.results)
        logFC[, i] = limma.fc
        rownames(top.table) = rownames(ebayes.results)
        limma.tops[[contr.names[i]]] = top.table
        de.genes = top.table[top.table[, "adj.P.Val"] <= fdr.cutoff, ]
        if (nrow(de.genes) == 0)
            cat(paste0("WARNING: it seems the contrast ",
                    contr.names[i],
                    " has no DE genes at the selected 'fdr.cutoff'.\n",
                    "The 'fdr.cutoff' was ignored in the calculations.\n"))
        if (nrow(de.genes) > 0){
            de.genes = de.genes[abs(de.genes[, "logFC"]) >= logFC.cutoff, ]
            if (nrow(de.genes) == 0)
                cat(paste0("WARNING: it seems the contrast ",
                    contr.names[i],
                    " has no DE genes at the selected 'logFC.cutoff'.\n",
                    "The 'logFC.cutoff' was ignored in the calculations.\n"))
        }
    }

    rownames(logFC) = rownames(ebayes.results)
    colnames(logFC) = contr.names

    return(list(logFC=logFC, limma.results=ebayes.results, limma.tops=limma.tops))
}

signifCal = function(combiPval, avgLFC){
  # combiPval : combined Pvalue after correction (correction  : BH, combination : fisher)
  #avgLFC : average logFC in geneset from edgerR/limma
  signif = -log10(combiPval)*abs(avgLFC)
  return(signif)
}

geneSetsLogFC = function(logFcMatrix, genesetCollection, condition){
  geneSetCollectionLogFC = c()
  for (geneset in names(genesetCollection)){
    setLogFC = c()
    for (gene in genesetCollection[[geneset]]){
      setLogFC = append(setLogFC, x = logFcMatrix$limma.tops[[condition]]$logFC[which(rownames(logFcMatrix$limma.tops[[condition]]) == gene)])
    }
    avgLogFC = mean(setLogFC)
    geneSetCollectionLogFC = append(geneSetCollectionLogFC, avgLogFC)
  }
  names(geneSetCollectionLogFC) = names(genesetCollection)
  return(geneSetCollectionLogFC)
}

adjustPval =  function(pvalTable, type){
  switch(type,
    holm = apply(pvalTable, 2, p.adjust, method = "holm"),
    hochberg = apply(pvalTable, 2, p.adjust, method = "hochberg"),
    hommel = apply(pvalTable, 2, p.adjust, method = "hommel"),
    bonferroni = apply(pvalTable, 2, p.adjust, method = "bonferroni"),
    BH = apply(pvalTable, 2, p.adjust, method = "BH"),
    BY = apply(pvalTable, 2, p.adjust, method = "BY"),
    fdr = apply(pvalTable, 2, p.adjust, method = "fdr"),
    none = apply(pvalTable, 2, p.adjust, method = "none"))
}

combinePval =  function(pvalAdjTable, type){
  if(!require(metap)){
    install.packages("metap")
    require(metap)
  }

  switch(type,
    sumz = apply(pvalAdjTable, 1, sumz),
    votep = apply(pvalAdjTable, 1, votep),
    minimump = apply(pvalAdjTable, 1, minimump),
    sumlog = apply(pvalAdjTable, 1, sumlog),
    sump = apply(pvalAdjTable, 1, sump),
    logitp = apply(pvalAdjTable, 1, logitp),
    meanp = apply(pvalAdjTable, 1, meanp),
    maximump = apply(pvalAdjTable, 1, maximump))
}

comparisonSummaryData = function(preparedDataResumPlot2){
  m = do.call(cbind, preparedDataResumPlot2)
  colnames(m) = rep(c("x.data","y.data","dir","rank","sig","id","gsSize"),length(names(preparedDataResumPlot2)))
  m = na.omit(m)
  m.id = base::abbreviate(rownames(m), minlength = 6, use.classes = FALSE, dot = FALSE, named = FALSE)
  m = m[,-grep("id", colnames(m))]
  colnames(m) = rep(c("x.data","y.data","dir","rank","sig","gsSize"),length(names(preparedDataResumPlot2)))
  nms = colnames(m)
  m2 = sapply(unique(nms), function(x) rowMeans (m[, nms == x]))
  m2 = as.data.frame(m2)
  m2 = cbind(m2,m.id)
  m2[,c("x.data","y.data","dir","rank","sig","gsSize")] = apply(m2[,c("x.data","y.data","dir","rank","sig","gsSize")],2,as.numeric)
  colnames(m2)[ncol(m2)] = "id"
  return(m2)
}

combineRanks = function(resTable){
    rankCol = grep("_Rank", colnames(resTable))
    meanRank = apply(resTable[,rankCol],1,mean)
    return(cbind(resTable, as.numeric(meanRank)))
}


prepareData = function(cGSEAcoreOutput, alpha = 0.05, directoryPath, pvalAdjMethod = "BH", pvalCombMethod = "sumlog", min.intersection.size = 1, shinyMode = FALSE){

  if(!require(SuperExactTest)){
    install.packages("SuperExactTest")
    require(SuperExactTest)
  }

  logFCTable = getlogFCFromLMFit(voom.results = cGSEAcoreOutput$Elist, contrast = cGSEAcoreOutput$contrast, logFC.cutoff = 0, fdr.cutoff = 1)

  time = unlist(cGSEAcoreOutput$time)
  names(time) = names(cGSEAcoreOutput$time)


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



  for (condi in colnames(cGSEAcoreOutput$contrast)){
    resTable = cGSEAOutputTable(cGSEAcoreOutput = cGSEAcoreOutput, contrastLevel = condi, geneSetLevels = names(cGSEAcoreOutput$collection))
    resTable = na.omit(resTable)
    # pvalue adjustment and combination
    pvalcol = grep("p.value", colnames(resTable))
    pvalTable = resTable[,pvalcol]
    #pvalAdjTable = apply(pvalTable, 2, p.adjust, method = "BH") #Benjamini Hochberg for pvalue adjustment
    pvalAdjTable = adjustPval(pvalTable, type = pvalAdjMethod)
    colnames(pvalAdjTable) = gsub("_p.value","_adj_p.value",colnames(pvalAdjTable))
    # pvalComb = apply(pvalAdjTable, 1, sumlog) #fisher method for pvalue combination
    pvalComb = combinePval(pvalAdjTable, type = pvalCombMethod)
    pvalComb = unlist(lapply(pvalComb, '[[',3))

    #result Table
    resTable2 = cbind(resTable, pvalAdjTable)
    resTable3 = cbind(resTable2, pvalComb)
    colnames(resTable3)[ncol(resTable3)] = "combined_p.value"
    result[[condi]] = na.omit(resTable3)
    result[[condi]] = combineRanks(result[[condi]])
    colnames(result[[condi]])[ncol(result[[condi]])] = "Avg_Rank"
    if (shinyMode == FALSE){
      write.csv(result[[condi]], file = paste(directoryPath,"/result_",condi,".csv", sep = ""))
    }
    #ranks
    rankCol = grep("Rank",colnames(resTable))
    rankTable = resTable[,rankCol]
    colnames(rankTable) = gsub("_Rank","", colnames(rankTable))

    #data for pVal binary heatmap
    heatmap[[condi]] = t(apply(pvalAdjTable,1, function(x) ifelse(x < alpha, 1, 0)))
    colnames(heatmap[[condi]]) = gsub("_adj_p.value","",colnames(heatmap[[condi]]))
    this.apply = apply(heatmap[[condi]], 1, sum)
    max_methods = max(this.apply)
    heatmap[[condi]] = heatmap[[condi]][c(which(this.apply == max_methods),which(this.apply == (max_methods-1))),]

    #data for PCA
    PCA[[condi]] = t(rankTable)

    #data for clustering
    clustering[[condi]] = hclust(dist(t(rankTable), method = "euclidean"), method = "ward.D2")

    #data for correlation plot
    correlation[[condi]] = cor(rankTable)

    #for SnailPlot
    gsFactor = as.factor(names(cGSEAcoreOutput$collection))
    listInput = apply(pvalAdjTable, 2, function(x) names(which(x < alpha)))
    cond = sapply(listInput, function(x) length(x) > 0)
    listInput = listInput[cond]
    listInput = sapply(listInput, function(x) factor(x, levels = levels(gsFactor)))
    names(listInput) = gsub("_adj_p.value","", names(listInput))
    snailPlot[[condi]] = listInput

    abrev_names_vector = base::abbreviate(gsub("HALLMARK_","H_",names(cGSEAcoreOutput$collection)), minlength = 6, use.classes = FALSE, dot = FALSE, named = FALSE, method = "left.kept")
    abrev_names_mat = cbind(cGSEAcoreOutput$collection,abrev_names_vector)
    colnames(abrev_names_mat) = c("original","abbreviation")
    abrev[[condi]] = t(as.data.frame(abrev_names_mat[,2]))
    abrev[[condi]] = cbind(rownames(abrev[[condi]]), abrev[[condi]])
    colnames(abrev[[condi]]) = c("Gene Set Name","abbreviation")
    if (shinyMode == FALSE){
      write.csv(abrev[[condi]], file = paste(directoryPath,"/abrreviations_",condi,".csv", sep = ""))
    }

    #for resumPlot1 - simple summary plot
    gsLogFC = geneSetsLogFC(logFcMatrix = logFCTable, genesetCollection = cGSEAcoreOutput$collection, condition = condi)
    resumPlot1[[condi]]$x = -log10(pvalComb)
    resumPlot1[[condi]]$y = gsLogFC
    resumPlot1[[condi]]$labels = abrev_names_vector

    #for resumPlot2 - advanced summary plot
    signifScore = signifCal(combiPval = pvalComb, avgLFC = gsLogFC)
    avgRank = apply(resTable[,rankCol],1,mean)

    resumPlot2[[condi]] = cbind(as.numeric(-log10(pvalComb)), abs(as.numeric(gsLogFC)), as.numeric(gsLogFC/abs(gsLogFC)), as.numeric(avgRank), as.numeric(signifScore), abrev_names_vector, unlist(lapply(cGSEAcoreOutput$collection, length)))
    colnames(resumPlot2[[condi]]) = c("x.data","y.data","dir","rank","sig","id", "gsSize")
    rownames(resumPlot2[[condi]]) = names(cGSEAcoreOutput$collection)
    resumPlot2[[condi]] = as.data.frame(resumPlot2[[condi]])
    resumPlot2[[condi]][,c("x.data","y.data","dir","rank","sig", "gsSize")] = apply(resumPlot2[[condi]][,c("x.data","y.data","dir","rank","sig", "gsSize")],2, as.numeric)

    # resumPlot2Table[[condi]] = cbind(as.numeric(-log10(pvalComb)), abs(as.numeric(gsLogFC)), as.numeric(gsLogFC/abs(gsLogFC)), as.numeric(avgRank), as.numeric(signifScore), abrev_names_vector, unlist(lapply(cGSEAcoreOutput$collection, length)))
    # colnames(resumPlot2Table[[condi]]) = c("x.data","y.data","dir","rank","sig","id", "gsSize")
    # rownames(resumPlot2Table[[condi]]) = names(cGSEAcoreOutput$collection)
    # resumPlot2Table[[condi]] = as.data.frame(resumPlot2[[condi]])
    # resumPlot2Table[[condi]][,c("x.data","y.data","dir","rank","sig", "gsSize")] = apply(resumPlot2[[condi]][,c("x.data","y.data","dir","rank","sig", "gsSize")],2, as.numeric)
  }
  if (length(colnames(cGSEAcoreOutput$contrast)) > 1){
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
