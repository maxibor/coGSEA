####################################################
################### camera #########################
####################################################

runcamera <- function(voom.results, contrast, genesetIdx,
                    num.workers=4, verbose = TRUE){

	if(!(require(limma))){
		source("http://www.bioconductor.org/biocLite.R")
		biocLite("limma")
		require(limma)
	}

    # run CAMERA and write out ranked 'gene sets' for each 'contrast'
    if (is.matrix(contrast)){
        contr.names = colnames(contrast)
        contr.num = ncol(contrast)
    }else{
        contr.names = names(contrast)
        contr.num = length(contrast)
    }
#    camera.results = vector("list", ncol(contrast))
    args.all = list()
    for(i in 1:contr.num){
        if (is.matrix(contrast))
            args.all[[i]] = list(voom.results = voom.results,
                contrast.name = contr.names[i],
                contrast = contrast[,i],
                genesetIdx = genesetIdx,
                verbose = verbose
            )
        else
            args.all[[i]] = list(voom.results = voom.results,
                contrast.name = contr.names[i],
                contrast = contrast[i],
                genesetIdx = genesetIdx,
                verbose = verbose
            )
#       print(args.all[[i]])
    }
    names(args.all) = contr.names
    if (Sys.info()['sysname'] == "Windows" || contr.num <= 1)
        camera.results = lapply(args.all, runcamera.contrast)
    else
        camera.results = mclapply(args.all, runcamera.contrast,
                mc.cores=num.workers)
#    names(camera.results) = colnames(contrast)
    return(camera.results)
}

runcamera.contrast <- function(args){
    if (args$verbose)
        print(paste0("   Running CAMERA for ", args$contrast.name))
    else
        cat(".")
    camera.results = camera(y=args$voom.results,
            index=args$genesetIdx,
            design=args$voom.results$design,
            contrast=args$contrast) # , allow.neg.cor=TRUE, inter.gene.cor=NA
#       print(head(camera.results[[i]]))

    camera.results =
            camera.results[order(camera.results[,"PValue"]),]

    camera.results = cbind(Rank=seq(1,
                    nrow(camera.results)), camera.results)
    colnames(camera.results)[which(colnames(camera.results) == "PValue")] = "p.value"
    return(camera.results)
}

####################################################
##################### fry ##########################
####################################################

runfry <- function(voom.results, contrast, genesetIdx,
        num.workers=4, verbose = TRUE){

	if(!(require(limma))){
		source("http://www.bioconductor.org/biocLite.R")
		biocLite("limma")
		require(limma)
	}
    # run fry and write out ranked 'gene sets' for each 'contrast'
    if (is.matrix(contrast)){
        contr.names = colnames(contrast)
        contr.num = ncol(contrast)
    }else{
        contr.names = names(contrast)
        contr.num = length(contrast)
    }

#    fry.results = vector("list", ncol(contrast))
    args.all = list()
    for(i in 1:contr.num){
        if (is.matrix(contrast))
            args.all[[i]] = list(voom.results = voom.results,
                    contrast.name = contr.names[i],
                    contrast = contrast[,i],
                    genesetIdx = genesetIdx,
                    verbose = verbose
            )
        else
            args.all[[i]] = list(voom.results = voom.results,
                    contrast.name = contr.names[i],
                    contrast = contrast[i],
                    genesetIdx = genesetIdx,
                    verbose = verbose
            )
    }
    names(args.all) = contr.names
    if (Sys.info()['sysname'] == "Windows" || contr.num <= 1)
        fry.results = lapply(args.all, runfry.contrast)
    else
        fry.results = mclapply(args.all, runfry.contrast,
                mc.cores=num.workers)
#    names(fry.results) = colnames(contrast)
    return(fry.results)
}

runfry.contrast <- function(args){
    if (args$verbose)
        print(paste0("   Running FRY for ", args$contrast.name))
    else
        cat(".")
    capture.output(fry.results <- fry(y=args$voom.results,
                    index=args$genesetIdx,
                    design=args$voom.results$design,
                    contrast=args$contrast))
    # returns PropDown/PropUp ==> proportion of genes that are
# down/up-regulated
    fry.results = fry.results[order(fry.results[,
                            "PValue"]),]
    fry.results = cbind(Rank=seq(1,
                    nrow(fry.results)), fry.results)
    colnames(fry.results)[which(colnames(fry.results) == "PValue")] = "p.value"
    return(fry.results)
}

####################################################
##################### gage #########################
####################################################

rungage <- function(voom.results, contrast, genesetIdx,
        num.workers=4, verbose = TRUE){
	if(!require(gage)){
		source("https://bioconductor.org/biocLite.R")
		biocLite("gage")
		require(gage)
	}
    # run GAGE and write out ranked 'gene sets' for each 'contrast'
    if (is.matrix(contrast)){
        contr.names = colnames(contrast)
        contr.num = ncol(contrast)
    }else{
        contr.names = names(contrast)
        contr.num = length(contrast)
    }

    groupData = prepareTwoGroupsData(voom.results, contrast, genesetIdx,
            min.samples = 3, verbose)
    #gage.results = vector("list", ncol(contrast))
    args.all = list()
    for(i in 1:contr.num){
        args.all[[i]] = list(contrast = contr.names[i],
                logCPM = groupData$data[[i]]$logCPM,
                group1 = groupData$data[[i]]$group1,
                group2 = groupData$data[[i]]$group2,
                gsets = groupData$gsets,
                genesetIdx = genesetIdx,
                verbose = verbose)
    }
    names(args.all) = contr.names
    if (Sys.info()['sysname'] == "Windows" || contr.num <= 1)
        gage.results = lapply(args.all, rungage.contrast)
    else
        gage.results = mclapply(args.all, rungage.contrast,
                mc.cores=num.workers)
    #stop("here")
    names(gage.results) = colnames(contrast)
    for (i in length(names(gage.results))){
      gage.results[[i]][which(gage.results[[i]]$set.size == 1),"p.value"] = 1
    }
    return(gage.results)
    # return(tmp)
}

rungage.contrast <- function(args){
    if (args$verbose)
        print(paste0("   Running GAGE for ", args$contrast))
    else
        cat(".")
    groupData = args$groupData
    # same.dir=FALSE ==> Two directional test
    gage.results = gage(exprs=args$logCPM,
            gsets=args$gsets,
            ref=args$group1,
            samp=args$group2, set.size = c(0,2000), same.dir=FALSE,
            compare = "unpaired")$greater[, 1:5]
    # returns PropDown/PropUp ==> proportion of genes that are
# down/up-regulated
    gage.results = gage.results[ order (
                    gage.results[,"p.val"]),]
    gage.results = cbind(Rank=seq(1, nrow(gage.results)),
            gage.results)
    colnames(gage.results)[which(colnames(gage.results) == "p.val")] = "p.value"
    return(as.data.frame(gage.results))
}

####################################################
################ globaltest ########################
####################################################

runglobaltest <- function(voom.results, contrast, genesetIdx,
        num.workers=4, verbose = TRUE){

	if(!require(globaltest)){
		source("https://bioconductor.org/biocLite.R")
		biocLite("globaltest")
		require(globaltest)
	}

    # run globaltest and write out ranked 'gene sets' for each 'contrast'
    if (is.matrix(contrast)){
        contr.names = colnames(contrast)
        contr.num = ncol(contrast)
    }else{
        contr.names = names(contrast)
        contr.num = length(contrast)
    }

    gt.options(transpose=TRUE)
    groupData = prepareTwoGroupsData(voom.results, contrast, genesetIdx,
            verbose = verbose)
#    globaltest.results = vector("list", ncol(contrast))
    args.all = list()
    for(i in 1:contr.num){
        args.all[[i]] = list(contrast = contr.names[i],
                logCPM = groupData$data[[i]]$logCPM,
                group1 = groupData$data[[i]]$group1,
                group2 = groupData$data[[i]]$group2,
                gsets = groupData$gsets,
                genesetIdx = genesetIdx,
                verbose = verbose)
    }
    names(args.all) = contr.names
    if (Sys.info()['sysname'] == "Windows" || contr.num <= 1)
        globaltest.results = lapply(args.all, runglobaltest.contrast)
    else
        globaltest.results = mclapply(args.all, runglobaltest.contrast,
                mc.cores=num.workers)
#    names(globaltest.results) = colnames(contrast)
    return(globaltest.results)
}

runglobaltest.contrast <- function(args){
    if (args$verbose)
        print(paste0("   Running GLOBALTEST for ",
                        args$contrast))
    else
        cat(".")

    groupIndx = c(args$group2,
            args$group1)
    data.log.sel = args$logCPM[, groupIndx]
    group = c(rep("d", length(args$group2)),
            rep("c", length(args$group1)))
    colnames(data.log.sel) = group
    # perform the globaltest
    globaltest.results = result(gt(factor(group),
                    data.log.sel, subsets=args$genesetIdx,
                    permutations=10000))

    globaltest.results = globaltest.results[
            order(globaltest.results[,"p-value"],
                    -globaltest.results[,"Statistic"]),]
    globaltest.results = cbind(
            Rank=seq(1, nrow(globaltest.results)),
            globaltest.results)
    colnames(globaltest.results)[which(
                    colnames(globaltest.results) == "p-value")] = "p.value"
    return(globaltest.results)
}

####################################################
################ GSVA ##############################
####################################################

rungsva <- function(method, voom.results, contrast, genesetIdx,
        num.workers=4, verbose = TRUE){
	if(!require(GSVA)){
		source("https://bioconductor.org/biocLite.R")
		biocLite("GSVA")
		require(GSVA)
	}

    # run gsva and write out ranked 'gene sets' for each 'contrast'
    if (is.matrix(contrast)){
        contr.names = colnames(contrast)
        contr.num = ncol(contrast)
    }else{
        contr.names = names(contrast)
        contr.num = length(contrast)
    }

    gsets = list()
    for (j in 1:length(genesetIdx)){
        gsets[[j]] = as.character(genesetIdx[[j]])
    }
    names(gsets) = names(genesetIdx)

    if (verbose)
        print(paste0("   Calculating gene set-level stats using ",
                        toupper(method)))
    else
        cat(".")
    # transform scores in gene set space using parallel computing
    data.log = voom.results$E
    rownames(data.log) = as.character(seq(1, nrow(data.log)))
    gs.es = calculateSetScores.parallel(data.log, gsets, method, num.workers)
    # fit the gene set scores and find DE gene sets
    gs.fit = lmFit(gs.es, design=voom.results$design)
    if (is.matrix(contrast)){
        gs.fit = contrasts.fit(gs.fit, contrast)
        coefs = 1:contr.num
    }else{
        coefs = contrast
    }
    gs.fit = eBayes(gs.fit)
    gsva.results = vector("list", contr.num)
    for(i in 1:contr.num){
        if (verbose)
            print(paste0("   Running ", toupper(method)," for ",
                            contr.names[i]))
        else
            cat(".")
        gsva.results[[i]] =  topTable(gs.fit, coef=coefs[i], number=Inf, sort.by="p",
                adjust.method="BH")
        gsva.results[[i]] = gsva.results[[i]][order(gsva.results[[i]][, "P.Value"],
                        -gsva.results[[i]][, "B"]),]
        gsva.results[[i]] = cbind(Rank=seq(1, nrow(gsva.results[[i]])),
                gsva.results[[i]])

        colnames(gsva.results[[i]])[which(
                        colnames(gsva.results[[i]]) == "P.Value")] = "p.value"
    }
    names(gsva.results) = contr.names
    return(gsva.results)
}

calculateSetScores.parallel <- function(data.log, gsets, method, num.workers){
    args.all = list()
    sets.per.task = 50
    total.tasks = ceiling(length(gsets) /  sets.per.task)
    for (i in 1:total.tasks){
        gsetsi = gsets[((i-1)*sets.per.task + 1):(i*sets.per.task)]
        gsetsi = gsetsi[!sapply(gsetsi, is.null)]
        args.all[[paste0("task", i)]] = list(
                data.log = data.log,
                gsets = gsetsi,
                method = method
        )
    }
    # parallelize the calculation of gene set scores to speed up the algorithms
    if (Sys.info()['sysname'] == "Windows")
        gs.es.all = lapply(args.all, rungsva.subcollection)
    else
        gs.es.all = mclapply(args.all, rungsva.subcollection,
                mc.cores=num.workers)
    # collect gene set scores from differet workers
    gs.es = c()
    for (i in 1:total.tasks)
        gs.es = rbind(gs.es, gs.es.all[[paste0("task", i)]])
    rownames(gs.es) = names(gsets)
    return(gs.es)
}

rungsva.subcollection <- function(args){
    set.seed(519863)
    gs.es = gsva(expr=args$data.log, gset.idx.list=args$gsets, mx.diff=TRUE,
            min.sz=1,
            method=args$method, parallel.sz=1,
            verbose=FALSE, rnaseq=FALSE)#$es.obs

    if (args$method == "gsva"){
        gs.es = gs.es$es.obs
    }

    return(gs.es)
}

####################################################
################ ora ###############################
####################################################


runora <- function(voom.results, contrast, genesetIdx, num.workers=4, verbose = TRUE){
    # run hypergeomteric test and write out ranked 'gene sets' for each
# 'contrast'
    # The p-value you want is the probability of getting 100 white balls in
    # a sample of size 400 from an urn with 3000 white balls and 12000
# black balls.
    # phyper(100, 3000, 12000, 400)
    #

	if(!(require(limma))){
		source("http://www.bioconductor.org/biocLite.R")
		biocLite("limma")
		require(limma)
	}

    if (is.matrix(contrast)){
        contr.names = colnames(contrast)
        contr.num = ncol(contrast)
    }else{
        contr.names = names(contrast)
        contr.num = length(contrast)
    }

    if (!is.null(voom.results$E)){
        vfit = lmFit(voom.results, voom.results$design)
        if (is.matrix(contrast)){
            vfit = contrasts.fit(vfit, contrast)
            coefs = 1:contr.num
        }else{
            coefs = contrast
        }
        vfit = eBayes(vfit)
        pvalue.cut=0.05
        logfc.cut=1
        universe = unlist(genesetIdx)
    }
    else if (!is.null(voom.results$featureIDs)){ # ORA Analysis
        universe = voom.results$featureIDs
    }

    ora.results = vector("list", contr.num)
    for(i in 1:contr.num){
        if (verbose)
            print(paste0("   Running ORA for ", contr.names[i]))
        else
            cat(".")
        if (!is.null(voom.results$E)){
            deGenes = rownames(topTable(vfit, coef=coefs[i], number=Inf,
                        p.value=pvalue.cut,
                        lfc=logfc.cut))
            if (length(deGenes) == 0){
                deGenes = rownames(topTable(vfit, coef=i, number=Inf,
                                p.value=pvalue.cut,
                                lfc=0))
                if (verbose)
                    print("ORA used a cut-off logFC = 0")
            }

        }else{
            deGenes = voom.results$ids
        }
        deGenes = deGenes[which(deGenes %in% universe)]
        ora.stats = runora.collection(deGenes, genesetIdx, universe)
        ora.results[[i]] = ora.stats
#       print(head(ora.results[[i]]))
        ora.results[[i]] =
ora.results[[i]][order(ora.results[[i]][,"p.value"]),]   # order by p.value
        ora.results[[i]] = cbind(Rank=seq(1, nrow(ora.results[[i]])),
ora.results[[i]])

    }
    names(ora.results) = colnames(contrast)
    return(ora.results)
}


runora.collection <- function(deGenes, genesetIdx, universe){
    tmp = rep(NA, length(genesetIdx))
    ora.stats = data.frame(p.value=tmp, p.adj = tmp) # , NumGenes=tmp
    totalDE = length(deGenes)
    n = length(universe) - totalDE
    for (j in 1:length(genesetIdx)){
        gset = genesetIdx[[j]]
        totalDEinS = length(intersect(gset, deGenes))
        totalSinUniverse = length(intersect(gset, universe))
        ora.stats[j, "p.value"] = phyper(q = totalDEinS- 0.5, m=
totalDE, n = n,
                k = totalSinUniverse, lower.tail = FALSE)
        #ora.stats[j, "NumGenes"] = totalDEinS
    }
    ora.stats[, "p.adj"] = p.adjust(ora.stats[, "p.value"], method = "BH")
    row.names(ora.stats) = names(genesetIdx)
    return ( ora.stats)
}

#####################################################
################ PADOG ##############################
#####################################################

runpadog <- function(voom.results, contrast, genesetIdx,
        num.workers=4, verbose = TRUE){
	if(!require(PADOG)){
		source("http://www.bioconductor.org/biocLite.R")
		biocLite("PADOG")
		require(PADOG)
	}

    # run padog and write out ranked 'gene sets' for each 'contrast'
    if (is.matrix(contrast)){
        contr.names = colnames(contrast)
        contr.num = ncol(contrast)
    }else{
        contr.names = names(contrast)
        contr.num = length(contrast)
    }

    groupData = prepareTwoGroupsData(voom.results, contrast, genesetIdx,
            min.samples = 3, verbose)
#    padog.results = vector("list", ncol(contrast))
    args.all = list()
    for(i in 1:contr.num){
        args.all[[i]] = list(contrast = contr.names[i],
                logCPM = groupData$data[[i]]$logCPM,
                group1 = groupData$data[[i]]$group1,
                group2 = groupData$data[[i]]$group2,
                gsets = groupData$gsets,
                genesetIdx = genesetIdx,
                verbose = verbose)
    }
    names(args.all) = contr.names
    if (Sys.info()['sysname'] == "Windows" || contr.num <= 1)
        padog.results = lapply(args.all, runpadog.contrast)
    else
        padog.results = mclapply(args.all, runpadog.contrast,
                mc.cores=num.workers)
#    names(padog.results) = colnames(contrast)
    return(padog.results)
}

runpadog.contrast <- function(args){
    if (args$verbose)
        print(paste0("   Running PADOG for ", args$contrast))
    else
        cat(".")

    group = c(rep("c", length(args$group1)),
            rep("d", length(args$group2)))

    padog.results = padog(esetm=args$logCPM,
            group=group, paired=FALSE,
            gslist=args$gsets, NI=100, Nmin = 1,
            verbose=FALSE)

    padog.results = padog.results[order(
                    padog.results[,"Ppadog"],
                    -padog.results[,"padog0"]),]
    padog.results = cbind(
            Rank = seq(1, nrow(padog.results)),
            padog.results)
    colnames(padog.results)[which(colnames(padog.results) == "Ppadog")] = "p.value"
    return(padog.results)
}

#####################################################
################ roast ###############################
#####################################################

runroast <- function(voom.results, contrast, genesetIdx,
        num.workers=4, nrot = 999, verbose = TRUE){
	if(!(require(limma))){
		source("http://www.bioconductor.org/biocLite.R")
		biocLite("limma")
		require(limma)
	}

    # run ROAST and write out ranked 'gene sets' for each 'contrast'
    if (is.matrix(contrast)){
        contr.names = colnames(contrast)
        contr.num = ncol(contrast)
    }else{
        contr.names = names(contrast)
        contr.num = length(contrast)
    }

#    roast.results = vector("list", ncol(contrast))
    args.all = list()
    for(i in 1:contr.num){
        if (is.matrix(contrast))
            args.all[[i]] = list(voom.results = voom.results,
                    contrast.name = contr.names[i],
                    contrast = contrast[,i],
                    genesetIdx = genesetIdx,
                    verbose = verbose
            )
        else
            args.all[[i]] = list(voom.results = voom.results,
                    contrast.name = contr.names[i],
                    contrast = contrast[i],
                    genesetIdx = genesetIdx,
                    verbose = verbose
            )
    }
    names(args.all) = contr.names
    if (Sys.info()['sysname'] == "Windows" || contr.num <= 1)
        roast.results = lapply(args.all, runroast.contrast)
    else
        roast.results = mclapply(args.all, runroast.contrast,
                mc.cores=num.workers)
#    names(roast.results) = colnames(contrast)
    return(roast.results)
}

runroast.contrast <- function(args){
    if (args$verbose)
        print(paste0("   Running ROAST for ", args$contrast.name))
    else
        cat(".")
    roast.results = mroast(y=args$voom.results,
            index=args$genesetIdx,
            design=args$voom.results$design,
            contrast=args$contrast, nrot=999)
    # returns PropDown/PropUp ==> proportion of genes that are
# down/up-regulated
    roast.results = roast.results[order(roast.results[,
                            "PValue"]),]
    roast.results = cbind(Rank=seq(1,
                    nrow(roast.results)), roast.results)
    colnames(roast.results)[which(colnames(roast.results) == "PValue")] = "p.value"
    return(roast.results)
}

#####################################################
################ safe ###############################
#####################################################


runsafe <- function(voom.results, contrast, genesetIdx,
        num.workers=4, verbose = TRUE){
	if(!(require(safe))){
		source("http://www.bioconductor.org/biocLite.R")
		biocLite("safe")
		require(safe)
	}

    # run safe and write out ranked 'gene sets' for each 'contrast'
    if (is.matrix(contrast)){
        contr.names = colnames(contrast)
        contr.num = ncol(contrast)
    }else{
        contr.names = names(contrast)
        contr.num = length(contrast)
    }

    groupData = prepareTwoGroupsData(voom.results, contrast, genesetIdx,
            min.samples = 3, verbose)
    capture.output(C.mat <- getCmatrix(keyword.list=groupData$gsets, min.size = 1,
                    present.genes=rownames(groupData$data[[1]]$logCPM)))
#    safe.results = vector("list", ncol(contrast))
    args.all = list()
    for(i in 1:contr.num){
        args.all[[i]] = list(contrast = contr.names[i],
                C.mat = C.mat,
                logCPM = groupData$data[[i]]$logCPM,
                group1 = groupData$data[[i]]$group1,
                group2 = groupData$data[[i]]$group2,
                gsets = groupData$gsets,
                genesetIdx = genesetIdx,
                verbose = verbose)
    }
    names(args.all) = contr.names
    if (Sys.info()['sysname'] == "Windows" || contr.num <= 1)
        safe.results = lapply(args.all, runsafe.contrast)
    else
        safe.results = mclapply(args.all, runsafe.contrast,
                mc.cores=num.workers)
#    names(safe.results) = colnames(contrast)
    return(safe.results)
}

runsafe.contrast <- function(args){
    if (args$verbose)
        print(paste0("   Running SAFE for ", args$contrast))
    else
        cat(".")
    group = c(rep("Ctr", length(args$group1)),
            rep("Trt", length(args$group2)))

    capture.output(temp <- safe(X.mat=args$logCPM,
                    y.vec=group, C.mat=args$C.mat,
                    print.it=FALSE, Pi.mat=100))#,
    safe.results = safe.toptable(temp,
            number=length(args$gsets),
            description = FALSE)
    # safe.results = na.omit(safe.results)
    # rownames(safe.results) = names(args$genesetIdx)
    safe.results[, "P.value"] = as.numeric(
            safe.results[, "P.value"])
    safe.results[, "Adj.p.value"] = as.numeric(
            safe.results[, "Adj.p.value"])
    rownames(safe.results) = safe.results[, "GenesetID"]
    safe.results = safe.results[, -c(1,6)]
    safe.results = safe.results[order(
                    safe.results[,"P.value"],
                    -safe.results[,"Statistic"]),]
    safe.results = cbind(Rank = seq(1,
                    nrow(safe.results)), safe.results)
    colnames(safe.results)[which(colnames(safe.results) == "P.value")] = "p.value"
    return(safe.results)
}

#####################################################
################ SetRank ############################
#####################################################

# runsetrank = function(voom.results, contrast, geneset = "H", specie = "Mus musculus", collection = NULL, num.workers=4, verbose = TRUE, lfcCuttof = 0, pCutoff = 1, fdrCutoff = 1, directory = getwd()){
#   if(!(require(SetRank))){
#     install.packages("SetRank")
#     require(SetRank)
#   }
#
#   if(!(require(edgeR))){
#     install.packages("edgeR")
#     require(edgeR)
#   }
#   collection_names = c("maxSetSize","referenceSet","sets","g","bigSets","intersection.p.cutoff","intersections","iMatrix")
#
#   options(mc.cores = num.workers)
#
#   if (specie == "Mus musculus" && is.null(collection) == TRUE){
#     if (geneset == "H"){
#     #   load(paste(directory,"/genesets/setRankMusCollection_H.RData", sep = ""))
#       theCollection = MusCollectionH
#     } else if (geneset == "C2_KEGG"){
#     #   load(paste(directory,"/genesets/setRankMusCollection_C2.RData", sep = ""))
#       theCollection = MusCollectionC2Kegg
#     } else if (geneset == "C2_REACTOME"){
#       theCollection = MusCollectionC2Reactome
#     } else {
#       cat(paste(geneset , "is not a valid MSigDB geneset for SetRank. Available genesets : \n - H \n - C2"))
#       return(NULL)
#     }
#   } else if (specie == "Homo sapiens"){
#     if (geneset == "H"){
# 	#   load(paste(directory,"genesets/setRankHomoCollection_H.RData", sep = ""))
#       theCollection = HomoCollectionH
#     } else if (geneset == "C2_KEGG"){
#     #   load(paste(directory,"genesets/setRankHomoCollection_C2.RData", sep = ""))
#       theCollection = HomoCollectionC2Kegg
#     } else if (geneset == "C2_REACTOME"){
#       theCollection = HomoCollectionC2Reactome
#     } else {
#        cat(paste(geneset , "is not a valid MSigDB geneset for SetRank. Available genesets : \n - H \n - C2"))
#       return(NULL)
#     }
#   } else if ( is.null(collection) == FALSE){
#       if (names(collection) == collection_names){
#         theCollection = collection
#       } else {
#         cat(paste(collection,"is not a valid SetRank collection object"))
#       }
#
#     }
#
#   else {
#     cat(paste(specie,"is not a valid specie for SetRank. Available species : \n - Mus musculus \n - Homo sapiens"))
#     return(NULL)
#   }
#
#   dir.create(file.path(getwd(), "SetRankResults"), showWarnings = FALSE)
#
#   voomFit = lmFit(voom.results, voom.results$design)
#   voomContrast = contrasts.fit(voomFit, contrasts = contrast)
#   voomTreat = treat(voomContrast, lfc = lfcCuttof)
#
#   setrank.results = list()
#
#   for (i in seq(1,length(colnames(contrast)))){
#     if (verbose == TRUE){
#       cat(paste("	Running SetRank for", colnames(contrast)[i], "contrast"))
#     }
#
#     network = setRankAnalysis(rownames(topTreat(voomTreat, coef = i, n = Inf)),setCollection = theCollection, use.ranks = TRUE, setPCutoff = pCutoff, fdrCutoff = fdrCutoff)
#     exportSingleResult(network = network, selectedGenes = rownames(topTreat(voomTreat, coef = i, n = Inf)), collection = theCollection, networkName = colnames(contrast)[i], outputPath = "./SetRankResults/")
#     fh = read.table(paste("./SetRankResults/",colnames(contrast)[i],"_pathways.txt", sep = ""), header = TRUE)
#     fh$Rank = seq(1,dim(fh)[1])
#
#     fh2 = cbind(fh$Rank,fh$size,fh$setRank,fh$pSetRank,fh$correctedPValue,fh$adjustedPValue)
#     rownames(fh2) = fh$description
#     colnames(fh2) = c("Rank","Size", "SetRank_Rank", "p.value", "Adj.p.value", "fdr")
#     fh = fh2
#     remove(fh2)
#
#     setrank.results[[colnames(contrast)[i]]] = fh
#   }
#   return(setrank.results)
# }
