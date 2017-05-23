####################################################
################### camera #########################
####################################################

runcamera <- function(voom.results, contrast, genesetIdx,
                    num.workers=4, verbose = TRUE){

    # run CAMERA and write out ranked 'gene sets' for each 'contrast'
    if (base::is.matrix(contrast)){
        contr.names = base::colnames(contrast)
        contr.num = base::ncol(contrast)
    }else{
        contr.names = base::names(contrast)
        contr.num = base::length(contrast)
    }
#    camera.results = vector("list", ncol(contrast))
    args.all = list()
    for(i in 1:contr.num){
        if (base::is.matrix(contrast))
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
    base::names(args.all) = contr.names
    if (base::Sys.info()['sysname'] == "Windows" || contr.num <= 1)
        camera.results = base::lapply(args.all, runcamera.contrast)
    else
        camera.results = parallel::mclapply(args.all, runcamera.contrast,
                mc.cores=num.workers)
#    names(camera.results) = colnames(contrast)
    return(camera.results)
}

runcamera.contrast <- function(args){
    if (args$verbose)
        print(paste0("   Running CAMERA for ", args$contrast.name))
    else
        cat(".")
    camera.results = limma::camera(y=args$voom.results,
            index=args$genesetIdx,
            design=args$voom.results$design,
            contrast=args$contrast) # , allow.neg.cor=TRUE, inter.gene.cor=NA
#       print(head(camera.results[[i]]))

    camera.results =
            camera.results[base::order(camera.results[,"PValue"]),]

    camera.results = base::cbind(Rank=base::seq(1,
                    base::nrow(camera.results)), camera.results)
    base::colnames(camera.results)[base::which(base::colnames(camera.results) == "PValue")] = "p.value"
    return(camera.results)
}

####################################################
##################### fry ##########################
####################################################

runfry <- function(voom.results, contrast, genesetIdx,
        num.workers=4, verbose = TRUE){

    # run fry and write out ranked 'gene sets' for each 'contrast'
    if (base::is.matrix(contrast)){
        contr.names = base::colnames(contrast)
        contr.num = base::ncol(contrast)
    }else{
        contr.names = base::names(contrast)
        contr.num = base::length(contrast)
    }

#    fry.results = vector("list", ncol(contrast))
    args.all = list()
    for(i in 1:contr.num){
        if (base::is.matrix(contrast))
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
    base::names(args.all) = contr.names
    if (base::Sys.info()['sysname'] == "Windows" || contr.num <= 1)
        fry.results = base::lapply(args.all, runfry.contrast)
    else
        fry.results = parallel::mclapply(args.all, runfry.contrast,
                mc.cores=num.workers)
#    names(fry.results) = colnames(contrast)
    return(fry.results)
}

runfry.contrast <- function(args){
    if (args$verbose)
        print(paste0("   Running FRY for ", args$contrast.name))
    else
        cat(".")
    capture.output(fry.results <- limma::fry(y=args$voom.results,
                    index=args$genesetIdx,
                    design=args$voom.results$design,
                    contrast=args$contrast))
    # returns PropDown/PropUp ==> proportion of genes that are
# down/up-regulated
    fry.results = fry.results[base::order(fry.results[,
                            "PValue"]),]
    fry.results = base::cbind(Rank=base::seq(1,
                    base::nrow(fry.results)), fry.results)
    base::colnames(fry.results)[base::which(base::colnames(fry.results) == "PValue")] = "p.value"
    return(fry.results)
}

####################################################
##################### gage #########################
####################################################

rungage <- function(voom.results, contrast, genesetIdx,
        num.workers=4, verbose = TRUE){

    # run GAGE and write out ranked 'gene sets' for each 'contrast'
    if (base::is.matrix(contrast)){
        contr.names = base::colnames(contrast)
        contr.num = base::ncol(contrast)
    }else{
        contr.names = base::names(contrast)
        contr.num = base::length(contrast)
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
    base::names(args.all) = contr.names
    if (Sys.info()['sysname'] == "Windows" || contr.num <= 1)
        gage.results = base::lapply(args.all, rungage.contrast)
    else
        gage.results = parallel::mclapply(args.all, rungage.contrast,
                mc.cores=num.workers)
    #stop("here")
    base::names(gage.results) = base::colnames(contrast)
    for (i in base::length(base::names(gage.results))){
      gage.results[[i]][base::which(gage.results[[i]]$set.size == 1),"p.value"] = 1
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
    gage.results = gage::gage(exprs=args$logCPM,
            gsets=args$gsets,
            ref=args$group1,
            samp=args$group2, set.size = c(0,2000), same.dir=FALSE,
            compare = "unpaired")$greater[, 1:5]
    # returns PropDown/PropUp ==> proportion of genes that are
# down/up-regulated
    gage.results = gage.results[ order (
                    gage.results[,"p.val"]),]
    gage.results = base::cbind(Rank=base::seq(1, base::nrow(gage.results)),
            gage.results)
    base::colnames(gage.results)[base::which(base::colnames(gage.results) == "p.val")] = "p.value"
    return(base::as.data.frame(gage.results))
}

####################################################
################ globaltest ########################
####################################################

runglobaltest <- function(voom.results, contrast, genesetIdx,
        num.workers=4, verbose = TRUE){


    # run globaltest and write out ranked 'gene sets' for each 'contrast'
    if (base::is.matrix(contrast)){
        contr.names = base::colnames(contrast)
        contr.num = base::ncol(contrast)
    }else{
        contr.names = base::names(contrast)
        contr.num = base::length(contrast)
    }

    globaltest::gt.options(transpose=TRUE)
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
    base::names(args.all) = contr.names
    if (base::Sys.info()['sysname'] == "Windows" || contr.num <= 1)
        globaltest.results = base::lapply(args.all, runglobaltest.contrast)
    else
        globaltest.results = parallel::mclapply(args.all, runglobaltest.contrast,
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

    groupIndx = base::c(args$group2,
            args$group1)
    data.log.sel = args$logCPM[, groupIndx]
    group = base::c(base::rep("d", base::length(args$group2)),
            base::rep("c", base::length(args$group1)))
    base::colnames(data.log.sel) = group
    # perform the globaltest
    globaltest.results = globaltest::result(globaltest::gt(factor(group),
                    data.log.sel, subsets=args$genesetIdx,
                    permutations=10000))

    globaltest.results = globaltest.results[
            base::order(globaltest.results[,"p-value"],
                    -globaltest.results[,"Statistic"]),]
    globaltest.results = base::cbind(
            Rank=base::seq(1, base::nrow(globaltest.results)),
            globaltest.results)
    base::colnames(globaltest.results)[base::which(
                    base::colnames(globaltest.results) == "p-value")] = "p.value"
    return(globaltest.results)
}

####################################################
################ GSVA ##############################
####################################################

rungsva <- function(method, voom.results, contrast, genesetIdx,
        num.workers=4, verbose = TRUE){


    # run gsva and write out ranked 'gene sets' for each 'contrast'
    if (base::is.matrix(contrast)){
        contr.names = base::colnames(contrast)
        contr.num = base::ncol(contrast)
    }else{
        contr.names = base::names(contrast)
        contr.num = base::length(contrast)
    }

    gsets = list()
    for (j in 1:base::length(genesetIdx)){
        gsets[[j]] = base::as.character(genesetIdx[[j]])
    }
    base::names(gsets) = base::names(genesetIdx)

    if (verbose)
        print(paste0("   Calculating gene set-level stats using ",
                        base::toupper(method)))
    else
        cat(".")
    # transform scores in gene set space using parallel computing
    data.log = voom.results$E
    base::rownames(data.log) = base::as.character(base::seq(1, base::nrow(data.log)))
    gs.es = calculateSetScores.parallel(data.log, gsets, method, num.workers)
    # fit the gene set scores and find DE gene sets
    gs.fit = limma::lmFit(gs.es, design=voom.results$design)
    if (base::is.matrix(contrast)){
        gs.fit = limma::contrasts.fit(gs.fit, contrast)
        coefs = 1:contr.num
    }else{
        coefs = contrast
    }
    gs.fit = limma::eBayes(gs.fit)
    gsva.results = base::vector("list", contr.num)
    for(i in 1:contr.num){
        if (verbose)
            print(paste0("   Running ", base::toupper(method)," for ",
                            contr.names[i]))
        else
            cat(".")
        gsva.results[[i]] =  limma::topTable(gs.fit, coef=coefs[i], number=Inf, sort.by="p",
                adjust.method="BH")
        gsva.results[[i]] = gsva.results[[i]][base::order(gsva.results[[i]][, "P.Value"],
                        -gsva.results[[i]][, "B"]),]
        gsva.results[[i]] = base::cbind(Rank=base::seq(1, base::nrow(gsva.results[[i]])),
                gsva.results[[i]])

        base::colnames(gsva.results[[i]])[base::which(
                        base::colnames(gsva.results[[i]]) == "P.Value")] = "p.value"
    }
    base::names(gsva.results) = contr.names
    return(gsva.results)
}

calculateSetScores.parallel <- function(data.log, gsets, method, num.workers){
    args.all = list()
    sets.per.task = 50
    total.tasks = base::ceiling(length(gsets) /  sets.per.task)
    for (i in 1:total.tasks){
        gsetsi = gsets[((i-1)*sets.per.task + 1):(i*sets.per.task)]
        gsetsi = gsetsi[!base::sapply(gsetsi, is.null)]
        args.all[[paste0("task", i)]] = list(
                data.log = data.log,
                gsets = gsetsi,
                method = method
        )
    }
    # parallelize the calculation of gene set scores to speed up the algorithms
    if (base::Sys.info()['sysname'] == "Windows")
        gs.es.all = base::lapply(args.all, rungsva.subcollection)
    else
        gs.es.all = parallel::mclapply(args.all, rungsva.subcollection,
                mc.cores=num.workers)
    # collect gene set scores from differet workers
    gs.es = c()
    for (i in 1:total.tasks)
        gs.es = base::rbind(gs.es, gs.es.all[[paste0("task", i)]])
    base::rownames(gs.es) = base::names(gsets)
    return(gs.es)
}

rungsva.subcollection <- function(args){
    base::set.seed(519863)
    gs.es = GSVA::gsva(expr=args$data.log, gset.idx.list=args$gsets, mx.diff=TRUE,
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

    if (base::is.matrix(contrast)){
        contr.names = base::colnames(contrast)
        contr.num = base::ncol(contrast)
    }else{
        contr.names = base::names(contrast)
        contr.num = base::length(contrast)
    }

    if (!base::is.null(voom.results$E)){
        vfit = limma::lmFit(voom.results, voom.results$design)
        if (base::is.matrix(contrast)){
            vfit = limma::contrasts.fit(vfit, contrast)
            coefs = 1:contr.num
        }else{
            coefs = contrast
        }
        vfit = limma::eBayes(vfit)
        pvalue.cut=0.05
        logfc.cut=1
        universe = base::unlist(genesetIdx)
    }
    else if (!base::is.null(voom.results$featureIDs)){ # ORA Analysis
        universe = voom.results$featureIDs
    }

    ora.results = base::vector("list", contr.num)
    for(i in 1:contr.num){
        if (verbose)
            print(paste0("   Running ORA for ", contr.names[i]))
        else
            cat(".")
        if (!base::is.null(voom.results$E)){
            deGenes = base::rownames(limma::topTable(vfit, coef=coefs[i], number=Inf,
                        p.value=pvalue.cut,
                        lfc=logfc.cut))
            if (base::length(deGenes) == 0){
                deGenes = base::rownames(limma::topTable(vfit, coef=i, number=Inf,
                                p.value=pvalue.cut,
                                lfc=0))
                if (verbose)
                    print("ORA used a cut-off logFC = 0")
            }

        }else{
            deGenes = voom.results$ids
        }
        deGenes = deGenes[base::which(deGenes %in% universe)]
        ora.stats = runora.collection(deGenes, genesetIdx, universe)
        ora.results[[i]] = ora.stats
#       print(head(ora.results[[i]]))
        ora.results[[i]] =
ora.results[[i]][base::order(ora.results[[i]][,"p.value"]),]   # order by p.value
        ora.results[[i]] = base::cbind(Rank=base::seq(1, base::nrow(ora.results[[i]])),
ora.results[[i]])

    }
    base::names(ora.results) = base::colnames(contrast)
    return(ora.results)
}


runora.collection <- function(deGenes, genesetIdx, universe){
    tmp = base::rep(NA, base::length(genesetIdx))
    ora.stats = data.frame(p.value=tmp, p.adj = tmp) # , NumGenes=tmp
    totalDE = base::length(deGenes)
    n = base::length(universe) - totalDE
    for (j in 1:base::length(genesetIdx)){
        gset = genesetIdx[[j]]
        totalDEinS = base::length(base::intersect(gset, deGenes))
        totalSinUniverse = base::length(base::intersect(gset, universe))
        ora.stats[j, "p.value"] = stats::phyper(q = totalDEinS- 0.5, m=
totalDE, n = n,
                k = totalSinUniverse, lower.tail = FALSE)
        #ora.stats[j, "NumGenes"] = totalDEinS
    }
    ora.stats[, "p.adj"] = stats::p.adjust(ora.stats[, "p.value"], method = "BH")
    base::row.names(ora.stats) = base::names(genesetIdx)
    return ( ora.stats)
}

#####################################################
################ PADOG ##############################
#####################################################

runpadog <- function(voom.results, contrast, genesetIdx,
        num.workers=4, verbose = TRUE){

    # run padog and write out ranked 'gene sets' for each 'contrast'
    if (base::is.matrix(contrast)){
        contr.names = base::colnames(contrast)
        contr.num = base::ncol(contrast)
    }else{
        contr.names = base::names(contrast)
        contr.num = base::length(contrast)
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
    base::names(args.all) = contr.names
    if (base::Sys.info()['sysname'] == "Windows" || contr.num <= 1)
        padog.results = base::lapply(args.all, runpadog.contrast)
    else
        padog.results = parallel::mclapply(args.all, runpadog.contrast,
                mc.cores=num.workers)
#    names(padog.results) = colnames(contrast)
    return(padog.results)
}

runpadog.contrast <- function(args){
    if (args$verbose)
        print(paste0("   Running PADOG for ", args$contrast))
    else
        cat(".")

    group = base::c(base::rep("c", base::length(args$group1)),
            base::rep("d", base::length(args$group2)))

    padog.results = PADOG::padog(esetm=args$logCPM,
            group=group, paired=FALSE,
            gslist=args$gsets, NI=100, Nmin = 1,
            verbose=FALSE)

    padog.results = padog.results[base::order(
                    padog.results[,"Ppadog"],
                    -padog.results[,"padog0"]),]
    padog.results = base::cbind(
            Rank = base::seq(1, base::nrow(padog.results)),
            padog.results)
    base::colnames(padog.results)[base::which(base::colnames(padog.results) == "Ppadog")] = "p.value"
    return(padog.results)
}

#####################################################
################ roast ###############################
#####################################################

runroast <- function(voom.results, contrast, genesetIdx,
        num.workers=4, nrot = 999, verbose = TRUE){


    # run ROAST and write out ranked 'gene sets' for each 'contrast'
    if (base::is.matrix(contrast)){
        contr.names = base::colnames(contrast)
        contr.num = base::ncol(contrast)
    }else{
        contr.names = base::names(contrast)
        contr.num = base::length(contrast)
    }

#    roast.results = vector("list", ncol(contrast))
    args.all = list()
    for(i in 1:contr.num){
        if (base::is.matrix(contrast))
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
    base::names(args.all) = contr.names
    if (base::Sys.info()['sysname'] == "Windows" || contr.num <= 1)
        roast.results = base::lapply(args.all, runroast.contrast)
    else
        roast.results = parallel::mclapply(args.all, runroast.contrast,
                mc.cores=num.workers)
#    names(roast.results) = colnames(contrast)
    return(roast.results)
}

runroast.contrast <- function(args){
    if (args$verbose)
        print(paste0("   Running ROAST for ", args$contrast.name))
    else
        cat(".")
    roast.results = limma::mroast(y=args$voom.results,
            index=args$genesetIdx,
            design=args$voom.results$design,
            contrast=args$contrast, nrot=999)
    # returns PropDown/PropUp ==> proportion of genes that are
# down/up-regulated
    roast.results = roast.results[base::order(roast.results[,
                            "PValue"]),]
    roast.results = base::cbind(Rank=base::seq(1,
                    base::nrow(roast.results)), roast.results)
    base::colnames(roast.results)[base::which(base::colnames(roast.results) == "PValue")] = "p.value"
    return(roast.results)
}

#####################################################
################ safe ###############################
#####################################################


runsafe <- function(voom.results, contrast, genesetIdx,
        num.workers=4, verbose = TRUE){

    # run safe and write out ranked 'gene sets' for each 'contrast'
    if (base::is.matrix(contrast)){
        contr.names = base::colnames(contrast)
        contr.num = base::ncol(contrast)
    }else{
        contr.names = base::names(contrast)
        contr.num = base::length(contrast)
    }

    groupData = prepareTwoGroupsData(voom.results, contrast, genesetIdx,
            min.samples = 3, verbose)
    utils::capture.output(C.mat <- safe::getCmatrix(keyword.list=groupData$gsets, min.size = 1,
                    present.genes=base::rownames(groupData$data[[1]]$logCPM)))
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
    base::names(args.all) = contr.names
    if (base::Sys.info()['sysname'] == "Windows" || contr.num <= 1)
        safe.results = base::lapply(args.all, runsafe.contrast)
    else
        safe.results = parallel::mclapply(args.all, runsafe.contrast,
                mc.cores=num.workers)
#    names(safe.results) = colnames(contrast)
    return(safe.results)
}

runsafe.contrast <- function(args){
    if (args$verbose)
        print(paste0("   Running SAFE for ", args$contrast))
    else
        cat(".")
    group = base::c(base::rep("Ctr", base::length(args$group1)),
            base::rep("Trt", base::length(args$group2)))

    utils::capture.output(temp <- safe::safe(X.mat=args$logCPM,
                    y.vec=group, C.mat=args$C.mat,
                    print.it=FALSE, Pi.mat=100))#,
    safe.results = safe::safe.toptable(temp,
            number=base::length(args$gsets),
            description = FALSE)
    # safe.results = na.omit(safe.results)
    # rownames(safe.results) = names(args$genesetIdx)
    safe.results[, "P.value"] = base::as.numeric(
            safe.results[, "P.value"])
    safe.results[, "Adj.p.value"] = base::as.numeric(
            safe.results[, "Adj.p.value"])
    base::rownames(safe.results) = safe.results[, "GenesetID"]
    safe.results = safe.results[, -c(1,6)]
    safe.results = safe.results[order(
                    safe.results[,"P.value"],
                    -safe.results[,"Statistic"]),]
    safe.results = base::cbind(Rank = base::seq(1,
                    base::nrow(safe.results)), safe.results)
    base::colnames(safe.results)[base::which(base::colnames(safe.results) == "P.value")] = "p.value"
    return(safe.results)
}

#####################################################
################ SetRank ############################
#####################################################

runsetrank = function(voom.results, contrast, geneset = "H", specie = "Mus musculus", collection = NULL, num.workers=4, verbose = TRUE, lfcCuttof = 0, pCutoff = 1, fdrCutoff = 1, directory = getwd()){

  collection_names = c("maxSetSize","referenceSet","sets","g","bigSets","intersection.p.cutoff","intersections","iMatrix")

  base::options(mc.cores = num.workers)

  if (specie == "Mus musculus" && base::is.null(collection) == TRUE){
    if (geneset == "H"){
    #   load(paste(directory,"/genesets/setRankMusCollection_H.RData", sep = ""))
      theCollection = MusCollectionH
    } else if (geneset == "C2_KEGG"){
    #   load(paste(directory,"/genesets/setRankMusCollection_C2.RData", sep = ""))
      theCollection = MusCollectionC2Kegg
    } else if (geneset == "C2_REACTOME"){
      theCollection = MusCollectionC2Reactome
    } else {
      cat(paste(geneset , "is not a valid MSigDB geneset for SetRank. Available genesets : \n - H \n - C2"))
      return(NULL)
    }
  } else if (specie == "Homo sapiens"){
    if (geneset == "H"){
	#   load(paste(directory,"genesets/setRankHomoCollection_H.RData", sep = ""))
      theCollection = HomoCollectionH
    } else if (geneset == "C2_KEGG"){
    #   load(paste(directory,"genesets/setRankHomoCollection_C2.RData", sep = ""))
      theCollection = HomoCollectionC2Kegg
    } else if (geneset == "C2_REACTOME"){
      theCollection = HomoCollectionC2Reactome
    } else {
       cat(paste(geneset , "is not a valid MSigDB geneset for SetRank. Available genesets : \n - H \n - C2"))
      return(NULL)
    }
  } else if ( base::is.null(collection) == FALSE){
      if (base::names(collection) == collection_names){
        theCollection = collection
      } else {
        cat(paste(collection,"is not a valid SetRank collection object"))
      }

    }

  else {
    cat(paste(specie,"is not a valid specie for SetRank. Available species : \n - Mus musculus \n - Homo sapiens"))
    return(NULL)
  }

  base::dir.create(file.path(base::getwd(), "SetRankResults"), showWarnings = FALSE)

  voomFit = limma::lmFit(voom.results, voom.results$design)
  voomContrast = limma::contrasts.fit(voomFit, contrasts = contrast)
  voomTreat = limma::treat(voomContrast, lfc = lfcCuttof)

  setrank.results = list()

  for (i in base::seq(1,base::length(base::colnames(contrast)))){
    if (verbose == TRUE){
      cat(paste("	Running SetRank for", base::colnames(contrast)[i], "contrast"))
    }

    network = SetRank::setRankAnalysis(base::rownames(limma::topTreat(voomTreat, coef = i, n = Inf)),setCollection = theCollection, use.ranks = TRUE, setPCutoff = pCutoff, fdrCutoff = fdrCutoff)
    SetRank::exportSingleResult(network = network, selectedGenes = base::rownames(limma::topTreat(voomTreat, coef = i, n = Inf)), collection = theCollection, networkName = colnames(contrast)[i], outputPath = "./SetRankResults/")
    fh = base::read.table(paste("./SetRankResults/",base::colnames(contrast)[i],"_pathways.txt", sep = ""), header = TRUE)
    fh$Rank = base::seq(1,base::dim(fh)[1])

    fh2 = base::cbind(fh$Rank,fh$size,fh$setRank,fh$pSetRank,fh$correctedPValue,fh$adjustedPValue)
    base::rownames(fh2) = fh$description
    base::colnames(fh2) = c("Rank","Size", "SetRank_Rank", "p.value", "Adj.p.value", "fdr")
    fh = fh2
    base::remove(fh2)

    setrank.results[[base::colnames(contrast)[i]]] = fh
  }
  return(setrank.results)
}
