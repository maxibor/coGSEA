if(!require(parallel)){
  install.packages("parallel")
  require(parallel)
}

cGSEA = function(ElistObject,
    contrastMatrix,
    geneSetCollection = "H",
    specie = "Mus musculus",
    ENTREZGenesIds,
    directoryPath,
    alpha = 0.05,
    pvalAdjMethod = "BH" ,
    pvalCombMethod = "sumlog",
    min.intersection.size = 1,
    GSEA.Methods = c("camera", "gage", "globaltest", "gsva", "ssgsea", "zscore", "plage", "ora", "padog", "roast", "safe"),
    num.workers = 4,
    shinyMode = FALSE){

  methods = method_checker(GSEA.Methods)
  coreResults = cGSEAcore(ElistObject = ElistObject, contrast = contrastMatrix, geneset = geneSetCollection, specie = specie, entrezGenesIds = ENTREZGenesIds, camera = methods$camera, gage = methods$gage, globaltest = methods$globaltest, gsva = methods$gsva, ssgsea = methods$ssgsea, zscore = methods$zscore, plage = methods$plage, ora = methods$ora, padog = methods$padog, roast = methods$roast, safe = methods$safe, setrank = methods$setrank, num.workers = num.workers)
  preparedData = prepareData(coreResults, alpha = alpha, directoryPath,  pvalAdjMethod = pvalAdjMethod , pvalCombMethod = pvalCombMethod , min.intersection.size = min.intersection.size)

  if (shinyMode == FALSE){
    cGSEAMakePlots(preparedData, directoryPath = directoryPath)
  } else if(shinyMode == TRUE){
    return(preparedData)
  }
}

cGSEAcore = function(ElistObject, contrast, geneset = "H", specie = "Mus musculus", entrezGenesIds, camera = TRUE, gage = TRUE, globaltest = TRUE, gsva = TRUE, ssgsea = TRUE, zscore = TRUE, plage = TRUE, ora = TRUE, padog = TRUE, roast = TRUE, safe = TRUE, setrank = FALSE, num.workers = 4, verbose = TRUE){



  genesetIdx = genesetToIdx(geneset = geneset, specie = specie, entrezGenesIds = entrezGenesIds)

  output = list()
  res = list()
  time = list()


  if(camera == TRUE){
    print("RUNNING CAMERA")
    try({
        camerastart = Sys.time()
        res$camera = runcamera(voom.results = ElistObject, contrast = contrast, genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
        camerastop = Sys.time()
        time$camera = difftime(camerastop, camerastart, units = c("secs"))

        })

  }

  if (gage == TRUE){
    print("RUNNING GAGE")
    try({
        gagestart = Sys.time()
        res$gage = rungage(voom.results = ElistObject, contrast = contrast, genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
        gagestop = Sys.time()
        time$gage = difftime(gagestop,gagestart, units = c("secs"))
        })

  }
  if (globaltest == TRUE){
    print("RUNNING GLOBALTEST")
    try({
        globalteststart = Sys.time()
        res$globaltest = runglobaltest(voom.results = ElistObject, contrast = contrast, genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
        globalteststop = Sys.time()
        time$globaltest = difftime(globalteststop,globalteststart, units = c("secs"))
        })

  }
  if (gsva == TRUE){
    print("RUNNING GSVA")
    try({
        gsvastart = Sys.time()
        res$gsva = rungsva(method = "gsva", voom.results = ElistObject, contrast = contrast,  genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
        gsvastop = Sys.time()
        time$gsva = difftime(gsvastop,gsvastart, units = c("secs"))
        })

  }
  if (ssgsea == TRUE){
    print("RUNNING SSGSEA")
    try({
        ssgseastart = Sys.time()
        res$ssgsea = rungsva(method = "ssgsea", voom.results = ElistObject, contrast = contrast,  genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
        ssgseastop = Sys.time()
        time$ssgsea = difftime(ssgseastop,ssgseastart, units = c("secs"))
        })

  }
  if (zscore == TRUE){
    print("RUNNING ZSCORE")
    try({
        zscorestart = Sys.time()
        res$zscore = rungsva(method = "zscore", voom.results = ElistObject, contrast = contrast,  genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
        zscorestop = Sys.time()
        time$zscore = difftime(zscorestop,zscorestart, units = c("secs"))
        })

  }
  if (plage == TRUE){
    print("RUNNING PLAGE ")
    try({
        plagestart = Sys.time()
        res$plage = rungsva(method = "plage", voom.results = ElistObject, contrast = contrast,  genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
        plagestop = Sys.time()
        time$plage = difftime(plagestop,plagestart, units = c("secs"))
        })

  }
  if (ora == TRUE){
    print("RUNNING ORA")
    try({
        orastart = Sys.time()
        res$ora = runora(voom.results = ElistObject, contrast = contrast,  genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
        orastop = Sys.time()
        time$ora = difftime(orastop, orastart, units = c("secs"))
        })

  }
  if (padog == TRUE){
    print("RUNNING PADOG")
    try({
        padogstart = Sys.time()
        res$padog = runpadog(voom.results = ElistObject, contrast = contrast,  genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
        padogstop = Sys.time()
        time$padog = difftime(padogstop, padogstart, units = c("secs"))
        })

  }
  if (roast == TRUE){
    print("RUNNING ROAST")
    try({
        roaststart = Sys.time()
        res$roast = runroast(voom.results = ElistObject, contrast = contrast,  genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
        roaststop = Sys.time()
        time$roast = difftime(roaststop,roaststart, units = c("secs"))
        })

  }
  if (safe == TRUE){
    print("RUNNING SAFE")
    try({
        safestart = Sys.time()
        res$safe = runsafe(voom.results = ElistObject, contrast = contrast,  genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
        safestop = Sys.time()
        time$safe = difftime(safestop, safestart, units = c("secs"))
        })

  }
  # if (setrank == TRUE){
  #   print("RUNNING SETRANK")
  #   setrankstart = Sys.time()
  #   res$setrank = runsetrank(voom.results = ElistObject, contrast = contrast,  geneset = geneset, specie = specie, num.workers = num.workers, verbose = verbose)
  #   setrankstop = Sys.time()
  #   time$setrank = difftime(setrankstop, setrankstart, units = c("secs"))
  # }
  output$time = time
  output$res = res
  output$contrast = contrast
  output$collection = genesetIdx
  output$Elist = ElistObject
  return(output)
}
