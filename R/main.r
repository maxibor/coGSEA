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
        camerastart = base::Sys.time()
        res$camera = runcamera(voom.results = ElistObject, contrast = contrast, genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
        camerastop = base::Sys.time()
        time$camera = base::difftime(camerastop, camerastart, units = c("secs"))

        })

  }

  if (gage == TRUE){
    print("RUNNING GAGE")
    try({
        gagestart = base::Sys.time()
        res$gage = rungage(voom.results = ElistObject, contrast = contrast, genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
        gagestop = base::Sys.time()
        time$gage = base::difftime(gagestop,gagestart, units = c("secs"))
        })

  }
  if (globaltest == TRUE){
    print("RUNNING GLOBALTEST")
    try({
        globalteststart = base::Sys.time()
        res$globaltest = runglobaltest(voom.results = ElistObject, contrast = contrast, genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
        globalteststop = base::Sys.time()
        time$globaltest = base::difftime(globalteststop,globalteststart, units = c("secs"))
        })

  }
  if (gsva == TRUE){
    print("RUNNING GSVA")
    try({
        gsvastart = base::Sys.time()
        res$gsva = rungsva(method = "gsva", voom.results = ElistObject, contrast = contrast,  genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
        gsvastop = base::Sys.time()
        time$gsva = base::difftime(gsvastop,gsvastart, units = c("secs"))
        })

  }
  if (ssgsea == TRUE){
    print("RUNNING SSGSEA")
    try({
        ssgseastart = base::Sys.time()
        res$ssgsea = rungsva(method = "ssgsea", voom.results = ElistObject, contrast = contrast,  genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
        ssgseastop = base::Sys.time()
        time$ssgsea = base::difftime(ssgseastop,ssgseastart, units = c("secs"))
        })

  }
  if (zscore == TRUE){
    print("RUNNING ZSCORE")
    try({
        zscorestart = base::Sys.time()
        res$zscore = rungsva(method = "zscore", voom.results = ElistObject, contrast = contrast,  genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
        zscorestop = base::Sys.time()
        time$zscore = base::difftime(zscorestop,zscorestart, units = c("secs"))
        })

  }
  if (plage == TRUE){
    print("RUNNING PLAGE ")
    try({
        plagestart = base::Sys.time()
        res$plage = rungsva(method = "plage", voom.results = ElistObject, contrast = contrast,  genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
        plagestop = base::Sys.time()
        time$plage = base::difftime(plagestop,plagestart, units = c("secs"))
        })

  }
  if (ora == TRUE){
    print("RUNNING ORA")
    try({
        orastart = base::Sys.time()
        res$ora = runora(voom.results = ElistObject, contrast = contrast,  genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
        orastop = base::Sys.time()
        time$ora = base::difftime(orastop, orastart, units = c("secs"))
        })

  }
  if (padog == TRUE){
    print("RUNNING PADOG")
    try({
        padogstart = base::Sys.time()
        res$padog = runpadog(voom.results = ElistObject, contrast = contrast,  genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
        padogstop = base::Sys.time()
        time$padog = base::difftime(padogstop, padogstart, units = c("secs"))
        })

  }
  if (roast == TRUE){
    print("RUNNING ROAST")
    try({
        roaststart = base::Sys.time()
        res$roast = runroast(voom.results = ElistObject, contrast = contrast,  genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
        roaststop = base::Sys.time()
        time$roast = base::difftime(roaststop,roaststart, units = c("secs"))
        })

  }
  if (safe == TRUE){
    print("RUNNING SAFE")
    try({
        safestart = base::Sys.time()
        res$safe = runsafe(voom.results = ElistObject, contrast = contrast,  genesetIdx = genesetIdx, num.workers = num.workers, verbose = verbose)
        safestop = base::Sys.time()
        time$safe = base::difftime(safestop, safestart, units = c("secs"))
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
