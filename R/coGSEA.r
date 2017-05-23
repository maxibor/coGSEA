
#' Wrapper function of coGSEA
#'
#' @param ElistObject  An Elist object.
#' @param contrastMatrix  A contrast matrix.
#' @param geneSetCollection  H, C2_KEGG, or C2_REACTOME, default to H (character). Can be custom of class list.
#' @param specie  Mus musculus or Homo sapiens, default to Mus Musculus (character).
#' @param ENTREZGenesIds  vector of gene ENTREZ identifiers.
#' @param directoryPath  path to existing directory where to save the results ("character").
#' @param alpha  alpha treshold, default to 0.05 (numeric).
#' @param pvalAdjMethod  p.value adjustment method (character).
#' @param pvalCombMethod  p.value combination method (character).
#' @param min.intersection.size  minimum intersection size for SnailPlot (integer)
#' @param GSEA.Methods  GSEA methods to run (vector)
#' @param num.workers  number of thread to use (integer)
#' @param shinyMode  Make plot or not (boolean)
#' @return A result list if shinyMode == TRUE, else nothing
#' @examples
#'\dontrun{
#'coGSEA(ElistObject = Elist, contrastMatrix = contrast, ENTREZGenesIds = Elist$genes$ENTREZID, directoryPath = "./")
#'}
#' @export


coGSEA = function(ElistObject,
    contrastMatrix,
    ENTREZGenesIds,
    geneSetCollection = "H",
    specie = "Mus musculus",
    directoryPath = "./",
    alpha = 0.05,
    pvalAdjMethod = "BH" ,
    pvalCombMethod = "sumlog",
    min.intersection.size = 1,
    GSEA.Methods = c("camera", "gage", "globaltest", "gsva", "ssgsea", "zscore", "plage", "ora", "padog", "roast", "safe"),
    num.workers = 4,
    shinyMode = FALSE){
  # loadLib()
  methods = method_checker(GSEA.Methods)
  coreResults = cGSEAcore(ElistObject = ElistObject, contrast = contrastMatrix, geneset = geneSetCollection, specie = specie, entrezGenesIds = ENTREZGenesIds, camera = methods$camera, gage = methods$gage, globaltest = methods$globaltest, gsva = methods$gsva, ssgsea = methods$ssgsea, zscore = methods$zscore, plage = methods$plage, ora = methods$ora, padog = methods$padog, roast = methods$roast, safe = methods$safe, setrank = methods$setrank, num.workers = num.workers)
  preparedData = prepareData(coreResults, alpha = alpha, directoryPath,  pvalAdjMethod = pvalAdjMethod , pvalCombMethod = pvalCombMethod , min.intersection.size = min.intersection.size, shinyMode = shinyMode)

  if (shinyMode == FALSE){
    cGSEAMakePlots(preparedData, directoryPath = directoryPath)
  } else if(shinyMode == TRUE){
    return(preparedData)
  }
}
