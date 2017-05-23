
#' Internal coGSEA function
#'
#' @param preparedData  output of prepareData
#' @param savePlot boolen
#' @param legend boolean
#' @param directoryPath character
#' @return nothing
#' @examples
#'\dontrun{
#'comparResumPlot(preparedData = preparedData, directoryPath = paste0(directoryPath,"/plots/"))
#'}
#' @export


comparResumPlot = function(preparedData, savePlot = TRUE, legend = TRUE, directoryPath = directoryPath){
    print("Plotting Summary comparison Plot for all conditions")
    generateSummaryPlots(preparedData$comparison, savePlot = savePlot, legend = legend, file.name = paste(directoryPath, "_comparison_sumplot", sep = ""), format = "pdf")

}
