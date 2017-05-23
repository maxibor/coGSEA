runTimePlot = function(preparedData, savePlot = TRUE, directoryPath = directoryPath){
  print("Plotting runtime")
  pdf(paste(directoryPath,"cGSEA_methods_runtime.pdf",sep = ""), width = 10, height = 7, useDingbats = FALSE)
  graphics::barplot(unlist(preparedData$time), las = 2, ylab = "second", main = "RunTime")
  dev.off()

}

clusteringPlot = function(preparedData, contrCondi, savePlot = TRUE, directoryPath = directoryPath){
  print(paste("Plotting cluster Plot for condition", contrCondi))
  pdf(paste(directoryPath, contrCondi, "_clustering.pdf", sep = ""), width = 10, height = 7, useDingbats = FALSE)
  graphics::plot(preparedData$clustering[[contrCondi]], hang = -0.1)
  dev.off()
}

heatmapPlot = function(preparedData, contrCondi, savePlot = TRUE, directoryPath = directoryPath){

    if(savePlot == TRUE){
        print(paste("Plotting heatmap for condition", contrCondi))
        pdf(paste(directoryPath, contrCondi, "_pvalue_heatmap.pdf", sep = ""), width = 10, height = 7, useDingbats = FALSE)
        pheatmap::pheatmap(preparedData$heatmap[[contrCondi]], fontsize_row = 5)
        dev.off()
    } else {
        pheatmap::pheatmap(preparedData$heatmap[[contrCondi]])
    }
}

pcaPlot = function(preparedData, contrCondi, savePlot = TRUE, directoryPath = directoryPath){

  res.pca = FactoMineR::PCA(preparedData$PCA[[contrCondi]], scale.unit = T, axes = c(1,2), graph = FALSE)
  eigen = res.pca$eig$`percentage of variance`
  base::names(eigen) = base::rownames(res.pca$eig)

  print(paste("Plotting Eigen values Plot for condition", contrCondi))
  pdf(paste(directoryPath, contrCondi, "_eigen_fall.pdf", sep = ""), width = 10, height = 7, useDingbats = FALSE)
  graphics::barplot(eigen, las = 2, ylab = "%")
  dev.off()

  print(paste("Plotting PCA Plot for condition", contrCondi))
  pdf(paste(directoryPath, contrCondi, "_pca.pdf", sep = ""), width = 10, height = 7, useDingbats = FALSE)
  graphics::plot(res.pca, choix = 'ind')
  dev.off()
}

corPlot = function(preparedData, contrCondi, savePlot = TRUE, directoryPath = directoryPath){

  print(paste("Plotting correlation Plot for condition", contrCondi))
  pdf(paste(directoryPath, contrCondi, "_correlation.pdf", sep = ""), width = 10, height = 7, useDingbats = FALSE)
  corrplot::corrplot(preparedData$correlation[[contrCondi]], method = "circle", type = "full", order = "hclust",mar=c(10, 4, 4, 2) + 0.1)
  dev.off()
}

upsetrPlot = function(preparedData, contrCondi, savePlot = TRUE, directoryPath = directoryPath){

  print(paste("Plotting UpsetR Plot for condition", contrCondi))
  pdf(paste(directoryPath, contrCondi, "_upsetr.pdf", sep = ""), width = 10, height = 7, useDingbats = FALSE)
  UpSetR::upset(fromList(preparedData$snailPlot[[contrCondi]]),order.by = "freq")
  dev.off()

}

snailPlot = function(preparedData, contrCondi, savePlot = TRUE, min.intersection.size = 1, directoryPath = directoryPath){

  snail = SuperExactTest::supertest(preparedData$snailPlot[[contrCondi]], n = length(levels(preparedData$snailPlot[[contrCondi]][[1]])))
  print(paste("Plotting Snail Plot for condition", contrCondi))
  if (savePlot == TRUE){
      pdf(paste(directoryPath, contrCondi, "_snailplot.pdf", sep = ""), width = 10, height = 7, useDingbats = FALSE)
      SuperExactTest::plot.msets(snail, sort.by="size", keep.empty.intersections=FALSE, min.intersection.size = min.intersection.size)
      dev.off()
  } else {
      SuperExactTest::plot.msets(snail, sort.by="size", keep.empty.intersections=FALSE, min.intersection.size = min.intersection.size)
  }

}

resumPlot1 = function(preparedData, contrCondi, savePlot = TRUE, directoryPath = directoryPath){

  print(paste("Plotting simple Summary Plot for condition", contrCondi))
  pdf(paste(directoryPath, contrCondi, "_simple_sumplot.pdf", sep = ""), width = 10, height = 7, useDingbats = FALSE)
  graphics::plot(preparedData$resumPlot1[[contrCondi]][["x"]], preparedData$resumPlot1[[contrCondi]][["y"]], xlab = '-log10(pVal)', ylab = "avg logFC")
  graphics::abline(h = 0)
  graphics::abline(v = 0)
  graphics::text(preparedData$resumPlot1[[contrCondi]][['x']], preparedData$resumPlot1[[contrCondi]][['y']], labels=preparedData$resumPlot1[[contrCondi]][['labels']], cex= 0.7)
  graphics::abline(v = -log10(0.05), col = "red")
  dev.off()
}

resumPlot2 = function(preparedData, contrCondi, savePlot = TRUE, directoryPath= directoryPath){
  print(paste("Plotting advanced Summary Plot for condition", contrCondi))
  generateSummaryPlots(preparedData$resumPlot2[[contrCondi]], savePlot = savePlot, file.name = paste(directoryPath, contrCondi, "_advanced_sumplot", sep = ""), format = "pdf")

}

cGSEAMakePlots = function(preparedData, directoryPath){
  base::dir.create(base::file.path(directoryPath,"/plots/"), showWarnings = FALSE)
  for (condi in base::names(preparedData$result)){
    base::dir.create(file.path(directoryPath,"/plots/", condi), showWarnings = FALSE)
    clusteringPlot(preparedData = preparedData, contrCondi = condi, directoryPath = paste0(directoryPath,"/plots/",condi,"/"))
    pcaPlot(preparedData = preparedData, contrCondi = condi, directoryPath = paste0(directoryPath,"/plots/",condi,"/"))
    corPlot(preparedData = preparedData, contrCondi = condi, directoryPath = paste0(directoryPath,"/plots/",condi,"/"))
    upsetrPlot(preparedData = preparedData, contrCondi = condi, directoryPath = paste0(directoryPath,"/plots/",condi,"/"))
    heatmapPlot(preparedData, contrCondi = condi, savePlot = TRUE, directoryPath =  paste0(directoryPath,"/plots/",condi,"/"))
    snailPlot(preparedData = preparedData, min.intersection.size = preparedData$minIntersectionSize, contrCondi = condi, directoryPath = paste0(directoryPath,"/plots/",condi,"/"))
    # resumPlot1(preparedData = preparedData, contrCondi = condi, directoryPath = paste0(directoryPath,"/plots/",condi,"/"))
    resumPlot2(preparedData = preparedData, contrCondi = condi, directoryPath = paste0(directoryPath,"/plots/",condi,"/"))
  }
  runTimePlot(preparedData = preparedData, directoryPath = paste0(directoryPath,"/plots/"))
  comparResumPlot(preparedData = preparedData, directoryPath = paste0(directoryPath,"/plots/"))
}
