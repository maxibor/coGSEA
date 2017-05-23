
#' Plotting function to make a eGSEA summary plot
#'
#' @param plot.data  A dataframe formatted for this function. Output of coGSEA$resumPlot2
#' @param file.name Name of  the  file to save. Not needed if savePlot = FALSE (character)
#' @param Xlab  label for X axis (character)
#' @param Ybal  label for Y axis (character)
#' @param firstN  N number of top N genes to highlight in blue (integer)
#' @param savePlot  Wheter to save the plot, or just display it (boolean)
#' @param legen whether to display or not the legend if savePlot = FALSE (boolean)

#' @return A plot if savePlot = TRUE, else nothing
#' @examples
#'\dontrun{
#' generateSummaryPlots(plot.data = coGSEAres$resumPlot2, file.name = "resumPlot", Xlab="-log10(p-value)",Ylab="Average Absolute logFC", format = "pdf", firstN = 10,  savePlot = TRUE)
#'}
#' @export

generateSummaryPlots <- function(plot.data, file.name = "resumPlot", Xlab="-log10(p-value)",
        Ylab="Average Absolute logFC", format = NULL, firstN = 10,  savePlot = TRUE, legend = TRUE){

    tryCatch({
        plot.data.sig = plot.data[plot.data[, "rank"] <= firstN, ]
        sig.cols = rep("black", nrow(plot.data.sig))
        if (min(plot.data[, "x.data"], na.rm=TRUE) > 0){
            xlimits = c(0.8 * min(plot.data[, "x.data"], na.rm=TRUE),
                max(plot.data[, "x.data"], na.rm=TRUE)*1.05)
        }else{
            xlimits = c(1.05 * min(plot.data[, "x.data"], na.rm=TRUE),
                    max(plot.data[, "x.data"], na.rm=TRUE)*0.8)
        }
        if (max(plot.data[, "y.data"], na.rm=TRUE) > 0){
            ylimits = c(min(plot.data[, "y.data"], na.rm=TRUE),
                    max(plot.data[, "y.data"], na.rm=TRUE) * 1.05)
        }else{
            ylimits = c(min(plot.data[, "y.data"], na.rm=TRUE),
                    max(plot.data[, "y.data"], na.rm=TRUE) * 0.9)
        }
    #   print(plot.data.sig)
    #       print(dim(plot.data))
        if (savePlot == TRUE){
            # plot rank-based coloured bubbles
            p = ggplot2::qplot(x.data, y.data, data=plot.data, size=gsSize,asp=1,
                    colour=rank,
                    xlab = Xlab, ylab = Ylab,
                    xlim=xlimits,
                    ylim=ylimits)
            # customize bubbles colour
            p = p + ggplot2::scale_colour_gradient(guide="colourbar", low="#56B1F7",
    high="#000000",
                    limits=c(1,100), na.value="black", name="Rank")
            # customize bubble size
            p = p + ggplot2::scale_size("Gene set size", range=c(2,20))
            if (is.null(format) || tolower(format) == "pdf"){
                pdf(paste0(file.name, ".rank.pdf"), width = 10, height = 7,
                        useDingbats = FALSE)

                # label the bubbles of the top 10 gene sets
                print(p + ggplot2::geom_text(size=5, mapping=ggplot2::aes(x=x.data, y=y.data,
                                label=id),
                                data=plot.data.sig,
                                colour=sig.cols, vjust=-1, hjust=1) )
                dev.off()
            }
            if (is.null(format) || tolower(format) == "png"){
                png(paste0(file.name, ".rank.png"), width = 800, height = 700)
                print(p + ggplot2::geom_text(size=5, mapping=ggplot2::aes(x=x.data, y=y.data,
        label=id),
                                data=plot.data.sig,
        colour=sig.cols, vjust=-1, hjust=1) )
                dev.off()
            }
        }



        # plot direction-based coloured bubbles
        top.10.ids = as.character(plot.data[plot.data[, "rank"] <= firstN,
"id"])
        sig.ids = ggplot2::setdiff(plot.data[rank(-plot.data[,"sig"], na.last =
TRUE) <= 5, "id"], top.10.ids)
        sig.cols = c(rep("black", length(top.10.ids)), rep("blue",
length(sig.ids)))
        plot.data.sig = plot.data[match(c(top.10.ids, sig.ids),
plot.data[, "id"]), ]
        p = qplot(x.data, y.data, data=plot.data, size=sig,asp=1,
                colour=dir,
                xlab = Xlab, ylab = Ylab,
                xlim=xlimits,
                ylim=ylimits)
        p = p + ggplot2::scale_colour_gradient(guide="colourbar", low="#56B1F7",
high="#E35F5F",
                limits=c(-1,1), na.value="black",
name="Regulation Direction") # low="#5FE377"
        p = p + ggplot2::scale_size("significance", range=c(2,20))

        if (savePlot == TRUE){
            if (is.null(format) || tolower(format) == "pdf"){
                pdf(paste0(file.name, "direction_.pdf"), width = 10, height = 7,
                        useDingbats = FALSE)

                print(p + ggplot2::geom_text(size=5, mapping=ggplot2::aes(x=x.data, y=y.data,
                                label=id),
                                data=plot.data.sig,
                                colour=sig.cols, vjust=-1, hjust=1) )
                dev.off()
            }
            if (is.null(format) || tolower(format) == "png"){
                png(paste0(file.name, "_direction.png"), width = 800, height = 700)
                print(p + ggplot2::geom_text(size=5, mapping=aes(x=x.data, y=y.data,
                                label=id),
                                data=plot.data.sig,
                                colour=sig.cols, vjust=-1, hjust=1) )
                dev.off()
            }
        } else if (savePlot == FALSE && legend == TRUE){
            print(p + ggplot2::geom_text(size=5, mapping=aes(x=x.data, y=y.data,
                            label=id),
                            data=plot.data.sig,
                            colour=sig.cols, vjust=-1, hjust=1))
        } else if (savePlot == FALSE && legend == FALSE){
            p
        }

    },
    error = function(e){
        print(paste0("WARNING: summary plots were not generated for ",
file.name))
    })
}
