% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/comparResumPlot.r
\name{comparResumPlot}
\alias{comparResumPlot}
\title{Internal coGSEA function}
\usage{
comparResumPlot(preparedData, savePlot = TRUE, legend = TRUE,
  directoryPath = directoryPath)
}
\arguments{
\item{preparedData}{output of prepareData}

\item{savePlot}{boolen}

\item{legend}{boolean}

\item{directoryPath}{character}
}
\value{
nothing
}
\description{
Internal coGSEA function
}
\examples{
\dontrun{
comparResumPlot(preparedData = preparedData, directoryPath = paste0(directoryPath,"/plots/"))
}
}
