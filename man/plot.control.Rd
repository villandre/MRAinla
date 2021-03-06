% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotFunctions.R
\name{plot.control}
\alias{plot.control}
\title{Control parameters for plot.INLAMRA}
\usage{
\method{plot}{control}(
  trim = 0,
  plotRaster = TRUE,
  graphicsEngine = jpeg,
  resolutionInMeters = NULL,
  timesToPlot = NULL,
  controlForScaleBar = NULL,
  controlForRasterPlot = NULL,
  controlForRasterLegend = NULL,
  controlForRasterColourScale = NULL
)
}
\arguments{
\item{trim}{Either 0, 1, or 2. Determines if predicted values should be trimmed before plotting to fit within the observation range. 0 (default): no trimming, 1: Values are trimmed to fall within the range of all observations, 2: Values are trimmed to fall within the range of observations on the corresponding day (if available).}

\item{plotRaster}{Logical. Should rasters (TRUE) or points (FALSE) be plotted? Ignored if type is "marginals". Points may be preferrable if data were not collected on a regular spatial grid.}

\item{graphicsEngine}{Function, e.g. jpeg, png, postscript. Corresponds to the graphics engine used to produce the plots,}

\item{resolutionInMeters}{Numeric value indicating the width, in meters, of each (square) cell in the raster. Ignored if plotRaster is FALSE.}

\item{timesToPlot}{Vector of dates in the same format as that used in the data used to fit INLAMRA. Only maps corresponding to those days will be plotted.}

\item{controlForScaleBar}{List of control parameters for the scale bars in the graphs, cf. ?mapmisc::scaleBar.}

\item{controlForRasterPlot}{List of control parameters for the raster plotting function, cf. ?raster::plot.}

\item{controlForRasterLegend}{List of control parameters for the raster legend, cf. ?mapmisc::legendBreaks.}

\item{controlForRasterColourScale}{List of control parameters for the raster colour scale, cf. ?mapmisc::colourScale.}
}
\value{
A list of control parameters.
}
\description{
Control parameters for the `control` argument of plot.
}
\details{
The function need not be called explicitly: it's just a convenient way to see/set control parameters.
}
\examples{
\dontrun{
}
}
