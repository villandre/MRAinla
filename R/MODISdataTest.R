#' MODIS Land Surface Temperature (LST) (Test set)
#'
#' A very small dataset containing covariate information for missing LST observations between May 19, 2012 and May 24, 2012 in the western part of the state of Maharashtra, India. It is meant to be used with MODISdataTraining.
#'
#' @format A spacetime::STIDF object with 385 observations, and four variables:
#' \describe{
#'   \item{landCover*}{Dummy variables indicating Land Cover Classification (UMD type II). Reference category is Land Cover = 0 (water).}
#'   \item{elevation}{elevation values in meters.}
#'   \item{time*}{Dummy variables indicating observation time. Reference day is May 18 (= time1).}
#'   \item{Aqua}{Binary variable indicating whether the reading was collected by the Aqua or Terra Satellite ('1' corresponds to Aqua).}
#' }
#' @source \url{https://lpdaac.usgs.gov/products/mod11a1v006/}, \url{https://lpdaac.usgs.gov/products/mcd12q1v006/}, \url{https://modis-land.gsfc.nasa.gov/pdf/MOD16UsersGuideV2.022019.pdf}, \url{https://asterweb.jpl.nasa.gov/gdem.asp}
#'
"MODISdataTest"
