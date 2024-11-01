# ---
# title: National Models 5.0 - occurrence binarization
# author: Mannfred Boehm
# created: October 30, 2024
# ---


#'@param aoi A raster of the "Area of Interest". Defined and created via Melina's function(s). 
#'
#'@param threshold A numeric between 0 and 1. Indicates the cutoff, above which all raster pixel values are assigned as "presence". Raster pixel values below the threshold are assigned "absent". 
#'
#'@param plot Default is `FALSE`. Set to `TRUE` to visualize occurrence patterns. Otherwise, see `terra::plot()` for customized plotting.
#'. 
#'@return A raster with pixels indicating presence (1) or absence (0) of the species in the input raster.
#'
#'@importFrom opticut lorenz
#'@importFrom opticut quantile
#'@importFrom terra rast
#'@importFrom terra values
#'
#'
#'
#'@examples ...tbd



library(dplyr)
library(opticut)
library(terra)

# not sure which folder to use for the V4 models, sorry if this is random!
raster <- terra::rast(x="G:/Shared drives/BAM_NationalModels4/NationalModels4.0/Feb2020/artifacts/CAWA/pred-CAWA-CAN-boot-7.tif")

occurrenceNM <- function(raster, threshold, plot=FALSE){
  
  # retrieve density per pixel
  dpp <- terra::values(raster)

  # remove NAs for `lorenz()`
  dpp_no_nas <- dpp[!is.na(dpp)]
  
  # generate the Lorenz curve
  lorenz_fit <- opticut::lorenz(dpp_no_nas)
  
  # `L` for ordered cumulative abundance quantiles (versus non-cumulative)
  # `threshold` partitions "1-threshold" proportion of values as presence (1) and the rest ("threshold") as absence (0)
  # e.g. for `threshold=0.8` the densest 20% of values are assigned as presence (1) and the rest as absence (0)
  current_threshold <- opticut:::quantile.lorenz(lorenz_fit, probs = threshold, type = "L")
  
  # assign pixels 1 or 0 based on the current threshold
  # preserve NA positions from the original raster to maintain raster range between input and output
  # Create a binary classification using the original dpp vector, 
  pixels_binary <- ifelse(dpp >= current_threshold, 1, ifelse(!is.na(dpp), 0, NA))
  
  # write the binary values into a raster object
  binary_raster <- terra::setValues(raster, pixels_binary)
  
  if (plot==TRUE){
    terra::plot(binary_raster, main=paste("species distribution with density threshold of", threshold))
  } 
  
  return(binary_raster)
}
