# ---
# title: National Models 5.0 - occurrence binarization
# author: Mannfred Boehm
# created: October 30, 2024
# ---


#'@param raster A raster of the "Area of Interest". Defined and created via Melina's function(s). 
#'
#'@param threshold A numeric between 0 and 1. Indicates the cutoff, above which all raster pixel values are assigned as "presence". Raster pixel values below the threshold are assigned "absent". 
#'
#'@param plot Default is `FALSE`. Set to `TRUE` to visualize occurrence patterns. Otherwise, see `terra::plot()` for customized plotting.
#'
#'@param ... Additional arguments to be passed to `opticut()`. 
#'
#'@return A raster with pixels indicating presence (1) or absence (0) of the species in the input raster.
#'
#'@importFrom opticut lorenz quantile opticut
#'@importFrom terra rast values
#'
#'
#'
#'@examples ...tbd



library(dplyr)
library(opticut)
library(terra)

library(parallel)
cl <- makeCluster(4)

# not sure which folder to use for the V4 models, sorry if this is random!
raster <- terra::rast(x="G:/Shared drives/BAM_NationalModels4/NationalModels4.0/bootpred/BBWA-BCR_ALL-boot-MeanRND.tif")

occurrenceNM <- function(raster, method=c("opticut", "lorenz"), quantile=NULL, plot=FALSE, ...){
  
  # choices are obtained from a default setting for the formal argument `arg` 
  # of the function from which match.arg was called.
  method <- match.arg(method)
  
  # retrieve density per pixel
  dpp <- terra::values(raster)

  # remove NAs for `lorenz()`
  dpp_no_nas <- dpp[!is.na(dpp)]
  
  
  # find optimal threshold based on user's specified method
  if (method == "opticut") {
   
    # subset predictions for faster computing (`sset=` argument in `opticut()` isn't working)
    sset <- sample(seq_along(dpp_no_nas), size = 15000)
    
    # create placeholder vector with an initial guess at the presence-absence split (median), 
    # this provides opticut() with an intial binary parition while it searches for a more accurate threshold.
    response <- as.numeric(dpp_no_nas[sset] > median(dpp_no_nas))
    
    # treat each pixel as a stratum (a column in `Z` used by .opticut1())
    strata <- as.factor(seq_along(sset))  
    
    # .opticut1() computes the log-likelihood ratio (LLR) for a given threshold compared to a null model (where no threshold is applied).
    # the model with the highest LLR is the best at separating high-density (presence) pixels from low-density (absence) pixels
    opticut_result <- opticut::opticut(response ~ dpp_no_nas[sset], strata = strata, dist = "binomial", cl = cl)
    
  
    # extract the threshold based on the optimal partition (maximum log likelihood ratio)
    max_llr <- which.max(opticut_result$species[[1]]$logLR)
    optimum_threshold <- opticut_result$species[[1]]$mu1[max_llr]
   
  
    # assign pixels 1 or 0 based on the current threshold
    # preserve NA positions from the original raster to maintain raster range between input and output
    # Create a binary classification using the original dpp vector, 
    pixels_binary <- ifelse(dpp >= optimum_threshold, 1, ifelse(!is.na(dpp), 0, NA))
    
    # write the binary values into a raster object
    binary_raster <- terra::setValues(raster, pixels_binary)
    
    print(paste("optimum density threshold estimated at", round(optimum_threshold, digits = 5)))
    
    if (plot==TRUE){
      terra::plot(binary_raster, main=paste("density threshold =", round(optimum_threshold, digits = 5)))
    } 
    
    
  } else if (method == "lorenz" & threshold=NULL) {
    
    # generate the Lorenz curve
    lorenz_fit <- opticut::lorenz(dpp_no_nas)
    
    # locate the pixel where the tangent intersect occurs (slope approaches 1:1)
    t_pixel <- summary(lorenz_fit)["t"]
    
    # "x" is actual bird density at the tangent intersection point
    # ("p" is the proportion of pixels, "L" is the proportion of birds)
    optimum_threshold <- lorenz_fit[t_pixel, "x"]
    
    # assign pixels 1 or 0 based on the current threshold
    # preserve NA positions from the original raster to maintain raster range between input and output
    # Create a binary classification using the original dpp vector, 
    pixels_binary <- ifelse(dpp >= optimum_threshold, 1, ifelse(!is.na(dpp), 0, NA))
    
    # write the binary values into a raster object
    binary_raster <- terra::setValues(raster, pixels_binary)
    
    print(paste("optimum density threshold estimated at", round(optimum_threshold, digits = 5)))
    
    if (plot==TRUE){
      terra::plot(binary_raster, main=paste("density threshold =", round(optimum_threshold, digits = 5)))
    } 
    
    
  } else if (method == "lorenz" & as.numeric(quantile)==TRUE)
 
    # `L` for ordered cumulative abundance quantiles (versus non-cumulative)
    # `threshold` partitions "1-threshold" proportion of values as presence (1) and the rest ("threshold") as absence (0)
    # e.g. for `threshold=0.8` the densest 20% of values are assigned as presence (1) and the rest as absence (0)
    user_quantile <- opticut:::quantile.lorenz(lorenz_fit, probs = quantile, type = "L")
  
    # assign pixels 1 or 0 based on the current threshold
    # preserve NA positions from the original raster to maintain raster range between input and output
    # Create a binary classification using the original dpp vector, 
    pixels_binary <- ifelse(dpp >= user_threshold, 1, ifelse(!is.na(dpp), 0, NA))
  
    # write the binary values into a raster object
    binary_raster <- terra::setValues(raster, pixels_binary)
  
    if (plot==TRUE){
    terra::plot(binary_raster, main=paste("quantile threshold =", user_threshold))
  } 
}
  
  
  
  
  return(binary_raster)
}
