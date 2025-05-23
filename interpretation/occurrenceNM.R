# ---
# title: National Models 5.0 - occurrence binarization
# author: Mannfred Boehm
# created: October 30, 2024
# ---


#'@param raster A raster of the "Area of Interest". Defined and created via Melina's function(s). 
#'
#'@param method One of "opticut" or "lorenz". The former  calls `opticut::opticut()` and finds the optimum 
#'threshold by fitting the binomial distribution to the various 2-partitions created by different density 
#'threshold cutoffs. The optimum threshold comes from the partitioning scheme with the highest 
#'log-likelihood ratio. By setting `method="lorenz"` with `quantile=NULL`, the  optimum threshold is estimated as 
#'the pixel density where the slope of the tangent of the Lorenz function is 1.   By setting `method="lorenz"` 
#'with e.g. `quantile=0.8`, 'the threshold is set by the user as the cumulative proportion of pixels >= 0.8. 
#'
#'@param subset Only used when `method="opticut"`. A numeric passed to `sample()` for subsetting the input raster. 
#'
#'@param quantile A numeric between 0 and 1. Indicates a cumulative proportion of pixels, above which all 
#'raster pixel values are assigned as "presence". Raster pixel values below the quantile are assigned "absent". 
#'E.g. By setting `quantile=0.8`, the threshold density for separating presence versus absence is whatever
#'pixel value accumulates 80% of the total pixels from the raster. 
#'
#'@param plot Default is `FALSE`. Set to `TRUE` to visualize occurrence patterns. 
#'Otherwise, see `terra::plot()` for customized plotting.
#'
#'@param ... Additional arguments to be passed to `opticut()`. 
#'
#'@return A list with two elements: (1) a raster with pixels indicating presence (1) or absence (0) of the species 
#'in the input raster. (2) A numeric indicating the estimated optimum threshold density value. 
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


# Set the robust log-likelihood option
options(ocoptions = list(
  collapse = "_",
  check_comb = TRUE,
  scale = "linear",
  try_error = TRUE,
  robust_loglik = TRUE  # Enables the safeguard for ill-defined likelihoods
))

# not sure which folder to use for the V4 models, sorry if this is random!
raster <- terra::rast(x="G:/Shared drives/BAM_NationalModels4/NationalModels4.0/May2020/pred250-ALFL-BCR_10-boot-1.tif")
raster <- terra::rast(x="G:/Shared drives/BAM_NationalModels4/NationalModels4.0/May2020/pred250-BTNW-BCR_4-boot-19.tif")






occurrenceNM <- function(raster, method=c("opticut", "lorenz"), subset=NULL, quantile=NULL, plot=FALSE, ...){
  
  # choices are obtained from a default setting for the formal argument `arg` 
  # of the function from which match.arg was called.
  method <- match.arg(method)
  
  # retrieve density per pixel
  dpp <- terra::values(raster)

  # remove NAs for `lorenz()`
  dpp_no_nas <- dpp[!is.na(dpp)]
  
  # find optimal threshold based on user's specified method
  if (method == "opticut") {
   
    # subset pixel-level predictions for faster computing (`sset=` argument in `opticut()` isn't working)
    if (is.null(subset)) {
      
      pixel_densities <- dpp_no_nas  
      strata <- as.factor(seq_along(pixel_densities)) # treat each pixel as a stratum (a column in `Z` used by .opticut1())
    
    } else if (is.numeric(subset)) {
      
      if (!is.null(subset) && sum(!is.na(dpp_no_nas)) < subset) {
        stop("subset size exceeds available non-NA values")
      }
      
      subset <- sample(seq_along(dpp_no_nas), size = subset)
      pixel_densities <- dpp_no_nas[subset] 
      strata <- as.factor(seq_along(subset)) # treat each pixel as a stratum (a column in `Z` used by .opticut1())
      
    }
    
    # .opticut1() computes the log-likelihood ratio (LLR) for a given threshold compared to a null model (where no threshold is applied).
    # the model with the highest LLR is the best at separating high-density (presence) pixels from low-density (absence) pixels
    opticut_result <- opticut::opticut(Y = pixel_densities, strata = strata, dist = "binomial", ...)
    
  
    # extract the threshold based on the optimal partition (maximum log likelihood ratio)
    # mu0 is the expected value of group 1 (absence), mu1 is the expected value of group 2 (presence)
    max_llr <- which.max(opticut_result$species[[1]]$logLR)
    optimum_threshold <- opticut_result$species[[1]]$mu1[max_llr]
   
    # assign pixels 1 or 0 based on the current threshold
    # preserve NA positions from the original raster to maintain raster range between input and output
    # Create a binary classification using the original dpp vector, 
    pixels_binary <- ifelse(dpp >= optimum_threshold, 1, ifelse(!is.na(dpp), 0, NA))
    
    
    } else if (method == "lorenz") {
    
    # generate the Lorenz curve
    lorenz_fit <- opticut::lorenz(dpp_no_nas)
    
      if (is.null(quantile)){
        
        # locate the pixel where the tangent intersect occurs (slope approaches 1:1)
        t_pixel <- summary(lorenz_fit)["t"]
    
        # "x" is actual bird density at the tangent intersection point
        # ("p" is the proportion of pixels, "L" is the proportion of birds)
        optimum_threshold <- lorenz_fit[t_pixel, "x"]
        
        } else {
        
        # `L` for ordered cumulative abundance quantiles (versus non-cumulative)
        # `threshold` partitions "1-threshold" proportion of values as presence (1) and the rest ("threshold") as absence (0)
        # e.g. for `threshold=0.8` the densest 20% of values are assigned as presence (1) and the rest as absence (0)
        optimum_threshold <- opticut:::quantile.lorenz(lorenz_fit, probs = quantile, type = "L")
  
        } 
    
    } # close all methods
      
  
  # binarize density raster, i.e.
  # assign pixels 1 or 0 based on the current threshold
  # note: need to preserve NA positions from the original raster to maintain raster range between input and output
  classify_pixels <- function(dpp, threshold) {
    ifelse(dpp >= threshold, 1, ifelse(!is.na(dpp), 0, NA))
  }
  
  pixels_binary <- classify_pixels(dpp, optimum_threshold)
  
  # write the binary values into a raster object
  binary_raster <- terra::setValues(raster, pixels_binary)
  
  # communicate threshold value
  print(paste("threshold = ", round(optimum_threshold, digits = 5)))
  
  # plot if asked for
  if (plot==TRUE){
    terra::plot(binary_raster, main=paste("density threshold =", round(optimum_threshold, digits = 5)))
  } 
  
  return(list(raster = binary_raster, threshold = optimum_threshold))
  
  
  
} # close function

  




# Extract log-likelihood values
log_likelihood_ratios <- opticut_result$species[[1]]$logLR

# Print a summary of the log-likelihoods to understand their range
summary(log_likelihood_ratios)

# Plot the subsetted raster
raster_df <- as.data.frame(raster, xy = TRUE)
colnames(raster_df)[3] <- "value" 
subset_raster_df <- raster_df |>  filter(row_number() %in% sset)

ggplot(raster_df,aes(x = x, y = y, fill = value)) +
  # Background layer: Full raster
  geom_raster()+
  scale_fill_viridis_c() +
  # Overlay layer: Subset of raster
  geom_tile(data = subset_raster_df, aes(x = x, y = y), fill = "red", width = 5000, height = 5000) +
  
  labs(title = "Subset Overlay on Original Raster", x = "Longitude", y = "Latitude") +
  theme_minimal()
