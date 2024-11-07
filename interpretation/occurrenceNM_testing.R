library(tidyverse)
library(opticut)
library(terra)


root <- "G:/Shared drives/BAM_NationalModels4/NationalModels4.0/Feb2020/artifacts"
priority_species <- read_csv(file="G:/Shared drives/BAM_NationalModels5/data/priority_spp_with_model_performance.csv")



testoccurrenceNM <- function(raster, method=c("opticut", "lorenz"), subset=NULL, quantile=NULL, plot=FALSE, ...){
  
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
    
    
  } else if (method == "lorenz") {
    
    # generate the Lorenz curve
    lorenz_fit <- opticut::lorenz(dpp_no_nas)
    
    if (is.null(quantile)){
      
      # locate the pixel where the tangent intersect occurs (slope approaches 1:1)
      t_pixel <- summary(lorenz_fit)["t"]
      
      # "x" is actual bird density at the tangent intersection point
      # ("p" is the proportion of pixels, "L" is the proportion of birds)
      optimum_threshold <- lorenz_fit[t_pixel, "x"]
      
    }
  }
  
  # original sum of all pixel densities before thresholding (for the x-axis)
  og_sum <- sum(dpp_no_nas)
  
  # retained pixel sum (for the y-axis)
  retained_sum <- sum(dpp_no_nas[dpp_no_nas >= optimum_threshold])
  
  return(data.frame(optimum_threshold = optimum_threshold, og_sum = og_sum, retained_sum = retained_sum))
}



flbc <- priority_species$species_code
list1 <- list()

for (i in 1:length(flbc)){
  
  current_path <- file.path(root, flbc[i], paste0("pred-", flbc[i], "-CAN-Mean.tif"))
  
  if (file.exists(current_path)) {
    
    raster <- terra::rast(x=file.path(root, flbc[i], paste0("pred-", flbc[i], "-CAN-Mean.tif")))
 
    list1[[i]] <- testoccurrenceNM(raster=raster, method = "lorenz")
  
  } else { list1[[i]] <- data.frame(optimum_threshold=NA, og_sum=NA, retained_sum=NA) }
  
  print(paste("iteration", i))

}

  
df <- 
  purrr::list_rbind(list1) |> 
  mutate(flbc = flbc) |> 
  as_tibble()

# plot with flbc as labels
ggplot(df, aes(x = og_sum, y = retained_sum, label = flbc)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +  # 1:1 line
  geom_text(size = 3) +  # avoid overlapping labels
  labs(
    x = "sum of original pixels",
    y = "sum of retained pixels"
  ) +
  scale_x_continuous(labels = scales::label_comma()) +  # format x-axis with commas
  scale_y_continuous(labels = scales::label_comma())+ 
  theme_minimal()
