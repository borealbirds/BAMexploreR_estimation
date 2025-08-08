
#1. attach packages----
print("* attaching packages on master *")
library(dplyr)
library(readr)
library(opticut)
library(terra)
library(parallel)



#2. define local or cluster
test <- FALSE
cc <- TRUE

#3. set number of tasks for local vs cluster----
if(cc){ n_tasks <- 32}
if(!cc | test){ n_tasks <- 4}


#4. create and register clusters----
# creates 64 copies of R running in parallel via 64 tasks, on one of the cluster's sockets (processors). 
# Belgua has ~965 nodes
print("* creating clusters *")
cl <- makePSOCKcluster(n_tasks, type="PSOCK")

# print number of tasks and host name for confirmation
cl


#5. set root path----
print("* setting root file path *")

if(!cc){root <- "G:/Shared drives/BAM_NationalModels4/NationalModels4.0/Feb2020/artifacts"}
if(cc){root <- "/home/mannfred/scratch/mean_rasters_v4"}

# print for confirmation
root

tmpcl <- clusterExport(cl, c("root"))

# prep some data and functions and export to cluster
priority_species <- read_csv(file.path(root, "priority_spp_with_model_performance.csv"))
flbc <- priority_species$species_code

# this function computes the optimum 1/0 threshold for a given species
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
    opticut_result <- opticut::opticut(Y = pixel_densities, strata = strata, dist = "gaussian", ...)
    
    
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



# this function loads in a raster, applies `testoccurrenceNM`, and organizes the output
process_species <- function(species_code) {
  current_path <- file.path(root, paste0("pred-", species_code, "-CAN-Mean.tif"))
  
  if (file.exists(current_path)) {
    print(paste("now working on", current_path))
    raster <- terra::rast(x = current_path)
    return(testoccurrenceNM(raster = raster, method = "opticut", subset=5000))
  } else {
    return(data.frame(optimum_threshold = NA, og_sum = NA, retained_sum = NA))
    print(paste("could not find", current_path))
  }
}




tempcl <- clusterExport(cl, c("root", "flbc", "testoccurrenceNM", "process_species"))


#6. attach packages on clusters----
# `clusterEvalQ` evaluates a literal expression on each cluster node. 
print("* Loading packages on workers *")

tmpcl <- clusterEvalQ(cl, library(opticut))
tmpcl <- clusterEvalQ(cl, library(terra))
tmpcl <- clusterEvalQ(cl, library(dplyr))


print("* starting parallel processing *")
list1 <- parLapply(cl=cl, X=flbc, fun=process_species)

stopCluster(cl)
print("* cluster stopped *")

print("*saving list1*")
saveRDS(list1, file=file.path(root, "list1.rds"))

if(cc){ q() }

# list1 <- readRDS(file.path("C:/Users/mannf/Proton Drive/mannfredboehm/My files/Drive/BAMexploreR_estimation/interpretation/list1_opticut.rds"))
# 
# df <-
#   purrr::list_rbind(list1) |>
#   mutate(flbc = flbc) |>
#   as_tibble()
# 
# #plot with flbc as labels
# ggplot(df, aes(x = og_sum, y = retained_sum, label = flbc)) +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +  # 1:1 line
#   geom_text(size = 3) +
#   geom_smooth(method = "lm", color = "pink", se =FALSE, linetype = "solid") +
#   labs(x = "sum of original pixels", y = "sum of retained pixels") +
#   scale_x_continuous(labels = scales::label_comma()) +  # format x-axis with commas
#   scale_y_continuous(labels = scales::label_comma())+
#   theme_minimal()
# 
# # get slope
# lmfit <- lm(retained_sum ~ og_sum, data = df)
# coef(lmfit)[2]
