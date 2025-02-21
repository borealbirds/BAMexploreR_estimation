# ---
# title: National Models 4.0 - estimating strength of 2 way interactions
# author: Mannfred Boehm
# created: February 20, 2025
# ---

# output: a dataframe with four columns:
# column 1 is cov 1, column 2 is cov 2
# column 3 is mean interaction strength, column 4 is SD

#1. attach packages----
library(gbm)
library(parallel)
library(tidyverse)


#2. setup local or cluster----
test <- TRUE
cc <- FALSE

# set number of tasks for local vs cluster
if(cc){ n_tasks <- 25}
if(!cc | test){ n_tasks <- 4}


# create and register clusters
print("* creating clusters *")
cl <- makePSOCKcluster(n_tasks, type="PSOCK")

# print number of tasks and host name for confirmation
cl

# set root path
print("* setting root file path *")

if(!cc){root <- "C:/Users/mannf/Proton Drive/mannfredboehm/My files/Drive"}
if(cc){root <- "/home/mannfred/scratch"}

print(paste("* currently working out of ", root, " *"))



#3. create indices needed for searching every 2-way interaction for available models

# create an index of every bcr x spp x 2-way interactions 
# loads as "bam_covariate_importance_v4"
load(file.path(root, "BAMexploreR", "data", "bam_covariate_importance_v4.rda"))

# `RcppAlgo::comboGrid` is like `expand.grid` and `tidyr::crossing` but avoids duplicates 
# e.g. for our purposes (var1=x, var2=y) is a duplicate of (var1=y, var2=x), so the Cartesian Product of comboGrid is smaller
boot_group_keys_i2 <- 
  RcppAlgos::comboGrid(unique(bam_covariate_importance_v4$var), 
                       unique(bam_covariate_importance_v4$var), 
                       bam_covariate_importance_v4$bcr, 
                       bam_covariate_importance_v4$spp, repetition =  FALSE) |> 
  tibble::as_tibble() |> 
  dplyr::rename(var_1=Var1, var_2=Var2, bcr=Var3, spp=Var4)


# create an index of all available bootstrap models
gbm_objs <- list.files(file.path(root, "v4_bootstraps"), pattern = "^gnmboot-.*\\.RData$", full.names = TRUE)
message("Found ", length(gbm_objs), " files: ", gbm_objs)




#4. load in V4 count data---- 
# needed for running `gbm::interact.gbm()`
try_load <- suppressWarnings(try(load(file.path(root, "BAMdb-GNMsubset-2020-01-08.RData")), silent=TRUE))


# check if `load` was successful and V4 data exists
if (inherits(try_load, "try-error") || !exists("dd")){
  message("could not V4 data from: ", root)
} else {
  message("successfully loaded V4 data from: ", root)
} #close else()

# FOR TESTING:
load(file.path("G:/Shared drives/BAM_NationalModels4/NationalModels4.0/Feb2020/out/boot/BANS/BCR_4/gnmboot-BANS-BCR_4-1.RData"))

#5. create function to estimate 2-way interactions for all covs in each gbm object----
process_gbm <- function(obj_path) {
  
  
  # load a bootstrap replicate (gbm object)
  try_load <- suppressWarnings(try(load(file.path(obj_path)), silent = TRUE))
  
  # check if `load` was successful and `out` exists
  if (inherits(try_load, "try-error") || !exists("out")) {
    message("could not load gbm_obj from: ", obj_path)
  } else {
    message("successfully loaded gbm_obj from: ", obj_path)
  } #close else()
  flush.console()
  
  
  # build `DAT` for the target species x bcr
  # define species and region
  spp <- attr(out, "__settings__")$species
  bcr <- attr(out, "__settings__")$region
  
  # create an index indicating where in `dd` we will find the current BCR's data
  # ss = subset
  ss <- dd[, bcr] == 1L
  
  # reconstruct DAT object (see: https://github.com/borealbirds/GNM/blob/master/R/04-boot.R)
  # `yy` contains the counts data, and is subset using `ss` and `spp` 
  DAT <- 
    data.frame(
      count = as.numeric(yy[ss, spp]),  
      offset = off[ss, spp],
      cyid = dd$cyid[ss],
      YEAR = dd$YEAR[ss],
      ARU = dd$ARU[ss],   
      dd2[ss, out$var.names[-c(1:2)]]  #-c(1:2) because we already have `YEAR` and `ARU` from `dd`
    )
  
  # keep just one observation (row) for each unique "cell x year" combination
  DAT <- DAT[sample.int(nrow(DAT)),]
  DAT <- DAT[!duplicated(DAT$cyid),]
  
  # ensure `DAT` has the covariates found in `out`
  required_vars <- out$var.names
  if (!all(required_vars %in% colnames(DAT))) {
    stop("one or more covariates in `out` are missing from `DAT`")
  } else {
    message("check: all covariates in `out` are in `DAT`")
  } 
  flush.console()
  
  # ensure `DAT` has the covariates found in `bam_covariate_importance_v4`, 
  # which only includes covariates with rel.inf > 1
  current_spp <- stringr::str_extract(basename(obj_path), "(?<=gnmboot-)[^-]+")
  current_bcr <- stringr::str_extract(basename(obj_path), "BCR_\\d+")
  
  influential_covs <- 
    bam_covariate_importance_v4 |> 
    dplyr::filter(spp == current_spp & bcr == current_bcr) |> 
    dplyr::pull(var)
  
  if (!all(required_vars2 %in% influential_covs)) {
    stop("one or more covariates in `bam_covariate_importance_v4` are missing from `DAT`")
  } else {
    message("check: all covariates in `bam_covariate_importance_v4` are in `DAT`")
  } 
  flush.console()
  
  
  # find evaluation points (`data.frame`) for every covariate permutation of degree 2 (indexed by i,j) 
  # NOTE: we limit covariates to those with a stand-alone rel.inf > 1, under the assumption that 
  # covariates with very low stand-alone rel.inf won't have strong 2-way interactions  # 
  pts <- list()
  n <- length(influential_covs); message("check: `out` has ", n, " covariates")
  interaction_index <- 1 # starts at 1 and increases for every covariate interaction computed. Resets at one when moving to the next bootstrap model.
  flush.console()
  
  # end at n-1 to avoid finding the interaction of variable n x variable n
  for (i in 1:(n-1)) {
    message("computing all interactions with ", influential_covs[i])
    flush.console()
    
    # start at i+1 to avoid finding the interaction of variable 1 x variable 1
    for (j in (i+1):n) {
      
      # set seed for reproducibility of the data subsampling procedure
      set.seed(interaction_index)
      
      # check if the response data is non-empty
      if (sum(DAT$count) < 1) {
        pts[[interaction_index]] <- NA  # assign NA if there are no occurrence records
        
      } else {
        
        # subsample 25% of the data using `slice_sample` to speed up `interact.gbm`
        DAT_sample <- 
          dplyr::slice_sample(DAT, prop = 0.25) |> 
          dplyr::select(all_of(influential_covs)) # select only the necessary covariates (i.e those in with high rel.inf)
        
        # check for variability in the selected covariates
        var_i_unique <- length(unique(DAT_sample[[influential_covs[i]]]))
        var_j_unique <- length(unique(DAT_sample[[influential_covs[j]]]))
        
        if (var_i_unique < 2 || var_j_unique < 2) {
          message("skipping interaction for ", influential_covs[i], " and ", influential_covs[j], " due to insufficient variability.")
          pts[[interaction_index]] <- NA
          flush.console()
        } else {
          
          
          pts[[interaction_index]] <- gbm::interact.gbm(x = out, data = DAT_sample, i.var = c(i, j))  # test the interaction between variable i and j
          
        } # close nested else()
        
      } # close top else ()
      
      # label the interaction
      names(pts)[interaction_index] <- paste(out$var.names[i], out$var.names[j], sep = ".")  # Label the interaction
      interaction_index <- interaction_index + 1
    } # close nested loop
    
  } # close top loop
  
  return(pts)  # return the list of interaction points
  
} # close function