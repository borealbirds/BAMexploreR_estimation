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
test <- FALSE
cc <- TRUE

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

if(!cc){root <- "G:/Shared drives/BAM_NationalModels4/NationalModels4.0"}
if(cc){root <- "/home/mannfred/scratch"}

print(paste("* currently working out of ", root, " *"))



#3. create indices needed for searching every 2-way interaction for available models

# create an index of every bcr x spp x 2-way interactions 
# loads as "bam_covariate_importance_v4"
load(file.path(root, "v4_bootstraps", "bam_covariate_importance_v4.rda"))


# create an index of all available bootstrap models
gbm_objs <- list.files(file.path(root, "v4_bootstraps"), pattern = "^gnmboot-.*\\.RData$", full.names = TRUE)
print(paste("* Found ", length(gbm_objs), " files: *"))



#4. load in V4 count data---- 
# needed for running `gbm::interact.gbm()`
try_load <- suppressWarnings(try(load(file.path(root, "v4_bootstraps", "BAMdb-GNMsubset-2020-01-08.RData")), silent=TRUE))


# check if `load` was successful and V4 data exists
if (inherits(try_load, "try-error") || !exists("dd")){
  print(paste("* could not V4 data from: ", root, " *"))
} else {
  print(paste("* successfully loaded BAMdb-GNMsubset data from: ", root, " *"))
} #close else()

# for testing:
# load("G:/Shared drives/BAM_NationalModels4/NationalModels4.0/Feb2020/out/boot/BANS/BCR_4/gnmboot-BANS-BCR_4-1.RData")

#5. create function to estimate 2-way interactions for all covs in each gbm object----
process_gbm <- function(obj_path) {
  
  
  # load a bootstrap replicate (gbm object)
  try_load <- suppressWarnings(try(load(file.path(obj_path)), silent = TRUE))
  
  # check if `load` was successful and `out` exists
  # if `load` isn't successful, return NULL. 
  if (inherits(try_load, "try-error") || !exists("out")) {
    message("could not load gbm_obj from: ", obj_path)
    return(NULL) 
  } else {
    message("successfully loaded gbm_obj from: ", obj_path)
  } #close else()
  flush.console()
  
  # check if `__settings__` attribute exists and is a list
  settings <- attr(out, "__settings__")
  if (!is.list(settings)) {
    message("missing `__settings__` attribute in out. skipping this bootstrap model." )
    return(NULL)
  }  
  
  # build `DAT` for the target species x bcr
  # define species and region
  spp <- settings$species
  bcr <- settings$region
  
  # check if species and bcr info exists
  if (is.null(spp) || is.null(bcr)) {
    message("species or bcr missing from `__settings__`. skipping this bootstrap model.")
    return(NULL)
  }
  

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

  # ensure `DAT` has the covariates found in `bam_covariate_importance_v4`, 
  # which only includes covariates with rel.inf > 1
  influential_covs <- 
    bam_covariate_importance_v4 |> 
    dplyr::filter(spp == !!spp & bcr == !!bcr) |> 
    dplyr::pull(var)
  
  # check that at least two high `rel.inf` covariates exist for this spp x bcr
  if (length(influential_covs) < 2) {
    print(paste("<2 influential covariates exist for ", spp, bcr, ": creating a NULL entry"))
    return(NULL)
  }

  # check that all high `rel.inf` covariates are in current gbm object 
  if (length(influential_covs[!influential_covs %in% out$var.names]) > 0) {
    message("warning: covariates are missing from out$var.names: ", paste(influential_covs[!influential_covs %in% out$var.names], collapse=", "))
  }

  # check that all high `rel.inf` covariates are in occurrence dataset
  if (!all(influential_covs %in% colnames(DAT))) {
    stop("one or more covariates in `bam_covariate_importance_v4` are missing from `DAT`")
  } else {
    message("check: all covariates in `bam_covariate_importance_v4` are in `DAT`")
  }   
  
  # find evaluation points (`data.frame`) for every covariate permutation of degree 2 (indexed by i,j) 
  # NOTE: we limit covariates to those with a stand-alone rel.inf > 1, under the assumption that 
  # covariates with very low stand-alone rel.inf won't have strong 2-way interactions  # 
  pts <- list()
  n <- length(influential_covs); message("check: `out` has ", n, " covariates with rel.inf > 1")
  interaction_index <- 1 # starts at 1 and increases for every covariate interaction computed. Resets at one when moving to the next bootstrap model.
  
  # identify indices in out$var.names that have high `rel.inf` 
  influential_indices <- match(influential_covs, out$var.names)
  
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
        
        # subsample 50% of the data using `slice_sample` to speed up `interact.gbm`
        DAT_sample <- 
          dplyr::slice_sample(DAT, prop = 0.50) |> 
          dplyr::select(all_of(influential_covs)) # select only the necessary covariates (i.e those in with high rel.inf)
        
        # check for variability in the selected covariates
        var_i_unique <- length(unique(DAT_sample[[influential_covs[i]]]))
        var_j_unique <- length(unique(DAT_sample[[influential_covs[j]]]))
        
        if (var_i_unique < 2 || var_j_unique < 2) {
          
          message("skipping interaction for ", influential_covs[i], " and ", influential_covs[j], " due to insufficient variability.")
          pts[[interaction_index]] <- NA
        } else if (!all(influential_covs[c(i, j)] %in% colnames(DAT_sample))) {

          message("Skipping interaction due to missing variables: ", influential_covs[i], " and ", influential_covs[j])
          pts[[interaction_index]] <- NA
        } else if (nrow(DAT_sample) == 0) {
          
          message("Skipping interaction due to empty DAT_sample")
          pts[[interaction_index]] <- NA  
        } else {  
          
          # tell `i.var=` to only use indices in `out$var.names` that have high rel.inf
          pts[[interaction_index]] <- gbm::interact.gbm(x = out, data = DAT_sample, i.var = influential_indices[c(i, j)])  # test the interaction between variable i and j
          
        } # close nested else()
        
      } # close top else ()
      
      # label the interaction
      names(pts)[interaction_index] <- paste(influential_covs[i], influential_covs[j], sep = ".")  # Label the interaction
      interaction_index <- interaction_index + 1
    } # close nested loop
    
  } # close top loop
  
  return(pts)  # return the list of interaction points
  
} # close function



#6. export the necessary variables and functions to the cluster----
print("* exporting gbm_objs and functions to cluster *")
clusterExport(cl, c("gbm_objs", "interact.gbm", "process_gbm", "dd", "dd2", "yy", "off", "bam_covariate_importance_v4"))



#7. run the function in parallel----
print("* running `process_gbm` in parallel *")
boot_pts_i2 <- parLapply(cl = cl, X = gbm_objs, fun = process_gbm)


#8. stop the cluster----
print("* stopping cluster :-)*")
stopCluster(cl)

#9. save list of 2-way interactions per bootstrap----
print("* saving boot_pts_i2_v4.rds  *")
saveRDS(boot_pts_i2, file=file.path(root, "boot_pts_i2_v4.rds"))



#10. sort 2-way interactions by species x bcr tuples----

# create an index from `gbm_objs` containing the species (FLBC), BCR, and bootstrap replicate
# NOTE: this has the same order and length as `gbm_objs`
sample_id <- 
  gbm_objs |> 
  sub("^.*gnmboot-", "", x = _) |> 
  sub("\\.RData", "", x = _) |>
  stringr::str_split_fixed(pattern="-", n=3) |> 
  tibble::as_tibble() |> 
  magrittr::set_colnames(c("spp", "bcr", "boot"))

# create an index of every possible permutation of var_1, var_2, bcr, and spp 
# `RcppAlgo::comboGrid` is like `expand.grid` and `tidyr::crossing` but avoids duplicates 
# e.g. for our purposes (var1=x, var2=y) is a duplicate of (var1=y, var2=x), so the Cartesian Product of comboGrid is smaller
boot_group_keys_i2 <- 
  RcppAlgos::comboGrid(unique(bam_covariate_importance_v4$var), 
                       unique(bam_covariate_importance_v4$var), 
                       bam_covariate_importance_v4$bcr, 
                       bam_covariate_importance_v4$spp, repetition =  FALSE) |> 
  tibble::as_tibble() |> 
  dplyr::rename(var_1=Var1, var_2=Var2, bcr=Var3, spp=Var4) |> 
  dplyr::semi_join(x=_, y=sample_id, by = c("spp", "bcr")) # keep only species x bcr tuples that exist in sample_id



# for every zth species x bcr x 2-way interaction tuple (rows in `boot_group_keys_i2`):
# gather the relevant bootstrap predictions from `boot_pts_i2`
# (each element of `boot_pts_i2` has many `y_means`; one for every 2-way interaction)
# and enter each wth `y_mean` (that matches the current 2-way interaction of interest) 
# as a sub-element of [[z]]
# output: a list of lists where the top level elements (n = nrow(boot_group_keys_i2)) 
# are species x bcr x 2-way interaction permutations, 
# the second-level elements are matrices of the associated average and sd 
# of the bootstrapped prediction spaces

# NOTE: `sample_id` must be created in `01_prepare_data.R`
boot_pts_sorted_i2 <- list()
for (z in 1:nrow(boot_group_keys_i2)){
  
  # zth bcr x species x 2-way interaction permutation
  key_z <- boot_group_keys_i2[z,]
  
  # problem: which elements in boot_pts_i2 match the species x bcr tuple defined in key_z?
  # approach: `sample_id` and `boot_pts_i2` have the same length and order (they both are derived from gbm_objs)
  # so we can identify elements in `boot_pts_i2` using info from `sample_id`
  # first, identify a unique species x bcr tuple, then gather covariate interactions into it
  boot_pts_index <- which(sample_id$bcr == key_z$bcr & sample_id$spp == key_z$spp)
  if (length(boot_pts_index) == 0) {
    next  # Skip this iteration if no matching species × BCR is found
  }
  
  # create a list of bootstrap predictions for zth species x bcr permutation 
  # this works because `sample_id` has the same length and order as 
  # the top-level elements of `boot_pts_i2`
  spp_bcr_list <- boot_pts_i2[boot_pts_index]
  
  # gather all relevant 2-way interaction bootstrap averages for the zth species x bcr permutation
  spp_bcr_var_i2 <- list()
  for (w in 1:length(spp_bcr_list)) {
    
    # get covariate names for the current spp x bcr permutation 
    # NOTE: this is admittedly wasteful because we're only interested in the interactions within `key_z`
    i2_names <- names(spp_bcr_list[[w]])
    
    # find which sub-element of spp_bcr_list[[w]] match the current 2-way interaction of interest
    # account for the possibility of being indexed as x*y or y*x
    if (paste(key_z$var_1, key_z$var_2, sep=".") %in% i2_names |
        paste(key_z$var_2, key_z$var_1, sep=".") %in% i2_names){
      
      matching_index <- 
        which(i2_names == paste(key_z$var_1, key_z$var_2, sep=".") | 
              i2_names == paste(key_z$var_2, key_z$var_1, sep="."))
      
      # assign current bcr x species x 2-way interaction x bootstrap tuple as a top-level element of `spp_bcr_var_i2[[w]]`
      # these bootstraps will eventually be gathered under a single bcr x species x 2-way interaction element in `boot_pts_sorted_i2`
      spp_bcr_var_i2[w] <- spp_bcr_list[[w]][matching_index]
      names(spp_bcr_var_i2)[w] <- paste("bootstrap replicate", w, sep="_")
      
    } #close if()
    
  } #close nested loop
  
  # take mean across bootstraps of the zth bcr x spp x 2-way interaction tuple
  if (purrr::is_empty(spp_bcr_var_i2) == FALSE){
    mean_z <- 
      spp_bcr_var_i2 |> 
      unlist() |> 
      tibble(value = _) |> 
      dplyr::summarise(mean = mean(value), sd = sd(value)) 
      
    
    if (!is.na(mean_z$mean)) {
      boot_pts_sorted_i2[[z]] <- mean_z
    } else {
      boot_pts_sorted_i2[[z]] <- list()
    } #close nested else
    
  } # close top level if()
  else {
    boot_pts_sorted_i2[[z]] <- list()
  } # close top level else
  
  names(boot_pts_sorted_i2)[z] <- paste(key_z$bcr, key_z$spp, key_z$var_1, key_z$var_2, sep=".")
  
  # print progress
  if (z %% 100000 == 0 || z == 1 || z == nrow(boot_group_keys_i2)) {
    message(paste("Iteration", z, "of", nrow(boot_group_keys_i2)))
  }
  
} # close top level loop



#11. tidy results into a tibble----

# import extraction lookup table to obtain covariate classes
# (for appending to covariate importance data)
nice_var_names <-
  readr::read_csv(file.path(root, "v4_bootstraps", "nice_var_names.csv")) |>
  dplyr::select(var_class, var)

boot_pts_reduced_i2 <- 
  boot_pts_sorted_i2 |> 
  purrr::discard(purrr::is_empty) |>  # remove empty elements
  purrr::imap_dfr(~ tibble(name = .y, mean = .x$mean, sd = .x$sd)) |> 
  tidyr::extract(name, into = c("bcr", "spp", "covariate_1", "covariate_2"), regex = "([^\\.]+)\\.([^\\.]+)\\.(.+)\\.(.+)$") |> 
  dplyr::arrange(desc(mean)) 
 # dplyr::left_join(nice_var_names, by = c("covariate_1" = "var")) |> 
 # dplyr::left_join(nice_var_names, by = c("covariate_2" = "var")) |> 
 # dplyr::select(var_nice.x, var_nice.y, mean, sd) 

print("* saving boot_pts_reduced_i2 *")  
saveRDS(boot_pts_reduced_i2, file=file.path(root, "bam_2way_interactions_v4.rds"))

