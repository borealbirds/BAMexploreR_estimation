# ---
# title: National Models 5.0 - preparing data for partial dependence plot function
# author: Mannfred Boehm
# created: September 13, 2024
# ---

#NOTES################################


#1. attach packages----
library(gbm)
library(parallel)
library(tidyverse)


#2. setup local or cluster----
test <- FALSE
cc <- TRUE

# set number of tasks for local vs cluster
if(cc){ n_tasks <- 32}
if(!cc | test){ n_tasks <- 4}

# create and register clusters
# creates `n_tasks` copies of R running in parallel via 32 tasks, on one of the cluster's sockets (processors). 
# Belgua has ~965 nodes
print("* creating clusters *")
cl <- makePSOCKcluster(n_tasks, type="PSOCK")

# print number of tasks and host name for confirmation
cl

# set root path
print("* setting root file path *")

if(!cc){root <- "G:/Shared drives/BAM_NationalModels5/output/06_bootstraps"}
if(cc){root <- "/home/mannfred/projects/def-bayne/NationalModels/06_bootstraps"}

print(paste("* currently working out of ", root, " *"))


# attach packages on clusters
# `clusterEvalQ` evaluates a literal expression on each cluster node. 
print("* Loading packages on workers *")
invisible(clusterEvalQ(cl, library(gbm)))
invisible(clusterEvalQ(cl, library(tidyverse)))



#3. access gbm objects and append spp/bcr/boot info----
gbm_objs <- list.files(file.path(root), pattern = "*\\.Rdata$", full.names = TRUE, recursive = TRUE)
if(!cc | test) {gbm_objs <- gbm_objs[1:2]}

print(paste("* Found", length(gbm_objs), " files *"))
print(head(gbm_objs))

# import covariate importance data from "PrepareDataCovariateImportanceV5.R"
# load(file="C:/Users/mannf/Proton Drive/mannfredboehm/My files/Drive/BAMexploreR/data-raw/bam_covariate_importance_v5.rda")
load(file="/home/mannfred/scratch/bam_covariate_importance_v5.rda")
print("* loading covariate importance data *")
head(bam_covariate_importance_v5)


# group_keys() finds the names of the unique species x bcr x var tuples in `bam_covariate_importnance`
boot_group_keys <- 
  bam_covariate_importance_v5 |> 
  dplyr::group_by(bcr, spp, var) |> 
  dplyr::group_keys()


# evaluate_grid_from_model takes a path to a .Rdata file containing a list of bootstrap models (`b.list`)
# and returns a nested list:
# top-level: each element is a species x bcr tuple 
# second-level: 32 elements, each element is a list representing a bootstrap
# third-level: each element is a `data.frame` containing a covariate's evaluation grid from `plot.gbm()`

evaluate_grid_from_model <- function(obj_path, target_vars) {
  
  # load list of 32 bootstraps for current species x bcr tuple
  load(obj_path)  
  
  # index current species x bcr tuple for searching in `bam_covariate_importance`
  spp <- stringr::str_extract(obj_path, "[A-Z]{4}")
  bcr <- stringr::str_extract(obj_path, "can[0-9]+")
  
  # match the current species x bcr in `bam_covariate_importance` and
  # filter out any covariates that aren't associated with the current tuple
  relevant_vars <- 
    target_vars |> 
    dplyr::filter(spp == !!spp, bcr == !!bcr) |> 
    dplyr::pull(var) |> 
    unique()
  
  
  
  # find evaluation points (`data.frame`) for every covariate (indexed by `r`)
  # this can sped up by reducing `continuous.resolution=100` in `plot.gbm()`
  boot_pts <- list()
  
  for (b in seq_along(b.list)) {
    model <- b.list[[b]]
    
    pts <- list()
    for (r in seq_along(model$var.names)) {
      pts[[r]] <- plot.gbm(x = model, return.grid = TRUE, i.var = r, type = "response")
    }
    
    names(pts) <- model$var.names
    boot_pts[[b]] <- pts
  }
  
  names(boot_pts) <- paste0("bootstrap_", seq_along(b.list))
  
  return(boot_pts)
}


#6. export the necessary variables and functions to the cluster----
print("* exporting gbm_objs and functions to cluster *")

clusterEvalQ(cl, {
  library(gbm)
  library(dplyr)
  library(stringr)
})

clusterExport(cl, c("gbm_objs", "evaluate_grid_from_model", "bam_covariate_importance_v5"))


#7. run the function in parallel----
print("* running `process_gbm` in parallel *")
boot_pts <- parLapply(cl = cl, X = gbm_objs, fun = evaluate_grid_from_model)


#8. stop the cluster----
print("* stopping cluster :-)*")
stopCluster(cl)


#9. save boot_pts----
print("* saving boot_pts_v5.rds  *")
saveRDS(boot_pts, file="/home/mannfred/scratch/boot_pts_v5.rds")



# create an index from `gbm_objs` containing the 
# species, BCR, and bootstrap replicate
sample_id <- 
  tibble::tibble(
    spp = str_extract(gbm_objs, "/([^/]+)/[^/]+\\.Rdata$") |> str_extract("[A-Z]{4}"),
    bcr = str_extract(gbm_objs, "[a-z]{3}[0-9]{2}")
  ) |> 
  dplyr::slice(rep(1:n(), each = 32)) |>  # 32 models per file
  dplyr::mutate(boot = rep(1:32, times = length(gbm_objs)))


# check that sample_id and boot_pts have same length
stopifnot(nrow(sample_id) == length(boot_pts)*32)


# for every zth bcr x species x var tuple (rows in `boot_group_keys`):
# search for the relevant bootstrap predictions in `boot_pts` by matching species and 
# bcr at the top level of `boot_pts` then matching `var` at the second level
# then, enter bootstrap dataframes for a given bcr x species x var tuple as a sub-element of [[z]]
# output: a list of lists where the top level elements are species x bcr x var tuples, 
# and the second-level elements are dataframes of that tuple's bootstrapped model predictions
# you can query every zth tuple for its bcr x species x var identity by: names(boot_pts_sorted)[z] 
boot_pts_sorted <- list()
for (z in 1:nrow(boot_group_keys)){
  
  # zth bcr x species x var tuple
  key_z <- boot_group_keys[z,]
  
  # `sample_id` and `boot_pts` have the same length and order (they are derived from `gbm_objs`)
  # so we can identify the species x bcr in the latter using info from the former
  # search for all bootstraps for a given bcr x sample tuple
  boot_pts_index <- which(sample_id$bcr == key_z$bcr & sample_id$spp == key_z$spp)
  
  # a list of bootstrap predictions for species x bcr tuple
  spp_bcr_list <- boot_pts[boot_pts_index]
  
  # for a given bcr x species tuple, search within every bootstrap model for
  spp_bcr_var <- list()
  for (w in 1:length(spp_bcr_list)) {
    
    # get all possible covariate names for the current bcr x species bootstrap
    var_names <- 
      lapply(spp_bcr_list[[w]], colnames) |> 
      lapply(X=_, `[[`, 1) |> #don't need the name of the y variable
      purrr::flatten_chr()
    
    # search the wth bootstrap to find the current covariate of interest
    spp_bcr_var[w] <- spp_bcr_list[[w]][which(var_names %in% key_z$var)]
    names(spp_bcr_var)[w] <- paste("bootstrap replicate", w, sep="_")
  }
  
  boot_pts_sorted[[z]] <- spp_bcr_var
  names(boot_pts_sorted)[z] <- paste(boot_group_keys[z, "bcr"], boot_group_keys[z, "spp"], boot_group_keys[z, "var"], sep = "_")
  
  # print progress
  cat(paste("\riteration", z))
  Sys.sleep(0.001)
}



boot_pts_sorted <- 
  
  # spp_bcr_entry inputs a list of 32 bootstraps for a species x bcr tuple (top level of boot_pts)
  # key_prefix is the name of species x bcr tuple
  purrr::imap(boot_pts, function(spp_bcr_bootstraps, spp_bcr_name) {
    
    # for every bootstrap (for a given spp x bcr tuple)
    # define a new function that requires a covariate dataframe
    purrr::imap(spp_bcr_bootstraps, function(vars_list, boot_name) {
    
      # for every covariate build a list:
      # the top element is the current species x bcr x covariate
      # the second-level elements are dataframes, one per 32 bootstraps
      purrr::imap(vars_list, function(df, var) {
      
        list(key = paste(spp_bcr_name, var, sep = "_"),
             boot = boot_name,
             df = df)
      
      # flatten() each level to turn nested lists into a single vector of lists
      }) 
  
    }) |> purrr::flatten()

  }) |> flatten()

# Then group by tuple key:

boot_pts_sorted <- 
  boot_pts_sorted |> 
  purrr::flatten() |> 
  split(x = _, f = purrr::map_chr(., "key")) |> 
  purrr::map(~ purrr::set_names(purrr::map(., "df"), nm = paste0("bootstrap_", seq_along(.))))

boot_pts_sorted <- split(map(boot_pts_sorted, "df"),
                         map_chr(boot_pts_sorted, "key"))


saveRDS(boot_pts_sorted, file=file.path(root, "boot_pts_sorted_v5.rds"))
