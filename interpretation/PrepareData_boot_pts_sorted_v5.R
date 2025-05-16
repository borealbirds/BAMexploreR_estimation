# ---
# title: National Models 4.0 - estimating and ranking 2-way covariate interactions for CONWA and CAWA
# author: Mannfred Boehm
# created: September 13, 2024
# ---

#NOTES################################


#1. attach packages----

print("* attaching packages on master *")
library(gbm)
library(parallel)
library(tidyverse)


#2. set up cluster----

# create and register clusters
print("* creating clusters *")
n_tasks <- 32
cl <- makePSOCKcluster(n_tasks, type="PSOCK")

# print number of tasks and host name for confirmation
cl

# set root path
root <- "/home/mannfred/scratch"
print(paste("* currently working out of ", root, " *"))


#3. stage data

gbm_objs <- list.files(path = root, pattern = "\\.Rdata$", full.names = TRUE, recursive = TRUE)

# for testing
gbm_objs <- gbm_objs[1:3]
gbm_objs

# import covariate importance data
load(file="bam_covariate_importance_v5.rda")
bam_covariate_importance_v5

# process_gbm takes a path to a .Rdata file containing a list of bootstrap models (`b.list`)
# and returns a nested list:
# top-level: each element is a bootstrap replicate (i.e., each model in `b.list`)
# second-level list: one `data.frame` per covariate, containing the evaluation grid from `plot.gbm()`
process_gbm <- function(obj_path) {
  
  load(obj_path)  # loads b.list
  
  boot_pts <- list()
  
  
  # find evaluation points (`data.frame`) for every covariate (indexed by `r`)
  # this can sped up by reducing `continuous.resolution=100` in `plot.gbm()`
  for (b in seq_along(b.list)) {
    model <- b.list[[b]]
    
    pts <- list()
    for (r in seq_along(model$var.names)) {
      pts[[r]] <- plot.gbm(x = model, return.grid = TRUE, i.var = r, type = "response")
    }
    
    boot_pts[[b]] <- pts
  }
  
  return(boot_pts)
}


#6. export the necessary variables and functions to the cluster----
print("* exporting gbm_objs and functions to cluster *")

clusterEvalQ(cl, {
  library(gbm)
  library(dplyr)
  library(stringr)
})

clusterExport(cl, c("gbm_objs", "process_gbm", "bam_covariate_importance_v5"))


#7. run the function in parallel----
print("* running `process_gbm` in parallel *")
boot_pts <- parLapply(cl = cl, X = gbm_objs, fun = process_gbm)


#8. stop the cluster----
print("* stopping cluster :-)*")
stopCluster(cl)


#9. save boot_pts----
print("* saving boot_pts_v5.rds  *")
saveRDS(boot_pts, file=file.path(root, "boot_pts_v5.rds"))


# from the family of sets (BCRs, species, covariates) 
# group_keys() finds the names of the unique tuples
boot_group_keys <- 
  bam_covariate_importance_v5 |> 
  dplyr::group_by(bcr, spp, var) |> 
  dplyr::group_keys()


# create an index from `gbm_objs` containing the 
# species (FLBC), BCR, and bootstrap replicate
sample_id <- 
  tibble::tibble(
    spp = str_extract(gbm_objs, "/([^/]+)/[^/]+\\.Rdata$") |> str_extract("[A-Z]{4}"),
    bcr = str_extract(gbm_objs, "[a-z]{3}[0-9]{2}")
  ) |> 
  dplyr::slice(rep(1:n(), each = 32)) |>  # 32 models per file
  dplyr::mutate(boot = rep(1:32, times = length(gbm_objs)))


# check that sample_id and boot_pts have same length
stopifnot(length(sample_id) == length(boot_pts))


# for every zth bcr x species x var tuple (rows in `boot_group_keys`):
# search for the relevant bootstrap predictions in `boot_pts` by matching species and bcr at the top level of `boot_pts` then matching `var` at the second level
# then, enter bootstrap dataframes for a given bcr x species x var tuple as a sub-element of [[z]]
# output: a list of lists where the top level elements are species x bcr x var tuples, 
# and the second-level elements are dataframes of that tuple's bootstrapped model predictions
# you can query every zth tuple for its bcr x species x var identity by: names(boot_pts_sorted)[z] 
boot_pts_sorted <- list()
for (z in 1:nrow(boot_group_keys)){
  
  # zth bcr x species x var tuple
  key_z <- boot_group_keys[z,]
  
  # `sample_id` and `boot_pts` have the same length and order (they are derived from `gbm_objs`)
  # so we can identify elements in the latter using info from the former
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

saveRDS(boot_pts_sorted, file=file.path(root, "boot_pts_sorted_v5.rds"))
