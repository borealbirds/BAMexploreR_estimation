# ---
# title: National Models 5.0 - data organization for analysing relative influence of covariates and species traits
# author: Mannfred Boehm
# created: January 13, 2025
# ---

#NOTES################################

# This script extracts covariate contributions to model predictions (VERSION 5).

# The outputs will be analysed for identifying covariates of importance across several metrics (e.g. BCR, ecology, etc.).

# This script generates the exported data for R package functions that summarise covariate performance. 


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
# creates `n_tasks` copies of R running in parallel via 32 tasks, on one of the cluster's sockets (processors). 
# Belgua has ~965 nodes
print("* creating clusters *")
cl <- makePSOCKcluster(n_tasks, type="PSOCK")

# print number of tasks and host name for confirmation
cl


# set root path
print("* setting root file path *")

if(!cc){root <- "G:/Shared drives/BAM_NationalModels5/output/bootstraps"}
if(cc){root <- "/home/mannfred/scratch"}

print(paste("* currently working out of ", root, " *"))

# attach packages on clusters
# `clusterEvalQ` evaluates a literal expression on each cluster node. 
print("* Loading packages on workers *")
invisible(clusterEvalQ(cl, library(gbm)))
invisible(clusterEvalQ(cl, library(tidyverse)))



#3. access gbm objects and append spp/bcr/boot info----
gbm_objs <- list.files(file.path(root), pattern = "*\\.R", full.names = TRUE, recursive = TRUE) |> 
(\(x) x[1:3])()
print(paste("* Found", length(gbm_objs), " files *"))
print(head(gbm_objs))




#4. extract the species (FLBC), BCR, and bootstrap replicate from `gbm_objs`----
sample_id <- 
  gbm_objs |> 
  basename() |> 
  stringr::str_split_fixed(pattern="_", n=3) |> 
  gsub("\\.R", "", x = _) |>
  tibble::as_tibble() |> 
  magrittr::set_colnames(c("spp", "bcr", "boot"))

if (exists("sample_id")) {
  print("* successfully constructed sample_id *")
} else {message("* sample_id not constructed *")}



#5. create master dataframe with indexing info----
# then, convert to list with each dataframe row as a list element 
# (easier to deal with a list for looping/`apply`ing)
gbm_data <- tibble(file_path = gbm_objs, spp = sample_id$spp, bcr = sample_id$bcr, boot = sample_id$boot) 
gbm_list <- split(x=gbm_data, f=seq(nrow(gbm_data)))




#6. define function that creates `bam_covariate_importance` data frame----


# this function works on the `gbm_data` data frame:
# for every gbm object (row) in `gbm_data` it creates a 
# dataframe describing relative influence per covariate for the given gbm object
# the output is a list with elements as rel. inf. dataframes 
compute_var_importance <- function(gbm_object) {
  
  
  tryCatch({
    
    # attempt to load the GBM object
    load(gbm_object[["file_path"]])
    
    # validate gbm object
    if (!exists("b.i") || is.null(b.i$n.trees) || b.i$n.trees <= 0) {
      warning(paste("Invalid GBM model in file:", gbm_object["file_path"]))
      return(NULL)
    } else {
      message("Successfully loaded: ", gbm_object["file_path"])
    }
    
    # compute variable importance
    # to reduce file size, filter for rel.inf >= 1
    importance_df <-
      b.i |>
      gbm::summary.gbm(plotit = FALSE) |>
      as_tibble() |>
      dplyr::bind_cols(x = _, c(gbm_object["spp"], gbm_object["bcr"], gbm_object["boot"]))
    
    message("Successfully extracted `rel.inf` from: ", gbm_object["file_path"])
    return(importance_df)
    
    # remove loaded gbm object to free up memory
    rm(out)
    gc()
    
    # give instructions to tryCatch if it runs into an error
  }, error = function(e) {
    warning(paste("Error processing", gbm_object["file_path"], ":", e$message))
    return(NULL) 
    
  }) # close tryCatch
  
} # close compute_var_importance()



#7. export the necessary variables and functions to the cluster----
print("* exporting objects and functions to cluster *")
clusterExport(cl, c("root", "gbm_list", "compute_var_importance"))


#8. run the covariate importance function in parallel----
print("* running `compute_var_importance` in parallel` *")
bam_covariate_importance_v5 <- parLapply(cl, gbm_list, compute_var_importance)
saveRDS(bam_covariate_importance_v5, file=file.path(root, "bam_covariate_importance_v5.rds"))


#9. stop the cluster----
print("* stopping cluster :-)*")
stopCluster(cl)



#10. post-cluster data cleanup (executed on local machine)----

# remove NULL (invalid model) entries
bam_covariate_importance_nonull_v5 <- bam_covariate_importance_v5[!sapply(bam_covariate_importance_v5, is.null)]

# reduce list of dataframes into a single dataframe
# took about 15 minutes with a list of 12808 tibbles
covariate_importance_merged <- suppressMessages(purrr::reduce(bam_covariate_importance_nonull_v5, full_join))

# THIS IS V4: GET V5!!!
# import extraction lookup table to obtain covariate classes----
# (for appending to covariate importance data)
# lookup table is missing "Year" and "Method", so manually adding here
nice_var_names <-
  readr::read_csv(file.path(root, "nice_var_names_v5.csv")) |>
  dplyr::select(var_class, var)

# check for missing bootstraps and fill in zeros where covariates 
# are missing from a bootstrap
covariate_importance_zeroed <-
  covariate_importance_merged |>
  dplyr::mutate(boot = as.integer(boot)) |> # change `chr` to `int`
  dplyr::group_by(spp, bcr, var) |>
  dplyr::summarise(mean_rel_inf = mean(rel.inf, na.rm = TRUE),  # calculate mean across all bootstraps (including 0s)
                   sd_rel_inf = sd(rel.inf, na.rm = TRUE),      # calculate standard deviation across all bootstraps (including 0s)
                   n_boots = sum(rel.inf > 0)) |>    # count the number of bootstraps where the variable had non-zero rel.inf
  dplyr::filter(mean_rel_inf >= 1) |> # filter out spp x bcr x var tuples with mean rel.inf <= 1
  dplyr::left_join(nice_var_names, by = "var") |> # append more interpretable variable names
  dplyr::arrange(spp, bcr, desc(mean_rel_inf))

saveRDS(covariate_importance_zeroed, file=file.path(root, "bam_covariate_importance_v5.rds"))


if(cc){ q() }
