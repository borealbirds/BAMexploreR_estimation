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



#4. extract the species (FLBC), BCR, and bootstrap replicate from `gbm_objs`----
sample_id <- 
  tibble::tibble(
    spp = str_extract(gbm_objs, "/([^/]+)/[^/]+\\.Rdata$") |> 
          str_extract("[A-Z]{4}"),
    bcr = str_extract(basename(gbm_objs), "[a-z]{3}[0-9]{1,2}")
  ) 

if (exists("sample_id")) {
  print("* successfully constructed sample_id *")
} else {message("* sample_id not constructed *")}

head(sample_id)


#5. create master dataframe with indexing info----
# then, convert to list with each dataframe row as a list element 
# (easier to deal with a list for looping/`apply`ing)
gbm_data <- tibble(file_path = gbm_objs, spp = sample_id$spp, bcr = sample_id$bcr)
gbm_list <- split(x=gbm_data, f=seq(nrow(gbm_data)))



#6. define function that creates `bam_covariate_importance` data frame----

# this function works on the `gbm_data` data frame:
# for every gbm object (row) in `gbm_data` it creates a 
# dataframe describing relative influence per covariate for the given gbm object
# the output is a list with elements as rel. inf. dataframes 
compute_var_importance <- function(gbm_list) {
  
    tryCatch({
    
    # attempt to load the GBM object
    load(gbm_list[["file_path"]])
    
    # validate gbm object
    if (!exists("b.list")) {
      warning(paste("Invalid b.list model in file:", gbm_list["file_path"]))
      return(NULL)
    } else {
      message("Successfully loaded: ", gbm_list["file_path"])
    }
    
    # define a function to compute variable importance from a list of gbm objects
    gbm_summary <- function(i) {
      
      model <- b.list[[i]]
      
      if (!inherits(model, "gbm") || is.null(model$n.trees) || model$n.trees <= 0) {
        return(NULL)
      }
      
      gbm::summary.gbm(model, plotit = FALSE) |>
        as_tibble() |>
        dplyr::bind_cols(x = _, c(gbm_list["spp"], gbm_list["bcr"], boot = i))
    }
    
    # estimate covariate importance from the current `b.list`
    out_list <- lapply(seq_along(b.list), gbm_summary)
    
    
    # return all importance data.frames from the current `b.list`
    out_list <- out_list[!sapply(out_list, is.null)]
    
    message("Successfully extracted `rel.inf` from: ", gbm_list["file_path"])
    return(out_list)
    
    # remove loaded gbm object to free up memory
    rm(b.list)
    gc()
    
    # give instructions to tryCatch if it runs into an error
  }, error = function(e) {
    warning(paste("Error processing", gbm_list["file_path"], ":", e$message))
    return(NULL) 
    
  }) # close tryCatch
  
} # close compute_var_importance()



#7. export the necessary variables and functions to the cluster----
print("* exporting objects and functions to cluster *")
clusterExport(cl, c("root", "gbm_list", "compute_var_importance"))


#8. run the covariate importance function in parallel----
print("* running `compute_var_importance` in parallel` *")
bam_covariate_importance_v5 <- parLapply(cl, gbm_list, compute_var_importance)
names(bam_covariate_importance_v5) <- paste(gbm_data$spp, gbm_data$bcr, sep="_")
saveRDS(bam_covariate_importance_v5, file=file.path("/home/mannfred/scratch/", "bam_covariate_importance_list_v5.rds"))


#9. stop the cluster----
print("* stopping cluster :-)*")
stopCluster(cl)



#10. post-cluster data cleanup (executed on local machine)----

# import extraction lookup table to obtain covariate classes----
# (for appending to covariate importance data)
# lookup table is missing "Year" and "Method", so manually adding here
nice_var_names <-
  readr::read_csv(file.path("/home/mannfred/scratch/nice_var_names_v5.csv")) |>
  dplyr::select(var_class, var) |> 
  dplyr::distinct(var, .keep_all=TRUE)

# combine all bootstraps per spp x bcr tuple
# imap applies a function to each element of a vector, and its index
merged_by_sppbcr <- 
  purrr::imap(bam_covariate_importance_v5, ~ bind_rows(.x)) |>  # merge 2nd level lists of tibbles for each spp x bcr tuple
  list_rbind()  # merge top level lists of tibbles into one master data frame

# get mean and SD covariate importance per species x bcr x var tuple
bam_covariate_importance_means_v5 <- 
  merged_by_sppbcr |> 
  dplyr::group_by(spp, bcr, var) |> 
  dplyr::summarise(
    mean_rel_inf = mean(rel.inf, na.rm = TRUE),
    sd_rel_inf   = sd(rel.inf, na.rm = TRUE),
    n_boots = sum(rel.inf > 0),
    .groups = "drop") |> 
  dplyr::left_join(nice_var_names, by = "var") |> # append more interpretable variable names
  dplyr::arrange(spp, bcr, desc(mean_rel_inf))



save(bam_covariate_importance_means_v5, file=file.path("/home/mannfred/scratch/bam_covariate_importance_v5.rda"))


if(cc){ q() }


# `bam_covariate_importance_v5` is a `data.frame` with rows as the mean relative influence (of 10 bootstraps) of a model covariate for a given species x BCR. 
# There are 7 columns: `spp` gives the Four-Letter Bird Code indicating the species, `bcr` is the Bird Conservation Region, and `var` is the covariate.
# `mean_rel_inf` and `sd_rel_inf` are the mean and standard deviation of the covariate relative influence across 10 bootstraps.
# Relative influence was reported only when the mean was greater than or equal to 1. 
# `n_boots` indicates how many bootstraps (for a given covariate) had a relative influence of > 1. 
# Even when relative influence was < 1 (for a single bootstrap), it still contributed to the mean. 
# `var_class` denotes a broad variable class to which the covariate belongs. 
# The code used to generate this dataset can be found in "PrepareDataCovariateImportanceV5.R"
