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

evaluate_grid_from_model <- function(gbm_obj){
  load(gbm_obj)  # loads `b.list`
  
  spp <- str_extract(gbm_obj, "[A-Z]{4}")
  bcr <- str_extract(gbm_obj, "(can|usa)[0-9]{1,2}")
  
  map_dfr(seq_along(b.list), function(boot_i) {
    
    model <- b.list[[boot_i]]
    vars <- model$var.names
    
    map_dfr(vars, function(var_i) {
      df <- plot.gbm(model, i.var = var_i, type = "response", return.grid = TRUE, continuous.resolution = 25)
      
      # convert all covariates to characters,
      # `bam_partial_dependence` will convert continuous back to numeric when needed
      # this allows all data types to be stored in one tibble
      covariate_values_chr <- as.character(df[[1]])
      covariate_type <- if (is.factor(df[[1]])) "categorical" else "continuous"
      
      tibble(
        species = spp,
        bcr = bcr,
        var = var_i,
        replicate = boot_i,
        covariate_value = covariate_values_chr,
        covariate_type = covariate_type,
        y = df$y
      )
    })
  })
}


#6. export the necessary variables and functions to the cluster----
print("* exporting gbm_objs and functions to cluster *")

clusterEvalQ(cl, {
  library(gbm)
  library(dplyr)
  library(stringr)
})

clusterExport(cl, c("gbm_objs", "evaluate_grid_from_model"))


#7. run the function in parallel----
print("* running `process_gbm` in parallel *")
boot_pts_sorted <- parLapply(cl = cl, X = gbm_objs
, fun = evaluate_grid_from_model)
boot_pts_df <- bind_rows(boot_pts_sorted)

#8. stop the cluster----
print("* stopping cluster :-)*")
stopCluster(cl)

#9. save data----
saveRDS(boot_pts_df, file="/home/mannfred/scratch/boot_pts_df_v5.rds")

