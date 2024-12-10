# ---
# title: National Models 4.0 - data organization for analysing relative influence of covariates and species traits
# author: Mannfred Boehm
# created: November 19, 2024
# ---

#NOTES################################

# This script extracts covariate contributions to model predictions (VERSION 4), and also synthesizes various trait databases.

# The outputs will be analysed for identifying covariates of importance across several metrics (e.g. BCR, ecology, etc.).

# This script generates the exported data for R package functions that summarise covariate performance. 


# attach packages----
library(gbm)
library(parallel)
library(tidyverse)

# define local or cluster----
test <- TRUE
cc <- FALSE

# set number of tasks for local vs cluster----
if(cc){ n_tasks <- 32}
if(!cc | test){ n_tasks <- 4}


# create and register clusters----
# creates 32 copies of R running in parallel via 32 tasks, on one of the cluster's sockets (processors). 
# Belgua has ~965 nodes
print("* creating clusters *")
cl <- makePSOCKcluster(n_tasks, type="PSOCK")

# print number of tasks and host name for confirmation
cl


# set root path----
print("* setting root file path *")

if(!cc){root <- "G:/My Drive/"}
if(cc){root <- "/home/mannfred/scratch/v4_bootstraps"}

tmpcl <- clusterExport(cl, c("root"))


# attach packages on clusters----
# `clusterEvalQ` evaluates a literal expression on each cluster node. 
print("* Loading packages on workers *")
tmpcl <- clusterEvalQ(cl, library(gbm))
tmpcl <- clusterEvalQ(cl, library(tidyverse))



# index gbm objects----
# access gbm objects 
gbm_objs <- 
  lapply(file.path(root, "v4_bootstraps"), list.files, full.names = TRUE) |> 
  unlist() |> 
  (\(x) x[1:3])()

print(paste("* Found", length(gbm_objs), " files *"))


# import extraction lookup table to obtain covariate classes----
# lookup table is missing "Year" and "Method", so manually adding here
nice_var_names <- 
  read_csv(file.path(root, "nice_var_names.csv")) |>
  dplyr::select(var_class, var) 

if (exists("nice_var_names")) {
  print("* Found nice_var_names.csv *")
  } else {print("* nice_var_names.csv not found *")}


# create an index from `gbm_objs` containing the species (FLBC), BCR, and bootstrap replicate----
sample_id <- 
  gbm_objs |> 
  sub("^.*gnmboot-", "", x = _) |> 
  sub("\\.RData", "", x = _) |>
  stringr::str_split_fixed(pattern="-", n=3) |> 
  tibble::as_tibble() |> 
  magrittr::set_colnames(c("spp", "bcr", "boot"))

if (exists("sample_id")) {
  print("* constructed sample_id *")
} else {message("* sample_id not constructed *")}



# define function that creates `bam_covariate_importance` data frame----

# create a list of dataframes containing relative influence per covariate
# every element of this list comes from a single bootstrap
compute_var_importance <- function(gbm_file, sample_row) {
  
  # load the GBM object
  load(gbm_file)
  
  # compute variable importance
  importance_df <- 
    out |> 
    gbm::summary.gbm(plotit = FALSE) |>
    as_tibble() |> 
    dplyr::cross_join(x = _, sample_row) |> 
    dplyr::mutate(file_name = gbm_file)
  
  return(importance_df)
}

# run the function in parallel----
covs <- parLapply(cl = cl, X = gbm_objs, fun = compute_var_importance, sample_id)


# merge list into dataframe
covariate_importance <- suppressMessages(purrr::reduce(covs, full_join))


# check for missing bootstraps and fill in zeros where covariates are missing from a bootstrap
covariate_importance_zeroed <- 
  covariate_importance |>  
  dplyr::mutate(boot = as.integer(boot)) |> # change `chr` to `int`
  tidyr::complete(spp, bcr, var, boot = 1:32, fill = list(rel.inf = 0)) |>  # create rows for every spp x bcr x var x bootstrap combo; fill missing rows with 0
  dplyr::group_by(spp, bcr, var) |>  
  dplyr::summarise(mean_rel_inf = mean(rel.inf, na.rm = TRUE),  # calculate mean across all bootstraps (including 0s)
                   sd_rel_inf = sd(rel.inf, na.rm = TRUE),      # calculate standard deviation across all bootstraps (including 0s)
                   n_boots = sum(rel.inf > 0)) |>    # count the number of bootstraps where the variable had non-zero rel.inf
  dplyr::filter(mean_rel_inf > 0) |> # filter out spp x bcr x var tuples that don't have rel.inf data
  dplyr::left_join(nice_var_names, by = "var") |> 
  dplyr::arrange(spp, bcr, desc(mean_rel_inf)) 

saveRDS(covariate_importance_zeroed, file=file.path(root, "covariate_importance_zeroed_v4.rds"))


# stop the cluster----
print("* stopping cluster :-)*")
stopCluster(cl)

if(cc){ q() }