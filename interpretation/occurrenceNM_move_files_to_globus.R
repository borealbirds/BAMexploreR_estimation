
root <- "G:/Shared drives/BAM_NationalModels4/NationalModels4.0/Feb2020/artifacts"
flbc <- 
  read_csv(file="G:/Shared drives/BAM_NationalModels5/data/priority_spp_with_model_performance.csv") |> 
  dplyr::select(species_code)

# List all directories in the root folder
all_species_dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)

# Filter directories that end with the bird codes
matching_dirs <- all_species_dirs[basename(all_species_dirs) %in% flbc$species_code]

# Now find all files containing 'mean' in each matching directory
mean_files <- unlist(sapply(matching_dirs, function(folder) {
  list.files(folder, pattern = "-Mean\\.tif$", full.names = TRUE)
}))

# Copy each file from mean_files to the new directory
file.copy(mean_files, "C:/Users/mannf/Downloads/mean_rasters_v4", overwrite = TRUE)
