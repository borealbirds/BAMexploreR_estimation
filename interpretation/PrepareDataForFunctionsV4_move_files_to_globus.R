library(tidyverse)


root <- "G:/Shared drives/BAM_NationalModels4/NationalModels4.0/Feb2020/out/boot/"

# ^[^0]*0[^0]*$ matches a string that:
# starts (^) with zero or more characters that are not 0 ([^0]*),
# followed by exactly one 0 (0),
# followed by zero or more characters that are not 0 ([^0]*),
# and ends ($) with this sequence.
# (this removes US BCRs)
all_files <- list.files(root, pattern="\\.RData$", recursive=FALSE, full.names = TRUE)
boot_files <- all_files[grep("^[^0]*0[^0]*$", basename(all_files))]


# Copy each file from mean_files to the new directory
file.copy(boot_files, "G:/My Drive/v4_bootstraps/")



# NEXT DAY: check which files didn't make it overnight, and continue transfer
donefiles <- list.files("G:/My Drive/v4_bootstraps/", pattern=".RData")
missed_files <- boot_files[!(basename(boot_files) %in% donefiles)]

file.copy(missed_files, "G:/My Drive/v4_bootstraps/")



