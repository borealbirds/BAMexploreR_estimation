# import extraction lookup table to obtain covariate classes----
# (for appending to covariate importance data)
# lookup table is missing "Year" and "Method", so manually adding here
nice_var_names <-
  readr::read_csv(file.path(root, "v4_bootstraps", "nice_var_names.csv")) |>
  dplyr::select(var_class, var)