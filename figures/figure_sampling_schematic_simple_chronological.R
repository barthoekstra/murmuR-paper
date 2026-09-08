# Chronological companion to figure_sampling_schematic_simple.R.
#
# The base script retains its original selection and output names by default.
# This entry point enables its consecutive-night search, orders the selected
# panels by date, and writes files with the `_chronological` suffix.

Sys.setenv(SAMPLING_SCHEMATIC_CHRONOLOGICAL = "true")
if (!nzchar(Sys.getenv("SAMPLING_SCHEMATIC_DATA_PATH"))) {
  Sys.setenv(
    SAMPLING_SCHEMATIC_DATA_PATH = "/data/birdcloudstorage-tvm/ibm-ml/data/"
  )
}
if (!nzchar(Sys.getenv("SAMPLING_SCHEMATIC_YEAR"))) {
  Sys.setenv(SAMPLING_SCHEMATIC_YEAR = "2022")
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1])))
} else {
  "figures"
}

source(file.path(script_dir, "figure_sampling_schematic_simple.R"))
