#' Load Chemical Fingerprint Data
#'
#' This function loads a file as a matrix. All explanatory data should be binary,
#' indicating presence or absence of a bit structure, the first column
#' should have the chemical name in order to join on additional physical data late.
#'
#' @param path Path to the input file
#' @param file_type type of file (xlsx, csv, or txt)
#' @return A matrix
#' @export
#'
load_chem_data <- function(path, file_type, filter = TRUE, filter_type = "Pubchem") {
  if (file_type == "xlsx") {
    data <- readxl::read_excel(here(path))
  }

  if (file_type == "csv") {
    data <- read.csv(here(path))
  }

  if (file_type == "txt") {
    data <- read.delim(here(path))
  }

  if (filter == TRUE) {
    data <- dplyr::select(
      data,
      "Name", contains("Label"), contains(filter_type)
    )
  }
  return(data)
}
