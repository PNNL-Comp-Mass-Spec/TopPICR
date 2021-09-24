#' Read in TopPIC files
#'
#' A wrapper function to the \code{read.delim} function that is specific to the
#' files output by TopPIC. Metadata is always included at the beginning of the
#' TopPIC file and a header is always provided. This function removes all lines
#' containing metadata and reads in the remaining lines of the file starting
#' with the header. Variables needed for the remainder of the TopPICR workflow
#' are added and variables not used by TopPICR are removed.
#'
#' @param file_path A character string specifying the path to the folder with
#'   the TopPIC output files.
#'
#' @param file_name A character vector of file names containing the files output
#'   from TopPIC. The output from these files will be combined into one data
#'   frame.
#'
#' @param ... Additional arguments to the \code{read.delim} function.
#'
#' @return A data.table where all the data from the input files are combined
#'   row wise. The following variables have been added/removed:
#'
#'   | Added             | Removed                    |
#'   | ----------------- | -------------------------- |
#'   | `mz`              | `Prsm ID`                  |
#'   | `Gene`            | `Spectrum ID`              |
#'   | `isDecoy`         | `Fragmentation`            |
#'   | `Dataset`         | `#peaks`                   |
#'   | `CV`              | `Proteoform ID`            |
#'   |                   | `Feature score`            |
#'   |                   | `First residue`            |
#'   |                   | `Last residue`             |
#'   |                   | `MIScore`                  |
#'   |                   | `#variable PTMs`           |
#'   |                   | `#matched peaks`           |
#'   |                   | `#matched fragment ions`   |
#'   |                   | `Spectrum-level Q-value`   |
#'   |                   | `Proteoform-level Q-value` |
#'
#' @export
#'
read_toppic <- function (file_path, file_name, ...) {

  # Read in the data from TopPIC files -----------------------------------------

  # Create a list that will be used to hold the data from each file.
  the_list <- vector(mode = "list",
                     length = length(file_name))

  # Loop through each file name and read the file into R.
  for (e in 1:length(file_name)) {

    # Find the number of lines preceding the header.
    the_lines <- system(paste("grep -n '^\\*' ",
                              file_path,
                              file_name[[e]],
                              sep = ""),
                        intern = TRUE)
    n_lines <- length(the_lines)
    n_prelim <- as.numeric(gsub("[^{0-9}]*", "", the_lines[[n_lines]]))

    # Read in the data skipping the lines preceding the header.
    the_list[[e]] <- read.delim(file = paste0(file_path, file_name[[e]]),
                                skip = n_prelim,
                                ...)

  }

  # Combine data from all files.
  x <- rbindlist(the_list)

  # Filter data and include additional variables -------------------------------

  # Create a character vector of variables that are not needed at any point in
  # the top down workflow. They will be removed with dplyr::any_of. The function
  # will not throw an error if any of the variables listed are not present in
  # the data set.
  useless <- c("Prsm ID", "Spectrum ID", "Fragmentation",
               "#peaks", "Proteoform ID", "Feature score", "MIScore",
               "#variable PTMs", "#matched peaks", "#matched fragment ions",
               "Q-value (spectral FDR)", "Proteoform FDR", "First residue",
               "Last residue", "Spectrum-level Q-value",
               "Proteoform-level Q-value", "First residue", "Last residue")

  x <- x %>%
    # Filter data based on protein accession and description variables.
    dplyr::filter(grepl("GN=", `Protein description`)) %>%
    # Add useful variables.
    dplyr::mutate(
      mz = (`Precursor mass` + Charge * 1.007276466621) / Charge,
      Gene = sub(".*GN=(\\S+).*","\\1",`Protein description`),
      isDecoy = grepl("^DECOY", `Protein accession`),
      Dataset = stringr::str_remove(`Data file name`, "_ms2.msalign"),
      # Extract the CV information (if it exists).
      CV = stringr::str_extract(Dataset, "[mp][0-9]*CV"),
      # Remove the CV information from Dataset (if it exists).
      Dataset = stringr::str_remove(Dataset, "[mp][0-9]*CV")
    ) %>%
  # Remove not so useful variables.
  dplyr::select(-dplyr::any_of(useless))

  return (x)

}
