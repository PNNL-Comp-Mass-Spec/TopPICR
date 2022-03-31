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
#'   from TopPIC. The output from these files will be combined into one
#'   \code{data.table}. Currently this function can read in data from either
#'   ms2_toppic_prsm.tsv files or ms1.feature files. When reading in data all
#'   files must be the same type (e.g., ms2_toppic_prsm.tsv).
#'
#' @param faims Logical. If true the data contain multiple FAIMS compensation
#'   voltages (CVs).
#'
#' @param ... Additional arguments to \code{\link[readr]{read_tsv}}.
#'
#' @return A data.table where all the data from the input files are combined
#'   row wise. When PrSM data is the input, the following variables are
#'   added/removed:
#'
#'   | Added             | Removed                    |
#'   | ----------------- | -------------------------- |
#'   | `Dataset`         | `Data file name`           |
#'   | `mz`              | `Prsm ID`                  |
#'   | `Gene`            | `Spectrum ID`              |
#'   | `isDecoy`         | `Fragmentation`            |
#'   | `CV`              | `#peaks`                   |
#'   |                   | `Proteoform ID`            |
#'   |                   | `Feature score`            |
#'   |                   | `First residue`            |
#'   |                   | `Last residue`             |
#'   |                   | `#variable PTMs`           |
#'   |                   | `#matched peaks`           |
#'   |                   | `#matched fragment ions`   |
#'   |                   | `Spectrum-level Q-value`   |
#'   |                   | `Proteoform-level Q-value` |
#'
#'   When unidentified (MS2) data is the input, the the following variables are
#'   added/removed:
#'
#'   | Added             | Removed                    |
#'   | ----------------- | -------------------------- |
#'   | `Dataset`         | `Sample_ID`                |
#'   | `CV`              | `Minimum_fraction_id`      |
#'   |                   | `Maximum_fraction_id`      |
#'
#' @md
#'
#' @author Evan A Martin
#'
#' @export
#'
read_toppic <- function (file_path, file_name, faims, ...) {

  # Create a list that will be used to hold the data from each file.
  the_list <- vector(mode = "list",
                     length = length(file_name))

  # Determine whether the data are MS1 or MS2 files and read the data in
  # accordingly.
  if (stringr::str_detect(file_name[[1]], "ms1.feature")) {

    # Loop through each MS1 file and read it into R.
    for (e in 1:length(file_name)) {

      # The ms1.feature files can be read directly into R with the read_table
      # function from the readr package. Add the Dataset variable which is the
      # name of the file that was just read in. It will also be modified
      # ("_ms1.feature" is removed) to match the format of Dataset from the PrSM
      # data.
      the_list[[e]] <- readr::read_table(
        file.path(file_path, file_name[[e]])
      ) %>%
        # Create the data set name from the file name.
        dplyr::mutate(
          Dataset = file_name[[e]],
          Dataset = stringr::str_remove(Dataset, "_ms1.feature")
        )

      # Rename all MS1 TopPIC variables that are used by TopPICR functions. We
      # take this step so only one instance of the TopPIC name has to be changed
      # if they change the name of any of these variables in a later version.
      the_list[[e]] <- the_list[[e]] %>%
        dplyr::mutate(
          Intensity = Intensity
        )

    }

  } else {

    # Loop through each PrSM (MS2) file, find and ignore all metadata rows, and
    # read the data rows into R.
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
      the_list[[e]] <- readr::read_tsv(
        file = file.path(file_path, file_name[[e]]),
        skip = n_prelim,
        ...
      )

      # Rename all MS2 TopPIC variables that are used by TopPICR functions. We
      # take this step so only one instance of the TopPIC name has to be changed
      # if they change the name of any of these variables in a later version.
      the_list[[e]] <- the_list[[e]] %>%
        dplyr::mutate(
          `Feature apex` = `Feature apex time`,
          `Feature intensity` = `Feature intensity`,
          `E-value` = `E-value`,
          Proteoform = Proteoform,
          `Scan(s)` = `Scan(s)`,
          `#unexpected modifications` = `#unexpected modifications`,
          `Precursor mass` = `Precursor mass`,
          `Adjusted precursor mass` = `Adjusted precursor mass`,
          MIScore = MIScore
        )

      # Create variables used by TopPICR from existing TopPIC variables.
      the_list[[e]] <- the_list[[e]] %>%
        dplyr::mutate(
          # If the TopPIC authorities change the names of the `Data file name`,
          # Charge, `Protein description`, or `Protein accession` variables they
          # will need to be adjusted accordingly in the following lines.
          Dataset = stringr::str_remove(`Data file name`, "_ms2.msalign"),
          Dataset = sapply(stringr::str_split(Dataset, "/"), tail, 1),
          mz = (`Precursor mass` + Charge * 1.007276466621) / Charge,
          Gene = sub(".*GN=(\\S+).*", "\\1", `Protein description`),
          isDecoy = grepl("^DECOY", `Protein accession`)
        ) %>%
        dplyr::filter(grepl("GN=", `Protein description`))

    }

  }

  # Combine data from all files.
  x <- rbindlist(the_list)

  # Include CV variable with data extracted from the Dataset variable if working
  # with FAIMS data.
  if (faims) {

    x <- x %>%
      dplyr::mutate(
        CV = sapply(stringr::str_split(Dataset, "_"), tail, 1),
        # Remove the CV information from the Dataset variable. If we don't do
        # this we will always be grouping by CV anytime we group by Dataset.
        Dataset = stringr::str_replace(Dataset, "_[0-9]*$", "")
      )

  } else {

    # Add a character string for the CV variable indicating there are no CVs.
    # The CV variable needs to exist even when there aren't any CVs because many
    # of the remaining functions will group by CV.
    x <- x %>%
      dplyr::mutate(CV = "noFAIMS")

  }

  # Create a character vector of variables that are not needed at any point in
  # the top down workflow. They will be removed with dplyr::any_of. The function
  # will not throw an error if any of the variables listed are not present in
  # the data set.
  useless <- c(
    # Identified feature variables:
    "Data file name", "Prsm ID", "Spectrum ID", "Fragmentation",
    "#peaks", "Proteoform ID", "Feature score",
    "#variable PTMs", "#matched peaks", "#matched fragment ions",
    "Q-value (spectral FDR)", "Proteoform FDR", "First residue",
    "Last residue", "Spectrum-level Q-value",
    "Proteoform-level Q-value", "First residue", "Last residue",
    # Unidentified feature variables:
    "Sample_ID", "Minimum_fraction_id", "Maximum_fraction_id"
  )

  return (dplyr::select(x, -dplyr::any_of(useless)))

}
