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
#'   row wise.
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
          Dataset = stringr::str_remove(Dataset, "_ms1.feature"),
          # Rename all MS1 TopPIC variables that are used by TopPICR functions.
          # We take this step so only one instance of the new TopPIC name has to
          # be updated if/when a variable name is changed in a later version.
          Mass = Mass,
          Intensity = Intensity,
          Time_apex = Time_apex
        ) %>%
        # Only keep the MS1 variables we use throughout TopPICR.
        dplyr::select(Dataset, Mass, Intensity, Time_apex)

    }

    # The following code reads in the MS2 data.
  } else {

    # Loop through each PrSM (MS2) file, find and ignore all metadata rows, and
    # read the data rows into R.
    for (e in 1:length(file_name)) {

      # Find the number of lines preceding the header.
      the_lines <- system2(command = "grep",
                           args = paste(" -n '^\\*' ",
                                        file_path,
                                        file_name[[e]],
                                        sep = ""),
                           stdout = TRUE)
      n_lines <- length(the_lines)
      n_prelim <- as.numeric(gsub("[^{0-9}]*", "", the_lines[[n_lines]]))

      # Read in the data skipping the lines preceding the header.
      the_list[[e]] <- readr::read_tsv(
        file = file.path(file_path, file_name[[e]]),
        skip = n_prelim,
        ...
      ) %>%
        # Rename all MS2 TopPIC variables that are used by TopPICR functions.
        # We take this step so only one instance of the new TopPIC name has to
        # be updated if/when a variable name is changed in a later version.
        dplyr::mutate(
          `Scan(s)` = `Scan(s)`,
          `Retention time` = `Retention time`,
          `Feature apex` = `Feature apex time`,
          Charge = Charge,
          `Precursor mass` = `Precursor mass`,
          `Adjusted precursor mass` = `Adjusted precursor mass`,
          `Feature intensity` = `Feature intensity`,
          `E-value` = `E-value`,
          `#unexpected modifications` = `#unexpected modifications`,
          AccMap = `#Protein hits`,
          MIScore = MIScore,
          `Special amino acids` = `Special amino acids`,
          Proteoform = Proteoform,
          UniProtAcc = `Protein accession`
        ) %>%
        # Only keep rows that have a gene in the protein description. Rows
        # without a gene aren't kept because ...
        dplyr::filter(grepl("GN=", `Protein description`)) %>%
        # Create variables used by TopPICR from existing TopPIC variables.
        #
        # If the TopPIC authorities change the names of the `Data file name`,
        # Charge, `Protein description`, or `Protein accession` variables they
        # will need to be adjusted accordingly in the following lines.
        dplyr::mutate(
          Dataset = stringr::str_remove(`Data file name`, "_ms2.msalign"),
          Dataset = purrr::map_chr(
            stringr::str_split(Dataset, "/"), dplyr::last
          ),
          mz = (`Precursor mass` + Charge * 1.007276466621) / Charge,
          Gene = sub(".*GN=(\\S+).*", "\\1", `Protein description`),
          isDecoy = grepl("^DECOY", `Protein accession`)
        ) %>%
        # Extract just the UniProt accession. This is the second element (when
        # splitting by |) of the `Protein accession` variable output by TopPIC.
        dplyr::mutate(
          UniProtAcc = purrr::map_chr(
            stringr::str_split(UniProtAcc, "\\|"), purrr::pluck, 2
          )
        ) %>%
        # Only keep variables we need throughout TopPICR.
        dplyr::select(
          Dataset,
          `Scan(s)`,
          `Retention time`,
          `Feature apex`,
          Charge,
          mz,
          `Precursor mass`,
          `Adjusted precursor mass`,
          `Feature intensity`,
          `E-value`,
          `#unexpected modifications`,
          AccMap,
          MIScore,
          `Special amino acids`,
          Proteoform,
          UniProtAcc,
          Gene,
          isDecoy
        ) %>%
        # Fill in missing values for all accessions that match a given sequence.
        # This step is necessary because TopPIC only includes information for
        # the `Protein accession` and `Protein description` variables when a
        # proteoform matches multiple sequences. All other variables are left
        # blank.
        tidyr::fill(`Scan(s)`:Proteoform, .direction = "down")

    }

  }

  # Combine data from all files.
  x <- rbindlist(the_list)

  # Include CV variable with data extracted from the Dataset variable if working
  # with FAIMS data.
  if (faims) {

    x <- x %>%
      dplyr::mutate(
        CV = purrr::map_chr(stringr::str_split(Dataset, "_"), dplyr::last),
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

  return (x)

}
