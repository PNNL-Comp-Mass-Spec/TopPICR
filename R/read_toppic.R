#' Read in TopPIC files
#'
#' A wrapper function to the \code{read.delim} function that is specific to the
#' files output by TopPIC. Metadata is always included at the beginning of the
#' TopPIC file and a header is always provided. This function removes all lines
#' containing metadata and reads in the remaining lines of the file starting
#' with the header.
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
#'   row wise.
#'
#' @export
#'
read_toppic <- function (file_path, file_name, ...) {

  # Create a list that will be used to hold the data from each file.
  the_list <- vector(mode = "list",
                     length = length(file_name))

  # Loop through each file name and read the file into R.
  for (e in length(file_name)) {

    # Find the number of lines preceeding the header ---------------

    the_lines <- system(paste("grep -n '^\\*' ",
                              file_path,
                              file_name[[e]],
                              sep = ""),
                        intern = TRUE)

    n_lines <- length(the_lines)

    n_prelim <- as.numeric(gsub("[^{0-9}]*", "", the_lines[[n_lines]]))

    # Read in the data skipping the lines preceeding the header ---------------

    the_list[[e]] <- read.delim(file = paste0(file_path, file_name[[e]]),
                                skip = n_prelim,
                                ...)

  }

  return (rbindlist(the_list))

}
