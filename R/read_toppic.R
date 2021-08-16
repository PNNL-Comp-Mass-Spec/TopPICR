#' Read in TopPIC files
#'
#' ...
#'
#' @param file_path ...
#'
#' @param file_name ...
#'
#' @param ... Additional arguments to the \code{read.delim} function.
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
