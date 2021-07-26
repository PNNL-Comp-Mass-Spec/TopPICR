#' ...
#'
#' ...
#'
#' @param x_cluster ...
#'
#' @param x ...
#'
#' @return ...
#'
#' @importFrom magrittr %>%
#'
#' @export
#'
choose_cv <- function (x_cluster, x) {

  # Join x_cluster and x. This is necessary to link the cluster with the variables
  # dropped before aligning the retention time, recalibrating the mass, and
  # performing hierarchical clustering. Then combine Gene and cluster_new to
  # be used as a new proteoform.
  x_gcc <- dplyr::inner_join(x, x_cluster) %>%
    dplyr:: mutate(gcc = paste(Gene, cluster_new, sep = "_"))

  # Remove rows with multiple Feature intensity values by only keeping rows with
  # the maximum feature intensity value.
  x_gcc <- x_gcc %>%
    dplyr::group_by(ProjID, CV, Fraction, gcc) %>%
    dplyr::summarise(`Feature intensity` = max(`Feature intensity`)) %>%
    dplyr::ungroup()

  # Fish out the unique proteoforms.
  unique_gcc <- unique(x_gcc$gcc)

  # Start a list to hold output for top_cv.
  gcc_list <- vector(mode = "list",
                     length = length(unique_gcc))

  # Loop through each proteoform.
  for (e in 1:length(unique_gcc)) {

    gcc_list[[e]] <- top_cv(x = x_gcc,
                            cur_gcc = unique_gcc[[e]])

  }

  # Combine all elements of the list into a data frame.
  x_msn_prep <- data.table::rbindlist(gcc_list)

  return (x_msn_prep)

}

# x_55 auxiliary functions -----------------------------------------------------

top_cv <- function (x, cur_gcc) {

  # Fabricate a character vector containing the CV values.
  volts <- c("-30", "-40", "-50")

  # Filter the entire data frame by gcc.
  x_frac <- x %>%
    dplyr::filter(gcc == cur_gcc)

  # Check the dimension of x_frac.
  if (dim(x_frac)[1] == 1) {

    return (x_frac %>%
              dplyr::select(ProjID, `Feature intensity`, gcc))

  }

  # Count the number of unique project IDs that have a feature intensity value
  # for each of the three CV values.
  counts <- c(x_frac %>%
                dplyr::filter(CV == "-30") %>%
                dplyr::distinct(ProjID) %>%
                nrow(),

              x_frac %>%
                dplyr::filter(CV == "-40") %>%
                dplyr::distinct(ProjID) %>%
                nrow(),

              x_frac %>%
                dplyr::filter(CV == "-50") %>%
                dplyr::distinct(ProjID) %>%
                nrow())

  # Determine which CV value has the most project IDs with a feature intensity.
  max_counts <- max(counts)

  # Determine where the max value occurs (which element contains the most
  # project IDs with a feature intensity).
  which_max <- which(counts == max_counts)

  # Check the number of matches for the highest number of project IDs per CV. If
  # this value is larger than one there is a tie between CVs for the highest
  # number of unique project IDs.
  if (length(which_max) > 1) {

    # Find the CV with the highest median feature intensity.
    which_cv <- max_fi(x = x_frac)

    # Subset the data with the CV corresponding to the highest median feature
    # feature intensity.
    x_max_fi <- x_frac %>%
      dplyr::filter(CV == which_cv)

    # Sum Across any ProteoForms that have multiple feature intensities per
    # project ID and CV.
    x_top <- top_fi(x = x_max_fi)

    return (x_top)

    # There is not a tie between CVs for the highest number of projects IDs.
  } else {

    # Subset the data by the CV with the highest number of unique project IDs.
    x_max_pids <- x_frac %>%
      dplyr::filter(CV %in% volts[which_max])

    # Sum Across any ProteoForms that have multiple feature intensities per
    # project ID and CV.
    x_top <- top_fi(x = x_max_pids)

    # Return the filtered data frame according to the CV value with the highest
    # project ID counts.
    return (x_top)

  }

}

# This function calculates the maximum of the feature intensity across multiple
# CV values for unique combinations of CV, project ID, and feature intensity.
top_fi <- function (x) {

  # Sum the feature intensities across CV values.
  x_top <- x %>%
    # Group by project ID because we will take the maximum feature intensity
    # value as a representative feature intensity for each project ID.
    dplyr::group_by(ProjID) %>%
    dplyr::summarize(`Feature intensity` = max(`Feature intensity`)) %>%
    # Add the proteoform column back in the data frame because the summarize
    # function removes all columns (except ProjID and Feature intensity).
    dplyr::mutate(gcc = unique(x$gcc))

  return (x_top)

}

# This function calculates the median feature intensity for each CV value. The
# CV corresponding to the highest median feature intensity will be used to
# subset the data matrix.
max_fi <- function (x) {

  # Compute the median feature intensity by CV.
  fi_by_cv <- x %>%
    dplyr::group_by(CV) %>%
    dplyr::summarize(median_fi = stats::median(`Feature intensity`))

  # Determine which CV the highest median feature intensity belongs to.
  which_max <- which.max(fi_by_cv$median_fi)

  # Extract the CV that corresponds to the highest median feature intensity.
  return (fi_by_cv$CV[[which_max]])

}
