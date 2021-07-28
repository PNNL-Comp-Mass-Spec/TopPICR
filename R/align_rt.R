#' ...
#'
#' ...
#'
#' @param x ...
#'
#' @param ref_ds A character string specifying the reference data set.
#'
#' @param ... Arguments for the \code{\link[stats]{loess}} function.
#'
#' @return ...
#'
#' @importFrom magrittr %>%
#'
#' @export
#'
align_rt <- function (x, ref_ds, ...) {

  # Find the apex or peak of the retention time.
  x_peak <- x %>%
    dplyr::group_by(Dataset, Proteoform) %>%
    # Calculate a retention time that is representative of each proteoform for
    # each Dataset.
    dplyr::summarize(RTest = find_peak(`P-value`, RTmin)) %>%
    dplyr::ungroup()

  # Create a data frame for the reference data set. This will be used to match
  # Proteoforms between the reference data set and every other data set in the
  # input data frame.
  x_ref <- x_peak %>%
    dplyr::filter(Dataset == ref_ds)

  # Align the retention times to the reference data set.
  x_align <- x_peak %>%
    dplyr::group_by(Dataset) %>%
    # Align retention times by Dataset to a reference Dataset.
    dplyr::mutate(RTalign = alignment(rt = RTest,
                                      pf = Proteoform,
                                      x_ref = x_ref,
                                      # The dots are for the loess function.
                                      ...)) %>%
    dplyr::ungroup()


  return (x_align)

}

# x_42 auxiliary functions -----------------------------------------------------

# Function for computing the median retention time depending on the number of
# observations for a given data set.  This function works within the mutate
# function.
find_peak <- function (pvalue, rtmin) {

  # Check for the number of p-values.
  if (length(pvalue) > 4) {

    # Order the rtmin values from smallest to largest. If there are more than 20
    # p-values only keep the first 20.
    series <- order(rtmin)[1:min(20, length(pvalue))]

    # Reorder the pvalue and rtmin vectors according to the order from above.
    pvalue <- pvalue[series]
    rtmin <- rtmin[series]

    # Find the p-value corresponding to the 25th percentile.
    quant <- stats::quantile(pvalue, 0.25)

    # Extract the indices of the p-values that are less than quant.
    idx <- which(pvalue <= quant)

    # Calculate the median RTmin for the smallest 25th percentile of p-values.
    rtmed <- stats::median(rtmin[idx])

    # Runs if the number of p-values is between 2 and 4.
  } else if (length(pvalue) > 1) {

    # Extract the indices of the two smallest p-values.
    idx <- order(pvalue)[1:2]

    # Calculate the median RTmin for the two smallest p-values.
    rtmed <- stats::median(rtmin[idx])

    # Runs if the number of p-values is equal to 1.
  } else {

    # Return the only RTmin value present in the data frame.
    rtmed <- rtmin[[1]]

  }

  # Return the median retention time.
  return (rtmed)

}

# Function for aligning retention times. This function works within the mutate
# function.
# rt: A vector of retention times extracted from the data matrix by the group_by
# function.
# pf: A vector of proteoforms extracted from the data matrix by the group_by
# function.
# x: the original data matrix (without any groupings).
# ref_ds: A character string. This is the name of the reference data set.
alignment <- function (rt, pf, x_ref, ...) {

  # Find the Proteoforms that occur in both the reference and current data sets.
  # The current data set comes from the group_by function (it will loop through
  # every unique data set in the data matrix).
  both <- dplyr::intersect(x_ref$Proteoform, pf)

  # Extract the indices where the intersecting Proteoforms occur in both the
  # reference data set and the current data set.
  idx_ref <- which(x_ref$Proteoform %in% both)
  idx_cur <- which(pf %in% both)

  # Run loess on the retention times from the reference and current data sets.
  model_loess <- stats::loess(x_ref$RTest[idx_ref] ~ rt[idx_cur], ...)

  # Return the predicted retention times.
  return (stats::predict(object = model_loess, newdata = rt))

}
