# align_rt preliminary functions -----------------------------------------------

#' Find reference data set
#'
#' \code{find_ref_ds} finds the data set that has the highest number of unique
#' proteoforms. This data set will be used as the reference data set when
#' aligning retention times and recalibrating the mass.
#'
#' @param x A \code{data.table} output from the \code{infer_pf} function.
#'
#' @return A character vector containing the name of the reference data set. The
#'   reference data set is the one that has the highest number of unique
#'   \code{Proteoforms}.
#'
#' @importFrom magrittr %>%
#'
#' @author Evan A Martin
#'
#' @export
#'
find_ref_ds <- function (x) {

  # Group the data by data set and count the number of unique proteoforms for
  # each data set. Select the data set that has the highest number of unique
  # proteoforms.
  ref_ds <- x %>%
    dplyr::group_by(Dataset) %>%
    dplyr::summarize(n_pfs = dplyr::n_distinct(Proteoform)) %>%
    dplyr::slice_max(n_pfs) %>%
    dplyr::pull(Dataset)

}

#' Use \code{loess} to create a model between data sets
#'
#' Run the \code{\link[stats]{loess}} function using the values from the
#' \code{Feature apex} variable between the reference data set and all other
#' data sets.
#'
#' @param x A \code{data.table} output from the \code{infer_pf} function.
#'
#' @param ref_ds A character string specifying the reference data set.
#'
#' @param ... Arguments for the \code{\link[stats]{loess}} function.
#'
#' @return A list containing the loess model between each data set and the
#'   reference data set.
#'
#' @importFrom magrittr %>%
#'
#' @author Evan A Martin
#'
#' @export
#'
form_model <- function (x, ref_ds, ...) {

  # Create a data frame for the reference data set. This will be used to match
  # Proteoforms between the reference data set and every other data set in the
  # input data frame. Only the distinct values of `Feature intensity`, `Feature
  # apex`, and Proteoform are kept. Then filter (keep) the rows with the maximum
  # feature intensity for each proteoform. This needs to be done to ensure there
  # is only one retention time (`Feature apex`) for each proteoform. If there
  # are multiple retention time values per proteoform the loess function will
  # not run.
  x_ref <- x %>%
    dplyr::filter(Dataset == ref_ds) %>%
    dplyr::distinct(`Feature intensity`, `Feature apex`, Proteoform) %>%
    dplyr::group_by(Proteoform) %>%
    dplyr::slice_max(`Feature intensity`) %>%
    # Remove rows where the max `Feature intensity` appears in more than one row
    # by selecting the min `Feature apex` value.
    dplyr::slice_min(`Feature apex`)

  # Create a vector of unique data set names. This will be used to create the
  # models between the reference and every other data set.
  ds <- unique(x$Dataset)

  x_model <- vector(mode = "list",
                    length = length(ds))

  # Loop through each data set and create a loess model between the current and
  # and reference data sets.
  for (e in 1:length(ds)) {

    # Filter the data to the rows corresponding to the current data set. See the
    # explanation above for explanation of why we use distinct, group_by, and
    # slice_max.
    x_cur <- x %>%
      dplyr::filter(Dataset == ds[[e]]) %>%
      dplyr::distinct(`Feature intensity`, `Feature apex`, Proteoform) %>%
      dplyr::group_by(Proteoform) %>%
      dplyr::slice_max(`Feature intensity`) %>%
      # Remove rows where the max `Feature intensity` appears in more than one
      # row by selecting the min `Feature apex` value.
      dplyr::slice_min(`Feature apex`)

    # Find Proteoforms that occur in both the reference and current data sets.
    both <- dplyr::intersect(x_ref$Proteoform,
                             x_cur$Proteoform)

    # Extract the indices where the intersecting Proteoforms occur in both the
    # reference data set and the current data set.
    idx_ref <- which(x_ref$Proteoform %in% both)
    idx_cur <- which(x_cur$Proteoform %in% both)

    # Run loess on the retention times from the reference and current data sets.
    model_loess <- stats::loess(
      x_ref$`Feature apex`[idx_ref] ~ x_cur$`Feature apex`[idx_cur],
      ...
    )

    x_model[[e]] <- model_loess

  }

  # Name the list with the Dataset names so the models can be applied correctly
  # when aligning the retention times.
  names(x_model) <- ds

  return (x_model)

}

# align_rt main function -------------------------------------------------------

#' Align retention times to a reference data set
#'
#' Aligns all retention times to a reference data set using the times from the
#' \code{Feature apex} column. The \code{loess} function from the \code{stats}
#' package is used for the alignment.
#'
#' @param x A \code{data.table} output from the \code{infer_pf} function.
#'
#' @param model A list containing the model output by the \code{loess} function.
#'   The elements of this list are the models between the reference data set and
#'   every other data set.
#'
#' @param var_name A character string. The name of the retention time variable
#'   that will be aligned with the reference data set (e.g., `Feature apex`).
#'
#' @return A \code{data.table} with the retention times aligned to a reference
#'   data set. The following variables have been added/removed:
#'
#'   | Added             | Removed                    |
#'   | ----------------- | -------------------------- |
#'   | `RTalign`         |                            |
#'
#' @md
#'
#' @importFrom magrittr %>%
#'
#' @author Evan A Martin
#'
#' @export
#'
align_rt <- function (x, model, var_name) {

  # Check if the data sets in x are also in the model list.
  if (!setequal(unique(x$Dataset), names(model))) {

    stop (paste("The data sets present in x and model do not match."))

  }

  # Align the retention times to the reference data set.
  x_align <- x %>%
    dplyr::group_by(Dataset) %>%
    # Align retention times by Dataset to a reference Dataset.
    dplyr::mutate(RTalign = alignment(rt = !!rlang::sym(var_name),
                                      ds = Dataset,
                                      mdl = model),
                  RTalign = round(RTalign, digits = 5)) %>%
    dplyr::ungroup()

  return (x_align)

}

# align_rt auxiliary functions -------------------------------------------------

# Function for aligning retention times. This function works within the mutate
# function.
# rt: A vector of retention times extracted from the data matrix by the group_by
# function.
# ds: A vector of data set names. This vector should only contain one data set
# name (repeated for the number of rows for the current data set).
# mdl: A list of loess models corresponding to each data set.
alignment <- function (rt, ds, mdl) {

  # Nab the model that corresponds to the current data set.
  idx <- which(names(mdl) == unique(ds))

  # Return the predicted retention times.
  return (stats::predict(object = mdl[[idx]], newdata = rt))

}
