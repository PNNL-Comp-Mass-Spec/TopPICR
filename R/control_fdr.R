# FDR control preliminary functions --------------------------------------------

#' Calculate p-value cutoff
#'
#' Calculates the p-value cutoff that keeps the FDR under the given threshold
#' within each annotation type.
#'
#' @param x A \code{data.table} output from the \code{augment_annotation}
#'   function.
#'
#' @param fdr_threshold A value between 0 and 1 indicating the desired FDR
#'   control level.
#'
#' @return A vector with three elements. These elements correspond to the
#'   E-value cutoff for the three annotation types (SwissProt, VarSplic, and
#'   TrEMBL). Note: the order of the annotation types depends on the order they
#'   appear in the AnnType column.
#'
#' @export
#'
pval_cutoff <- function (x, fdr_threshold) {

  # Create a character vector of annotation types.
  anns <- unique(x$AnnType)

  # Compute the p-value cutoff for each of the annotation types.
  e_vals <- sapply(anns, find_cutoff, x = x, fdr_threshold = fdr_threshold)

  return (e_vals)

}

# FDR control main function ----------------------------------------------------

#' Filter the data by p-value
#'
#' The \code{control_fdr()} function filters the data by E-value for each
#' annotation type. If the E-values from the \code{eval_cutoff()} function are
#' used the overall FDR will fall below the given threshold.
#'
#' @param x A \code{data.table} output from the \code{augment_annotation}
#'   function.
#'
#' @param e_vals A vector of three numeric values to be used as E-value cutoffs.
#'   The E-value cutoffs must be in the order of unique(x$AnnType). For example,
#'   if unique(x$AnnType) = "SwissProt", "TrEMBL", "VarSplic" then the E-value
#'   cutoff for SwissProt must be the first element of e_vals, the E-value
#'   cutoff for TrEMBL must be the second element of e_vals, and the E-value
#'   cutoff for VarSplic must be the third element of e_vals.
#'
#' @return A \code{data.table} with the rows removed that will produce an FDR
#'   below the given threshold at the gene level.
#'
#' @importFrom magrittr %>%
#'
#' @export
#'
control_fdr <- function (x, e_vals) {

  # Create a character vector of annotation types.
  anns <- unique(x$AnnType)

  # Only keep the rows with a p-value below the calculated FDR threshold for
  # each of the three annotation types (SwissProt, VarSplic, TrEMBL). Remove the
  # column that specifies if a protein is scrambled.
  x <- x %>%
    dplyr::filter((AnnType == anns[[1]] & `E-value` <= e_vals[[1]]) |
                    (AnnType == anns[[2]] & `E-value` <= e_vals[[2]]) |
                    (AnnType == anns[[3]] & `E-value` <= e_vals[[3]])) %>%
    dplyr::filter(!isDecoy) %>%
    dplyr::select(-isDecoy)

  return (x)

}

# FDR control auxiliary functions ----------------------------------------------

find_cutoff <- function (x, ann_type, fdr_threshold) {

  # Create a copy of x filtered by annotation type. This subsetted data frame
  # will be used multiple times throughout this function.
  x_filtered <- x %>%
    dplyr::filter(AnnType == ann_type) %>%
    dplyr::select(`E-value`, Gene, isDecoy)

  # Initialize the interval vector to the minimum and maximum values of the
  # `E-value` vector. This will be used to start the while loop.
  pval_interval <- c(min(x_filtered$`E-value`),
                     max(x_filtered$`E-value`))

  # Reduce the `E-value` vector until the interval consists on only one value.
  # This will happen when all p-values in the `E-value` vector produce and FDR
  # less than the threshold.
  while (length(pval_interval) != 1) {

    # Only keep the rows with p-values within the p-value interval.
    filtad <- x_filtered %>%
      # Only keep rows whose p-value is larger than the lower limit AND smaller
      # than the upper limit.
      dplyr::filter(`E-value` > pval_interval[[1]] &
                      `E-value` < pval_interval[[2]])

    # Generate a new p-value interval based on the filtered data set.
    pval_interval <- make_interval(x = filtad,
                                   fdr_threshold = fdr_threshold)

  }

  return (pval_interval)

}

make_interval <- function (x, fdr_threshold) {

  # Create a summary of the p-values for the given annotation type. This summary
  # will be used to find a smaller range of p-values to search in order to find
  # the largest p-value cutoff that produces an FDR below the threshold.
  summary_e_vals <- x %>%
    dplyr::pull(`E-value`) %>%
    quantile() %>%
    as.numeric()

  # Compute the FDR for each p-value in the summary.
  fdrs <- sapply(summary_e_vals, compute_fdr, x = x)

  # Find which FDR values is above the FDR threshold. This vector will be used
  # to find the p-values that produce FDR values that are just above and just
  # below the threshold. The largest value in this vector corresponds to the
  # index of the p-value in summary_e_vals that still has an FDR below the
  # threshold. This p-value and the following p-value (in the summary_e_vals
  # vector) will be used later for the fine tune search.
  which_fdrs <- which(fdrs < fdr_threshold)

  # Compare the length of the which_fdrs and summary_e_vals vectors. If they are
  # the same length the largest p-value for the given annotation type still
  # produces an FDR below the threshold and this p-value can be returned for
  # later filtering.
  if (length(summary_e_vals) == length(which_fdrs)) {

    # Return the maximum p-value.
    return (summary_e_vals[[5]])

    # Runs when the length of which_fdrs is smaller than summary_e_vals. This
    # occurs when there is at least one p-value in summary_e_vals that creates an
    # FDR less than the threshold and at least one p-value in summary_e_vals that
    # creates an FDR greater than the threshold.
  } else {

    # Return the p-value interval that will be searched.
    return (c(summary_e_vals[[max(which_fdrs)]],
              summary_e_vals[[(max(which_fdrs) + 1)]]))

  }

}

compute_fdr <- function (pval_cutoff, x) {

  fdr <- x %>%
    dplyr::filter(`E-value` <= pval_cutoff) %>%
    dplyr::distinct(Gene, isDecoy) %>%
    dplyr::pull(isDecoy) %>%
    # Mean of a logical vector gives the proportion of TRUE values.
    mean()

}
