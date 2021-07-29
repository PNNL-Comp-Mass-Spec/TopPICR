# FDR control preliminary functions --------------------------------------------

#' Calculate p-value cutoff
#'
#' Calculates the p-value cutoff that keeps the FDR under the given threshold
#' within each annotation type.
#'
#' @param x ...
#'
#' @param fdr_threshold ...
#'
#' @return ...
#'
#' @export
#'
pval_cutoff <- function (x, fdr_threshold) {

  # Create a character vector of annotation types.
  anns <- unique(x$AnnType)

  # Compute the p-value cutoff for each of the annotation types.
  pvals <- sapply(anns, find_cutoff, x = x, fdr_threshold = fdr_threshold)

  return (pvals)

}

# FDR control main function ----------------------------------------------------

#' Filter the data by p-value
#'
#' The \code{control_fdr()} function filters the data by p-value for each
#' annotation type. If the p-values from the \code{pval_cutoff()} function are
#' used the overall FDR will fall below the given threshold.
#'
#' @param x ...
#'
#' @param pvals A vector of three numeric values to be used as p-value cutoffs.
#'   The p-value cutoffs must be in the order of unique(x$AnnType). For example,
#'   if unique(x$AnnType) = "SwissProt", "TrEMBL", "VarSplic" then the p-value
#'   cutoff for SwissProt must be the first element of pvals, the p-value cutoff
#'   for TrEMBL must be the second element of pvals, and the p-value cutoff for
#'   VarSplic must be the third element of pvals.
#'
#' @return ...
#'
#' @importFrom magrittr %>%
#'
#' @export
#'
control_fdr <- function (x, pvals) {

  # Create a character vector of annotation types.
  anns <- unique(x$AnnType)

  # Only keep the rows with a p-value below the calculated FDR threshold for
  # each of the three annotation types (SwissProt, VarSplic, TrEMBL). Remove the
  # column that specifies if a protein is scrambled.
  x <- x %>%
    dplyr::filter((AnnType == anns[[1]] & `P-value` <= pvals[[1]]) |
                    (AnnType == anns[[2]] & `P-value` <= pvals[[2]]) |
                    (AnnType == anns[[3]] & `P-value` <= pvals[[3]])) %>%
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
    dplyr::select(`P-value`, Gene, isDecoy)

  # Initialize the interval vector to the minimum and maximum values of the
  # `P-value` vector. This will be used to start the while loop.
  pval_interval <- c(min(x_filtered$`P-value`),
                     max(x_filtered$`P-value`))

  # Reduce the `P-value` vector until the interval consists on only one value.
  # This will happen when all p-values in the `P-value` vector produce and FDR
  # less than the threshold.
  while (length(pval_interval) != 1) {

    # Only keep the rows with p-values within the p-value interval.
    filtad <- x_filtered %>%
      # Only keep rows whose p-value is larger than the lower limit AND smaller
      # than the upper limit.
      dplyr::filter(`P-value` > pval_interval[[1]] &
                      `P-value` < pval_interval[[2]])

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
  summary_pvals <- x %>%
    dplyr::pull(`P-value`) %>%
    quantile() %>%
    as.numeric()

  # Compute the FDR for each p-value in the summary.
  fdrs <- sapply(summary_pvals, compute_fdr, x = x)

  # Find which FDR values is above the FDR threshold. This vector will be used
  # to find the p-values that produce FDR values that are just above and just
  # below the threshold. The largest value in this vector corresponds to the
  # index of the p-value in summary_pvals that still has an FDR below the
  # threshold. This p-value and the following p-value (in the summary_pvals
  # vector) will be used later for the fine tune search.
  which_fdrs <- which(fdrs < fdr_threshold)

  # Compare the length of the which_fdrs and summary_pvals vectors. If they are
  # the same length the largest p-value for the given annotation type still
  # produces an FDR below the threshold and this p-value can be returned for
  # later filtering.
  if (length(summary_pvals) == length(which_fdrs)) {

    # Return the maximum p-value.
    return (summary_pvals[[5]])

    # Runs when the length of which_fdrs is smaller than summary_pvals. This
    # occurs when there is at least one p-value in summary_pvals that creates an
    # FDR less than the threshold and at least one p-value in summary_pvals that
    # creates an FDR greater than the threshold.
  } else {

    # Return the p-value interval that will be searched.
    return (c(summary_pvals[[max(which_fdrs)]],
              summary_pvals[[(max(which_fdrs) + 1)]]))

  }

}

compute_fdr <- function (pval_cutoff, x) {

  fdr <- x %>%
    dplyr::filter(`P-value` <= pval_cutoff) %>%
    dplyr::distinct(Gene, isDecoy) %>%
    dplyr::pull(isDecoy) %>%
    # Mean of a logical vector gives the proportion of TRUE values.
    mean()

}
