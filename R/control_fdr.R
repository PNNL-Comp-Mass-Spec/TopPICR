# FDR control preliminary functions --------------------------------------------

#' Calculate E-value cutoff
#'
#' Calculates the E-value cutoff that keeps the FDR under the given threshold
#' within each annotation type.
#'
#' @param x A \code{data.table} output from either the \code{augment_annotation}
#'   or \code{map_proteoform} function.
#'
#' @param fdr_threshold A value between 0 and 1 indicating the desired FDR
#'   level.
#'
#' @return A vector with three elements. These elements correspond to the
#'   E-value cutoff for the three annotation types (SwissProt, VarSplic, and
#'   TrEMBL). Note: the order of the annotation types depends on the order they
#'   appear in the AnnType column.
#'
#' @author Evan A Martin
#'
#' @export
#'
find_evalue_cutoff <- function (x, fdr_threshold) {

  # Create a character vector of annotation types.
  anns <- unique(x$AnnType)

  # Compute the E-value cutoff for each of the annotation types.
  e_vals <- sapply(anns, find_cutoff, x = x, fdr_threshold = fdr_threshold)

  return (e_vals)

}

#' Compute the FDR
#'
#' Computes the FDR at the PrSM, sequence, proteoform, accession, and gene
#' levels.
#'
#' @param x A \code{data.table} that contains decoy proteins.
#'
#' @return A list containing the FDR and counts at the PrSM, sequence,
#'   proteoform, accession, and gene levels. The FDR values and counts for each
#'   level are also printed when the function is called. All FDR values are
#'   reported as percentages.
#'
#' @author Evan A Martin
#'
#' @export
#'
compute_fdr <- function (x) {

  prsm <- x %>%
    dplyr::distinct(Dataset, `Scan(s)`, isDecoy) %>%
    dplyr::summarize(fdr = mean(isDecoy),
                     count = sum(!isDecoy))

  proteoform <- x %>%
    dplyr::distinct(Proteoform, isDecoy) %>%
    dplyr::summarize(fdr = mean(isDecoy),
                     count = sum(!isDecoy))

  cleanseq <- x %>%
    dplyr::distinct(cleanSeq, isDecoy) %>%
    dplyr::summarize(fdr = mean(isDecoy),
                     count = sum(!isDecoy))

  accession <- x %>%
    dplyr::distinct(UniProtAcc, isDecoy) %>%
    dplyr::summarize(fdr = mean(isDecoy),
                     count = sum(!isDecoy))

  gene <- x %>%
    dplyr::distinct(Gene, isDecoy) %>%
    dplyr::summarize(fdr = mean(isDecoy),
                     count = sum(!isDecoy))

  # Print the output then return a vector containing the FDR values.
  cat("PrSM level FDR: ", round(prsm$fdr * 100, 3), "% \n",
      "Number of unique PrSMs: ", prsm$count, "\n\n",

      "Proteoform level FDR: ", round(proteoform$fdr * 100, 3), "% \n",
      "Number of unique Proteoforms: ", proteoform$count, "\n\n",

      "Sequence level FDR: ", round(cleanseq$fdr * 100, 3), "% \n",
      "Number of unique Sequences: ", cleanseq$count, "\n\n",

      "UniProt accession level FDR: ", round(accession$fdr * 100, 3), "% \n",
      "Number of unique UniProt accessions: ", accession$count, "\n\n",

      "Gene level FDR: ", round(gene$fdr * 100, 3), "% \n",
      "Number of unique genes: ", gene$count, "\n",

      sep = "")

  x <- list(PrSM = c(round(prsm$fdr * 100, 3), prsm$count),
            Proteoform = c(round(proteoform$fdr * 100, 3), proteoform$count),
            cleanSeq = c(round(cleanseq$fdr * 100, 3), cleanseq$count),
            UniProtAcc = c(round(accession$fdr * 100, 3), accession$count),
            Gene = c(round(gene$fdr * 100, 3), gene$count))

}

# FDR control main function ----------------------------------------------------

#' Filter the data by E-value
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
#'   if \code{unique(x$AnnType) = "SwissProt", "TrEMBL", "VarSplic"} then the
#'   E-value cutoff for SwissProt must be the first element of e_vals, the
#'   E-value cutoff for TrEMBL must be the second element of e_vals, and the
#'   E-value cutoff for VarSplic must be the third element of e_vals.
#'
#' @return A \code{data.table} with the rows removed that fall below the E-value
#'   cutoff.
#'
#' @md
#'
#' @author Evan A Martin
#'
#' @importFrom magrittr %>%
#'
#' @export
#'
apply_evalue_cutoff <- function (x, e_vals) {

  # Create a character vector of annotation types.
  anns <- unique(x$AnnType)

  # Only keep the rows with a E-value below the calculated FDR threshold for
  # each of the three annotation types (SwissProt, VarSplic, TrEMBL). Remove the
  # column that specifies if a protein is scrambled.
  x <- x %>%
    dplyr::filter((AnnType == anns[[1]] & `E-value` <= e_vals[[1]]) |
                    (AnnType == anns[[2]] & `E-value` <= e_vals[[2]]) |
                    (AnnType == anns[[3]] & `E-value` <= e_vals[[3]]))

  return (x)

}

# FDR control auxiliary functions ----------------------------------------------

# @author Evan A Martin
find_cutoff <- function (x, ann_type, fdr_threshold) {

  # Create a copy of x filtered by annotation type. This subsetted data frame
  # will be used multiple times throughout this function.
  x_filtered <- x %>%
    dplyr::filter(AnnType == ann_type) %>%
    dplyr::select(`E-value`, Gene, isDecoy)

  # Extract the minimum and maximum e-values for the given annotation type.
  # These values will be used throughout the function.
  mini <- min(x_filtered$`E-value`)
  maxi <- max(x_filtered$`E-value`)

  # Compute the FDR for the minimum and maximum values of the E-value vector. If
  # the FDR for the minimum value is larger than the threshold throw an error
  # because the FDR can never be below the threshold. If the FDR for the maximum
  # value is already below the threshold no filtering needs to occur.
  min_fdr <- calc_fdr(eval_cutoff = mini, x = x_filtered)
  max_fdr <- calc_fdr(eval_cutoff = maxi, x = x_filtered)


  # Check the FDR for the minimum E-value. If it is larger than the cutoff
  # provided give a warning that the FDR for the current annotation type cannot
  # be below the FDR threshold.
  if (min_fdr > fdr_threshold) {

    warning (paste("For the given data the ",
                   ann_type,
                   " FDR cannot be below the specified threshold.",
                   sep = ""))

    return (mini)

  }

  # Check the FDR for the maximum E-value. If it is smaller than the cutoff
  # provided no further steps need to be taken in order to keep the FDR below
  # the cutoff.
  if (max_fdr < fdr_threshold) {

    return (maxi)

  }

  # Order the E-values from smallest to largest. These values will be used to
  # calculate the FDR until the FDR passes the threshold.
  ordered <- x_filtered %>%
    dplyr::select(`E-value`) %>%
    dplyr::arrange(`E-value`) %>%
    dplyr::pull(`E-value`)

  # Divide the ordered E-values into equal sections. The E-values from the
  # selected indices will be used to calculate the FDR in search of the E-value
  # that produces and FDR as close as possible to the given FDR threshold.
  inds <- floor(quantile(1:dim(x_filtered)[[1]]))

  # Create a variable that is the distance between the first value of inds and
  # the last value of inds. This difference is the number of E-values between
  # the min and max of the index vector (or number of E-values in the vector).
  differ <- inds[[5]] - inds[[1]]

  while (differ >= 2) {

    # Extract the E-values from the ordered E-value vector that correspond to
    # the indices from the index summary.
    summary_evals <- ordered[inds]

    # Compute the FDR for each E-value in the summary. This vector will be used
    # to further subset the E-value vector until it finds the largest E-value
    # that produces an FDR below the specified threshold.
    fdrs <- sapply(summary_evals, calc_fdr, x = x_filtered)

    # Find which FDR values are above the FDR threshold. The smallest value in
    # this vector corresponds to the index of the E-value in summary_evals that
    # has an FDR above the threshold. This E-value will be used to subset and
    # order the `E-value` vector and the FDR will be calculated for each
    # subsequent E-value until the FDR falls below the threshold.
    which_fdrs <- which(fdrs < fdr_threshold)

    # Nab the largest E-value that produces an FDR below the threshold. This is
    # the E-value that will be returned when either of the conditions that break
    # the while loop are met.
    maximin <- ordered[[inds[[max(which_fdrs)]]]]

    # Create a new inds vector for the next iteration of the while loop.
    # Reduce the space of E-values searched by only considering the indices
    # between the index for the E-value just above and just below the FDR
    # threshold.
    inds <- floor(quantile(
      inds[[max(which_fdrs)]]:inds[[(max(which_fdrs) + 1)]]
    ))

    # Update the difference in indices.
    differ <- inds[[5]] - inds[[1]]

  }

  # Return the largest E-value that still produces an FDR below the threshold.
  return (maximin)

}

# @author Evan A Martin
calc_fdr <- function (eval_cutoff, x) {

  fdr <- x %>%
    dplyr::filter(`E-value` <= eval_cutoff) %>%
    dplyr::distinct(Gene, isDecoy) %>%
    dplyr::pull(isDecoy) %>%
    # Mean of a logical vector gives the proportion of TRUE values.
    mean()

}

# Post FDR control functions ---------------------------------------------------

#' Remove decoy proteins
#'
#' This function filters out the decoy proteins and also removes the
#' \code{isDecoy} variable.
#'
#' @param x A \code{data.table} that contains decoy proteins. This function
#'   should be called after \code{apply_evalue_cutoff} and before the
#'   \code{align_rt} function.
#'
#' @return A \code{data.table} with the rows containing decoy proteins and the
#'   \code{isDecoy} column removed.
#'
#' @author Evan A Martin
#'
#' @export
#'
remove_decoys <- function (x) {

  return (x %>%
            dplyr::filter(!isDecoy) %>%
            dplyr::select(-isDecoy))

}
