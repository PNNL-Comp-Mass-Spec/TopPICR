#' Include additional variables and UniProt accessions
#'
#' All UniProt accessions that match each amino acid sequence will be added to
#' the \code{data.table}. Currently the \code{data.table} only consists of the
#' highest scoring UniProt accession for each amino acid sequence. Additional
#' variables are also created from the current variables.
#'
#' @param x A \code{data.table} output from the \code{read_toppic} function.
#'
#' @param fst_path A character string specifying the path to the protein
#'   database (.fasta) file.
#'
#' @param fst_name A character string containing the name of the .fasta file.
#'
#' @return A \code{data.table} with all protein accessions that match each amino
#'   acid sequence, not just the highest scoring accession. The number of rows
#'   in the output could be much larger than the number of rows in the input
#'   because all accessions for each sequence are included. The following
#'   variables have been added/removed:
#'
#'   | Added             | Removed                    |
#'   | ----------------- | -------------------------- |
#'   | `cleanSeq`        |                            |
#'   | `AccMap`          |                            |
#'   | `UniProtAcc`      |                            |
#'   | `protLength`      |                            |
#'   | `firstAA`         |                            |
#'   | `lastAA`          |                            |
#'   | `AnnType`         |                            |
#'   | `PercentCoverage` |                            |
#'
#' @md
#'
#' @importFrom magrittr %>%
#'
#' @author Evan A Martin
#'
#' @export
#'
augment_annotation <- function (x,
                                fst_path,
                                fst_name) {

  # Assemble AAStringSet objects -----------------------------------------------

  # Produce an AAStringSet object from the combined fasta file. This will be
  # separated into normal and decoy AAStringSet objects later.
  fst_both <- Biostrings::readAAStringSet(file.path(fst_path,
                                                    fst_name)) %>%
    Biostrings::chartr("X", "A", .) %>%
    Biostrings::chartr("B", "D", .) %>%
    Biostrings::chartr("Z", "E", .) %>%
    Biostrings::chartr("J", "I", .)

  # Augment annotation: normal -------------------------------------------------

  # Combine information from the fasta file with the input data to create the
  # Gene, UniProtAcc, AnnType, firstAA, lastAA, protLength, and percentCoverage
  # variables.
  x_norm <- enhance_annotation(
    x = dplyr::filter(x, !isDecoy),
    fst = ss_to_df(fst_both[grepl("^(sp|tr)", fst_both@ranges@NAMES)])
  )

  # Augment annotation: decoy --------------------------------------------------

  # Combine information from the fasta file with the input data to create the
  # Gene, UniProtAcc, AnnType, firstAA, lastAA, protLength, and percentCoverage
  # variables.
  x_decoy <- enhance_annotation(
    x = dplyr::filter(x, isDecoy),
    fst = ss_to_df(fst_both[grepl("^DECOY", fst_both@ranges@NAMES)])
  )

  # Combine the normal and decoy data matrices.
  return (rbind(x_norm, x_decoy))

}

# augment_annotation auxiliary functions ---------------------------------------

# Convert an AAStringSet to a data frame.
#
# @author Evan A Martin
ss_to_df <- function (stringset) {

  return (
    data.frame(
      accession = names(stringset),
      sequence = stringset,
      protLength = BiocGenerics::width(stringset),
      row.names = NULL,
      check.names = FALSE
    ) %>%
      dplyr::mutate(
        UniProtAcc = sub(".*[|]([^|]+)[|].*", "\\1", accession),
        AnnType = dplyr::case_when(grepl("^(DECOY_)?sp\\|[^-]*-\\d+\\|.*",
                                         accession) ~ "VarSplic",
                                   grepl("^(DECOY_)?tr.*",
                                         accession) ~ "TrEMBL",
                                   grepl("^(DECOY_)?sp\\|[^-]*\\|.*",
                                         accession) ~ "SwissProt",
                                   TRUE ~ NA_character_)
      ) %>%
      dplyr::select(-accession)
  )

}

# Add the Gene, UniProtAcc, AnnType, protLength, fisrtAA, lastAA, and
# percentCoverage variables to the data.
#
# @author Evan A Martin
enhance_annotation <- function (x, fst) {

  return (
    x %>%
      dplyr::mutate(
        cleanSeq = gsub("\\[.+?\\]", "", Proteoform),
        cleanSeq = gsub("\\(|\\)", "", cleanSeq),
        cleanSeq = sub("^[A-Z]?\\.(.*)\\.[A-Z]?", "\\1", cleanSeq)
      ) %>%
      dplyr::inner_join(fst) %>%
      dplyr::relocate(AnnType, .after = UniProtAcc) %>%
      dplyr::mutate(
        fl = purrr::map2(sequence, cleanSeq, stringr::str_locate),
        firstAA = purrr::map_int(fl, purrr::pluck, 1),
        lastAA = purrr::map_int(fl, purrr::pluck, 2),
        percentCoverage = signif(100 * (lastAA - firstAA + 1) / protLength, 2)
      ) %>%
      dplyr::select(-fl, -sequence)
  )

}
