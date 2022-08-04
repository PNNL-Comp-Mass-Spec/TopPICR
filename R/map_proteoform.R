#' Map proteoforms to proteins
#'
#' Add variables for a given proteoform in relation to the protein(s) it maps
#' to. For example, the location of the first and last amino acid where the
#' proteoform matches the protein, percentage of the protein the proteoform
#' covers, and length of the protein.
#'
#' @param x A \code{data.table} output from the \code{read_toppic} function.
#'
#' @param fst_path A character string specifying the path to the protein
#'   database (.fasta) file.
#'
#' @param fst_name A character string containing the name of the .fasta file.
#'
#' @return A \code{data.table} with additional variables containing information
#'   on the sequence in relation to its parent protein. For example, the
#'   location of the first and last amino acid, percent coverage of the parent
#'   protein, and protein length.
#'
#' @md
#'
#' @importFrom magrittr %>%
#'
#' @author Evan A Martin
#'
#' @export
#'
map_proteoform <- function (x,
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
  x_norm <- create_mapping(
    x = dplyr::filter(x, !isDecoy),
    fst = ss_to_df(fst_both[grepl("^(sp|tr)", fst_both@ranges@NAMES)])
  )

  # Augment annotation: decoy --------------------------------------------------

  # Combine information from the fasta file with the input data to create the
  # Gene, UniProtAcc, AnnType, firstAA, lastAA, protLength, and percentCoverage
  # variables.
  x_decoy <- create_mapping(
    x = dplyr::filter(x, isDecoy),
    fst = ss_to_df(fst_both[grepl("^DECOY", fst_both@ranges@NAMES)])
  )

  # Combine the normal and decoy data matrices.
  return (rbind(x_norm, x_decoy))

}

# map_proteoform auxiliary functions ---------------------------------------

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
create_mapping <- function (x, fst) {

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
