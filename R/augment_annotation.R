#' Include additional variables and UniProt accessions
#'
#' All UniProt accessions that match each amino acid sequence will be added to
#' the \code{data.table}. Currently the \code{data.table} only consits of the
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

  # Extract all entries for normal proteins (proteins whose names start with
  # either sp or tr).
  fst_norm <- fst_both[grepl("^(sp|tr)", fst_both@ranges@NAMES)]

  # Augment annotation: normal -------------------------------------------------

  # Augment the normal data frame. This will search all three normal data bases
  # and keep all matches to each clean sequence (each value of cleanSeq).
  x_norm <- x %>%
    dplyr::filter(!isDecoy) %>%
    dplyr::mutate(
      cleanSeq = gsub("\\[.+?\\]", "", Proteoform),
      cleanSeq = gsub("\\(|\\)", "", cleanSeq),
      cleanSeq = sub("^[A-Z]?\\.(.*)\\.[A-Z]?", "\\1", cleanSeq)
    )

  # Extract all unique clean sequences. This vector will be used to find all
  # matching UniProt accessions for each sequence.
  peptides <- unique(x_norm$cleanSeq)

  # Find all matching UniProt accessions for each clean sequence.
  augmented <- map_peptides_to_fasta(peptides,
                                     fst_norm) %>%
    dplyr::mutate(
      AnnType2 = dplyr::case_when(
        grepl("sp\\|", accession) ~ "SwissProt",
        grepl("tr\\|", accession) ~ "TrEMBL"
      )
    ) %>%
    dplyr::group_by(cleanSeq, AnnType2) %>%
    dplyr::add_count(name = "AccMap") %>%
    dplyr::ungroup() %>%
    dplyr::select(-AnnType2)

  x_norm <- dplyr::inner_join(x_norm, augmented) %>%
    dplyr::mutate(
      UniProtAcc = sub(".*[|]([^|]+)[|].*", "\\1", accession)
    ) %>%
    dplyr::select(-accession)

  # Create the protLength, firstAA, and lastAA variables.
  acc_info <- add_acc_info(data = x_norm,
                           fst = fst_norm)

  # Add the variables created in the acc_info object to the x_norm object.
  x_norm <- dplyr::inner_join(x = x_norm,
                              y = acc_info)

  # Add AnnType and percentCoverage (the percentage each observed amino
  # acid sequence matches the reference amino acid sequence) to x_norm.
  x_norm <- x_norm %>%
    dplyr::mutate(
      AnnType = dplyr::case_when(grepl("^(DECOY_)?sp\\|[^-]*-\\d+\\|.*",
                                       `Protein accession`) ~ "VarSplic",
                                 grepl("^(DECOY_)?tr.*",
                                       `Protein accession`) ~ "TrEMBL",
                                 grepl("^(DECOY_)?sp\\|[^-]*\\|.*",
                                       `Protein accession`) ~ "SwissProt",
                                 TRUE ~ NA_character_),
      percentCoverage = signif(100 * (lastAA - firstAA + 1) / protLength, 2)
    )

  # Augment annotation: decoy/scrambled ----------------------------------------

  # Augment the decoy data frame. Most of the variables included in the normal
  # data frame are not necessary in the decoy data frame because the decoy rows
  # will be deleted in the FDR control step.
  x_decoy <- x %>%
    dplyr::filter(isDecoy) %>%
    dplyr::mutate(
      cleanSeq = NA,
      AccMap = NA,
      UniProtAcc = NA,
      protLength = NA,
      firstAA = NA,
      lastAA = NA,
      AnnType = dplyr::case_when(grepl("^(DECOY_)?sp\\|[^-]*-\\d+\\|.*",
                                       `Protein accession`) ~ "VarSplic",
                                 grepl("^(DECOY_)?tr.*",
                                       `Protein accession`) ~ "TrEMBL",
                                 grepl("^(DECOY_)?sp\\|[^-]*\\|.*",
                                       `Protein accession`) ~ "SwissProt",
                                 TRUE ~ NA_character_),
      percentCoverage = NA
    )

  # Combine the normal and decoy proteins.
  x <- rbind(x_norm, x_decoy)

  return (x)

}

# augment_annotation auxiliary functions ---------------------------------------

map_peptides_to_fasta <- function (peptides,
                                   fst,
                                   numCores = parallel::detectCores() - 1) {

  fasta.df <- data.frame(accession = names(fst),
                         sequence = fst,
                         row.names = NULL,
                         check.names = FALSE)

  cl <- parallel::makeCluster(numCores)

  # Add the pipe and the find_matching_accessions function so they can be used
  # in parallel.
  `%>%` <- magrittr::`%>%`
  parallel::clusterExport(cl, "%>%", envir = 1)
  find_matching_accessions <- TopPICR:::find_matching_accessions
  parallel::clusterExport(cl, "find_matching_accessions", envir = 1)

  fasta.df.split <- lapply(parallel::clusterSplit(cl, 1:nrow(fasta.df)),
                           function(idx) {
                             fasta.df[idx, , drop = FALSE]
                           })

  out <- parallel::clusterApply(
    cl,
    fasta.df.split,
    find_matching_accessions,
    peptides = peptides
  )

  parallel::stopCluster(cl)

  out <- unlist(out, recursive = FALSE)
  out <- rbindlist(out)

}

find_matching_accessions <- function (peptides, fasta.df) {
  lapply(peptides, function (peptide) {
    fasta.df %>%
      dplyr::mutate(is_match = grepl(peptide, sequence)) %>%
      dplyr::filter(is_match) %>%
      dplyr::mutate(cleanSeq = peptide) %>%
      dplyr::select(cleanSeq, accession, -sequence, -is_match)
  })
}

add_acc_info <- function (data, fst) {

  # Only keep the unique combinations of UniProtAcc and cleanSeq to speed up the
  # for loop (only compute values once).
  uniqueD <- data %>%
    dplyr::distinct(UniProtAcc, cleanSeq)

  # Initialize vectors to full length.
  prot_len <- vector(length = dim(uniqueD)[[1]])
  first_aa <- vector(length = dim(uniqueD)[[1]])
  last_aa <- vector(length = dim(uniqueD)[[1]])

  # Loop through each accession and observed sequence to create the protein
  # length, first amino acid position, and last amino acid position variables.
  for (e in 1:length(prot_len)) {

    # Find the current accession in the fasta file.
    idx <- grep(uniqueD$UniProtAcc[[e]], names(fst))[[1]]

    # Convert the protein sequence to a string so it can be matched to the
    # observed amino acid sequence.
    real <- Biostrings::toString(fst[idx])

    # Find the first and last position of the observed sequence in the reference
    # sequence.
    fl <- stringr::str_locate(real, uniqueD$cleanSeq[[e]])

    # Assign computed values to their corresponding vector.
    prot_len[[e]] <- BiocGenerics::width(fst[idx])
    first_aa[[e]] <- fl[[1]]
    last_aa[[e]] <- fl[[2]]

  }

  # Return the first and last indices of the amino acid sequence and the
  # observed protein length.
  return (data.frame(UniProtAcc = uniqueD$UniProtAcc,
                     cleanSeq = uniqueD$cleanSeq,
                     protLength = prot_len,
                     firstAA = first_aa,
                     lastAA = last_aa))

}
