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
#'   because all accessions for each sequence are included.
#'
#' @md
#'
#' @importFrom magrittr %>%
#'
#' @author Vlad Petyuk, Evan A Martin
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

  # Extract all entries for scrambled proteins (proteins whose names start with
  # DECOY).
  fst_decoy <- fst_both[grepl("^DECOY", fst_both@ranges@NAMES)]

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
    dplyr::group_by(cleanSeq) %>%
    dplyr::add_count(name = "AccMap") %>%
    dplyr::ungroup()

  # Add the UniProtAcc and AnnType variables using the accession variable.
  x_norm <- dplyr::inner_join(x_norm, augmented) %>%
    dplyr::mutate(
      Gene = stringr::str_extract(accession, "GN=\\S+"),
      Gene = stringr::str_replace(Gene, "GN=", ""),
      Gene = dplyr::case_when(is.na(Gene) ~ "",
                              TRUE ~ Gene),
      UniProtAcc = sub(".*[|]([^|]+)[|].*", "\\1", accession),
      AnnType = dplyr::case_when(grepl("^(DECOY_)?sp\\|[^-]*-\\d+\\|.*",
                                       accession) ~ "VarSplic",
                                 grepl("^(DECOY_)?tr.*",
                                       accession) ~ "TrEMBL",
                                 grepl("^(DECOY_)?sp\\|[^-]*\\|.*",
                                       accession) ~ "SwissProt",
                                 TRUE ~ NA_character_),
      percentCoverage = signif(100 * (lastAA - firstAA + 1) / protLength, 2)
    ) %>%
    dplyr::select(-accession)

  # Augment annotation: decoy --------------------------------------------------

  # Augment the decoy data frame.
  x_decoy <- x %>%
    dplyr::filter(isDecoy) %>%
    dplyr::mutate(
      cleanSeq = gsub("\\[.+?\\]", "", Proteoform),
      cleanSeq = gsub("\\(|\\)", "", cleanSeq),
      cleanSeq = sub("^[A-Z]?\\.(.*)\\.[A-Z]?", "\\1", cleanSeq)
    )

  # Extract all unique clean sequences. This vector will be used to find all
  # matching UniProt accessions for each sequence.
  peptides <- unique(x_decoy$cleanSeq)

  # Find all matching UniProt accessions for each clean sequence.
  augmented <- map_peptides_to_fasta(peptides,
                                     fst_decoy) %>%
    dplyr::group_by(cleanSeq) %>%
    dplyr::add_count(name = "AccMap") %>%
    dplyr::ungroup()

  # Add the UniProtAcc and AnnType variables using the accession variable.
  x_decoy <- dplyr::inner_join(x_decoy, augmented) %>%
    dplyr::mutate(
      Gene = stringr::str_extract(accession, "GN=\\S+"),
      Gene = stringr::str_replace(Gene, "GN=", ""),
      Gene = dplyr::case_when(is.na(Gene) ~ "",
                              TRUE ~ Gene),
      UniProtAcc = sub(".*[|]([^|]+)[|].*", "\\1", accession),
      AnnType = dplyr::case_when(grepl("^(DECOY_)?sp\\|[^-]*-\\d+\\|.*",
                                       accession) ~ "VarSplic",
                                 grepl("^(DECOY_)?tr.*",
                                       accession) ~ "TrEMBL",
                                 grepl("^(DECOY_)?sp\\|[^-]*\\|.*",
                                       accession) ~ "SwissProt",
                                 TRUE ~ NA_character_),
      percentCoverage = signif(100 * (lastAA - firstAA + 1) / protLength, 2)
    ) %>%
    dplyr::select(-accession)

  # Combine the normal and decoy data tables and return the combined table.
  return (rbind(x_norm, x_decoy))

}

# augment_annotation auxiliary functions ---------------------------------------

# @author Michael Nestor
map_peptides_to_fasta <- function (peptides,
                                   fst,
                                   numCores = parallel::detectCores() - 1) {

  fasta.df <- data.frame(accession = names(fst),
                         sequence = fst,
                         width = BiocGenerics::width(fst),
                         row.names = NULL,
                         check.names = FALSE)

  cl <- parallel::makeCluster(numCores)

  # Add the pipe and the find_matching_accessions function so they can be used
  # in parallel.
  parallel::clusterCall(cl, assign,
                        "find_matching_accessions",
                        find_matching_accessions,
                        envir = .GlobalEnv)
  parallel::clusterCall(cl, assign,
                        "%>%",
                        magrittr::`%>%`,
                        envir = .GlobalEnv)

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

  return (rbindlist(out))

}

# @author Michael Nestor
#
find_matching_accessions <- function (peptides, fasta.df) {
  lapply(peptides, function (peptide) {
    fasta.df %>%
      dplyr::mutate(is_match = grepl(peptide, sequence)) %>%
      dplyr::filter(is_match) %>%
      dplyr::mutate(cleanSeq = peptide,
                    protLength = width,
                    firstAA = stringr::str_locate(sequence, cleanSeq)[, 1],
                    lastAA = stringr::str_locate(sequence, cleanSeq)[, 2]) %>%
      dplyr::select(cleanSeq, accession, protLength, firstAA, lastAA)
  })
}

#'
#' @importFrom foreach `%dopar%`
#'
#' @author Evan A Martin
#'

add_acc_info <- function (data, fst) {

  # Only keep the unique combinations of UniProtAcc and cleanSeq to speed up the
  # for loop (only compute values once).
  uniqueD <- data %>%
    dplyr::distinct(UniProtAcc, cleanSeq)

  # Prepare to run in parallel.
  cores <- parallel::detectCores() - 1
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)

  # Loop through each accession and observed sequence to create the protein
  # length, first amino acid position, and last amino acid position variables.
  posi <- foreach::foreach(e = 1:dim(uniqueD)[[1]]) %dopar% {

    # Find the current accession in the fasta file.
    idx <- grep(uniqueD$UniProtAcc[[e]], names(fst))[[1]]

    # Convert the protein sequence to a string so it can be matched to the
    # observed amino acid sequence.
    real <- Biostrings::toString(fst[idx])

    # Find the first and last position of the observed sequence in the reference
    # sequence.
    fl <- stringr::str_locate(real, uniqueD$cleanSeq[[e]])

    # Assign computed values to a data frame.
    positions <- data.frame(
      prot_len = BiocGenerics::width(fst[idx]),
      first_aa = fl[[1]],
      last_aa = fl[[2]]
    )

    # Using the foreach function with %dopar% will assign the last element
    # within the curly brackets to the object when foreach is called. In this
    # case positions will be assigned to the eth element of posi.
    positions

  }

  parallel::stopCluster(cl)

  # Combine the list containing prot_len, first_aa, and last_aa into a
  # data.table.
  posi <- data.table::rbindlist(posi)

  # Return the first and last indices of the amino acid sequence and the
  # observed protein length.
  return (data.frame(UniProtAcc = uniqueD$UniProtAcc,
                     cleanSeq = uniqueD$cleanSeq,
                     protLength = posi$prot_len,
                     firstAA = posi$first_aa,
                     lastAA = posi$last_aa))

}
