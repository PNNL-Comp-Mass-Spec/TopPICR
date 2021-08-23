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
#'   | Added             | Removed                  |
#'   | ----------------- | ------------------------ |
#'   | `mz`              | `Prsm ID`                |
#'   | `Gene`            | `Spectrum ID`            |
#'   | `isDecoy`         | `Fragmentation`          |
#'   | `RTmin`           | `Retention time`         |
#'   | `ProtLen`         | `#peaks`                 |
#'   | `PercentCoverage` | `Proteoform ID`          |
#'   |                   | `Feature score`          |
#'   |                   | `MIScore`                |
#'   |                   | `#variable PTMs`         |
#'   |                   | `#matched peaks`         |
#'   |                   | `#matched fragment ions` |
#'   |                   | `Q-value (spectral FDR)` |
#'   |                   | `Proteoform FDR`         |
#'   |                   | `First residue`          |
#'   |                   | `Last residue`           |
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

  #!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!
  # Need to decide how the Data file name and CV variables will be treated for
  # the rest of the workflow.
  #!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!

  # Filter data and include additional variables -------------------------------

  # Create a character vector of variables that are not needed at any point in
  # the top down workflow. They will be removed with dplyr::any_of. The function
  # will not throw an error if any of the variables listed are not present in
  # the data set.
  useless <- c("Prsm ID", "Spectrum ID", "Fragmentation", "Retention time",
               "#peaks", "Proteoform ID", "Feature score", "MIScore",
               "#variable PTMs", "#matched peaks", "#matched fragment ions",
               "Q-value (spectral FDR)", "Proteoform FDR", "First residue",
               "Last residue")

  x <- x %>%
    # Filter data based on protein accession and description variables.
    dplyr::filter(grepl("GN=", `Protein description`)) %>%
    # Add useful variables.
    dplyr::mutate(
      mz = (`Precursor mass` + Charge * 1.007276466621) / Charge
    ) %>%
    dplyr::mutate(Gene = sub(".*GN=(\\S+).*","\\1",`Protein description`)) %>%
    dplyr::mutate(isDecoy = grepl("^DECOY", `Protein accession`)) %>%
    dplyr::mutate(RTmin = `Retention time` / 60) %>%
    # Remove not so useful variables.
    dplyr::select(-dplyr::any_of(useless))

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
    augment(fst_norm)

  # Shorten the names of the names attribute of the AAStringSet to just the
  # UniProt accession number.
  names(fst_norm) <- sub(".*[|]([^|]+)[|].*",
                         "\\1",
                         names(fst_norm))

  # Create a data frame with just the accession number and the length of the
  # amino acid sequence.
  prot_len_norm <- data.frame(UniProtAcc = names(fst_norm),
                              ProtLen = BiocGenerics::width(fst_norm),
                              stringsAsFactors = FALSE)

  # Add the ProtLen column to x_norm by matching to the UniProtAcc column.
  x_norm <- dplyr::inner_join(x_norm, prot_len_norm)

  # Add PercentCoverage to x_norm. This is the percentage each observed amino
  # acid sequence matches the reference protein's amino acid sequence.
  x_norm <- x_norm %>%
    dplyr::mutate(PercentCoverage = signif(100 * (Last_AA - First_AA + 1) /
                                             ProtLen, 2))

  # Augment annotation: decoy/scrambled ----------------------------------------

  # Augment the decoy data frame. This will search all three scrambled data
  # bases and keep all matches to each clean sequence (each value of cleanSeq).
  x_decoy <- x %>%
    dplyr::filter(isDecoy) %>%
    augment(fst_decoy)

  # Shorten the names of the names attribute of the AAStringSet to just the
  # UniProt accession number.
  names(fst_decoy) <- sub(".*[|]([^|]+)[|].*",
                          "\\1",
                          names(fst_decoy))

  # Create a data frame with just the accession number and the length of the
  # amino acid sequence.
  prot_len_decoy <- data.frame(UniProtAcc = names(fst_decoy),
                               ProtLen = BiocGenerics::width(fst_decoy),
                               stringsAsFactors = FALSE)

  # Add the ProtLen column to x_norm by matching to the UniProtAcc column.
  x_decoy <- dplyr::inner_join(x_decoy, prot_len_decoy)

  # Add PercentCoverage to x_decoy. This is the percentage each observed amino
  # acid sequence matches the reference protein's amino acid sequence.
  x_decoy <- x_decoy %>%
    dplyr::mutate(PercentCoverage = signif(100 * (Last_AA - First_AA + 1) /
                                             ProtLen, 2))

  # Combine the normal and decoy proteins.
  x <- rbind(x_norm, x_decoy)

  # Add a column for the ProteoForm variable.
  x <- x %>%
    dplyr::mutate(ProteoForm = paste(Gene, `First_AA`, `Last_AA`,
                                     round(`Adjusted precursor mass`,
                                           digits = -1),
                                     sep = "_"))

  return (x)

}

# augment_annotation auxiliary functions ---------------------------------------

augment <- function (x,
                     fst) {

  # Create a new data frame from the input data frame x.
  gene_sequence_pairs <-  x %>%
    # Create a new variable (or column) from the original Proteoform column.
    # This variable is a character string that only contains uppercase letters
    # from the Proteoform column. All other characters are removed.
    dplyr::mutate(
      cleanSeq = gsub("\\[.+?\\]", "", Proteoform),
      cleanSeq = gsub("\\(|\\)", "", cleanSeq),
      cleanSeq = sub("^[A-Z]?\\.(.*)\\.[A-Z]?", "\\1", cleanSeq)
    ) %>%
    # Only keep the newly added cleanSeq variable and the Gene column.
    dplyr::select(cleanSeq, Gene) %>%
    # Keep all unique combinations of cleanSeq and Gene. Each unique cleanSeq
    # or Gene can appear in multiple rows if they are associated with multiple
    # unique values of the other variable.
    dplyr::distinct()

  # Loop through the data frame, gene_sequence_pairs, row-wise and extract
  # information from the reference data bases that match the given cleanSeq.
  # This will increase the number of rows for each cleanSeq by the number of
  # times each cleanSeq matches an entry in the reference data bases.
  y <- apply(gene_sequence_pairs, 1,
             function(x) get_matching_uniprotacc(x["Gene"],
                                                 x["cleanSeq"],
                                                 fst = fst))

  # Combine the augmented data frames row-wise.
  y <- Reduce(rbind, y)

  # Create a new data frame from x.
  res <- x %>%
    # In the new data frame create the cleanSeq variable (to match with the
    # augmented data frame created above).
    dplyr::mutate(
      cleanSeq = gsub("\\[.+?\\]", "", Proteoform),
      cleanSeq = gsub("\\(|\\)", "", cleanSeq),
      cleanSeq = sub("^[A-Z]?\\.(.*)\\.[A-Z]?", "\\1", cleanSeq)
    ) %>%
    # Match rows in res and y based on the keys (cleanSeq, UniProtAcc, and
    # AnnType) and keep all columns from both res and y.
    dplyr::inner_join(y)

  return (res)

}

# Finds rows in fst whose gene names appear in x (the data frame created from
# the x1, x2, ... data frames). The fst argument is the object of class
# AAStringSet. This comes from the Biostrings package and is the reference
# database.
get_matching_uniprotacc <- function (gene,
                                     sequence,
                                     fst) {

  # fst_i is an AAStringSet object whose names match the given gene. Only the
  # rows with unique width and seq values are kept.
  fst_i <- fst[grep(paste0("GN=", gene, "( |$)"), names(fst))] %>%
    # Only keep the unique rows of fst whose names contain the string "gene".
    unique()

  # The "sequence" comes from the augment function. The cleanSeq string (which
  # is the "sequence") is created there. names_i is the name of the sequence
  # from fst_i that matches the sequence in the input of get_matching_uniprotacc
  # function. The full name of the sequence is the output of this line
  names_i <- sapply(fst_i, grepl, pattern = sequence, USE.NAMES = TRUE) %>%
    # Finds the indices of fst_i where the matches occur.
    which() %>%
    # Extracts the names from the fst object. These will be used for subsetting
    # and string matching later in this function.
    names()

  # acc_names_i is just the uniprot accession number. The other information from
  # the matching fst name has been removed.
  acc_names_i <- names_i %>%
    sub(".*[|]([^|]+)[|].*","\\1",.)

  # Finds the elements in fst_i that match the value in names_i and compares the
  # sequence from the reference data base (fst_i) to the given sequence. The
  # first element of matches_i is the index of the reference sequence where the
  # observed sequence begins to match the reference sequence. For example, if
  # the reference sequence is MPKRKVSSAEGAAKEE and the observed sequence is
  # PKRKVSSAEGAAKEE the first element in matches_i will be 2 because the
  # observed sequence starts to match the reference sequence at the second
  # position.
  matches_i <- lapply(fst_i[names_i], regexpr, pattern = sequence)

  # Takes the number returned from the previous line and coerces it to a numeric
  # value.
  first_aa <- as.numeric(unlist(matches_i))
  last_aa <- first_aa + nchar(sequence) - 1

  # Combine the output from names_i, acc_names_i, first_aa, and last_aa into a
  # data frame. This data frame will have multiple rows depending on the number
  # of times a sequence has a match in the reference data frames.
  out <- data.frame(Gene = gene,
                    names_i = names_i,
                    cleanSeq = sequence,
                    UniProtAcc = acc_names_i,
                    First_AA = first_aa,
                    Last_AA = last_aa,
                    stringsAsFactors = FALSE) %>%
    # Create the AnnType variable and determine which data base the sequence
    # comes from.
    dplyr::mutate(
      AnnType = dplyr::case_when(grepl("^(DECOY_)?sp\\|[^-]*-\\d+\\|.*",
                                       names_i) ~ "VarSplic",
                                 grepl("^(DECOY_)?tr.*",
                                       names_i) ~ "TrEMBL",
                                 grepl("^(DECOY_)?sp\\|[^-]*\\|.*",
                                       names_i) ~ "SwissProt",
                                 TRUE ~ NA_character_)
    ) %>%
    # Remove the names_i variable.
    dplyr::select(-names_i)

  return (out)

}
