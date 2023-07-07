#' Read in TopPIC results from DMS.
#'
#' Works only within PNNL's DMS.
#'
#' @param data_package_num DMS data package number that contain TopPIC search
#'    results. Note, it may not work correctly if multipe TopPIC searches
#'    were done with different parameter or FASTA files.
#'
#' @return A list with MS/MS identification and MS features.
#'
#' @md
#'
#' @author Vlad Petyuk
#'
#' @export
#'
read_TopPIC_DMS <- function(data_package_num){

  jobRecords <- get_job_records_by_dataset_package(data_package_num)
  toppic_output <- get_results_for_multiple_jobs.dt(jobRecords, expected_multiple_files = TRUE)

  ids <- bind_rows(toppic_output$`_TopPIC_PrSMs.txt`)

  ids <- ids %>%
    dplyr::mutate(
      `Feature apex` = `Feature apex time`,
      firstAA = `First residue`,
      lastAA = `Last residue`,
      AccMap = `#Protein hits`,
      mz = (`Precursor mass` + Charge * 1.007276466621) / Charge,
      Gene = sub(".*GN=(\\S+).*", "\\1", `Protein description`),
      UniProtAcc = `Protein accession`,
      isDecoy = grepl("^DECOY", `Protein accession`))

  # annotation type. Note this is UniProt specific
  ids <- ids %>%
    dplyr::mutate(
      AnnType = dplyr::case_when(grepl("^(DECOY_)?sp\\|[^-]*-\\d+\\|.*",
                                       `Protein accession`) ~ "VarSplic",
                                 grepl("^(DECOY_)?tr.*",
                                       `Protein accession`) ~ "TrEMBL",
                                 grepl("^(DECOY_)?sp\\|[^-]*\\|.*",
                                       `Protein accession`) ~ "SwissProt",
                                 TRUE ~ NA_character_))

  # adding cleanSeq
  ids <- ids %>%
    dplyr::mutate(
      cleanSeq = gsub("\\[.+?\\]", "", Proteoform),
      cleanSeq = gsub("\\(|\\)", "", cleanSeq),
      cleanSeq = sub("^[A-Z]?\\.(.*)\\.[A-Z]?", "\\1", cleanSeq))


  # If there is not a special amino acid change the NA to -. This is
  # necessary otherwise special amino acids get incorrectly copied in the
  # next step.
  ids <- ids %>%
    dplyr::mutate(
      `Special amino acids` = dplyr::case_when(
        is.na(`Special amino acids`) ~ "-",
        TRUE ~ as.character(`Special amino acids`)))


  # Extract just the UniProt accession. This is the second element (when
  # splitting by |) of the `Protein accession` variable output by TopPIC.
  ids <- ids %>%
    dplyr::mutate(
      UniProtAcc = case_when(grepl("\\|", UniProtAcc) ~
                               purrr::map_chr(stringr::str_split(UniProtAcc, "\\|"), purrr::pluck, 2),
                             TRUE ~ UniProtAcc))


  # Only keep variables we need throughout TopPICR.
  ids <- ids %>%
    dplyr::select(
      Dataset,
      `Scan(s)`,
      `Retention time`,
      `Feature apex`,
      Charge,
      mz,
      `Precursor mass`,
      `Adjusted precursor mass`,
      `Feature intensity`,
      `E-value`,
      `#unexpected modifications`,
      AccMap,
      MIScore,
      Proteoform,
      `Special amino acids`,
      `Protein accession`,
      `Protein description`,
      # `First residue`,
      # `Last residue`,
      firstAA,
      lastAA,
      `Proteoform mass`,
      `Proteoform-level Q-value`,
      UniProtAcc,
      Gene,
      AnnType,
      cleanSeq,
      isDecoy)

  # Fill in missing values for all accessions that match a given sequence.
  # This step is necessary because TopPIC only includes information for
  # the `Protein accession` and `Protein description` variables when a
  # proteoform matches multiple sequences. All other variables are left
  # blank.


  cols_to_fill <- c("Scan(s)", "Retention time", "Feature apex", "Charge",
                    "mz", "Precursor mass", "Adjusted precursor mass",
                    "Feature intensity", "E-value",
                    "#unexpected modifications",
                    "AccMap", "MIScore", "Proteoform", "Proteoform mass",
                    "Proteoform-level Q-value", "cleanSeq")

  ids <- ids %>%
    tidyr::fill(all_of(cols_to_fill), .direction = "down")

  # ids <- ids %>%
  #    tidyr::fill(`Scan(s)`:Proteoform, .direction = "down") %>%
  #    tidyr::fill(`Proteoform mass`, .direction = "down") %>%
  #    tidyr::fill(`Proteoform-level Q-value`, .direction = "down") %>%
  #    tidyr::fill(`cleanSeq`, .direction = "down")

  # Only keep the MS1 variables we use throughout TopPICR.
  feat <- bind_rows(toppic_output$`_ms1.feature`) %>%
    dplyr::select(Dataset, Mass, Intensity, Time_apex)

  return(list(ms2identifications = ids, ms1features = feat))
}


