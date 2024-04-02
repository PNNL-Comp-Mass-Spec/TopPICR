#' A non-exported function. This function determines the version of the TopPIC
#' given the column names in the identification and features files.
#'
#'

determine_toppic_version <- function (x){


   if(all(colnames(x$`_TopPIC_PrSMs.txt`) == id_cols_1.5.4) &
      all(colnames(x$`_ms1.feature`) == features_cols_1.5.4)){
      return("1.5.4")
   }

   if(all(colnames(x$`_TopPIC_PrSMs.txt`) == id_cols_1.7.0) &
      all(colnames(x$`_ms1.feature`) == features_cols_1.7.0)){
      return("1.7.0")
   }

   warning("TopPIC version is unknown!")

   return("unknown")
}




id_cols_1.5.4 <- c("Dataset", "Data file name", "Prsm ID", "Spectrum ID",
                   "Fragmentation", "Scan(s)", "Retention time", "#peaks",
                   "Charge", "Precursor mass", "Adjusted precursor mass",
                   "Proteoform ID", "Feature intensity", "Feature score",
                   "Feature apex time", "#Protein hits", "Protein accession",
                   "Protein description", "First residue", "Last residue",
                   "Special amino acids", "Proteoform", "Proteoform mass",
                   "Protein N-terminal form", "#unexpected modifications",
                   "#variable PTMs", "MIScore", "#matched peaks",
                   "#matched fragment ions", "E-value",
                   "Spectrum-level Q-value", "Proteoform-level Q-value")

id_cols_1.7.0 <- c("Dataset", "Data file name", "Prsm ID", "Spectrum ID",
                   "Fragmentation", "Scan(s)", "Retention time", "#peaks",
                   "Charge", "Precursor mass", "Adjusted precursor mass",
                   "Proteoform ID", "Feature ID", "Feature intensity",
                   "Feature score", "Feature apex time", "#Protein hits",
                   "Protein accession", "Protein description", "First residue",
                   "Last residue", "Special amino acids",
                   "Database protein sequence", "Proteoform", "Proteoform mass",
                   "Protein N-terminal form", "Fixed PTMs",
                   "#unexpected modifications", "unexpected modifications",
                   "#variable PTMs", "variable PTMs", "MIScore",
                   "#matched peaks", "#matched fragment ions", "E-value",
                   "Spectrum-level Q-value", "Proteoform-level Q-value")

features_cols_1.5.4 <- c("Dataset", "Sample_ID", "ID", "Mass", "Intensity",
                         "Time_begin", "Time_end", "Time_apex",
                         "Minimum_charge_state", "Maximum_charge_state",
                         "Minimum_fraction_id", "Maximum_fraction_id")

features_cols_1.7.0 <- c("Dataset", "Sample_ID", "Feature_ID",
                         "Monoisotopic_mass", "Intensity", "Min_time",
                         "Max_time", "Min_scan", "Max_scan", "Min_charge",
                         "Max_charge", "Apex_time", "Apex_scan",
                         "Apex_intensity", "Rep_charge", "Rep_average_mz",
                         "Envelope_number", "EC_score", "Min_fraction_ID",
                         "Max_fraction_ID", "Elution_length")

