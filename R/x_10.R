# For the first function x is a list of data frames instead of a single data 
# frame as the rest of the functions will be.
x_10 <- function (x = list()) {
  
  x <- data.table::rbindlist(x) %>%
    # Remove data.table class.
    # tibble() %>%
    # Add useful variables.
    filter(!grepl("Contaminant", `Protein accession`)) %>%
    filter(grepl("GN=", `Protein description`)) %>%
    mutate(ProjID = str_extract(Dataset, "\\d{8}")) %>%
    mutate(LC_Column = str_sub(Dataset, start = -8)) %>%
    mutate(mz = (`Precursor mass` + Charge * 1.007276466621) / Charge) %>%
    mutate(Gene = sub(".*GN=(\\S+).*","\\1",`Protein description`)) %>%
    mutate(isDecoy = grepl("^XXX", `Protein accession`)) %>%
    mutate(RTmin = `Retention time` / 60) %>%
    # Remove not so useful variables. In other words, remove variables that will
    # never be used or thought of again.
    select(-c(`Data file name`, `Prsm ID`, `Spectrum ID`, Fragmentation,
              `Retention time`, `#peaks`, `Proteoform ID`, `Feature score`,
              MIScore, `#variable PTMs`, `#matched peaks`,
              `#matched fragment ions`, `Q-value (spectral FDR)`,
              `Proteoform FDR`))
  
  return (x)
  
}
