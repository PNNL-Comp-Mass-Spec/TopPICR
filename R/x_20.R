prime <- function (cutoff, threshold, x, type) {
  
  # Create a vector of ordered (smallest to largest) p-values. This vector will
  # be used to search for the p-value cutoff that is the closest to the
  # threshold value.
  pvals <- x$`P-value`[order(x$`P-value`)]
  
  fdr <- x %>%
    filter(AnnType == type) %>%
    mutate(lpval = log10(`P-value`)) %>%
    filter(lpval < log_cutoff) %>%
    distinct(Gene, isDecoy) %>%
    pull(isDecoy) %>%
    table() %>%
    prop.table()
  
  # Check if only one value is returned by prop.table. This will occur if there
  # are no counts for either TRUE or FALSE.
  if (length(dimnames(fdr)$.) == 1) {
    
    # Check if the proportion returned is TRUE. If it is this is the value that
    # needs to be returned. If it FALSE then return 0.
    if (dimnames(fdr)$.[1]) {
      
      # Return the first element returned by prop.table because the only counts
      # are TRUE.
      # return (if (fdr[[1]] < threshold) fdr[[1]] else fdr[[1]] + 1e6)
      return (fdr[[1]])
      
      # Runs when only a value for FALSE is returned by prop.table.
    } else {
      
      # Return 0 because all counts are FALSE.
      return (0)
      
    }
    
    # Runs if there is a proportion for both TRUE and FALSE.
  } else {
    
    # Return the proportion of counts for TRUE.
    # return (if (fdr[[2]] < threshold) fdr[[2]] else fdr[[2]] + 1e6)
    return (fdr[[2]])
    
  }
  
}

optimize(prime, interval = c(0, 1),
         x = x, type = "TrEMBL", threshold = 0.01)

optim(par = 1e-10, fn = prime, x = x, type = "TrEMBL")
