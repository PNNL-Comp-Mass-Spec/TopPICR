x_47 <- function (x, x_align, ref_ds) {
  
  # Compute the robust median of the precursor mass within Gene, Dataset, and
  # Proteoform. This variable will be used to calculate the error (in ppm)
  # between the reference and all other data sets.
  x <- x %>%
    select(c(Dataset, Proteoform, Gene,
             `Precursor mass`, `Adjusted precursor mass`)) %>%
    group_by(Dataset, Proteoform, Gene) %>%
    summarize(Mass = robust_median(mass = `Precursor mass`,
                                   adj_mass = `Adjusted precursor mass`))
  
  # Combine the data frames with the robust median mass and the aligned
  # retention times.
  x_align <- x_align %>%
    ungroup() %>%
    inner_join(x)
  
  # Extract the reference data set.
  x_ref <- x_align %>%
    filter(Dataset == ref_ds) %>%
    select(-Dataset, -Gene)
  
  # Calculate the ppm and retention time errors.
  x_error <- x_align %>%
    group_by(Dataset) %>%
    summarize(calc_error(mass = Mass,
                         rt = RTalign,
                         pf = Proteoform,
                         x_ref = x_ref)) %>%
    ungroup()
  
  # Recalibrate the mass.
  x_recal <- inner_join(x = x_align,
                        y = x_error,
                        by = "Dataset") %>%
    mutate(RecalMass = Mass * (1 - ppm_median / 1e6))
  
  # Compute the mass and retention time errors.
  errors <- x_error %>%
    filter(Dataset != ref_ds) %>%
    summarize(ppm_error = median(ppm_mad),
              rt_error = median(rt_mad))
  
  # Return the x_recal data frame that contains the RecalMass variable along
  # with the ppm and retention time errors.
  return (list(x = x_recal,
               ppm_error = errors$ppm_error[[1]],
               rt_error = errors$rt_error[[1]]))
  
}

# x_47 auxiliary functions -----------------------------------------------------

# mass: A vector of masses extracted from the data matrix by the group_by
# function.
# rt: A vector of retention times extracted from the data matrix by the group_by
# function.
# pf: A vector of proteoforms extracted from the data matrix by the group_by
# function.
# x_ref: the x data matrix filtered by the reference data set.
calc_error <- function (mass, rt, pf, x_ref) {
  
  # Combine the reference data set along with the data subsetted by the group_by
  # function by Gene and Proteoform.
  error <- inner_join(x = x_ref[, c("Mass", "RTalign", "Proteoform")],
                      y = tibble(Mass = mass,
                                 RTalign = rt,
                                 Proteoform = pf),
                      by = "Proteoform") %>%
    # Calculate the ppm and retention time errors.
    mutate(ppm_error = 1e6 * (Mass.y - Mass.x) / Mass.x,
           rt_error = RTalign.y - RTalign.x)
  
  # Only keep ppm values between plus/minus 20 for computing the median. The
  # ppm median is used for recalibrating the mass later. We only keep the values
  # within 20 ppm for a more accurate representation of the group.
  reduced <- error %>%
    filter(ppm_error > -20 & ppm_error < 20)
  
  # Calculate the median absolute deviation, median and standard deviation of
  # the remaining ppm values.
  mms <- data.frame(ppm_mad = mad(error$ppm_error),
                    ppm_median = median(reduced$ppm_error),
                    rt_mad = mad(error$rt_error))
  
  return (mms)
  
}

# Calculates the median Precursor mass of the mass values that are within 0.5
# Dalton of the Adjusted precursor mass. There are cases when all of the
# Precursor masses are outside this 0.5 Dalton range of the Adjusted precursor
# mass. If this is the case the Precursor mass that is closest to the Adjusted
# precursor mass is returned.
robust_median <- function (mass, adj_mass) {
  
  # Calculate the absolute difference between the Precursor mass (mass) and the
  # Adjusted precursor mass (adj_mass). This vector will be used to determine
  # which Precursor masses will be used to calculate the representative mass for
  # the Dataset, Proteoform, Gene combination.
  diffs <- abs(adj_mass - mass)
  
  # Only keep mass (Precursor mass) values that are within 0.5 Dalton of
  # adj_mass (Adjusted precursor mass).
  goodies <- which(diffs < 0.5)
  
  # Check the length of goodies. It is possible that all of the differences are
  # outside the 0.5 threshold. If this is the case we will use the mass that is
  # the closest to the adj_mass.
  if (length(goodies) == 0) {
    
    # Plus n isotope case.
    idx <- mass - adj_mass > 0.75
    
    # Subtract carbon12/13 difference from mass.
    mass[idx] <- mass[idx] - 1.003355
    
    idx <- mass - adj_mass < -0.75
    
    # Will add an amino group to the mass.
    mass[idx] <- mass[idx] + 0.984016
    
    # Runs if there is at least one element in goodies.
  } else {
    
    # Only keep the masses are within the 0.5 Dalton threshold.
    mass <- mass[goodies]
    
  }
  
  # Take the median of the mass values that are within 0.5 Dalton of the
  # adj_mass to use as the representative mass.
  return (median(mass, na.rm = TRUE))
  
  
}
