#' Create feature metadata
#'
#' This function uses the information from the identified feature data to create
#' a feature metadata data frame. This data frame will be used at the end of the
#' `TopPICR` pipeline when creating an `MSnSet` object.
#'
#' @param x A \code{data.table} output from the \code{create_pcg} function.
#'
#' @param errors A \code{list} output from the \code{calc_error} function. The
#'   first element of the list contains the standard deviation and median of the
#'   mass measurement error for each data set. The second element is the
#'   standard deviation of the mass measurement error across all data sets. The
#'   third element is the standard deviation of the retention time in seconds
#'   across all data sets.
#'
#' @param n_mme_sd A numeric value indicating the number of standard deviations
#'   in ppm. This value is used as the threshold for determining if two clusters
#'   (from two different genes) have the same centroid in mass space.
#'
#' @param n_rt_sd A numeric value representing the number of standard
#'   deviations in seconds. This value is used as the threshold for determining
#'   if two clusters (from two different genes) have the same centroid in
#'   retention time space.
#'
#' @return A \code{data.table} with a row for each unique combination of `Gene`,
#'   `pcGroup`, and `Proteoform`. The `collision` variable indicates which
#'   clusters have the same centroid across different genes. The
#'   \code{data.table} also contains other proteoform information such as
#'   protein length, UniProt accession, and first and last amino acid among
#'   others.
#'
#' @importFrom magrittr %>%
#'
#' @author Evan A Martin
#'
#' @export
#'
create_mdata <- function (x, errors, n_mme_sd, n_rt_sd) {

  # Convert n_mme_sd to ppm. This is necessary so the remainder of the
  # function does not have to be altered to reflect the change from using a ppm
  # cutoff to a cutoff in standard deviations.
  ppm_cutoff <- n_mme_sd * errors$ppm_sd

  # Convert n_rt_sd to seconds. This is necessary so the remainder of the
  # function does not have to be altered to reflect the change from using a
  # retention time cutoff in seconds to a cutoff in standard deviations.
  rt_cutoff <- n_rt_sd * errors$rt_sd

  # First create the collision variable.
  x_collision <- find_collider(x = x,
                               ppm_cutoff = ppm_cutoff,
                               rt_cutoff = rt_cutoff) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      # First remove anything in square brackets that starts with a letter.
      # The string in the square brackets may contain numbers, special
      # characters, and letters (in that order) after the initial letter. For
      # example, this will remove [Acetyl], [Plus1Oxy], and [NH3_Loss].
      ModMass = stringr::str_remove_all(
        Proteoform, "\\[[a-zA-Z]+[0-9]*.*?[a-zA-Z]*\\]"
      ) %>%
        # Remove all upper and lower case letters.
        stringr::str_remove_all("[a-zA-Z]") %>%
        # Remove beginning and ending periods.
        stringr::str_remove_all("^\\.|\\.$") %>%
        # Remove all parentheses and square brackets.
        stringr::str_remove_all("\\(|\\)|\\[|\\]"),
      ModMass = suppressWarnings(
        round(as.numeric(ModMass), digits = 0)
      )
    ) %>%
    dplyr::group_by(Gene, pcGroup, cleanSeq, ModMass) %>%
    dplyr:: mutate(
      spectralCount = dplyr::n(),
      rt = stats::median(RTalign, na.rm = TRUE),
      mass = stats::median(RecalMass, na.rm = TRUE)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Gene, pcGroup) %>%
    dplyr::slice_max(spectralCount) %>%
    dplyr::slice_min(`E-value`) %>%
    dplyr::ungroup() %>%
    dplyr::select(Gene, pcGroup, collision, UniProtAcc, mass, rt,
                  firstAA, lastAA, protLength, Proteoform, AnnType,
                  `#unexpected modifications`)

  # Return the metadata data frame.
  return (
    x_collision
  )

}

# Finds clusters with the same centroid (or are within the specified
# mass/retention time threshold) but the genes of the two clusters are different
# and marks them as colliding genes.
# @author Evan A Martin
find_collider <- function (x, ppm_cutoff, rt_cutoff) {

  # Keep only unique Gene/cluster combinations (excluding the 0 cluster).
  # Calculate the centroid of each cluster. The centroid is the median mass and
  # median retention time of all points in the cluster.
  unique_gccs <- x %>%
    dplyr::filter(cluster != 0) %>%
    dplyr::group_by(Gene, cluster) %>%
    dplyr::mutate(
      gcc = paste(Gene, cluster, sep = "_"),
      collision = "-",
      cntr_mass = stats::median(RecalMass, na.rm = TRUE),
      cntr_rt = stats::median(RTalign, na.rm = TRUE)
    ) %>%
    dplyr::distinct(Gene, cluster, collision, pcGroup,
                    gcc, cntr_mass, cntr_rt) %>%
    dplyr::ungroup()

  # Nab the number of unique gccs. This will be used to keep track of which gccs
  # need to be combined into one gcc.
  n_gccs <- nrow(unique_gccs)

  # Loop through each gcc and compare with all other gccs to determine if they
  # share the same centroid (within the tolerance given).
  for (e in 1:n_gccs) {

    for (v in 1:n_gccs) {

      # Only perform computations on the upper triangular matrix.
      if (e < v) {

        # Check if the retention time falls within the threshold. If it does
        # proceed with calculating the ppm difference between the two gccs.
        if (abs(unique_gccs$cntr_rt[[e]] -
                unique_gccs$cntr_rt[[v]]) < rt_cutoff) {

          # Compute the ppm difference between the two gccs.
          diff_ppm <- abs(
            (unique_gccs$cntr_mass[[e]] - unique_gccs$cntr_mass[[v]]) /
              unique_gccs$cntr_mass[[e]] * 1e6
          )

          # Check if the vth gcc falls within the threshold of the eth gcc and
          # if the genes for the two gccs are different. If the genes are
          # different the two gccs will be combined.
          if (diff_ppm < ppm_cutoff &&
              unique_gccs$Gene[[e]] != unique_gccs$Gene[[v]]) {

            # Mark each gene with the gene it collides with. For example, if
            # APP_3 and MBP_5 have the same centroid the collision column for
            # rows containing APP_3 will be marked with "MBP" and the rows
            # containing MBP_5 will will be marked with "APP".
            unique_gccs$collision[[e]] <- paste(
              unique_gccs$Gene[[v]],
              unique_gccs$pcGroup[[v]],
              sep = "_"
            )
            unique_gccs$collision[[v]] <- paste(
              unique_gccs$Gene[[e]],
              unique_gccs$pcGroup[[e]],
              sep = "_"
            )

          } # End check if ppm falls within the threshold.

        } # End check if retention time falls within the threshold.

      } # End check if e is less than v.

    } # End loop across columns (v)

  } # End loop down rows (e)

  return (dplyr::inner_join(x,
                            dplyr::select(unique_gccs,
                                          Gene,
                                          cluster,
                                          pcGroup,
                                          collision)))

}
