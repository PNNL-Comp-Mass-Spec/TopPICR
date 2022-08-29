# The main function will create a new column pcg for Proteoform Cluster Group
# (similar to cluster_new from the merge_clusters function). This function will
# start with the most abundant cluster's centroid and compare all remaining
# cluster centroids to it and group them if the centroids fall within the
# mass/rt envelope of the main cluster. If clusters are grouped with the main
# cluster their centroids will be removed from the list and will not be
# considered in creating their own group with other clusters.

#' Combines clusters based on mass and retention time
#'
#' Compares the centroid of each cluster within a gene and combines the two
#' clusters if the centroids fall within the specified threshold. The clusters
#' are combined by creating a new cluster membership variable \code{pcGroup}.
#'
#' @param x A \code{data.table} output from the \code{cluster} function.
#'
#' @param errors A \code{list} output from the \code{calc_error} function. The
#'   first element of the list contains the standard deviation and median of the
#'   mass measurement error for each data set. The second element is the
#'   standard deviation of the mass measurement error across all data sets. The
#'   third element is the standard deviation of the retention time in seconds
#'   across all data sets.
#'
#' @param n_mme_sd The number of standard deviations used to create an envelope
#'   around a cluster's centroid in mass space. The standard deviation is in
#'   ppm. This threshold will be used to determine if two clusters will be
#'   combined.
#'
#' @param n_Da The number of Daltons used to create an envelope around a
#'   cluster's centroid in mass space. Any cluster whose centroid falls within
#'   this envelope will be included in the pool of clusters that could be
#'   combined. This value is used to account for isotopic error and is separate
#'   from \code{n_mme_sd}.
#'
#' @param n_rt_sd The number of standard deviations used to create an envelope
#'   around a cluster's centroid in retention time space. Any cluster whose
#'   centroid falls within this envelope will be included in the pool of
#'   clusters that could be combined.
#'
#' @return A \code{data.table} with the updated cluster assignment for each
#'   gene/cluster combination.
#'
#' @md
#'
#' @importFrom magrittr %>%
#'
#' @author Evan A Martin
#'
#' @export
#'
create_pcg <- function (x, errors, n_mme_sd, n_Da, n_rt_sd) {

  # Convert n_mme_sd to ppm. This is necessary so the remainder of the
  # function does not have to be altered to reflect the change from using a ppm
  # cutoff to a cutoff in standard deviations.
  ppm_cutoff <- n_mme_sd * errors$ppm_sd

  x_group <- x %>%
    dplyr::filter(cluster != 0) %>%
    dplyr::group_by(Gene) %>%
    dplyr::mutate(
      pcGroup = groupate(x = dplyr::across(),
                         ppm_cutoff = ppm_cutoff,
                         n_Da = n_Da,
                         rt_sd = errors$rt_sd,
                         n_rt_sd = n_rt_sd)
    ) %>%
    dplyr::ungroup()

}

# @author Evan A Martin
# This function will carry out the bulk of the calculations (sorting, grouping,
# creating group numbers, ...)
groupate <- function (x, ppm_cutoff, n_Da, rt_sd, n_rt_sd) {

  c12c13 <- 1.0033548378

  # Naming key/explanation:
  #
  # ridley - The data frame containing each cluster centroid. This data frame
  # will be looped through until all clusters have been grouped. It is possible
  # that a cluster will not be grouped with any others.
  #
  # samus - A data frame that will contain the proteoform/cluster groups. The
  # proteoform/cluster group variable in this data frame will be returned at the
  # end of the function after being joined with the input data frame to ensure
  # it is the proper length.
  #
  # metroid - A temporary single row data frame. This data frame contains the
  # centroid of the current cluster. All remaining centroids in the ridley data
  # frame will be compared to the centroid in metroid. If any centroid in ridley
  # falls within the mas/rt envelope and the ppm threshold it will be grouped
  # with the centroid in metroid.
  #
  # metroid_group - Another temporary data frame that contains the rows of
  # ridley that fall within the mass/rt envelope and below the ppm threshold.
  # These rows represent clusters that can be combined with the current cluster
  # (in metroid).

  # Determine the size of each cluster. This will be used to combine clusters
  # (smaller clusters are added to larger clusters). Calculate the mass/rt
  # centroid.
  ridley <- x %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarize(n_members = dplyr::n(),
                     cntr_mass = median(RecalMass, na.rm = TRUE),
                     cntr_rt = median(RTalign, na.rm = TRUE)) %>%
    dplyr::arrange(-n_members) %>%
    dplyr::ungroup()

  # Check if there is only one cluster for the given gene. If there is only one
  # cluster the rest of the algorithm does not need to be run.
  if (dim(ridley)[[1]] == 1) {

    # Return the original cluster vector because nothing needs to be done to it.
    return (ridley$cluster)

  }

  # Create an empty tibble that will be added two as the rows of ridley are
  # looped through.
  samus <- tibble::tibble(NULL)

  # Loop through each cluster from largest to smallest and group clusters based
  # on the mass/retention time criteria and ppm cutoff.
  while (TRUE) {

    # Nab the first row of ridley. The mass/rt centroid of this temporary data
    # frame will be used to compare with other cluster centroids.
    metroid <- ridley[1, ]

    # Remove the first row from ridley. This data frame will be used to compare
    # centroids to the current centroid (found in metroid).
    ridley <- ridley[-1, ]

    # Check if there are zero rows in ridley. If there are then we have made it
    # through each of the clusters.
    if (dim(ridley)[[1]] == 0) {

      # Set the proteoform cluster group to the cluster number of the last
      # cluster.
      samus <- metroid %>%
        dplyr::mutate(pcg = cluster) %>%
        dplyr::bind_rows(samus, .)

      # Exit the while loop because all clusters have been combined into groups.
      break

    }

    # Only keep the centroids that fall within the mass/rt envelope and below
    # the ppm threshold.
    metroid_group <- ridley %>%
      dplyr::mutate(diff_mass = abs(metroid$cntr_mass - cntr_mass),
                    diff_rt = abs(metroid$cntr_rt - cntr_rt)) %>%
      dplyr::filter(diff_mass < n_Da * c12c13,
                    diff_rt < n_rt_sd * rt_sd) %>%
      dplyr::mutate(
        diff_ppm = abs(
          (metroid$cntr_mass - cntr_mass -
             (c12c13 * round(metroid$cntr_mass - cntr_mass,
                             0))) / metroid$cntr_mass * 1e6
        )
      ) %>%
      dplyr::filter(diff_ppm < ppm_cutoff)

    # Check if there are any clusters that meet the criteria for being grouped
    # with the current cluster. If there aren't then the current cluster cannot
    # be grouped with the others and the pcg will be set to the cluster number.
    if (dim(metroid_group)[[1]] == 0) {

      # Add the current cluster to the data frame that will be the final output.
      samus <- metroid %>%
        dplyr::mutate(pcg = cluster) %>%
        dplyr::bind_rows(samus, .)

      # Runs when there is one or more rows in metroid_group.
    } else {

      # Combine clusters in a group and add these rows to samus.
      samus <- metroid_group %>%
        dplyr::select(-diff_mass, -diff_rt, -diff_ppm) %>%
        # Add the current cluster data frame (metroid) to the data frame
        # (metroid_group) of the other clusters that will form a group with the
        # current cluster.
        dplyr::bind_rows(metroid, .) %>%
        dplyr::mutate(pcg = metroid$cluster) %>%
        dplyr::bind_rows(samus, .)

      # Remove the grouped clusters from the ridley data frame so they are not
      # considered again.
      ridley <- ridley[-(which(ridley$cluster %in% metroid_group$cluster)), ]

    }

  } # end while loop

  # Join the data input data frame with samus so the output will be the same
  # length as the input.
  suppressMessages(
    da_end <- dplyr::inner_join(x = x,
                                y = dplyr::select(samus, cluster, pcg))

  )

  return (da_end$pcg)

}
