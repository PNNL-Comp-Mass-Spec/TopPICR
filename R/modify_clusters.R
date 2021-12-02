# Create gene/cluster combinations. These could contain clusters from two
# different genes.

#' Combine clusters across genes
#'
#' Combine two clusters if they have the same centroid (or are within the
#' specified mass/retention time threshold) but the genes of the two clusters
#' are different.
#'
#' @param x ...
#'
#' @param cutoff_ppm ...
#'
#' @param cutoff_rt ...
#'
#' @return ...
#'
#' @importFrom magrittr %>%
#'
#' @author Evan A Martin
#'
#' @export
#'
create_gcc <- function (x, cutoff_ppm, cutoff_rt) {

  # Keep only unique Gene/cluster combinations (excluding the 0 cluster).
  # Calculate the centroid of each cluster. The centroid is the median mass and
  # median retention time of all points in the cluster.
  unique_gccs <- x %>%
    filter(cluster != 0) %>%
    group_by(Gene, cluster) %>%
    mutate(gcc = paste(Gene, cluster, sep = "_"),
           cntr_mass = median(RecalMass, na.rm = TRUE),
           cntr_rt = median(RTalign, na.rm = TRUE)) %>%
    distinct(Gene, cluster, gcc, cntr_mass, cntr_rt)

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
                unique_gccs$cntr_rt[[v]]) < cutoff_rt) {

          # Compute the ppm difference between the two gccs.
          diff_ppm <- abs(
            (unique_gccs$cntr_mass[[e]] - unique_gccs$cntr_mass[[v]]) /
              unique_gccs$cntr_mass[[e]] * 1e6
          )

          # Check if the vth gcc falls within the threshold of the eth gcc and
          # if the genes for the two gccs are different. If the genes are
          # different the two gccs will be combined.
          if (diff_ppm < cutoff_ppm &&
              unique_gccs$Gene[[e]] != unique_gccs$Gene[[v]]) {

            # Combine the two gccs with the same centroid into one gcc. This new
            # gcc will span two different genes.
            unique_gccs$gcc[[e]] <- paste(
              paste(unique_gccs$Gene[[e]],
                    unique_gccs$cluster[[e]],
                    sep = "_"),
              paste(unique_gccs$Gene[[v]],
                    unique_gccs$cluster[[v]],
                    sep = "_"),
              sep = ";"
            )
            unique_gccs$gcc[[v]] <- paste(
              paste(unique_gccs$Gene[[e]],
                    unique_gccs$cluster[[e]],
                    sep = "_"),
              paste(unique_gccs$Gene[[v]],
                    unique_gccs$cluster[[v]],
                    sep = "_"),
              sep = ";"
            )


          }

        }

      }

    }

  }

  return (dplyr::inner_join(x,
                            dplyr::select(unique_gccs,
                                          Gene,
                                          cluster,
                                          gcc)))

}

# The main function will create a new column pcg for Proteoform Cluster Group
# (similar to cluster_new from the merge_clusters function). This function will
# start with the most abundant cluster's centroid and compare all remaining
# cluster centroids to it and group them if the centroids fall within the
# mass/rt envelope of the main cluster. If clusters are grouped with the main
# cluster their centroids will be removed from the list and will not be
# considered in creating their own group with other clusters.

#' ...
#'
#' ...
#'
#' @param x ...
#'
#' @param errors ...
#'
#' @param ppm_cutoff ...
#'
#' @param n_Da ...
#'
#' @param n_rt_sd ...
#'
#' @importFrom magrittr %>%
#'
#' @author Evan A Martin
#'
#' @export
#'
gen_pcg <- function (x, errors, ppm_cutoff, n_Da, n_rt_sd) {

  x_group <- x %>%
    dplyr::filter(cluster != 0) %>%
    dplyr::group_by(Gene) %>%
    dplyr::mutate(
      pcGroup = groupate(x = .,
                         gene = Gene[[1]],
                         ppm_cutoff = ppm_cutoff,
                         n_Da = n_Da,
                         rt_sd = errors$rt_sd,
                         n_rt_sd = n_rt_sd)
    )

}

# @author Evan A Martin
# This function will carry out the bulk of the calculations (sorting, grouping,
# creating group numbers, ...)
groupate <- function (x, gene, ppm_cutoff, n_Da, rt_sd, n_rt_sd) {

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

  # Extract all rows corresponding to a given gene.
  ridley <- x %>%
    dplyr::filter(Gene == gene) %>%
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
    da_end <- x %>%
      dplyr::filter(Gene == gene) %>%
      dplyr::inner_join(x = .,
                        y = dplyr::select(samus, cluster, pcg))
  )

  return (da_end$pcg)

}
