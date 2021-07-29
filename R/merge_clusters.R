#' Merge clusters using at mass/retention time envelope
#'
#' ...
#'
#' @param x ...
#'
#' @param min_size ...
#'
#' @param ppm_cutoff ...
#'
#' @param ppm_range ...
#'
#' @param rt_range ...
#'
#' @return ...
#'
#' @importFrom magrittr %>%
#'
#' @export
#'
merge_clusters <- function (x, min_size, ppm_cutoff, ppm_range, rt_range) {

  x_cluster <- x[[1]] %>%
    dplyr::group_by(Gene) %>%
    dplyr::mutate(
      cluster_new = merge_mrt(x = .,
                              gene_name = unique(Gene),
                              rt_error = x[[3]],
                              ppm_cutoff = ppm_cutoff,
                              ppm_range = ppm_range,
                              rt_range = rt_range)
    )

  # Find indices of all rows with a cluster size smaller than 10. These rows
  # will be changed to 0 (the noise cluster).
  noise <- x_cluster %>%
    dplyr::group_by(Gene, cluster_new) %>%
    # Create a variable where TRUE if the number of points in a cluster is less
    # than 10 and FALSE if number of points is 10 or more.
    dplyr::mutate(noise = dplyr::n() < min_size) %>%
    dplyr::pull(noise)

  # Change the cluster membership for all clusters with fewer than 10 members to
  # 0 (the noise cluster).
  x_cluster[noise, "cluster_new"] <- 0

  return (x_cluster)

}

# x_53 auxiliary functions -----------------------------------------------------

# The correct cluster function determines if a point outside a cluster should be
# combined with the given cluster based on the distance in ppm between the
# cluster median and the corrected isotopic mass of the point. This function
# considers points within a mass/retention time envelope around a given cluster.
correction <- function (med_mass, min_rt, max_rt, clstr,
                        mass, rt, top_clstr, rt_error,
                        ppm_cutoff, ppm_range, rt_range) {

  cutoff_mass <- ppm_range * 1.0033548378
  cutoff_rt <- rt_range * rt_error

  # Find indices that fall within mass and rt bounds.
  idx_mass <- which(mass < med_mass + cutoff_mass &
                      mass > med_mass - cutoff_mass)
  idx_rt <- which(rt < max_rt + cutoff_rt &
                    rt > min_rt - cutoff_rt)

  # Find the indices that are in BOTH idx_mass and idx_rt
  idx1 <- dplyr::intersect(idx_mass, idx_rt)

  # If there are not any points that fall within this envelope return the
  # original cluster vector.
  if (length(idx1) == 0) {

    return (clstr)

  }

  # Subset the mass vector with the points that fall within the envelope.
  rm <- mass[idx1]

  # Calculate the corrected isotopic mass in ppm.
  im <- ((med_mass - rm - (1.0033548378 * round(med_mass - rm, 0))) /
           sapply(rm, function (x) mean(c(med_mass, x))) * 1e6)

  # Find indices whose absolute value of the corrected isotopic mass falls below
  # the ppm cutoff.
  idx2 <- which(abs(im) < ppm_cutoff)

  # Change the points that fall within the mass/rt envelope and below the ppm
  # threshold to the current top cluster. The first subset extracts the points
  # that are within the mass/rt envelope and the second subset extracts the
  # points that fall below the ppm threshold.
  clstr[idx1][idx2] <- top_clstr

  return (clstr)

}

# The merge_mrt function (merge clusters based on a Mass/Retention Time
# envelope around each cluster) changes cluster membership based on the
# following algorithm:
# 1. Group the x_cluster data frame by gene
# 2. Find the cluster (excluding 0) with the most elements
# 3. Compute the median/mean RecalMass for most abundant cluster
# 4. For each other point in the graph compute the difference betwen the
#    RecalMass of each point and the median/mean of the cluster.
# 5. Calculate the value from MSnID to correct for isotopic error.
# 6. If the value from 5 is smaller than 3 ppm convert the point to the cluster
#    being compared to otherwise leave the cluster as is
# 7. For efficiency, remove the points belonging to the most abundant cluster
# 8. Repeat steps 3-7 until all clusters have been iterated through
merge_mrt <- function (x, gene_name, rt_error, ppm_cutoff,
                       ppm_range, rt_range) {

  # Extract all rows corresponding to a given gene.
  x_gn_clstr <- x %>%
    dplyr::filter(Gene == gene_name) %>%
    dplyr::select(Gene, RecalMass, RTalign, cluster)

  # Compute the number of unique clusters for the current gene. This will be
  # used to determine if the algorithm to change cluster membership needs to be
  # run. For example, if there is only one cluster for a gene then nothing needs
  # to be done.
  n_clstrs <- length(unique(x_gn_clstr$cluster))

  # Check if there is only one cluster for the given gene. This is used to
  # determine if the change cluster membership algorithm needs to be run.
  if (n_clstrs == 1) {

    # Return the original cluster vector because nothing needs to be done to it.
    return (x_gn_clstr$cluster)

  }

  # Start the used_clusters vector at zero. This vector is used to keep track of
  # the clusters that have been considered by the change cluster algorithm. It
  # starts with 0 because points cannot be added to the 0 (noise) cluster based
  # on their mass or retention time. A point can only be changed to a noise
  # point (changed to the 0 cluster) if it belongs to a cluster with fewer than
  # 10 points.
  used_clusters <- c(0)

  # Use mass/rt envelope and ppm to determine membership ---------------

  # Loop through each cluster and correct isotopic error.
  while (TRUE) {

    # Determine most numerous cluster (excluding 0).
    top_cluster <- x_gn_clstr %>%
      dplyr::filter(!(cluster %in% used_clusters)) %>%
      dplyr::group_by(cluster) %>%
      dplyr::tally() %>%
      dplyr::arrange(-n) %>%
      dplyr::slice(1) %>%
      dplyr::pull(cluster)

    # Find median for mass and rt for the top cluster.
    med_mass <- x_gn_clstr %>%
      dplyr::filter(cluster == top_cluster) %>%
      dplyr::summarize(med_mass = stats::median(RecalMass, na.rm = TRUE)) %>%
      dplyr::pull(med_mass)
    min_rt <- x_gn_clstr %>%
      dplyr::filter(cluster == top_cluster) %>%
      dplyr::summarize(min_rt = min(RTalign, na.rm = TRUE)) %>%
      dplyr::pull(min_rt)
    max_rt <- x_gn_clstr %>%
      dplyr::filter(cluster == top_cluster) %>%
      dplyr::summarize(max_rt = max(RTalign, na.rm = TRUE)) %>%
      dplyr::pull(max_rt)

    used_clusters <- append(used_clusters, top_cluster)

    # Find indices of x_gn_clstr that don't correspond to past clusters (keep
    # 0).
    idx <- which(!(x_gn_clstr$cluster %in% used_clusters[-1]))

    # Check the length of idx. If it is greater than zero proceed with changing
    # cluster membership for points within the mass/rt envelope.
    if (length(idx) > 0) {

      # Change cluster membership according to the isotopic mass.
      x_gn_clstr[idx, ] <- x_gn_clstr %>%
        dplyr::filter(!(cluster %in% used_clusters[-1])) %>%
        dplyr::mutate(
          cluster = correction(med_mass = med_mass,
                               min_rt = min_rt,
                               max_rt = max_rt,
                               clstr = cluster,
                               mass = RecalMass,
                               rt = RTalign,
                               top_clstr = top_cluster,
                               rt_error = rt_error,
                               ppm_cutoff = ppm_cutoff,
                               ppm_range = ppm_range,
                               rt_range = rt_range)
        )

    }

    # Check if all clusters have been considered. It is possible that some
    # clusters that existed before have all been assigned to new clusters. This
    # if statement checks if the used_clusters vector is greater than or equal
    # the current number of unique clusters because it is possible that a gene
    # does not have any points classified as noise. If this is the case the
    # used_clusters vector could be larger than the unique clusters because this
    # vector always includes the 0 or noise cluster.
    if (length(used_clusters) >= length(unique(x_gn_clstr$cluster))) {

      break

    }

  }

  return (x_gn_clstr$cluster)

}
