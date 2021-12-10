#' Merge clusters using a mass/retention time envelope
#'
#' Merges points with clusters based on their mass and retention time. Starting
#' with the largest cluster, each point is compared to a mass/retention time
#' envelope around the centroid of the cluster. If a point falls within this
#' envelope the point is merged with the cluster. This process is repeated for
#' each cluster (except the noise cluster).
#'
#' @param x A \code{data.table} output from the \code{cluster} function.
#'
#' @param errors A \code{list} output from the \code{calc_error} function.
#'
#' @param min_size An integer value indicating the minimum number of points a
#'   cluster must have. All clusters with fewer members than min_size will be
#'   reclassified as "noise" points.
#'
#' @param ppm_cutoff The threshold used to determine if a point is close enough
#'   in mass to be merged with a new cluster.
#'
#' @param n_Da The number of Daltons (multiplied by the difference in mass
#'   between Carbon 12 and Carbon 13) to use when creating the mass portion of
#'   the mass/retention time envelope.
#'
#' @param n_rt_sd The number of standard deviations to use when creating the
#'   retention time portion of the mass/retention time envelope.
#'
#' @return A \code{data.table} with like clusters combined if the points in the
#'   clusters fall within the mass/retention time envelope of the larger
#'   cluster. The following variables have been added/removed:
#'
#'   | Added             | Removed                    |
#'   | ----------------- | -------------------------- |
#'   | `cluster_new`     |                            |
#'
#' @md
#'
#' @importFrom magrittr %>%
#'
#' @author Evan A Martin
#'
#' @export
#'
merge_clusters <- function (x, errors, min_size, ppm_cutoff, n_Da, n_rt_sd) {

  x_cluster <- x %>%
    dplyr::group_by(Gene) %>%
    dplyr::mutate(
      cluster_new = merge_mrt(x = dplyr::across(),
                              rt_error = errors$rt_sd,
                              ppm_cutoff = ppm_cutoff,
                              n_Da = n_Da,
                              n_rt_sd = n_rt_sd)
    )

  # Find indices of all rows with a cluster size smaller than min_size. These
  # rows will be changed to 0 (the noise cluster).
  noise <- x_cluster %>%
    dplyr::group_by(Gene, cluster_new) %>%
    # Create a variable where TRUE if the number of points in a cluster is less
    # than min_size and FALSE if number of points is min_size or more.
    dplyr::mutate(noise = dplyr::n() < min_size) %>%
    dplyr::pull(noise)

  # Change the cluster membership for all clusters with fewer than min_size
  # members to 0 (the noise cluster).
  x_cluster[noise, "cluster_new"] <- 0

  return (x_cluster)

}

# merge_clusters auxiliary functions -------------------------------------------

# @author Evan A Martin
# The merge_mrt function (merge clusters based on a Mass/Retention Time
# envelope around each cluster) changes cluster membership based on the
# following algorithm:
# 1. Group the x_cluster data frame by gene
# 2. Find the cluster (excluding 0) with the most elements
# 3. Compute the median/mean RecalMass or repMass for most abundant cluster
# 4. For each other point in the graph compute the difference between the
#    RecalMass/repMass of each point and the median/mean of the cluster.
# 5. Calculate the value from MSnID to correct for isotopic error.
# 6. If the value from 5 is smaller than 3 ppm convert the point to the cluster
#    being compared to otherwise leave the cluster as is
# 7. For efficiency, remove the points belonging to the most abundant cluster
# 8. Repeat steps 3-7 until all clusters have been iterated through
merge_mrt <- function (x, rt_error, ppm_cutoff,
                       n_Da, n_rt_sd) {

  # Compute the number of unique clusters for the current gene. This will be
  # used to determine if the algorithm to change cluster membership needs to be
  # run. For example, if there is only one cluster for a gene then nothing needs
  # to be done.
  n_clstrs <- dplyr::n_distinct(x$cluster)

  # Check if there is only one cluster for the given gene. This is used to
  # determine if the change cluster membership algorithm needs to be run.
  if (n_clstrs == 1) {

    # Return the original cluster vector because nothing needs to be done to it.
    return (x$cluster)

  }

  # Start the used_clusters vector at zero. This vector is used to keep track of
  # the clusters that have been considered by the change cluster algorithm. It
  # starts with 0 because points cannot be added to the 0 (noise) cluster based
  # on their mass or retention time. A point can only be changed to a noise
  # point (changed to the 0 cluster) if it belongs to a cluster with fewer than
  # min_size points.
  used_clusters <- c(0)

  # Use mass/rt envelope and ppm to determine membership ---------------

  # Loop through each cluster and correct isotopic error.
  while (TRUE) {

    # Determine most numerous cluster (excluding 0).
    top_cluster <- x %>%
      dplyr::filter(!(cluster %in% used_clusters)) %>%
      dplyr::group_by(cluster) %>%
      dplyr::tally() %>%
      dplyr::arrange(-n) %>%
      dplyr::slice(1) %>%
      dplyr::pull(cluster)

    # Find median mass and min/max rt for the top cluster.
    moss_piglet <- x %>%
      dplyr::filter(cluster == top_cluster) %>%
      dplyr::summarize(
        med_mass = stats::median(RecalMass, na.rm = TRUE),
        min_rt = min(RTalign, na.rm = TRUE),
        max_rt = max(RTalign, na.rm = TRUE)
      )

    used_clusters <- append(used_clusters, top_cluster)

    # Find indices of x that don't correspond to past clusters (keep
    # 0).
    idx <- which(!(x$cluster %in% used_clusters[-1]))

    # Check the length of idx. If it is greater than zero proceed with changing
    # cluster membership for points within the mass/rt envelope.
    if (length(idx) > 0) {

      # Change cluster membership according to the isotopic mass.
      x[idx, ] <- x %>%
        dplyr::filter(!(cluster %in% used_clusters[-1])) %>%
        dplyr::mutate(
          cluster = correction(med_mass = dplyr::pull(moss_piglet, med_mass),
                               min_rt = dplyr::pull(moss_piglet, min_rt),
                               max_rt = dplyr::pull(moss_piglet, max_rt),
                               clstr = cluster,
                               mass = RecalMass,
                               rt = RTalign,
                               top_clstr = top_cluster,
                               rt_error = rt_error,
                               ppm_cutoff = ppm_cutoff,
                               n_Da = n_Da,
                               n_rt_sd = n_rt_sd)
        )

    }

    # Check if all clusters have been considered. It is possible that some
    # clusters that existed before have all been assigned to new clusters. This
    # if statement checks if the used_clusters vector is greater than or equal
    # the current number of unique clusters because it is possible that a gene
    # does not have any points classified as noise. If this is the case the
    # used_clusters vector could be larger than the unique clusters because this
    # vector always includes the 0 or noise cluster.
    if (length(used_clusters) >= dplyr::n_distinct(x$cluster)) {

      break

    }

  }

  return (x$cluster)

}

# @author Evan A Martin
# The correct cluster function determines if a point outside a cluster should be
# combined with the given cluster based on the distance in ppm between the
# cluster median and the corrected isotopic mass of the point. This function
# considers points within a mass/retention time envelope around a given cluster.
correction <- function (med_mass, min_rt, max_rt, clstr,
                        mass, rt, top_clstr, rt_error,
                        ppm_cutoff, n_Da, n_rt_sd) {

  cutoff_mass <- n_Da * 1.0033548378
  cutoff_rt <- n_rt_sd * rt_error

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

#' Merge clusters using the overlap coefficient
#'
#' ...
#'
#' @param x ...
#'
#' @param oc_cutoff ...
#'
#' @return ...
#'
#' @importFrom magrittr %>%
#'
#' @author Vlad Petyuk
#'
#' @export
#'
merge_clusters_oc <- function (x, oc_cutoff) {

  output <- x[[1]] %>%
    dplyr::group_by(Gene) %>%
    dplyr::do(out = merge_clusters_wrapper(., oc_cutoff = oc_cutoff)) %>%
    tidyr::unnest(cols=c(out))

  return (output)

}

# @author Vlad Petyuk
oc <- function (xi, xj) {

  ci <- xi[, "Proteoform", drop = TRUE]
  cj <- xj[, "Proteoform", drop = TRUE]

  min(sum(ci %in% cj), sum(cj %in% ci))/min(length(ci), length(cj))

}

# @author Vlad Petyuk
get_mergeable_clusters <- function(x, oc_cutoff){
  clusters <- unique(x$cluster)
  for(i in clusters){
    xi <- x[x$cluster == i,]
    for(j in clusters){
      if(i != j){
        xj <- x[x$cluster == j,]
        over_coef_val <- oc(xi, xj)
        # condition of a mergeable cluster
        if(over_coef_val >= oc_cutoff)
          return(c(i,j))
      }
    }
  }
  return(NULL)
}

# @author Vlad Petyuk
# Recursion
# not an efficient algorithm, especially given the updates of cluster membership
# after every iteration. But it works.
merge_clusters2 <- function(clusters_df, oc_cutoff){
  ij <- get_mergeable_clusters(clusters_df, oc_cutoff)
  if(is.null(ij)){
    return(clusters_df)
  }else{
    # assumes that "i" (ij[1]) is larger cluster than "j" (ij[2]). This can be
    # achieved by pre-ordering.
    clusters_df[clusters_df$cluster == ij[2], "cluster"] <- ij[1]
    merge_clusters2(clusters_df, oc_cutoff)
  }
}

# @author Vlad Petyuk
merge_clusters_wrapper <- function (x, oc_cutoff) {
  idx0 <- x$cluster == 0
  x$cluster_new <- 0
  y <- merge_clusters2(x[!idx0, c("Proteoform","cluster")], oc_cutoff)
  x$cluster_new[!idx0] <- y$cluster
  x$Gene <- NULL # removes duplicate gene due to the group/do functions.
  return (x)
}
