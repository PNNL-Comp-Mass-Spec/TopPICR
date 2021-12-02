#' Cluster data based on mass and retention time
#'
#' Cluster the data using hierarchical clustering with the
#' \code{\link[stats]{hclust}} function. The data are clustered with the
#' normalized recalibrated mass and the normalized aligned retention times.
#'
#' @param x A \code{data.table} output from the \code{recalibrate_mass}
#'   function.
#'
#' @param errors A \code{list} output from the \code{calc_error} function.
#'
#' @param repMass Logical. If TRUE the representative mass will be used when
#'   clustering.
#'
#' @param method A character string indicating what agglomeration method should
#'   be used in the hclust function. See \code{\link[stats]{hclust}} for more
#'   details.
#'
#' @param height An number specifying the height at which the tree created by
#'   the hclust function should be cut. See \code{\link[stats]{hclust}} for more
#'   details.
#'
#' @param min_size An integer value indicating the minimum number of points a
#'   cluster must have. All clusters with fewer members than min_size will be
#'   reclassified as "noise" points.
#'
#' @return A \code{data.table} with the cluster assignment for each observation.
#'   The following variables have been added/removed:
#'
#'   | Added             | Removed                    |
#'   | ----------------- | -------------------------- |
#'   | `cluster`         |                            |
#'
#' @md
#'
#' @importFrom magrittr %>%
#'
#' @author Evan A Martin
#'
#' @export
#'
cluster <- function (x, errors, repMass = TRUE, method, height, min_size) {

  # Set the mass object according to the repMass input and copy the x object
  # input to join with after clustering.
  if (repMass) {

    # Create a variable name that points to the correct column name based on the
    # input of the repMass argument. We will either cluster on the RecalMass or
    # repMass variable.
    mass <- "repMass"

    # Check if repMass is present in the x data frame.
    if (!("repMass" %in% names(x))) {

      stop ("The variable repMass is not present in x.")

    }

    # Save the input to merge with after clustering. This will retain the same
    # number of rows in the output as the input.
    x_input <- x

    # Remove redundant rows. This is necessary for determining the noise points
    # later in the clustering step. If there is an RTalign/repMass point that is
    # repeated more times than the value of min_size the point will be
    # classified as a cluster even though there is only one point in the
    # cluster.
    x <- x %>%
      dplyr::distinct(Dataset, Gene, RTalign, repMass)

  } else {

    mass <- "RecalMass"

  }

  x_cluster <- x %>%
    dplyr::group_by(Gene) %>%
    # Remove rows corresponding to Genes with only one observation because
    # hclust will throw an error unless there are at least two observations.
    dplyr::add_count(name = "obs") %>%
    dplyr::filter(obs > 1) %>%
    dplyr::select(-obs) %>%
    dplyr::ungroup() %>%
    # Normalize the mass according to the ppm error that was computed in the
    # calc_error function. The normalized recalibrated mass will be used for
    # clustering. This means the h argument in the cutree function will
    # correspond to the standard deviation.
    dplyr::mutate(
      NormRecalMass = log10(!!rlang::sym(mass)) / log10(1 + errors$ppm_sd / 1e6)
    ) %>%
    # Normalize the aligned rt according to the rt error that was computed in
    # the calc_error function. The normalized rt will be used for clustering.
    # This means the h argument in the cutree function will correspond to the
    # standard deviation.
    dplyr::mutate(NormRTalign = RTalign / errors$rt_sd) %>%
    dplyr::ungroup() %>%
    dplyr::nest_by(Gene) %>%
    dplyr::mutate(
      cluster = list(
        stats::cutree(
          stats::hclust(stats::dist(dplyr::select(data,
                                                  NormRTalign,
                                                  NormRecalMass)),
                        method = method),
          h = height
        )
      )
    ) %>%
    # Find all clusters that have fewer than min_size members. In the next step
    # these points will be reclassified as "noise" points.
    dplyr::mutate(noise = list(tibble::tibble(clust = cluster) %>%
                                 dplyr::group_by(clust) %>%
                                 dplyr::tally() %>%
                                 dplyr::filter(n < min_size) %>%
                                 dplyr::pull(clust))) %>%
    # Convert all clusters with fewer than min_size members to cluster 0.
    dplyr::mutate(cluster = list(replace(cluster, cluster %in% noise, 0))) %>%
    dplyr::select(Gene, data, cluster) %>%
    tidyr::unnest(cols = c(data, cluster)) %>%
    dplyr::select(-NormRecalMass, -NormRTalign) %>%
    dplyr::ungroup()

  # If repMass is TRUE the x_cluster object needs to be joined with the original
  # input to preserve all the rows and columns in the original input.
  if (repMass) {

    # Combine the clustered data with the input data.
    x_cluster <- dplyr::inner_join(x_input, x_cluster)

  }

  # Return the cluster data frame.
  return (x_cluster)

}

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
