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
cluster <- function (x, errors, method, height, min_size) {

  return (
    x %>%
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
        NormRecalMass = log10(RecalMass) / log10(1 + errors$ppm_sd / 1e6)
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
  )

}
