#' ...
#'
#' ...
#'
#' @param x A list with three elements output from the x_47 function. The first
#'   element is the the data matrix with the aligned retention times and
#'   recalibrated mass. The second and third elements are the median of the ppm
#'   and retention time errors respectively. These values are calculated by
#'   Dataset.
#'
#' @param method A character string indicating what agglomeration method should
#'   be used in the hclust function.
#'
#' @param height An number specifying the height at which the tree created by
#'   the hclust function should be cut.
#'
#' @param min_size An integer value indicating the minimum number of points a
#'   cluster must have. All clusters with fewer members than min_size will be
#'   reclassified as "noise" points.
#'
#' @return ...
#'
#' @importFrom magrittr %>%
#'
#' @export
#'
cluster <- function (x, method, height, min_size) {

  x_cluster <- x[[1]] %>%
    group_by(Gene) %>%
    # Normalize the mass according to the ppm error that was computed in the
    # x_47 function. The normalized recalibrated mass will be used for
    # clustering. This means the h argument in the cutree function will
    # correspond to the standard deviation.
    mutate(NormRecalMass = log10(RecalMass) / log10(1 + x[[2]] / 1e6)) %>%
    # Normalize the aligned rt according to the rt error that was computed in
    # the x_47 function. The normalized rt will be used for clustering. This
    # means the h argument in the cutree function will correspond to the
    # standard deviation.
    mutate(NormRTalign = RTalign / x[[3]]) %>%
    ungroup() %>%
    nest_by(Gene) %>%
    mutate(hcluster = list(hclust(dist(select(data,
                                              NormRTalign,
                                              NormRecalMass)),
                                  method = method))) %>%
    mutate(cluster = list(cutree(hcluster, h = height))) %>%
    # Find all clusters that have fewer than min_size members. In the next step
    # these points will be reclassified as "noise" points.
    mutate(noise = list(tibble(clust = cluster) %>%
                          group_by(clust) %>%
                          tally() %>%
                          filter(n < min_size) %>%
                          pull(clust))) %>%
    # Convert all clusters with fewer than min_size members to cluster 0.
    mutate(cluster = list(replace(cluster, cluster %in% noise, 0))) %>%
    select(Gene, data, cluster) %>%
    unnest(cols = c(data, cluster))

  # Return the cluster data frame along with the ppm and rt errors.
  return (list(x = x_cluster,
               ppm_error = x[[2]],
               rt_error = x[[3]]))

}
