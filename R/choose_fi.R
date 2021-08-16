#' Choose a representative feature intensity
#'
#' ...
#'
#' @param x ...
#'
#' @param x_cluster ...
#'
#' @return ...
#'
#' @importFrom magrittr %>%
#'
#' @export
#'
choose_fi <- function (x, x_cluster) {

  # Join x and x_cluster. This is necessary to link the cluster with the
  # variables dropped before aligning the retention time, recalibrating the
  # mass, and performing hierarchical clustering. Then combine Gene and
  # cluster_new to be used as a new proteoform.
  xC <- dplyr::inner_join(x, x_cluster) %>%
    dplyr::filter(cluster_new != 0) %>%
    dplyr::mutate(gcc = paste(Gene, cluster_new, sep = "_")) %>%
    dplyr::select(gcc, ProjID, CV, Fraction, `Feature intensity`)

  chosen_fi <- xC %>%
    dplyr::group_by(gcc, ProjID, CV, Fraction) %>%
    # Remove multiple feature intensity values for each combination of the
    # grouping variables (if they exist).
    dplyr::summarize(`Feature intensity` = max(`Feature intensity`)) %>%
    dplyr::group_by(gcc, CV) %>%
    dplyr::mutate(counts = dplyr::n_distinct(ProjID)) %>%
    dplyr::group_by(gcc) %>%
    dplyr::slice_max(counts) %>%
    dplyr::group_by(gcc, CV) %>%
    dplyr::mutate(median_fi = stats::median(`Feature intensity`)) %>%
    dplyr::group_by(gcc) %>%
    dplyr::slice_max(median_fi) %>%
    dplyr::group_by(gcc, ProjID) %>%
    dplyr::summarize(`Feature intensity` = max(`Feature intensity`))

  return (chosen_fi)


}
