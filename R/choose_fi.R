#' Choose a representative feature intensity
#'
#' ...
#'
#' @param x ...
#'
#' @param x_cluster ...
#'
#' @param by A character string corresponding to the variables that will be used
#'   to select a representative feature intensity value. Can be either "cv" or
#'   "cv_fraction".
#'
#' @return ...
#'
#' @importFrom magrittr %>%
#'
#' @export
#'
choose_fi <- function (x, x_cluster, by) {

  # Make sure by is either "cv" or "cv_fraction".
  if (!(by %in% c("cv", "cv_fraction")))
    stop ("by must be either 'cv' or 'cv_fraction'")

  # Choose the representative feature intensity values based on the input to by.
  if (by == "cv") {

    cv_only(x = x, x_cluster = x_cluster)

  } else if (by == "cv_fraction") {

    cv_frac(x = x, x_cluster = x_cluster)

  }


}

# choose_fi auxiliary functions ------------------------------------------------

cv_only <- function (x, x_cluster) {

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

cv_frac <- function (x, x_cluster) {

  output <- x_cluster %>%
    dplyr::filter(cluster_new != 0) %>%
    dplyr::select(-cluster) %>%
    dplyr::rename(cluster = cluster_new) %>%
    dplyr::select(Gene, Dataset, Proteoform, cluster)

  x2 <- x %>%
    dplyr::select(Dataset, Gene, `Feature intensity`, ProjID, CV, Fraction) %>%
    dplyr::inner_join(output) %>%
    dplyr::mutate(gcc = paste(Gene, cluster, sep = "_")) %>%
    dplyr::select(-Gene, -cluster)

  rep_cv_fraction <- x2 %>%
    dplyr::group_by(ProjID, gcc, CV, Fraction) %>%
    dplyr::summarize(max_int = max(`Feature intensity`)) %>%
    dplyr::group_by(gcc, CV, Fraction) %>%
    dplyr::summarize(n_obs = dplyr::n(),
                     int = sum(max_int)) %>%
    dplyr::group_by(gcc) %>%
    dplyr::slice(which.max(n_obs))

  x_rep <- x2 %>%
    dplyr::group_by(ProjID, gcc, CV, Fraction) %>%
    dplyr::summarize(max_int = max(`Feature intensity`)) %>%
    dplyr::semi_join(rep_cv_fraction) %>%
    dplyr::ungroup() %>%
    dplyr::select(ProjID, gcc, max_int) %>%
    dplyr::rename(`Feature intensity` = max_int)

  return (x_rep)

}
