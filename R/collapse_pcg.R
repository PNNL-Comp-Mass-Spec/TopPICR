#' Choose a representative feature intensity
#'
#' ...
#'
#' @param x A \code{data.table} output from the \code{gen_pcg} function. This
#'   \code{data.table} contains the \code{pcGroup} variable which will be used
#'   when choosing a representatvie feature intensity value.
#'
#' @param feature A \code{data.table} containing the unidentified feature data
#'   output from the \code{impute_fi} function.
#'
#' @return A \code{data.table} with the representative feature intensity values.
#'
#' @importFrom magrittr %>%
#'
#' @export
#'
collapse_pcg <- function (x, feature) {

  # Generate the Fraction variable within the unidentified data.
  feature <- feature %>%
    dplyr::mutate(Fraction = grepl("_WO_", Dataset))

  fused <- dplyr::inner_join(
    feature,
    dplyr::distinct(x, Gene, cluster, pcGroup)
  ) %>%
    # The following code keeps only the highest intensity value when a given
    # Dataset has multiple Intensity values within a pcGroup. This computation
    # is carried out within each Gene, CV, and Fraction.
    dplyr::group_by(Dataset, CV, Gene, Fraction, pcGroup) %>%
    dplyr::summarize(Intensity = max(Intensity, na.rm = TRUE)) %>%
    dplyr::ungroup()

  # Keep the Intensity value corresponding to the best CV/Fraction. The best is
  # defined as the maximum Intensity value within the CV/Fraction combination.
  chosen <- fused %>%
    # First, select the "best" Fraction for each Gene and pcGroup. This is the
    # Fraction with the highest number of data sets. If there are two fractions
    # with the same Dataset counts choose the fraction with the highest median
    # Intensity.
    dplyr::group_by(Fraction, Gene, pcGroup) %>%
    dplyr::mutate(fCounts = dplyr::n_distinct(Dataset)) %>%
    dplyr::group_by(Gene, pcGroup) %>%
    dplyr::slice_max(fCounts) %>%
    dplyr::group_by(Fraction, Gene, pcGroup) %>%
    dplyr::mutate(medianFi = stats::median(Intensity)) %>%
    dplyr::group_by(Gene, pcGroup) %>%
    dplyr::slice_max(medianFi) %>%
    # Second, select the "best" CV for each Gene and pcGroup. This process is
    # the same as above except with using CV instead of Fraction when grouping
    # the data.
    dplyr::group_by(CV, Gene, pcGroup) %>%
    dplyr::mutate(counts = dplyr::n_distinct(Dataset)) %>%
    dplyr::group_by(Gene, pcGroup) %>%
    dplyr::slice_max(counts) %>%
    dplyr::group_by(CV, Gene, pcGroup) %>%
    dplyr::mutate(median_fi = stats::median(Intensity)) %>%
    dplyr::group_by(Gene, pcGroup) %>%
    dplyr::slice_max(median_fi) %>%
    dplyr::group_by(Dataset, Gene, pcGroup) %>%
    dplyr::summarize(Intensity = max(Intensity))

  return (chosen)

}
