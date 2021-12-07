#' Choose feature intensity from the best CV/Fraction
#'
#' ...
#'
#' @param x A \code{data.table} output from the \code{collapse_pcg} function.
#'
#' @return A \code{data.table} ...
#'
#' @importFrom magrittr %>%
#'
#' @author Evan A Martin, Vlad Petyuk
#'
#' @export
#'
choose_fi <- function (x) {

  # Generate the Fraction variable within the unidentified data.
  x <- x %>%
    dplyr::mutate(Fraction = grepl("_WO_", Dataset)) %>%
    # Keep the Intensity value corresponding to the best CV/Fraction. The best
    # is defined as the maximum Intensity value within the CV/Fraction
    # combination. First, select the "best" Fraction for each Gene and pcGroup.
    # If there are two CV/Fraction combinations with the same counts (number of
    # unique Datasets) choose the CV/Fraction with the highest median Intensity.
    dplyr::group_by(Fraction, CV, Gene, pcGroup) %>%
    dplyr::mutate(fCounts = dplyr::n_distinct(Dataset)) %>%
    dplyr::group_by(Gene, pcGroup) %>%
    dplyr::slice_max(fCounts) %>%
    dplyr::group_by(Fraction, CV, Gene, pcGroup) %>%
    dplyr::mutate(medianFi = stats::median(Intensity)) %>%
    dplyr::group_by(Gene, pcGroup) %>%
    dplyr::slice_max(medianFi) %>%
    dplyr::select(-fCounts, -medianFi)

  return (x)

}
