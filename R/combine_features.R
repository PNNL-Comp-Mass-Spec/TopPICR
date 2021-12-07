#' Combine identified and unidentified feature data
#'
#' Combines the identified and unidentified feature data. The maximum feature
#' intensity is selected for each combination of Dataset, CV, Gene, and pcGroup
#' after joining the two data tables.
#'
#' @param x A \code{data.table} output from the \code{gen_pcg} function.
#'
#' @param feature A \code{data.table} containing the unidentified feature data
#'   output from the \code{impute_fi} function.
#'
#' @return A \code{data.table} ...
#'
#' @importFrom magrittr %>%
#'
#' @author Vlad Petyuk
#'
#' @export
#'
combine_features <- function (x, feature) {

  fused <- dplyr::inner_join(
    feature,
    dplyr::distinct(x, Gene, cluster, pcGroup)
  ) %>%
    # The following code keeps only the highest intensity value when a given
    # Dataset has multiple Intensity values within a pcGroup. This computation
    # is carried out within each Gene, CV, and Fraction.
    dplyr::group_by(Dataset, CV, Gene, pcGroup) %>%
    dplyr::summarize(Intensity = max(Intensity, na.rm = TRUE),
                     RTalign = mean(RTalign, na.rm = TRUE),
                     RTapex = mean(Time_apex, na.rm = TRUE)) %>%
    dplyr::ungroup()

  return (fused)

}
