#' Combine identified and unidentified feature data
#'
#' Combines the identified and unidentified feature data. The maximum feature
#' intensity is selected for each combination of Dataset, CV, Gene, and pcGroup
#' after joining the two data tables.
#'
#' @param ms2 A \code{data.table} containing the identified feature data.
#'
#' @param ms1 A \code{data.table} containing the unidentified feature data.
#'
#' @param summary_fn A character string specifying the function that will be
#'   used to select a representative feature intensity. The current functions
#'   allowed are "max" or "sum". The default input is "max".
#'
#' @return A \code{data.table} containing the combined ms1 and ms2 feature
#'   intensity values.
#'
#' @importFrom magrittr %>%
#'
#' @author Vlad Petyuk
#'
#' @export
#'
combine_features <- function (ms2, ms1, summary_fn = "max") {

  # Check if summary_fn is either "max" or "sum". For now we aren't allowing any
  # other function to be used to select a representative feature intensity.
  if (!(summary_fn %in% c("max", "sum"))) {

    stop ("The input for summary_fn must either be 'max' or 'sum'.")

  }

  # Create an object that points to the function that will be used to create a
  # representative feature intensity value.
  the_summary <- match.fun(summary_fn)

  fused <- dplyr::inner_join(
    ms1,
    dplyr::distinct(ms2, Gene, cluster, pcGroup)
  ) %>%
    # The following code keeps only the highest intensity value when a given
    # Dataset has multiple Intensity values within a pcGroup. This computation
    # is carried out within each Gene, CV, and Fraction.
    dplyr::group_by(Dataset, CV, Gene, pcGroup) %>%
    dplyr::summarize(Intensity = the_summary(Intensity, na.rm = TRUE),
                     RTalign = mean(RTalign, na.rm = TRUE),
                     RTapex = mean(Time_apex, na.rm = TRUE)) %>%
    dplyr::ungroup()

  return (fused)

}
