#' Convert TopPICR output to an MSnSet object
#'
#' Takes the combined data (output from the match between runs step) and the
#' meta data object and converts them to an MSnSet object.
#'
#' @param combined A \code{data.table} output from the \code{combine_features}
#'   function.
#'
#' @param meta A \code{data.table} output from the \code{create_mdata} function.
#'
#' @return An MSnSet object
#'
# ' @export
#'
#' @author Evan A Martin
#'
deprecated_to_msnset <- function (combined, meta) {

  # Check if FAIMS data. If not there is no need to have CV in the feature name.
  if (combined$CV[[1]] == "noFAIMS") {

    x_msn <- combined %>%
      dplyr::mutate(
        feature_name = stringr::str_c(Gene, pcGroup, sep = "_")
      )

    # Code runs if there are multiple FAIMS CVs.
  } else {

    x_msn <- combined %>%
      dplyr::mutate(
        feature_name = stringr::str_c(Gene, pcGroup, CV, sep = "_")
      )

  }

  # Create MSnSet expression object.
  x_expr <- x_msn %>%
    tidyr::pivot_wider(
      id_cols = "feature_name",
      names_from = "Dataset",
      values_from = "Intensity"
    ) %>%
    data.frame(check.names = FALSE) %>%
    `rownames<-`(.$feature_name) %>%
    dplyr::select(-feature_name) %>%
    data.matrix()

  # Create MSnSet feature object.
  x_feat <- x_msn %>%
    dplyr::group_by(Gene, pcGroup, CV, feature_name) %>%
    dplyr::summarize(median_intensity = median(Intensity),
                     count = dplyr::n(),
                     .groups = "keep") %>%
    dplyr::ungroup() %>%
    dplyr::left_join(meta, by = c("Gene", "pcGroup")) %>%
    dplyr::mutate(proteoform_id = paste(Gene, pcGroup, sep = "_")) %>%
    data.frame() %>%
    # For some reason there is one Gene pcGroup that occurs more than once. We
    # will use distinct to select just one of them.
    dplyr::distinct(Gene, pcGroup, .keep_all = TRUE) %>%
    `rownames<-`(.$feature_name)

  return (
    MSnbase::MSnSet(x_expr, x_feat[rownames(x_expr), ])
  )

}
