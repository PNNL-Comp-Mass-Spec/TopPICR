#' Select one gene per feature
#'
#' A given feature can occur in multiple genes. This function selects a single
#' gene for each feature. The gene with the minimum E-value is kept and all
#' other genes are removed.
#'
#' @param x A \code{data.table} output from the \code{read_toppic} function.
#'
#' @return A \code{data.table}.
#'
#' @importFrom magrittr %>%
#'
#' @author Vlad Petyuk
#'
#' @export
#'
rm_false_gene <- function (x) {

  y <- x %>%
    dplyr::select(Dataset, CV, `Feature apex`, `Feature intensity`,
                  Gene, `E-value`) %>%
    dplyr::distinct() %>%
    dplyr::group_by(Dataset, CV, `Feature apex`, `Feature intensity`) %>%
    dplyr::slice_min(`E-value`) %>%
    dplyr::ungroup() %>%
    dplyr::select(Dataset, CV, `Feature apex`, `Feature intensity`, Gene)

  return (dplyr::semi_join(x, y))

}
