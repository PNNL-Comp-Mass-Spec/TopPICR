#' Filter UniProt accessions by count
#'
#' Filters out UniProt accessions that occur in fewer than \code{threshold}
#' project IDs. Accessions are counted by the \code{count_by} argument. For
#' example, if \code{threshold} is 10 and \code{count_by} is Gene all rows
#' containing UniProt accessions that occur in fewer than 10 project IDs will be
#' removed. UniProt accessions are counted within Gene.
#'
#' @param x ...
#'
#' @param count_by A Character string indicating what variable observations
#'   should be counted by (e.g., Gene, ProteoForm, ...).
#'
#' @param threshold An integer indicating the number of subjects (or project
#'   IDs) a UniProt accession must occur in to be kept. The default is 10.
#'
#' @return ...
#'
#' @importFrom magrittr %>%
#'
#' @export
#'
filter_accession <- function (x, count_by, threshold = 10) {

  # Filter by the number of samples each ProteoForm occurs in ------------------

  # Remove from x any rows that contain values in the ProteoForm column that do
  # not occur in at least 10 subjects (ProjID).
  x <- x %>%
    dplyr::distinct(ProjID, !!rlang::sym(count_by), UniProtAcc) %>%
    dplyr::group_by(!!rlang::sym(count_by), UniProtAcc) %>%
    dplyr::tally() %>%
    dplyr::filter(n >= threshold) %>%
    dplyr::select(-n) %>%
    dplyr::semi_join(x, .) %>%
    dplyr::ungroup()

  return (x)

}
