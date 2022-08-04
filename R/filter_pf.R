#' Filter low count proteoforms
#'
#' Counts the number of samples each proteoform occurs in and filters out any
#' proteoforms that occur in fewer samples than the threshold.
#'
#' @param x A \code{data.table} output from either \code{augment_annotation} or
#'   \code{map_proteoform}.
#'
#' @param count_within A character vector specifying the variables that will be
#'   grouped before counting the number of occurrences of the \code{count}
#'   variable.
#'
#' @param count A character string indicating the variable that will be counted.
#'
#' @param threshold An integer. The value that will be used to remove rows whose
#'   count falls below the threshold.
#'
#' @return A \code{data.table} where all rows containing proteoforms with counts
#'   below the threshold are filtered out.
#'
#' @importFrom magrittr %>%
#'
#' @export
#'
filter_by_count <- function (x, count_within, count, threshold) {

  return (
    x %>%
      dplyr::select(dplyr::all_of(count_within), !!rlang::sym(count)) %>%
      dplyr::distinct() %>%
      dplyr::group_by(!!rlang::sym(count)) %>%
      dplyr::tally() %>%
      dplyr::filter(n >= threshold) %>%
      dplyr::select(-n) %>%
      dplyr::semi_join(x, .) %>%
      dplyr::ungroup()
  )

}
