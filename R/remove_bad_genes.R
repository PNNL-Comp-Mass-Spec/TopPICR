#' ...
#'
#' ...
#'
#' @param x ...
#'
#' @return ...
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
    dplyr::slice(which.min(`E-value`)) %>%
    dplyr::ungroup() %>%
    dplyr::select(Dataset, CV, `Feature apex`, `Feature intensity`, Gene)

  x2 <- x %>%
    dplyr::semi_join(y)

  return (x2)

}
