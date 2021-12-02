#' Create feature metadata
#'
#' This function uses the information from the identified feature data to create
#' a feature metadata data frame. This data frame will be used at the end of the
#' `TopPICR` pipeline when creating an `MSnSet` object.
#'
#' @param x A \code{data.table} output from the \code{mk_pcg} function.
#'
#' @importFrom magrittr %>%
#'
#' @author Evan A Martin
#'
#' @export
#'
mk_metadata <- function (x) {

  x_meta <- x %>%
    dplyr::ungroup()

  # Return the metadata data frame.
  return (x_meta)

}

# x_meta <- x_grp %>%
#   ungroup() %>%
#   group_by(Gene, pcGroup, Proteoform) %>%
#   mutate(n_pf = n(),
#          rt = median(RTalign, na.rm = TRUE),
#          mass = median(RecalMass, na.rm = TRUE)) %>%
#   ungroup() %>%
#   group_by(Gene, pcGroup) %>%
#   slice_max(n_pf) %>%
#   slice_max(`Feature intensity`) %>%
#   slice_min(`E-value`) %>%
#   ungroup() %>%
#   select(Gene, pcGroup, UniProtAcc, mass, rt,
#          firstAA, lastAA, protLength, Proteoform)

# Test gene and pc group:
# Gene - ACBD7; pcGroup - 21
