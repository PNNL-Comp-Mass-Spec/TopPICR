#' Determine the proteoform level for each \code{Proteoform}
#'
#' Assigns each \code{Proteoform} a level based on the classification system
#' proposed by Smith et al. 2019.
#'
#' @param x A \code{data.table} output from the \code{augment_annotation} (or
#'   another downstream) function.
#'
#' @return A \code{data.table} with the \code{Proteoform Level} variable added.
#'
#' @references Smith, L.M., Thomas, P.M., Shortreed, M.R. et al. A five-level
#'   classification system for proteoform identifications. Nat Methods 16,
#'   939–940 (2019). https://doi.org/10.1038/s41592-019-0573-x
#'
#' @author James Fulcher
#'
#' @export
#'
set_pf_level <- function (x) {

  x <- x %>%
    dplyr::mutate(
      PTMScore = stringr::str_replace_all(
        Proteoform, "\\.\\(([A-Z])\\)\\[Acetyl\\]", ""
      ),
      PTMScore = stringr::str_replace_all(
        PTMScore,
        "\\(\\(C\\)\\[Carbamidomethyl\\]\\(C\\)\\[Carbamidomethyl\\]\\)",
        "\\(C\\)"
      ),
      PTMScore = stringr::str_replace_all(
        PTMScore, "\\(\\(C\\)\\[Carbamidomethyl\\]\\)", "\\(C\\)"
      ),
      PTMScore = stringr::str_replace_all(
        PTMScore, "\\(C\\)\\[Carbamidomethyl\\]", ""
      ),
      PTMScore = stringr::str_replace_all(PTMScore, "(?=\\[).*?(?<=\\])", ""),
      PTMScore = stringr::str_replace_all(PTMScore, "\\(\\(M", "\\(M"),
      PTMScore = stringr::str_replace_all(PTMScore, "M\\)\\)", "M\\)"),
      PTMScore = gsub("[[:digit:]]+", "", PTMScore),
      PTMScore = gsub( "[[\\.]]","", PTMScore),
      PTMScore = gsub( "\\[","", PTMScore),
      PTMScore = gsub( "\\-","", PTMScore),
      PTMScore = gsub("\\.","", PTMScore),
      PTMScore = stringr::str_replace_all(PTMScore, "(?<=\\)).*(?=\\()", ""),
      PTMScore = stringr::str_replace_all(PTMScore, "\\)\\(", ""),
      PTMScore = stringr::str_extract(PTMScore, "(?<=\\().*(?=\\))"),
      PTMScore = stringr::str_replace_all(PTMScore, "[[:punct:]]", "")
    ) %>%
    dplyr::mutate(
      `Proteoform Level` = dplyr::case_when(
        ## No ambiguity
        `#unexpected modifications` == 0 & AccMap == 1 ~ "1",
        ## TopPIC error assignment
        stringr::str_detect(MIScore, "\\[\\]") ~ "5",
        ## No ambiguity
        (MIScore != "-" & `#unexpected modifications` == 1 & AccMap == 1 &
           nchar(PTMScore) == 1) ~ "1",
        ## No ambiguity
        (MIScore != "-" & `#unexpected modifications` == 2 & AccMap == 1 &
           nchar(PTMScore) == 2) ~ "1",
        ## PTM localization ambiguity
        (MIScore != "-" & `#unexpected modifications` == 1 & AccMap == 1 &
           nchar(PTMScore) > 1) ~ "2A",
        ## PTM localization ambiguity
        (MIScore != "-" & `#unexpected modifications` == 2 & AccMap == 1 &
           nchar(PTMScore) > 2) ~ "2A",
        ## PTM identification ambiguity
        (MIScore == "-" & `#unexpected modifications` == 1 & AccMap == 1 &
           nchar(PTMScore) == 1) ~ "2B",
        ### Gene of origin ambiguity
        `#unexpected modifications` == 0 & AccMap > 1 ~ "2D",
        ## Gene ambiguity
        (MIScore != "-" & `#unexpected modifications` == 1 & AccMap > 1 &
           nchar(PTMScore) == 1) ~ "2D",
        ## Gene ambiguity
        (MIScore != "-" & `#unexpected modifications` == 1 & AccMap > 1 &
           is.na(PTMScore)) ~ "2D",
        ## Gene ambiguity
        (MIScore != "-" & `#unexpected modifications` == 2 & AccMap > 1 &
           nchar(PTMScore) == 2) ~ "2D",
        ## PTM identification and gene ambiguity
        (MIScore == "-" & `#unexpected modifications` == 1 & AccMap > 1 &
           nchar(PTMScore) == 1) ~ "3E",
        ## PTM identification and localization ambiguity
        (MIScore == "-" & `#unexpected modifications` == 1 & AccMap == 1 &
           nchar(PTMScore) > 1) ~ "3A",
        ## Gene ambiguity and localization of PTM
        (MIScore != "-" & `#unexpected modifications` == 1 & AccMap > 1 &
           nchar(PTMScore) > 1) ~ "3D",
        ## Gene ambiguity and localization of PTM
        (MIScore != "-" & `#unexpected modifications` == 2 & AccMap > 1 &
           nchar(PTMScore) > 2) ~ "3D",
        ## Gene ambiguity, PTM ambiguity, and localization ambiguity
        (MIScore == "-" & `#unexpected modifications` == 1 & AccMap > 1 &
           nchar(PTMScore) > 1) ~ "4B"
      )
    ) %>%
    dplyr::select(-PTMScore)

}
