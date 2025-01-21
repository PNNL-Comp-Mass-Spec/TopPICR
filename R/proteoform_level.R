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
#'   939-940 (2019). https://doi.org/10.1038/s41592-019-0573-x
#'
#' @author James Fulcher
#'
#' @export
#'
set_pf_level <- function (x) {

  x$`Special amino acids`[x$`Special amino acids` == "-"] <- NA
  
  x  <- x %>% dplyr::mutate(
    PTMScore = stringr::str_replace_all(
      Proteoform, "\\.\\(([A-J,L-Z])\\)\\[Acetyl\\]", ""
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
    PTMScore = gsub("[[\\.]]", "", PTMScore),
    PTMScore = gsub("\\[", "", PTMScore),
    PTMScore = gsub("\\-", "", PTMScore),
    PTMScore = gsub("\\.", "", PTMScore),
    PTMScore = stringr::str_replace_all(PTMScore, "(?<=\\)).*(?=\\()", ""),
    PTMScore = stringr::str_replace_all(PTMScore, "\\)\\(", ""),
    PTMScore = stringr::str_extract(PTMScore, "(?<=\\().*(?=\\))"),
    PTMScore = stringr::str_replace_all(PTMScore, "[[:punct:]]", "")
  ) %>%
    dplyr::mutate(
      `Proteoform Level` = dplyr::case_when(
        MIScore == "-" &
          `#unexpected modifications` ==
          0 &
          AccMap == 1  &
          is.na(PTMScore) &
          is.na(`Special amino acids`) ~ "1",
        stringr::str_detect(MIScore, "\\[\\]") ~ "5",
        (
          MIScore != "-" & `#unexpected modifications` ==
            0 &
            AccMap == 1 & nchar(PTMScore) == 1 &
            is.na(`Special amino acids`)
        ) ~ "1",
        (
          MIScore !=
            "-" &
            `#unexpected modifications` == 0 & AccMap == 1 &
            nchar(PTMScore) == 2 &
            is.na(`Special amino acids`)
        ) ~ "1",
        (
          MIScore != "-" & `#unexpected modifications` ==
            0 &
            AccMap == 1 & nchar(PTMScore) > 1 &
            is.na(`Special amino acids`)
        ) ~ "2A",
        (
          MIScore !=
            "-" &
            `#unexpected modifications` == 1 & AccMap == 1 &
            nchar(PTMScore) > 2 &
            is.na(`Special amino acids`)
        ) ~ "2A",
        (
          MIScore == "-" & `#unexpected modifications` ==
            1 &
            AccMap == 1 &
            nchar(PTMScore) == 1 &
            is.na(`Special amino acids`)
        ) ~ "2B",
        (
          MIScore == "-" & `#unexpected modifications` ==
            2 &
            AccMap == 1 &
            nchar(PTMScore) == 2 &
            is.na(`Special amino acids`)
        ) ~ "2B",
        (MIScore !=
           "-" &
           `#unexpected modifications` ==
           0 &
           AccMap == 1 &
           !is.na(`Special amino acids`)
        ) ~ "2C",
        ( MIScore !=
            "-" &
            `#unexpected modifications` ==
            1 &
            AccMap == 1 &
            nchar(PTMScore) == 2 &
            !is.na(`Special amino acids`)
        ) ~ "2C",
        ( MIScore ==
            "-" &
            `#unexpected modifications` ==
            0 &
            AccMap == 1 &
            is.na(PTMScore) &
            !is.na(`Special amino acids`)
        ) ~ "2C",
        (  MIScore != "-" & `#unexpected modifications` ==
             0 &
             AccMap > 1 & nchar(PTMScore) == 1 &
             is.na(`Special amino acids`)
        ) ~ "2D",
        (
          MIScore !=
            "-" &
            `#unexpected modifications` == 0 & AccMap > 1 &
            is.na(PTMScore) &
            is.na(`Special amino acids`)
        ) ~ "2D",
        (
          MIScore != "-" & `#unexpected modifications` ==
            1 &
            AccMap > 1 & nchar(PTMScore) == 2 &
            is.na(`Special amino acids`)
        ) ~ "2D",
        (
          MIScore == "-" & `#unexpected modifications` ==
            0 &
            AccMap > 1 & is.na(PTMScore) &
            is.na(`Special amino acids`)
        ) ~ "2D",
        (
          MIScore == "-" & `#unexpected modifications` ==
            1 &
            AccMap == 1 & nchar(PTMScore) > 1 &
            is.na(`Special amino acids`)
        ) ~ "3A",
        (
          !MIScore == "-" & `#unexpected modifications` ==
            0 &
            AccMap == 1 & nchar(PTMScore) > 1 &
            !is.na(`Special amino acids`)
        ) ~ "3B",
        (
          !MIScore == "-" & `#unexpected modifications` ==
            1 &
            AccMap == 1 & nchar(PTMScore) > 2 &
            !is.na(`Special amino acids`)
        ) ~ "3B",
        (
          MIScore == "-" & `#unexpected modifications` ==
            1 &
            AccMap == 1 & nchar(PTMScore) == 1 &
            !is.na(`Special amino acids`)
        ) ~ "3C",
        (
          MIScore == "-" & `#unexpected modifications` ==
            2 &
            AccMap == 1 & nchar(PTMScore) == 2 &
            !is.na(`Special amino acids`)
        ) ~ "3C",
        (
          MIScore !=
            "-" &
            `#unexpected modifications` == 0 & AccMap > 1 &
            nchar(PTMScore) > 1 &
            is.na(`Special amino acids`)
        ) ~ "3D",
        (
          MIScore != "-" & `#unexpected modifications` ==
            1 &
            AccMap > 1 & nchar(PTMScore) > 2 &
            is.na(`Special amino acids`)
        ) ~ "3D",
        (
          MIScore ==
            "-" &
            `#unexpected modifications` == 1 & AccMap > 1 &
            nchar(PTMScore) == 1 &
            is.na(`Special amino acids`)
        ) ~ "3E",
        (
          MIScore ==
            "-" &
            `#unexpected modifications` == 2 & AccMap > 1 &
            nchar(PTMScore) == 2 &
            is.na(`Special amino acids`)
        ) ~ "3E",
        (
          !MIScore ==
            "-" &
            `#unexpected modifications` == 0 & AccMap > 1 &
            nchar(PTMScore) == 1 &
            !is.na(`Special amino acids`)
        ) ~ "3F",
        (
          !MIScore ==
            "-" &
            `#unexpected modifications` == 1 & AccMap > 1 &
            nchar(PTMScore) == 2 &
            !is.na(`Special amino acids`)
        ) ~ "3F",
        (
          MIScore ==
            "-" &
            `#unexpected modifications` == 0 & AccMap > 1 &
            is.na(PTMScore) &
            !is.na(`Special amino acids`)
        ) ~ "3F",
        (
          MIScore ==
            "-" &
            `#unexpected modifications` == 1 & AccMap == 1 &
            nchar(PTMScore) > 1 &
            !is.na(`Special amino acids`)
        ) ~ "4A",
        (
          MIScore ==
            "-" &
            `#unexpected modifications` == 2 & AccMap == 1 &
            nchar(PTMScore) > 2 &
            !is.na(`Special amino acids`)
        ) ~ "4A",
        (
          MIScore ==
            "-" &
            `#unexpected modifications` == 1 & AccMap > 1 &
            nchar(PTMScore) > 1 &
            is.na(`Special amino acids`)
        ) ~ "4B",
        (
          MIScore ==
            "-" &
            `#unexpected modifications` == 2 & AccMap > 1 &
            nchar(PTMScore) > 2 &
            is.na(`Special amino acids`)
        ) ~ "4B",
        (
          !MIScore ==
            "-" &
            `#unexpected modifications` == 0 & AccMap > 1 &
            nchar(PTMScore) > 1 &
            !is.na(`Special amino acids`)
        ) ~ "4C",
        (
          !MIScore ==
            "-" &
            `#unexpected modifications` == 1 & AccMap > 1 &
            nchar(PTMScore) > 2 &
            !is.na(`Special amino acids`)
        ) ~ "4C",
        (
          MIScore ==
            "-" &
            `#unexpected modifications` == 1 & AccMap > 1 &
            nchar(PTMScore) == 1 &
            !is.na(`Special amino acids`)
        ) ~ "4D",
        (
          MIScore ==
            "-" &
            `#unexpected modifications` == 2 & AccMap > 1 &
            nchar(PTMScore) == 2 &
            !is.na(`Special amino acids`)
        ) ~ "4D",
        (
          MIScore ==
            "-" &
            `#unexpected modifications` == 1 & AccMap > 1 &
            nchar(PTMScore) > 1 &
            !is.na(`Special amino acids`)
        ) ~ "5",
      )
    ) %>%
    dplyr::select(-PTMScore)

}
