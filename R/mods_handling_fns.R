#' Report the modifications found on a protein
#'
#' @param x A \code{data.frame} fData of an `MSnSet` object.
#'
#' @param accession A character string specifying which UniProt accession will
#'   be plotted.
#'
#' @export
#'
#'
get_mods_counts <- function (x, accession) {
  y <- filter(x, UniProtAcc == accession) %>%
    mutate(stuff = map(mods, `$`, mod_names)) %>%
    pull(stuff) %>%
    unlist() %>%
    table() %>%
    sort() %>%
    rev()
  return(y)
}

#' @describeIn get_mods_counts
#'
#' @export
#'
get_distinct_mods <- function(x, accession){
  y <- filter(x, UniProtAcc == accession) %>%
    mutate(stuff = map(mods, `$`, mods)) %>%
    pull(stuff) %>%
    unlist() %>%
    unique() %>%
    sort() %>%
    rev()
  return(y)
}


NULL



#' Functions for removing modifications
#'
#' This type of preprocessing makes sense prior plotting
#'
#' @param x A \code{data.frame} fData of an `MSnSet` object.
#'
#' @param mod_name A character string specifying the modification to retain
#'
#' @export
#'
#'
retain_one_mod <- function(x, mod_name){

  for(i in 1:nrow(x)){
    y <- x[i,"mods"][[1]]
    if(length(y$mods) > 0){
      if(mod_name %in% y$mod_names){

        idx_to_retain <- y$mod_names == mod_name
        for(e in setdiff(names(y),"mods_str")){
          y[[e]] <- y[[e]][idx_to_retain]
        }
        # update mod_str
        posi_str <- map2_chr(y$mods_left_border, y$mods_right_border, paste, sep="-")
        y$mods_str <- paste(map2_chr(y$mods, posi_str, paste, sep="@"), collapse=", ")
        x[i,"mods"][[1]] <- list(y)
      }else{
        y <- list(mods=character(0), mods_masses=numeric(0),
                  mods_left_border=numeric(0), mods_right_border=numeric(0),
                  mods_str="", AAs=character(0), mod_names=character(0))
        x[i,"mods"][[1]] <- list(y)
      }
    }
  }
  return(x)
}


#' @describeIn retain_one_mod
#'
#' @export
#'
remove_a_mod <- function(x, mod_name){

  for(i in 1:nrow(x)){
    y <- x[i,"mods"][[1]]
    if(length(y$mods) > 0){
      if(mod_name %in% y$mod_names){

        idx_to_retain <- y$mod_names != mod_name
        for(e in setdiff(names(y),"mods_str")){
          y[[e]] <- y[[e]][idx_to_retain]
        }
        # update mod_str
        posi_str <- map2_chr(y$mods_left_border, y$mods_right_border, paste, sep="-")
        y$mods_str <- paste(map2_chr(y$mods, posi_str, paste, sep="@"), collapse=", ")
        x[i,"mods"][[1]] <- list(y)
      }
    }
  }
  return(x)

}


NULL

