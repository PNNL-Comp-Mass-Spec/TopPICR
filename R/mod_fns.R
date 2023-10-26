# Create mod data --------------------------------------------------------------

#' @export
create_mod_data <- function (mod_file, mod_path, use_unimod = FALSE) {

  # Create a data frame for carbon 13.
  isotopic_errors <- tibble::tribble(
    ~name, ~mass,
    "+1C13", 1.003355,
    "+2C13", 2.00671,
    "+3C13", 3.010065,
    "-1C13", -1.003355,
    "-2C13", -2.00671,
    "-3C13", -3.010065
  ) %>%
    dplyr::select(mass, name)

  # Determine what to return based on use_unimod argument.
  if (use_unimod) {

    # Read in mod information output by TopPIC.
    toppic_mods <- readr::read_csv(
      file.path(mod_path, mod_file),
      comment = "#",
      col_names = FALSE
    ) %>%
      dplyr::select(2:1) %>%
      `colnames<-`(c("mass", "name"))

    is_present <- rep(FALSE, nrow(unimods))

    # Determine what mods from TopPIC are also in unimod.
    for(i in 1:nrow(unimods)){
      mass_diff <- abs(unimods$mass[i] - toppic_mods$mass)
      if(any(mass_diff < 0.001)){
        is_present[i] <- TRUE
      }
    }

    # Only keep mods not found in TopPIC mods.
    unimods <- unimods[!is_present, ]

    # leave only one out of duplicated masses
    unimods <- unimods %>%
      dplyr::group_by(mass) %>%
      dplyr::slice_min(name)
    # What does the min function do on character vectors? Is this the behavior
    # you want for slice_min?

    unimods <- dplyr::bind_rows(unimods, isotopic_errors)

    return (unimods)

  } else {

    return (isotopic_errors)

  }

}

# Update mod data and combine with metadata ------------------------------------

#' @export
add_mods <- function (x,
                      mods,
                      nterm_tol = 3,
                      acetyl_id = "Acetyl",
                      centroid_tol = 0.3,
                      matching_tol = 0.1,
                      matching_tol_self = .Machine$double.eps) {

  # Extract mods from the Proteoform column. The output from this step will be a
  # list. The information in this list will be used/modified by the remaining
  # function calls.
  x2 <- x %>%
    dplyr::mutate(mods = purrr::map(Proteoform, extract_mods))

  # Convert acetylation to N-terminal acetylation if the following conditions
  # are met: 1. if near N-terminus of the protein AND 2. on the N-terminus of
  # the peptide--rename N-Acetyl.
  x2 <- annotate_Nterm_acetyls(x = x2,
                               nterm_tol = nterm_tol,
                               acetyl_id = acetyl_id)

  # Cluster masses of mods and match to unimod names.
  ann_table <- get_mass_annotation_table(x = x2,
                                         unimods = mods,
                                         centroid_tol = centroid_tol,
                                         matching_tol = matching_tol)

  # If the mass matches a unimod mass use the unimod name. Otherwise, use the
  # rounded centroid mass as the name.
  x2 <- x2 %>%
    dplyr::mutate(mods = purrr::map(mods,
                                    annotate_masses,
                                    mass_annotation_table = ann_table,
                                    matching_tol_self = matching_tol_self))

  return (x2)

}

extract_mods <- function (pform) {

  # add N- and C-term flanking character to make them the same
  pform <- sub("(^\\.)(.*)", "-\\1\\2", pform)
  pform <- sub("(.*)(\\.$)", "\\1\\2-", pform)

  # extract mod notation. Whatever is in []
  mods_masses <- stringr::str_extract_all(
    pform, "(?<=\\[)[^\\]\\[]*(?=\\])", simplify = TRUE
  ) # extracting what is in [...]
  mods_masses <- as.character(mods_masses)
  mods <- mods_masses

  suppressWarnings(
    mods_masses <- as.numeric(mods_masses)
  )

  # extracting modification location borders
  s <- pform
  s <- gsub("\\[.+?\\]", "", s)
  s <- sub(".\\.(.+?)\\..", "\\1", s) # note, this requires flanking AAs

  # the core loop of extracting borders
  # a horrible piece of code that needs to be replaced on a simpler one
  s_chrs <- substring(s, 1:nchar(s), 1:nchar(s))
  stk <- c()
  count_openings <- 0
  aa_idx <- 1
  # left_borders <- c()
  # right_borders <- c()
  # closing_order <- c()
  left_borders <- numeric(0)
  right_borders <- numeric(0)
  closing_order <- numeric(0)
  for(i in 1:length(s_chrs)){
    # sort out parenthesis vs AA
    if(s_chrs[i] %in% c("(",")")){
      if(s_chrs[i] == "("){
        count_openings <- count_openings + 1 # this counter goes one way
        stk <- c(stk, count_openings)
        left_borders <- c(left_borders, aa_idx)
      }else{
        closing_order <- c(closing_order, stk[length(stk)])
        stk <- stk[-length(stk)]
        right_borders <- c(right_borders, aa_idx - 1)
      }
    }else{
      aa_idx <- aa_idx + 1
    }
  }
  left_borders <- left_borders[closing_order]

  # extract AAs
  s <- gsub("[()]","",s)
  aas <- character(0)
  if(length(left_borders) > 0)
    aas <- substring(s, left_borders, right_borders)

  posi_str <- purrr::map2_chr(left_borders, right_borders, paste, sep = "-")
  mods_str <- paste(
    purrr::map2_chr(mods, posi_str, paste, sep = "@"), collapse = ", "
  )

  out <- list(mods, mods_masses, left_borders, right_borders, mods_str, aas)
  names(out) <- c("mods", "mods_masses", "mods_left_border",
                  "mods_right_border", "mods_str", "AAs")

  return (out)

}

# @export
#
# This piece of code isn't working. Will delete with the next clean-up.
# Leaving it for "just in case".
#
# annotate_Nterm_acetyls <- function (x,
#                                     nterm_tol = 3,
#                                     acetyl_id = "Acetyl") {
#
#   for (i in 1:nrow(x)) {
#
#     # Determine if an acetyl could be an n-terminus acetyl.
#     if (x[i, "firstAA"] <= nterm_tol) {
#
#
#       # Check if any mods are present.
#       if (length(x[i, "mods"][[1]][[1]]$mods) > 0) {
#
#         # NOTE: For the rest of the if statement we will work with just the
#         # first element of mods, mods_left_border, and mods_str because there
#         # could be multiple mods. We are only interested in the first mod if it
#         # is Acetyl.
#
#         # Check if the mod is acetyl and if it is on the first amino acid.
#         if (x[i, "mods"][[1]][[1]]$mods[1] == acetyl_id &&
#             x[i, "mods"][[1]][[1]]$mods_left_border[1] == 1) {
#
#           # Change from Acetyl to N-Acetyl because the current acetyl occurs on
#           # the n-terminus.
#           x[i, "mods"][[1]][[1]]$mods[1] <- paste("N-", acetyl_id, sep="")
#
#           # Update the name of the mod in the mods_str vector. This will become
#           # N-Acetyl@<start position>-<end position>
#           x[i, "mods"][[1]][[1]]$mods_str[1] <- paste(
#             "N-", x[i, "mods"][[1]][[1]]$mods_str[1], sep = ""
#           )
#
#         }
#
#       }
#
#     }
#
#   }
#
#   return (x)
#
# }



#' @export
annotate_Nterm_acetyls <- function(x, nterm_tol = 3, acetyl_id = "Acetyl"){

  temp_mod_names <- "mods"
  if("mod_names" %in% names(x$mods[[1]]))
    temp_mod_names <- "mod_names"

  for(i in 1:nrow(x)){
    if(x[i,"firstAA"] <= nterm_tol){
      y <- x[i,"mods"][[1]]
      if(length(y$mods) > 0){
        if(y$mods[1] == acetyl_id & y$mods_left_border[1] == 1){
          y$mods[1] <- paste("N-", acetyl_id, sep="")
          y[[temp_mod_names]][1] <- y$mods[1]
          posi_str <- map2_chr(y$mods_left_border, y$mods_right_border, paste, sep="-")
          y$mods_str <- paste(map2_chr(y[[temp_mod_names]], posi_str, paste, sep="@"), collapse=", ")
          x[i,"mods"][[1]] <- list(y)
        }
      }
    }
  }
  return(x)
}







#' @export
get_mass_annotation_table <- function (x,
                                       unimods,
                                       centroid_tol = 0.3,
                                       matching_tol = 0.1) {
  suppressWarnings(
    mods_masses <- purrr::map(x$mods, ~ .x$mods_masses) %>%
      unlist() %>%
      as.numeric() %>%
      purrr::keep(~!is.na(.x))
  )

  d <- stats::density(mods_masses, n = 2^17, bw = 0.0002)
  dd <- data.frame(x=d$x, y=d$y)
  out <- tryCatch (
    .centroiding(dd, centroid_tol, numeric()),
    error = function (e) NULL
  )

  if (is.null(out)) {

    stop ("Try increasing the value of centroid_tol.")

  }

  # group PTMs
  mod_mass_centroid <- rep(NA, length(mods_masses))
  for(i in seq_along(out)){
    idx <- abs(mods_masses - out[i]) < centroid_tol
    mod_mass_centroid[idx] <- out[i]
  }

  # match to UniMod
  mod_name <- purrr::map_chr(mod_mass_centroid,
                             .get_matching_mod,
                             unimods,
                             matching_tol)

  return(data.frame(mod_mass = mods_masses,
                    mod_mass_centroid,
                    mod_name))
}

# # recursive centroiding
# .centroiding <- function(x, centroid_tol, cs){
#   if(nrow(x) == 0){
#     return(cs)
#   }else{
#     i <- which.max(x$y)
#     ci <- x$x[i] # centroid
#     idxs <- abs(x$x - ci) < centroid_tol
#     x <- x[!idxs,]
#     .centroiding(x, centroid_tol, c(cs, ci))
#   }
# }
#

# iterative centroiding
.centroiding <- function(x, centroid_tol, cs){
  # cs <- numeric(0)
  while(nrow(x) > 0){
    i <- which.max(x$y)
    ci <- x$x[i] # centroid
    idxs <- abs(x$x - ci) < centroid_tol
    x <- x[!idxs,]
    cs <- c(cs, ci)
  }
  return(cs)
}



.get_matching_mod <- function(x, reference_table, matching_tol){
  d <- reference_table$mass - x
  i <- which.min(abs(d))
  if(abs(d[i]) < matching_tol)
    return(reference_table$name[i])
  else
    return(NA_character_)
}

#' @export
annotate_masses <- function (x,
                             mass_annotation_table,
                             matching_tol_self = .Machine$double.eps) {

  for (i in seq_along(x)) {

    names <- purrr::map_chr(x$mods_masses,
                            .annotate_one_mass,
                            mass_annotation_table = mass_annotation_table,
                            matching_tol_self = matching_tol_self)
    chr_idx <- is.na(names)
    names[chr_idx] <- x$mods[chr_idx]
    x$mod_names <- names

  }

  posi_str <- purrr::map2_chr(x$mods_left_border,
                              x$mods_right_border,
                              paste,
                              sep = "-")

  x$mods_str <- paste(purrr::map2_chr(x$mod_names,
                                      posi_str,
                                      paste,
                                      sep = "@"),
                      collapse = ", ")

  return (x)

}

# applied to x_meta/fData
.annotate_one_mass <- function (x,
                                mass_annotation_table,
                                matching_tol_self) {
  name <- NA_character_
  if(!is.na(x)){
    d <- mass_annotation_table$mod_mass - x
    i <- which.min(abs(d))
    if(abs(d[i]) < matching_tol_self)
      name <- ifelse(
        is.na(mass_annotation_table$mod_name[i]),
        as.character(round(mass_annotation_table$mod_mass_centroid[i], 3)),
        mass_annotation_table$mod_name[i]
      )
    else
      message("missed in the annotation table")
  }
  return(name)
}

# Modify metadata object -------------------------------------------------------

# @author Vlad Petyuk
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
        posi_str <- purrr::map2_chr(y$mods_left_border,
                                    y$mods_right_border,
                                    paste,
                                    sep="-")
        y$mods_str <- paste(purrr::map2_chr(y$mods,
                                            posi_str,
                                            paste,
                                            sep="@"),
                            collapse=", ")
        x[i,"mods"][[1]] <- list(y)
      }
    }
  }
  return(x)

}

# @author Vlad Petyuk
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
        posi_str <- purrr::map2_chr(y$mods_left_border,
                                    y$mods_right_border,
                                    paste,
                                    sep="-")
        y$mods_str <- paste(purrr::map2_chr(y$mods,
                                            posi_str,
                                            paste,
                                            sep="@"),
                            collapse=", ")
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
