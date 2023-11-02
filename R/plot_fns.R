#' Plot all fragments for a given UniProt accession
#'
#' Plots all fragments associated with a given accession. If a gene name is
#' provided a separate plot will be produced for each accession within the
#' specified gene.
#'
#' @param x A \code{data.table} output from the \code{create_mdata} function.
#'
#' @param gene A character string. A plot will be produced for each UniProt
#'   accession mapping to the specified gene. When a \code{gene} is specified
#'   \code{accession} does not also need to be specified.
#'
#' @param accession A character string specifying which UniProt accession will
#'   be plotted.
#'
#' @param turbo A logical value indicating whether the Turbo colormap will be
#'   used.
#'
#' @param viridis A character string indicating which colormap option will be
#'   used. The available options are: "magma" (or "A"), "inferno" (or "B"),
#'   "plasma" (or "C"), "viridis" (or "D") and "cividis" (or "E").
#'
#' @param size An integer specifying the size of the border around each
#'   fragment. The default size is 1.
#'
#' @param color A character string indicating the color of the border around
#'   each fragment. The default is white.
#'
#' @param save_plot A logical value. If TRUE the plot will be saved as a PNG
#'   file. When a plot is saved to a file it is not also displayed in R/RStudio.
#'
#' @param file_path The path to the folder where each plot will be saved. This
#'   must be specified if \code{save_plot} is TRUE.
#'
#' @return A plot of the fragments for a given accession. When \code{save_plot}
#'   is TRUE a PNG file is created for each plot.
#'
#' @export
#'
#' @author Evan A Martin
#'
plot_fragments <- function (x, gene = NULL, accession = NULL,
                            turbo = TRUE, viridis = NULL,
                            size = 1, color = "white",
                            save_plot = FALSE, file_path = NULL) {

  # Throw an error if both gene and accession are NULL.
  if (is.null(gene) && is.null(accession)) {

    stop ("Either gene or accession must be specified.")

  }

  # Make sure either turbo or viridis is specified.
  if (!turbo && is.null(viridis)) {

    stop (paste("Either turbo needs to be TRUE or viridis needs to be",
                "specified with one of the available options.",
                sep = " "))

  }

  # Check if viridis is one of the approved options.
  if (!is.null(viridis)) {

    if (!(viridis %in% c("magma", "inferno", "plasma", "viridis", "cividis",
                         "A", "B", "C", "D", "E"))) {

      stop ("The input to viridis is not one of the available options.")

    }

  }

  # Throw an error if save_plot is TRUE and file_path is not provided.
  if (save_plot && is.null(file_path)) {

    stop ("If save_plot = TRUE file_path must also be specified.")

  }

  # Plot by the accession provided if gene is NULL.
  if (is.null(gene)) {

    plot_accession(x = x,
                   accession = accession,
                   turbo = turbo,
                   viridis = viridis,
                   size = size,
                   color = color,
                   save_plot = save_plot,
                   file_path = file_path)

    # Plot all accessions for a given gene.
  } else {

    unique_acc <- x %>%
      dplyr::filter(Gene == gene) %>%
      dplyr::pull(UniProtAcc) %>%
      unique()

    # Loop through each accession corresponding to the given gene and plot the
    # associated fragments.
    for (e in 1:length(unique_acc)) {

      p <- plot_accession(x = x,
                          accession = unique_acc[[e]],
                          turbo = turbo,
                          viridis = viridis,
                          size = size,
                          color = color,
                          save_plot = save_plot,
                          file_path = file_path)

      if (!save_plot) {

        print (p)

      }

    }

  }

}

# @author Vlad Petyuk
plot_accession <- function (x, accession, turbo, viridis, size,
                            color, save_plot, file_path) {

  prot_len <- x %>%
    dplyr::filter(UniProtAcc == accession) %>%
    dplyr::slice(1) %>%
    dplyr::pull(protLength) %>%
    as.numeric()

  prot_name <- x %>%
    dplyr::filter(UniProtAcc == accession) %>%
    dplyr::slice(1) %>%
    dplyr::select(Gene, UniProtAcc, AnnType) %>%
    paste0(collapse = ", ")

  prot <- x %>%
    dplyr::filter(UniProtAcc == accession) %>%
    dplyr::group_by(firstAA, lastAA) %>%
    dplyr::summarize(n = sum(spectralCount)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Length = lastAA - firstAA + 1) %>%
    dplyr::arrange(firstAA, -Length, -n)


  # setting staggered ymin
  min_y <- 0.1
  step_y <- 0.033
  width_y <- 0.025

  prot$ymin <- min_y

  if (nrow(prot) > 1) {

    for (i in 2:nrow(prot)) {

      current_y <- min_y

      while (TRUE) {

        # is there a conflict
        max_last_residue <- prot %>%
          dplyr::slice(1:(i-1)) %>%
          dplyr::filter(ymin == current_y) %>%
          dplyr::pull(lastAA) %>%
          max()

        if (max_last_residue + 0 >= prot[i, "firstAA", drop = TRUE]) {

          current_y <- current_y + step_y

        } else {

          break()

        }

      }

      prot[i,"ymin"] <- current_y

    }

  }

  prot$ymax <- prot$ymin + width_y

  p <- ggplot2::ggplot(data = prot) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0,
                                    xmax = prot_len + 2,
                                    ymin = -0.04,
                                    ymax = 0.04)) +
    ggplot2::geom_rect(ggplot2::aes(xmin = firstAA - 0.5,
                                    xmax = lastAA + 0.5,
                                    ymin = ymin,
                                    ymax = ymax,
                                    fill = n),
                       size = size,
                       color = color)

  # Use turbo colors.
  if (turbo) {

    n_colors <- dplyr::n_distinct(prot$n)

    # use gradient if there is only one color.
    if (n_colors == 1) {

      p <- p +
        ggplot2::scale_fill_gradient(low = turbo(2)[[1]],
                                     high = turbo(2)[[2]])

      # Use gradientn when there is more than one color.
    } else if (n_colors > 1){

      p <- p +
        ggplot2::scale_fill_gradientn(colors = turbo(n_colors))

    }

    # Use viridis colors if turbo is FALSE.
  } else if (!turbo && !is.null(viridis)) {

    p <- p +
      ggplot2::scale_fill_viridis_c(option = viridis)

  }

  p <- p +
    ggplot2::ylab(NULL) +
    ggplot2::xlab("residue") +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.ticks.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.line.y = ggplot2::element_blank(),
                   axis.line.x = ggplot2::element_blank()) +
    ggplot2::scale_x_continuous(breaks = seq(0, prot_len, 20)) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 16)
    ) +
    ggplot2::ggtitle(prot_name)

  if (max(prot$ymax) < 0.5)
    p <- p + ggplot2::ylim(-0.04, 0.5)

  file_name <- gsub(", ", "_", prot_name)

  if (save_plot) {

    ggplot2::ggsave(filename = paste0(file_name,".png"),
                    plot = p,
                    path = file_path)

  } else {

    return (p)

  }

}

#' ...
#'
#' ...
#'
#' @param x A \code{data.table} output from the \code{create_mdata} function.
#'
#' @param accession ...
#'
#' @param fill_by ...
#'
#' @param name_from ...
#'
#' @param mods_to_name ...
#'
#' @param border_color ...
#'
#' @param border_size ...
#'
#' @param aa_step ...
#'
#' @param na_symbol ...
#'
#' @param save_plot ...
#'
#' @param file_path ...
#'
#' @return ...
#'
#' @export
#'
#' @author Vlad Petyuk
#'
plot_accession_ptm <- function (x,
                                accession,
                                fill_by = "spectralCount",
                                turbo = TRUE,
                                viridis = NULL,
                                name_from = c("Gene", "UniProtAcc"),
                                mods_to_name = "all",
                                border_color = "white",
                                border_size = NULL,
                                aa_step = 20,
                                na_symbol = " ",
                                save_plot = FALSE,
                                file_path = NULL) {

  # Make sure fill_by is a variable in the data.
  if (!fill_by %in% names(x)) {

    stop ("The variable specified in fill_by is not present in x.")

  }

  # Combine the Gene and pcGroup into one variable. This is used later when
  # producing the plot with the PTMs.
  x <- x %>%
    dplyr::mutate(feature_name = paste(Gene, pcGroup, sep = "_"))

  # sort out mod_to_name
  if(length(mods_to_name) == 1 && grepl("top\\d+", mods_to_name)){
    mods <- get_mods_counts(x, accession)
    top_n_mods <- as.numeric(sub("top", "", mods_to_name))
    top_n_mods <- min(top_n_mods, 9) # we'll restrict with 9
    mods_to_name <- names(mods)[seq_len(min(length(mods), top_n_mods))]
  }

  prot_len <- dplyr::filter(x, UniProtAcc == accession) %>%
    dplyr::distinct(protLength) %>%
    as.numeric()

  if(is.null(border_size)){
    border_size <- 1.5/log10(prot_len)
  }

  prot_name_cols <- x %>%
    dplyr::filter(UniProtAcc == accession) %>%
    dplyr::distinct(!!! rlang::syms(name_from))

  prot_name <- paste0(prot_name_cols, collapse = ", ")

  # if(nrow(prot_name_cols) == 1){
  #   prot_name <- paste0(prot_name_cols, collapse = ", ")
  # }else{
  #   stop ("name_from isn't unique!")
  # }

  prot <- x %>%
    dplyr::filter(UniProtAcc == accession) %>%
    dplyr::mutate(Length = lastAA - firstAA + 1) %>%
    dplyr::arrange(firstAA, -Length)

  # setting staggered ymin
  min_y <- 0.05
  step_y <- 0.033
  width_y <- 0.025

  prot$ymin <- min_y
  if(nrow(prot) > 1){
    for(i in 2:nrow(prot)){
      current_y <- min_y
      while(TRUE){
        # is there a conflict
        max_last_residue <- prot %>%
          dplyr::slice(1:(i-1)) %>%
          dplyr::filter(ymin == current_y) %>%
          dplyr::pull(lastAA) %>%
          max()
        if(max_last_residue + 0 >= prot[i, "firstAA", drop = TRUE]){
          current_y <- current_y + step_y
        }else{
          break()
        }
      }
      prot[i, "ymin"] <- current_y
    }
  }

  prot$ymax <- prot$ymin + width_y

  y_height <- prot$ymax[[1]] - prot$ymin[[1]]

  p <- ggplot2::ggplot(data = prot) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 1 - 0.5,
                                    xmax = prot_len + 0.5,
                                    ymin = 0,
                                    ymax = y_height)) +
    ggplot2::geom_rect(ggplot2::aes(xmin = firstAA - 0.5,
                                    xmax = lastAA + 0.5,
                                    ymin = ymin,
                                    ymax = ymax,
                                    fill = !!rlang::sym(fill_by)),
                       color = border_color,
                       size = border_size)

  # Use turbo colors.
  if (turbo) {

    n_colors <- dplyr::n_distinct(prot[, fill_by])

    # use gradient if there is only one color.
    if (n_colors == 1) {

      p <- p +
        ggplot2::scale_fill_gradient(low = turbo(2)[[1]],
                                     high = turbo(2)[[2]])

      # Use gradientn when there is more than one color.
    } else if (n_colors > 1){

      p <- p +
        ggplot2::scale_fill_gradientn(colors = turbo(n_colors))

    }

    # Use viridis colors if turbo is FALSE.
  } else if (!turbo && !is.null(viridis)) {

    p <- p +
      ggplot2::scale_fill_viridis_c(option = viridis)

  }

  prot_mods <- list()
  for(i in 1:nrow(prot)){
    mod_i <- as.data.frame(prot$mods[[i]][c("mod_names",
                                            "mods_left_border",
                                            "mods_right_border")])
    if(nrow(mod_i) > 0)
      mod_i$feature_name <- prot$feature_name[i]
    prot_mods <- c(prot_mods, list(mod_i))
  }

  prot_mods <- dplyr::bind_rows(prot_mods) %>%
    dplyr::mutate(label = mod_names)

  # plotting mods piece
  if (!(identical(mods_to_name, "none") || nrow(prot_mods) == 0)) {

    prot_mods_full <- dplyr::inner_join(prot,
                                        prot_mods,
                                        by = "feature_name") %>%
      dplyr::mutate(
        x = (mods_right_border + mods_left_border) / 2 + firstAA - 1
      ) %>%
      dplyr::mutate(y = (ymin + ymax) / 2)

    if (identical(mods_to_name, "all")) {
      p <- p +
        ggplot2::geom_label(
          ggplot2::aes(x = x,
                       y = y,
                       label = label),
          data = prot_mods_full,
          fill = ggplot2::alpha("white", alpha = 0.9),
          show.legend = FALSE
        )
    }else{
      prot_mods_full <- prot_mods_full %>%
        dplyr::mutate(
          label_selected = dplyr::case_when(
            label %in% mods_to_name ~ label,
            TRUE ~ "Other"
          )
        )

      # updates mods_to_name
      mods_to_name <- seq_along(mods_to_name) %>%
        stats::setNames(mods_to_name)

      # update label here
      prot_mods_full <- prot_mods_full %>%
        dplyr::mutate(label_short = mods_to_name[mod_names]) %>%
        dplyr::mutate(label_short = as.character(label_short)) %>%
        dplyr::mutate(
          label_short = dplyr::case_when(
            is.na(label_short) ~ na_symbol,
            TRUE ~ label_short
          )
        )

      namer <- prot_mods_full %>%
        dplyr::distinct(label_short, label_selected) %>%
        {temp <- .$label_short; names(temp) <- .$label_selected;temp}
      namer <- sort(namer, na.last = NA)
      if (namer[1] == na_symbol) {
        namer <- c(namer[-1], namer[1]) # putting blank space last
      }

      p <- p +
        ggplot2::geom_label(ggplot2::aes(x = x, y = y),
                            label = " ",
                            data = prot_mods_full,
                            fill = ggplot2::alpha("white", alpha = 0.9),
                            show.legend = FALSE) +
        ggplot2::geom_point(
          ggplot2::aes(x = x, y = y, shape = label_selected),
          size = 4,
          data = prot_mods_full
        ) +
        ggplot2::scale_shape_manual(values = namer, name = "modification")
    }
  }

  p <- p +
    ggplot2::ylab(NULL) +
    ggplot2::xlab("residue") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_blank(),
      legend.key = ggplot2::element_rect(color = "black")
    ) +
    ggplot2::scale_x_continuous(
      breaks = c(1, seq(aa_step, prot_len, aa_step)),
      # breaks = floor(seq(from = 1, to = prot_len, length.out = aa_step)),
      expand = expansion(mult = c(0.03, 0)),
      limits = c(1 - 0.5, prot_len + 0.5)
    ) +
    ggplot2::scale_y_continuous(
      expand = c(0, 0),
      limits = if (max(prot$ymax) < 0.5) c(0, 0.5) else NULL
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 16)
    ) +
    ggplot2::ggtitle(prot_name)

  file_name <- gsub(", ", "_", prot_name)

  if (save_plot) {

    ggplot2::ggsave(filename = paste0(file_name,".png"),
                    plot = p,
                    path = file_path)

  } else {

    return (p)

  }

}




#' @export
#'
#' @importFrom ggplot2 scale_fill_manual
#'
plot_pform_sig <- function(x, accession, fill_by, ...){

  cuts <- c(0,0.0001, 0.001, 0.01, 0.05, 1)
  x[[fill_by]] <- cut(x[[fill_by]], cuts, ordered_result=T)

  fillings <- scale_fill_manual(
    breaks = levels(x[[fill_by]]),
    values = c("#E53935", "#FDD835", "#43A047", "#1E88E5", "#BDBDBD"),
    limits = levels(x[[fill_by]]),
    labels = c("< 1e-04", "(1e-04, 1e-03]", "(1e-03, 0.01]",
               "(0.01, 0.05]", "> 0.05"),
  )

  p <- plot_accession_ptm(x, accession, fill_by, ...)
  p <- p + fillings
  return(p)
}



# library(scales)
#' @export
#'
#' @importFrom scales trans_new pretty_breaks log_breaks
#' @importFrom ggplot2 scale_fill_gradientn
#'
plot_pform_sig_soft <- function(x, accession, fill_by, ...){

  p <- plot_accession_ptm(x, accession, fill_by, ...)

  log10_rev_trans <- trans_new("log10_rev",
                               transform = function(x) -log10(x),
                               inverse = function(x) 10 ^ (-x),
                               breaks = function(x) rev(log_breaks(10)(x)))

  breaks <- min(p$data[[fill_by]], na.rm = TRUE) %>%
    {signif(10 ^ (pretty_breaks()(0 : floor(log10(.) * 2)) / 2), 1)}

  labels <- ifelse(breaks > 1e-3,
                   format(breaks, scientific = FALSE, drop0trailing = TRUE),
                   format(breaks))

  p <- p + scale_fill_gradientn(trans=log10_rev_trans, colours=MSnSet.utils::jet2.colors(21), # hack!
                                breaks=breaks, labels=breaks, limits=c(1, min(breaks)))
  return(p)
}










# @author Vlad Petyuk
get_mods_counts <- function(x, acc){
  y <- dplyr::filter(x, UniProtAcc == acc) %>%
    dplyr::mutate(stuff = purrr::map(mods, `$`, mod_names)) %>%
    dplyr::pull(stuff) %>%
    unlist() %>%
    table() %>%
    sort() %>%
    rev()
  return(y)
}

# Spectral count by accession --------------------------------------------------

#' ...
#'
#' ...
#'
#' @param x The x_meta data table (the output from create_mdata)
#'
#' @param gene A character string indicating the accession to plot.
#'
#' @param ptm_annotation A data frame that contains the mass of each PTM along
#'   with its name. In the case when a PTM does not have a name the PTM's mass,
#'   rounded to three decimal places, is used.
#'
#' @export
#'
data_pf <- function (x, accession, ptm_annotation) {

  # Filter x by the input to type ----------------------------------------------

  # Filter the data by accession.
  x <- x %>%
    dplyr::ungroup() %>%
    dplyr::filter(UniProtAcc == accession)


  # Extract PTMs ---------------------------------------------------------------

  unique_seq <- x %>%
    dplyr::distinct(Proteoform) %>%
    dplyr::pull(Proteoform)

  # Create a list to hold the data frame output by id_ptm. The data frames will
  # be combined into one and used for plotting the PTMs present.
  the_list <- vector(mode = "list",
                     length = length(unique_seq))
  ptm_info <- vector(mode = "list",
                     length = length(unique_seq))
  prot_info <- vector(mode = "list",
                      length = length(unique_seq))

  for (e in 1:length(unique_seq)) {

    ptm_info[[e]] <- id_ptm(unique_seq[[e]],
                            ptm_annotation = ptm_annotation)
    prot_info[[e]] <- id_prot(x = x, pf = unique_seq[[e]])
    the_list[[e]] <- cbind(ptm_info[[e]], prot_info[[e]])

  }

  the_df <- rbindlist(the_list)

  # Nab start/end points for protein and fragments -----------------------------

  # Determine the starting point of the horizontal line that will represent the
  # protein.
  if (min(the_df$firstAA) <= 5) {
    the_df$startProt <- 0
  } else {
    the_df$startProt <- min(the_df$firstAA) - 5
  }

  # Determine the ending point of the horizontal line that will represent the
  # protein.
  if (max(the_df$protLength) - max(x$lastAA) <= 5) {
    the_df$endProt <- max(x$protLength)
  } else {
    the_df$endProt <- max(x$lastAA) + 5
  }

  # Remove any rows that contain NA in the ptm_id column. These rows should not
  # be plotted.
  the_df <- the_df %>%
    dplyr::filter(!is.na(ptm_id))

  # Sum spectral counts within ptm_id at the same location. Each PTM (at the
  # same location) should only be plotted once.
  the_df <- the_df %>%
    dplyr::group_by(ptm_id, location) %>%
    dplyr::mutate(spectralCount = sum(spectralCount))

  # Add the accession to the data frame as an attribute. This will be used by
  # the plotting function to add the main title.
  attr(the_df, "accession") <- accession

  return (the_df)

}

#' Plot PTMs with associated spectral counts
#'
#' ...
#'
#' @param x A data frame output by \code{data_pf}.
#'
#' @export
plot_acc <- function (x) {

  # Plot PTMs ------------------------------------------------------------------

  p <- ggplot2::ggplot(data = x) +
    ggplot2::geom_point(ggplot2::aes(x = location + firstAA - 1,
                                     y = spectralCount,
                                     color = ptm_id),
                        shape = 17,
                        size = 3) +
    ggplot2::geom_segment(ggplot2::aes(x = startProt - 0.5,
                                       y = 0,
                                       xend = endProt + 0.5,
                                       yend = 0),
                          size = 1) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.line = ggplot2::element_line(colour = "black"),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 16)
    ) +
    ggplot2::scale_x_continuous(breaks = seq(min(x$firstAA),
                                             max(x$protLength),
                                             10)) +
    ggplot2::coord_cartesian(ylim = c(0, max(x$spectralCount))) +
    ggplot2::xlab("Location") +
    ggplot2::ylab("Spectral Count") +
    ggplot2::ggtitle(attr(x, "accession"))

  return (p)

}

# pf: A character string from the Proteoform column.
id_ptm <- function (pf, ptm_annotation) {

  # Extract all PTMs in the current proteoform. The identified or named PTMs
  # will be converted to a mass and the unidentified PTMs will be left alone.
  ptm <- pf %>%
    stringr::str_extract_all("\\[.+?\\]", simplify = TRUE) %>%
    stringr::str_remove_all("\\[|\\]")

  # If there are no PTMs return a data frame of NAs.
  if (length(ptm) == 0) {

    return (data.frame(ptm_id = NA,
                       location = NA,
                       length = NA))

  }

  # Remove any amino acids at the beginning and end of the sequence along with
  # the period separating them from the rest of the sequence. This is needed to
  # determine the location in the sequence where the first ( occurs.
  scs <- pf %>%
    stringr::str_remove("^.??\\.") %>%
    stringr::str_remove("\\..??$")

  # Locate the first (. This will be used to determine if the acetylation PTM is
  # on the N-terminus
  paren_lctn <- stringr::str_locate(scs, "\\(")[[1]]

  # and remove any
  # parentheses.
  scs <- scs %>%
    stringr::str_remove_all("\\(|\\)")

  ptm_name <- vector(length = length(ptm))

  # Loop through all PTMs and create symbol IDs based on the PTM's name or mass.
  for (e in 1:length(ptm)) {

    # Check if the PTM is named.
    if (stringr::str_detect(ptm[[e]], "[:alpha:]")) {

      # Check if the PTM is acetyl. Its name will change if the PTM position
      # includes the first amino acid.
      if (e == 1 && ptm[[e]] == "Acetyl" && paren_lctn == 1) {

        ptm_name[[e]] <- "N-Acetyl"

      } else {

        ptm_name[[e]] <- ptm[[e]]

      }

      # The following runs for unidentified modifications (there is just a
      # mass in the square brackets).
    } else {

      temp_mass <- round(as.numeric(ptm[[e]]), 3)

      ptm_name[[e]] <- to_name(mass = temp_mass,
                               ptm_annotation = ptm_annotation)

    }

  }

  ptm_lctn <- vector(mode = "numeric",
                     length = length(ptm))

  # Loop through all PTMs and find their location in the sequence.
  for (e in 1:length(ptm)) {

    # Find the location of the first PTM. Subtract one from the output of
    # str_locate because we find the PTM by searching for the first occurrence
    # of "[". The amino acid where the PTM occurs is one place before the "[".
    ptm_lctn[[e]] <- stringr::str_locate(scs, "\\[")[[1]] - 1

    # Remove PTM. The PTM needs to be removed from the proteoform so the next
    # PTM in line will show up at the correct location in the amino acid
    # sequence.
    scs <- stringr::str_remove(scs, "\\[.+?\\]")

  }

  return (
    data.frame(ptm_id = ptm_name,
               location = ptm_lctn,
               length = stringr::str_length(scs))
  )

}

# Converts mass to a name using the file from Vlad.
to_name <- function (mass, ptm_annotation) {

  temp_id <- ptm_annotation %>%
    dplyr::filter(ptm_weight == mass) %>%
    dplyr::slice(1) %>%
    dplyr::pull(mod_id)

  if (length(temp_id) == 0) temp_id <- "0"

  return (temp_id)

}

id_prot <- function (x, pf) {

  return (
    x %>%
      dplyr::ungroup() %>%
      dplyr::filter(Proteoform == pf) %>%
      dplyr::select(firstAA, lastAA, protLength, Gene, pcGroup, spectralCount)
  )

}

# Turbo color map --------------------------------------------------------------

# The following functions are from the GitHub Gist jlmelville/turbo_colormap.R
# The matrix turbo_colormap_data also comes from this gist.

# Better, shorter, more idiomatic, vectorized versions of the interpolate and
# turbo functions thanks to @onesandzeroes
interpolate <- function (colormap, x) {
  x <- pmax(0, pmin(1, x))
  a <- floor(x * 255)
  b <- pmin(255, a + 1)
  f <- x * 255 - a
  a <- a + 1
  b <- b + 1

  colormap[a, ] + (colormap[b, ] - colormap[a, ]) * f
}

turbo <- function (n, start = 0, end = 1) {
  xs <- seq.int(from = start, to = end, length.out = n)
  grDevices::rgb(interpolate(turbo_colormap_data, xs))
}
