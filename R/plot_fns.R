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

  } else if (!is.null(viridis)) {

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
