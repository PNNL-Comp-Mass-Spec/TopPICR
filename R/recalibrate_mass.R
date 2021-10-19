# recalibrate_mass preliminary functions ---------------------------------------

#' Calculate ppm and retention time error
#'
#' Computes the error in ppm and retention time between each data set and a
#' reference data set. The default reference data set (returned by the
#' \code{find_ref_ds} function) is the one with the most unique proteoforms.
#'
#' @param x A \code{data.table} output from the \code{align_rt} function.
#'
#' @param ref_ds A character string specifying the name of the reference data
#'   set.
#'
#' @return A list with three elements. The first element is the data matrix that
#'   contains the error terms for the mass and retention times. These terms will
#'   be used to both recalibrate the mass and cluster the data. The second
#'   element is the median absolute deviation (MAD) of the error in mass between
#'   the reference and non-reference data sets. This is reported in ppm. The
#'   third element is the MAD of the retention time error between the reference
#'   data set and all other data sets.
#'
#' @importFrom magrittr %>%
#'
#' @export
#'
calc_error <- function (x, ref_ds) {

  # Compute the mass error in ppm by data set. This is the difference between
  # the `Precursor mass` and the `Adjusted precursor mass` converted to parts
  # per million. Also computes the median and standard deviation of the ppm
  # errors (by Dataset).
  ppm_stats <- x %>%
    dplyr::filter(`#unexpected modifications` == 0) %>%
    dplyr::mutate(ppm_error = (1e6 * (`Precursor mass` -
                                       `Adjusted precursor mass`) /
                                 `Adjusted precursor mass`)) %>%
    dplyr::group_by(Dataset) %>%
    dplyr::summarize(ppm_sd = stats::mad(ppm_error),
                     ppm_median = stats::median(ppm_error)) %>%
    dplyr::ungroup()

  # Extract the reference data set.
  x_ref <- x %>%
    dplyr::filter(Dataset == ref_ds) %>%
    dplyr::select(RTalign, Proteoform, `Feature intensity`, `Feature apex`)

  # Compute the median absolute deviation of the retention time errors by data
  # set then take the median of the mad values.
  rt_stats <- x %>%
    dplyr::group_by(Dataset) %>%
    dplyr::summarize(rt_sd = calc_rt_stats(rt = RTalign,
                                           pf = Proteoform,
                                           fi = `Feature intensity`,
                                           fa = `Feature apex`,
                                           x_ref = x_ref)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(Dataset != ref_ds) %>%
    dplyr::summarize(the_sd = stats::median(rt_sd)) %>%
    dplyr::pull(the_sd)

  # Return the x_error data frame that contains the repMass variable along
  # with the ppm and retention time errors.
  return (list(
    ppm_error = ppm_stats,
    ppm_sd = stats::median(ppm_stats$ppm_sd),
    rt_sd = rt_stats
  ))

}

# recalibrate_mass main function -----------------------------------------------

#' Recalibrate the mass
#'
#' Calibrates the mass to a reference data set using the error terms returned by
#' the \code{calc_error} function.
#'
#' @param x The \code{data.table} output from the \code{align_rt} function.
#'
#' @param errors A list output from the \code{calc_error} function.
#'
#' @param repMass Logical. If TRUE a representative mass will be calculated from
#'   the \code{RecalMass} variable. Rows will be grouped by \code{Proteoform}.
#'
#' @return A \code{data.table} containing the recalibrated mass variable. The
#'   following variables were added/removed:
#'
#'   | Added             | Removed                    |
#'   | ----------------- | -------------------------- |
#'   | `RecalMass`       |                            |
#'   | `repMass` (when repMass = TRUE) |              |
#'
#' @md
#'
#' @importFrom magrittr %>%
#'
#' @export
#'
recalibrate_mass <- function (x, errors, repMass = TRUE, feature = FALSE) {

  # Check to make sure the data sets in x are also in errors$ppm_stats. If the
  # data sets in each data frame are not identical the mass cannot be
  # recalibrated. Here x is the data.table output by the align_rt function and
  # errors$ppm_stats is the data.table output by the calc_error function.
  if (!setequal(unique(x$Dataset), unique(errors$ppm_error$Dataset))) {

    stop (paste("The data sets present in x and errors$ppm_error do not match."))

  }

  if (feature) {

    mass <- "Mass"

  } else {

    mass <- "Precursor mass"

  }

  # Add the ppm_median variable to the data. This variable will be used to carry
  # out the actual recalibration.
  x_recal <-  dplyr::inner_join(x = x,
                                y = errors$ppm_error,
                                by = "Dataset") %>%
    # Recalibrate the mass.
    dplyr::mutate(
      RecalMass = !!rlang::sym(mass) * (1 - ppm_median / 1e6)
    ) %>%
    # Remove intermediate variables used to recalibrate the mass. These
    # variables are in the output of the calc_error function and can be accessed
    # there if they are needed.
    dplyr::select(-ppm_sd, -ppm_median)

  # Check if the representative mass will be calculated.
  if (repMass) {

    # Compute the robust median of the precursor mass within Dataset and
    # Proteoform. This variable will be used to calculate the error (in ppm)
    # between the reference and all other data sets.
    x_abridged <- x_recal %>%
      dplyr::select(c(Dataset, Proteoform, RecalMass,
                      `Adjusted precursor mass`)) %>%
      dplyr::group_by(Dataset, Proteoform) %>%
      # Calculate a representative mass (repMass) within each Dataset and
      # Proteoform. This mass will be used to recalibrate the mass in the
      # recalibrate_mass function.
      dplyr::summarize(
        repMass = robust_median(mass = RecalMass,
                                adj_mass = `Adjusted precursor mass`)
      )

    # Combine the data frames with the robust median mass and the aligned
    # retention times.
    x_recal <- x_recal %>%
      dplyr::ungroup() %>%
      dplyr::inner_join(x_abridged)

  }

  # Return the x_recal data frame that contains the RecalMass and repMass (when
  # applicable) variables.
  return (x_recal)

}

# calc_error auxiliary functions -----------------------------------------------

# Calculates the median Precursor mass of the mass values that are within 0.5
# Dalton of the Adjusted precursor mass. There are cases when all of the
# Precursor masses are outside this 0.5 Dalton range of the Adjusted precursor
# mass. If this is the case the Precursor mass that is closest to the Adjusted
# precursor mass is returned.
robust_median <- function (mass, adj_mass) {

  # Calculate the absolute difference between the Precursor mass (mass) and the
  # Adjusted precursor mass (adj_mass). This vector will be used to determine
  # which Precursor masses will be used to calculate the representative mass for
  # the Dataset and Proteoform combination.
  diffs <- abs(adj_mass - mass)

  # Only keep mass (Precursor mass) values that are within 0.5 Dalton of
  # adj_mass (Adjusted precursor mass).
  goodies <- which(diffs < 0.5)

  # Check the length of goodies. It is possible that all of the differences are
  # outside the 0.5 threshold. If this is the case we will use the mass that is
  # the closest to the adj_mass.
  if (length(goodies) == 0) {

    # Plus n isotope case.
    idx <- mass - adj_mass > 0.75

    # Subtract carbon12/13 difference from mass.
    mass[idx] <- mass[idx] - 1.003355

    idx <- mass - adj_mass < -0.75

    # Will add an amino group to the mass.
    mass[idx] <- mass[idx] + 0.984016

    # Runs if there is at least one element in goodies.
  } else {

    # Only keep the masses are within the 0.5 Dalton threshold.
    mass <- mass[goodies]

  }

  # Take the median of the mass values that are within 0.5 Dalton of the
  # adj_mass to use as the representative mass.
  return (stats::median(mass, na.rm = TRUE))

}

# rt: A vector of retention times grouped by Dataset.
# pf: A vector of proteoforms grouped by Dataset.
# fi: A vector of feature intensity values. These values are grouped by Dataset.
# fa: A vector of feature apex values. These values are grouped by Dataset.
# x_ref: the x data matrix filtered by the reference data set.
calc_rt_stats <- function (rt, pf, fi, fa, x_ref) {

  # Combine the reference data set along with the data subsetted by the group_by
  # function by Dataset.
  rt_error <- dplyr::inner_join(
    x = x_ref %>%
      dplyr::group_by(Proteoform) %>%
      dplyr::slice_max(`Feature intensity`) %>%
      # Remove rows where the max `Feature intensity` appears in more than one
      # row by selecting the min `Feature apex` value.
      dplyr::slice_min(`Feature apex`) %>%
      dplyr::select(RTalign, Proteoform) %>%
      # We only want to keep one row with the repMass value for each Proteoform.
      # When joining with the non-reference data set the number of rows is
      # inflated because the join function creates a row for each combination of
      # values that differ between the two data sets.
      dplyr::distinct(),
    y = tibble::tibble(RTalign = rt,
                       Proteoform = pf,
                       `Feature intensity` = fi,
                       `Feature apex` = fa) %>%
      dplyr::group_by(Proteoform) %>%
      dplyr::slice_max(`Feature intensity`) %>%
      # Remove rows where the max `Feature intensity` appears in more than one
      # row by selecting the min `Feature apex` value.
      dplyr::slice_min(`Feature apex`) %>%
      dplyr::select(RTalign, Proteoform) %>%
      # We only want to keep one row with the repMass value for each Proteoform.
      # When joining with the non-reference data set the number of rows is
      # inflated because the join function creates a row for each combination of
      # values that differ between the two data sets.
      dplyr::distinct(),
    by = "Proteoform"
  ) %>%
    # Calculate the median absolute deviation of the rt errors.
    dplyr::mutate(error = RTalign.y - RTalign.x) %>%
    dplyr::pull(error)

  return (stats::mad(rt_error))

}
