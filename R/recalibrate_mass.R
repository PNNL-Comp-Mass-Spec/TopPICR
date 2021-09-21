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

  # Compute the robust median of the precursor mass within Gene, Dataset, and
  # Proteoform. This variable will be used to calculate the error (in ppm)
  # between the reference and all other data sets.
  x_abridged <- x %>%
    dplyr::select(c(Dataset, Proteoform, Gene,
                    `Precursor mass`, `Adjusted precursor mass`)) %>%
    dplyr::group_by(Dataset, Proteoform, Gene) %>%
    dplyr::summarize(Mass = robust_median(mass = `Precursor mass`,
                                          adj_mass = `Adjusted precursor mass`))

  # Combine the data frames with the robust median mass and the aligned
  # retention times.
  x_mass <- x %>%
    dplyr::ungroup() %>%
    dplyr::inner_join(x_abridged)

  # Extract the reference data set.
  x_ref <- x_mass %>%
    dplyr::filter(Dataset == ref_ds) %>%
    dplyr::select(-Dataset, -Gene)

  # Compute the mass and retention time errors.
  x_error <- x_mass %>%
    dplyr::group_by(Dataset) %>%
    dplyr::summarize(ppm_rt_error(mass = Mass,
                                  rt = RTalign,
                                  pf = Proteoform,
                                  x_ref = x_ref)) %>%
    dplyr::ungroup()

  # Compute the median of the mass and retention time errors across all data
  # sets (excluding the reference data set). These values will be used in the
  # cluster functions to determine a mass/retention time envelope around each
  # cluster.
  errors <- x_error %>%
    dplyr::filter(Dataset != ref_ds) %>%
    dplyr::summarize(ppm_error = stats::median(ppm_mad),
                     rt_error = stats::median(rt_mad))

  # Return the x_error data frame that contains the RecalMass variable along
  # with the ppm and retention time errors.
  return (list(x = x_error,
               ppm_error = errors$ppm_error[[1]],
               rt_error = errors$rt_error[[1]]))

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
#' @return A \code{data.table} containing the recalibrated mass variable. The
#'   following variables were added/removed:
#'
#'   | Added             | Removed                    |
#'   | ----------------- | -------------------------- |
#'   | `RecalMass`       |                            |
#'
#' @importFrom magrittr %>%
#'
#' @export
#'
recalibrate_mass <- function (x, errors) {

  # Check to make sure the data sets in x are also in errors$x. If the data sets
  # in each data frame are not identical the mass cannot be recalibrated. Here x
  # is the data.table output by the align_rt function and errors$x is the
  # data.table output by the calc_error function.
  if (!setequal(unique(x$Dataset), unique(errors$x$Dataset))) {

    stop (paste("The data sets present in x and errors$x do not match."))

  }

  #!#!#!#!#!#!#! START #!#!#!#!#!#!#!

  # Test code to make sure the changes produce the same output as the original
  # function.

  x_abridged <- x %>%
    dplyr::group_by(Dataset, Proteoform, Gene) %>%
    dplyr::summarize(
      Mass = robust_median(mass = `Precursor mass`,
                           adj_mass = `Adjusted precursor mass`)
    )

  x2 <- x %>%
    dplyr::ungroup() %>%
    dplyr::inner_join(x_abridged)

  # NOTE: When this code is removed the variables in the calls to the inner_join
  # and mutate functions will need to be changed.

  #!#!#!#!#!#!#! END #!#!#!#!#!#!#!

  # Recalibrate the mass.
  x_recal <- dplyr::inner_join(x = x2,
                               y = errors$x,
                               by = "Dataset") %>%
    dplyr::mutate(RecalMass = Mass * (1 - ppm_median / 1e6)) %>%
    # Remove intermediate variables used to recalibrate the mass. These
    # variables are in the output of the calc_error function and can be accessed
    # there if they are needed.
    dplyr::select(-ppm_mad, -ppm_median, -rt_mad)

  # Return the x_recal data frame that contains the RecalMass variable.
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
  # the Dataset, Proteoform, Gene combination.
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

# mass: A vector of masses extracted from the data matrix by the group_by
# function.
# rt: A vector of retention times extracted from the data matrix by the group_by
# function.
# pf: A vector of proteoforms extracted from the data matrix by the group_by
# function.
# x_ref: the x data matrix filtered by the reference data set.
ppm_rt_error <- function (mass, rt, pf, x_ref) {

  # Combine the reference data set along with the data subsetted by the group_by
  # function by Gene and Proteoform.
  error <- dplyr::inner_join(x = x_ref[, c("Mass", "RTalign", "Proteoform")],
                             y = tibble::tibble(Mass = mass,
                                                RTalign = rt,
                                                Proteoform = pf),
                             by = "Proteoform") %>%
    # Calculate the ppm and retention time errors.
    dplyr::mutate(ppm_error = 1e6 * (Mass.y - Mass.x) / Mass.x,
                  rt_error = RTalign.y - RTalign.x)

  # Only keep ppm values between plus/minus 20 for computing the median. The
  # ppm median is used for recalibrating the mass later. We only keep the values
  # within 20 ppm for a more accurate representation of the group.
  reduced <- error %>%
    dplyr::filter(ppm_error > -20 & ppm_error < 20)

  # Calculate the median absolute deviation, median and standard deviation of
  # the remaining ppm values.
  mms <- data.frame(ppm_mad = stats::mad(error$ppm_error),
                    ppm_median = stats::median(reduced$ppm_error),
                    rt_mad = stats::mad(error$rt_error))

  return (mms)

}
