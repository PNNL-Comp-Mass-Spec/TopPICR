#' Extract feature intensity values from unidentified feature data
#'
#' Matches unidentified feature intensities with identified feature clusters.
#' The mass and retention time of each unidentified feature is compared to the
#' centroid of each identified cluster. If the unidentified feature falls within
#' the specified ppm threshold it is added to the cluster.
#'
#' @param ms2 A \code{data.table} containing the identified feature data.
#'
#' @param ms1 A \code{data.table} containing the unidentified feature data.
#'   This data is created from the files output by TopPIC ending with
#'   "ms1.feature" and must contain the \code{RTalign} and \code{RecalMass}
#'   variables. These variables are created with the unidentified feature data
#'   in the \code{align_rt} and \code{recalibrate_mass} functions.
#'
#' @param errors A \code{list} output from the \code{calc_error} function.
#'
#' @param n_mme_sd A numeric value indicating the number of standard
#'   deviations to use when creating a cutoff in the mass dimension. This
#'   threshold is used for determining whether an unidentified feature is close
#'   enough in mass to be considered part of an identified feature cluster. The
#'   mean \code{RecalMass} of all points in the cluster is used for the
#'   comparison.
#'
#' @param n_rt_sd A numeric value representing the number of standard deviations
#'   to use when creating a retention time cutoff. This value is the threshold
#'   used to determine if an unidentified feature is close enough in retention
#'   time to be considered part of an identified feature cluster. The mean
#'   retention time of all points in the cluster is used for the comparison.
#'
#' @param summary_fn A character string specifying the function to use when
#'   summarizing the feature intensity. The function must contain the
#'   \code{na.rm} argument. Some examples of functions that are allowed are
#'   \code{max}, \code{sum}, or \code{median}.
#'
#' @return A \code{data.table} containing all unidentified features that fall
#'   within the threshold of an identified feature gene/cluster combination.
#'
#' @importFrom foreach `%dopar%`
#' @importFrom magrittr %>%
#'
#' @author Evan A Martin
#'
#' @export
#'
match_features <- function(ms2, ms1, errors, n_mme_sd, n_rt_sd,
                           summary_fn = "max") {

  # Convert n_mme_sd to ppm. This is necessary so the remainder of the
  # function does not have to be altered to reflect the change from using a ppm
  # cutoff to a cutoff in standard deviations.
  ppm_cutoff <- n_mme_sd * errors$ppm_sd

  # Convert n_rt_sd to seconds. This is necessary so the remainder of the
  # function does not have to be altered to reflect the change from using a
  # retention time cutoff in seconds to a cutoff in standard deviations.
  rt_cutoff <- n_rt_sd * errors$rt_sd

  # Match between runs ---------------------------------------------------------

  centroids_mass <- ms2 %>%
    # Remove the noise cluster. We don't want to impute feature intensity values
    # for it.
    dplyr::filter(cluster != 0) %>%
    # Group by Gene and cluster--the combination of these two variables is a
    # proteoform. We will calculate the centroid of the mass and retention time
    # for each proteoform. The centroids will be used to extract feature
    # intensity values from the unidentified feature data.
    dplyr::group_by(Gene, cluster, pcGroup) %>%
    # Calculate the median RecalMass which will be used as the centroid in
    # finding like unidentified features.
    dplyr::summarize(mediMass = stats::median(RecalMass, na.rm = TRUE))

  centroids_rt <- ms2 %>%
    # Remove the noise cluster. We don't want to impute feature intensity values
    # for it.
    dplyr::filter(cluster != 0) %>%
    # Group by Gene and pcGroup to calculate the centroid for the retention
    # time. The centroids will be used to extract feature intensity values from
    # the unidentified feature data.
    dplyr::group_by(Gene, pcGroup) %>%
    # Calculate the median RTalign which will be used as the centroid in finding
    # like unidentified features.
    dplyr::summarize(mediRT = stats::median(RTalign, na.rm = TRUE))

  centroids <- dplyr::inner_join(x = centroids_mass,
                                 y = centroids_rt)

  # Prepare to run in parallel.
  cores <- parallel::detectCores() - 1
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)

  # Loop through each centroid and find corresponding unidentified features.
  unidentified <- foreach::foreach(e = 1:dim(centroids)[[1]],
                                   .export = "%>%") %dopar% {

    features <- ms1 %>%
      dplyr::filter(
        (1e6 * abs(RecalMass - centroids$mediMass[[e]]) /
           centroids$mediMass[[e]]) < ppm_cutoff
      ) %>%
      dplyr::filter(
        abs(RTalign - centroids$mediRT[[e]]) < rt_cutoff
      ) %>%
      dplyr::select(Dataset, CV, Intensity) %>%
      dplyr::mutate(
        Gene = centroids$Gene[[e]],
        cluster = centroids$cluster[[e]],
        pcGroup = centroids$pcGroup[[e]]
      )

    # Using the foreach function with %dopar% will assign the last element
    # within the curly brackets to the object when foreach is called. In this
    # case features will be assigned to the eth element of unidentified.
    features

  }

  parallel::stopCluster(cl)

  # Create an object that points to the function that will be used to create a
  # representative feature intensity value.
  the_summary <- match.fun(summary_fn)

  # Combine the list containing the Dataset, Intensity, CV, Gene, and cluster
  # variables for the unidentified features that are within the threshold of the
  # centroid.
  output <- data.table::rbindlist(unidentified) %>%
    dplyr::group_by(Dataset, CV, Gene, cluster, pcGroup) %>%
    dplyr::summarize(
      Intensity = the_summary(Intensity, na.rm = TRUE)
    )

  # Combine MS1 and MS2 features -----------------------------------------------

  ms2 <- ms2 %>%
    dplyr::group_by(Dataset, CV, Gene, cluster, pcGroup) %>%
    # Just want to keep the maximum feature intensity per cluster.
    dplyr::summarize(
      `Feature intensity` = max(`Feature intensity`, na.rm = TRUE)
    )

  fused <- dplyr::full_join(output, ms2) %>%
    dplyr::mutate(
      Intensity = pmax(Intensity, `Feature intensity`, na.rm = TRUE)
    ) %>%
    dplyr::select(-`Feature intensity`)

  fused <- fused %>%
    # The following code keeps only the highest intensity value when a given
    # Dataset has multiple Intensity values within a pcGroup. This computation
    # is carried out within each Gene, CV, and Fraction.
    dplyr::group_by(Dataset, CV, Gene, pcGroup) %>%
    dplyr::summarize(Intensity = the_summary(Intensity, na.rm = TRUE)) %>%
    dplyr::ungroup()

  return (fused)

}
