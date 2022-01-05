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
#' @param errors A \code{list} output from the \code{calc_error} function. The
#'   ppm and retention time standard deviations (across all data sets) are the
#'   last two elements of this list. They will be used along with the
#'   \code{n_ppm_sd} and \code{n_rt_sd} arguments for determining whether or not
#'   an unidentified feature belongs to a given cluster.
#'
#' @param n_ppm_sd The number of standard deviations used to create a ppm cutoff
#'   when determining if an unidentified feature belongs to a given cluster. A
#'   feature belongs to a cluster if its recalibrated mass and aligned retention
#'   time fall within the mass/retention time envelope around the cluster's
#'   centroid.
#'
#' @param n_rt_sd The number of standard deviations used to create a retention
#'   time cutoff when determining if an unidentified feature belongs to a given
#'   cluster. A feature belongs to a cluster if its recalibrated mass and
#'   aligned retention time fall within the mass/retention time envelope around
#'   the cluster's centroid.
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
#'
#' @author Evan A Martin
#'
#' @export
#'
match_features <- function(ms2, ms1, errors, n_ppm_sd, n_rt_sd,
                           summary_fn = "max") {

  # Create the mass and retention time cutoffs based on the respective standard
  # deviations and number of standard deviations.
  ppm_cutoff <- errors$ppm_sd * n_ppm_sd
  rt_cutoff <- errors$rt_sd * n_rt_sd

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

  output <- data.table::rbindlist(unidentified) %>%
    dplyr::group_by(Dataset, CV, Gene, cluster, pcGroup) %>%
    dplyr::summarize(
      Intensity = the_summary(Intensity, na.rm = TRUE)
    )

  # Combine the list containing the Dataset, Intensity, CV, Gene, and cluster
  # variables for the unidentified features that are within the threshold of the
  # centroid.
  return (output)

}
