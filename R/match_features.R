#' Extract feature intensity values from unidentified feature data
#'
#' ...
#'
#' @param x A \code{data.table} output from the \code{merge_clusters} function.
#'
#' @param feature A \code{data.table} containing the unidentified feature data.
#'   This data is created from the files output by TopPIC ending with
#'   "ms1.feature" and must contain the \code{RTalign} and \code{RecalMass}
#'   variables. These variables are created with the unidentified feature data
#'   in the \code{align_rt} and \code{recalibrate_mass} functions.
#'
#' @param ppm_cutoff An integer representing the retention time cutoff in
#'   seconds. This value is the threshold used to determine if an unidentified
#'   feature is close enough in retention time to be considered part of an
#'   identified feature cluster. The mean retention time of all points in the
#'   cluster is used for the comparison.
#'
#' @param rt_cutoff An integer indicating the mass cutoff in ppm. This value is
#'   used as the threshold for determining if an unidentified feature is close
#'   enough in mass to be considered part of an identified feature cluster. The
#'   mean \code{RecalMass} of all points in the cluster is used for the
#'   comparison.
#'
#' @return A \code{data.table} containing all unidentified features that fall
#'   within the threshold of an identified feature gene/cluster combination.
#'
#' @author Evan A Martin
#'
#' @export
#'
impute_fi <- function(x, feature, ppm_cutoff, rt_cutoff) {

  centroids <- x %>%
    # Remove the noise cluster. We don't want to impute feature intensity values
    # for it.
    dplyr::filter(cluster != 0) %>%
    # Group by Gene and cluster--the combination of these two variables is a
    # proteoform. We will calculate the centroid of the mass and retention time
    # for each proteoform. The centroids will be used to extract feature
    # intensity values from the unidentified feature data.
    dplyr::group_by(Gene, cluster) %>%
    # Calculate the mean RecalMass and mean RTalign which will be used as the
    # centroid in finding like unidentified features.
    dplyr::summarize(meanMass = mean(RecalMass, na.rm = TRUE),
                     meanRT = mean(RTalign, na.rm = TRUE))

  # Prepare to run in parallel.
  cores <- parallel::detectCores() - 1
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)

  # Loop through each centroid and find corresponding unidentified features.
  unidentified <- foreach::foreach(e = 1:dim(centroids)[[1]],
                                   .export = "%>%") %dopar% {

    features <- feature %>%
      dplyr::filter(
        (1e6 * abs(RecalMass - centroids$meanMass[[e]]) /
           centroids$meanMass[[e]]) < ppm_cutoff
      ) %>%
      dplyr::filter(
        abs(RTalign - centroids$meanRT[[e]]) < rt_cutoff
      ) %>%
      dplyr::select(Dataset, Intensity, CV, Time_apex,
                    RTalign, RecalMass) %>%
      dplyr::mutate(
        Gene = centroids$Gene[[e]],
        cluster = centroids$cluster[[e]]
      )

    # Using the foreach function with %dopar% will assign the last element
    # within the curly brackets to the object when foreach is called. In this
    # case features will be assigned to the eth element of unidentified.
    features

  }

  parallel::stopCluster(cl)

  # Combine the list containing the Dataset, Intensity, CV, Gene, and cluster
  # variables for the unidentified features that are within the threshold of the
  # centroid.
  return (data.table::rbindlist(unidentified))

}
