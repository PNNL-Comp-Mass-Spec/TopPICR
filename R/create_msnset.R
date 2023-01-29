# create_msnset

#' Convert TopPICR output to an MSnSet object
#'
#' Takes the combined data (output from the match between runs step),
#' dataset description and the proteoform metadata
#' and creates an MSnSet object.
#'
#' @param x A \code{data.table} output from the \code{match_features}
#'   function.
#'
#' @param meta A \code{data.table} output from the \code{create_mdata} function.
#'
#' @param dataset_feature_data A \code{data.frame} optional argument
#' that links `Dataset` with fractions if sample were fractionated. This is
#' important because observation of proteoform in a fraction defined as separate
#' features.
#'
#' @param dataset_sample_data A \code{data.frame} optional argument
#' that links `Dataset` with `sample_name`. Again, if the sample was
#' fractionated, multiple datsets correspond to a single sample. Additional
#' information in this table can be any other technical (e.g. batch, replicate)
#' or phenotypic information about samples.
#'
#' @param proteoform_def A \code{character vector} describing how proteoform ID
#' is constructed. By default it is "Gene" and "pcGroup" (proteform cluster group).
#'
#' @param partitioning A \code{character vector} describing how proteoform ID
#' is partitioned into features. For example, if the sample is fractionated then
#' `fraction` (defined in `dataset_feature_data` table) is part of the feature
#' definition. Also, if the dataset was acquired on an instrument with FAIMS
#' capability, then `CV` value (present in the first argument) is also part of
#' the feature definition. That is if the data was aquired with fractionation
#' and CV stepping then a proteoform is split into multiple features
#' (spread across fractions and CV values). Default is NULL, implying that there
#' was no partitioning of the proteoforms.
#'
#'
#' @return An MSnSet object
#'
#' @export
#'
#' @author Vlad Petyuk
#'


create_msnset <- function(x,
                          meta,
                          dataset_feature_data,
                          dataset_sample_data,
                          proteoform_def = c("Gene", "pcGroup"),
                          partitioning = NULL){

   # add info encoded in dataset names
   if(!missing(dataset_feature_data))
      x <- inner_join(x, dataset_feature_data, by = "Dataset")
   if(!missing(dataset_sample_data)){
      x <- inner_join(x, dataset_sample_data, by = "Dataset")
   }else{
      x <- mutate(x, sample_name = Dataset)
      dataset_sample_data <- distinct(x, Dataset, sample_name)
   }

   # # define features
   feature_def <- c(proteoform_def, partitioning)
   x <- x %>%
      tidyr::unite("feature_id", all_of(feature_def), sep = "_", remove = FALSE) %>%
      tidyr::unite("proteoform_id", all_of(proteoform_def), sep = "_", remove = FALSE)

   # EXPRESSION
   y <- select(x, feature_id, sample_name)
   if(any(duplicated(y))){
      stop("Features and samples do not result in unique combinations.
   Most likely there are fractions or CV values that aren't accounted in the features definition.",
           call. = FALSE)
   }
   x_expr <- x %>%
      pivot_wider(id_cols = "feature_id",
                  names_from = "sample_name",
                  values_from = "Intensity") %>%
      as.data.frame() %>%
      # {rownames(.) <- .$feature_id;.} %>%
      `rownames<-`(.$feature_id) %>%
      select(-feature_id) %>%
      as.matrix()

   # FEATURES
   # summarize intensities
   x_feat <- x %>%
      group_by(feature_id, proteoform_id, across(all_of(feature_def))) %>%
      summarize(median_intensity = median(Intensity),
                count = n(),
                .groups = "keep") %>%
      ungroup()

   # add proteoform meta data
   if(!missing(meta)){
      x_feat <- left_join(x_feat, meta, by = proteoform_def)
   }
   x_feat <- x_feat %>%
      as.data.frame() %>%
      # {rownames(.) <- .$feature_id;.}
      `rownames<-`(.$feature_id)


   # PHENO
   x_pheno <- dataset_sample_data %>%
      select(-Dataset) %>%
      distinct() %>%
      as.data.frame() %>%
      # {rownames(.) <- .$sample_name;.}
      `rownames<-`(.$sample_name)

   m <- MSnSet(x_expr,
               x_feat[rownames(x_expr),,drop=FALSE],
               x_pheno[colnames(x_expr),,drop=FALSE])
   return(m)
}



