# This script contains the functions used to prepare the unimod data to be
# combined with the modifications for a specific study from TopPIC. The unimods
# data frame is an internal object. It will be combined with output from TopPIC
# to create the modification data frames used by TopPICR plotting functions.

unimods <- preprocess_unimod("~/Downloads")

# Functions for processing unimod data -----------------------------------------

preprocess_unimod <- function (path_to_unimod) {

  x <- XML::xmlParse(
    file.path(path_to_unimod, "unimod.xml")
  )
  y <- XML::xmlToList(x)

  titles <- vector(length = length(y$modifications))
  mono_masses <- vector(length = length(y$modifications))
  classification <- vector(length = length(y$modifications))

  for (i in seq_along(y$modifications)) {

    titles[[i]] <- as.character(y$modifications[[i]]$.attrs["title"])
    mono_masses[[i]] <- as.numeric(
      y$modifications[[i]]$delta$.attrs["mono_mass"]
    )
    classification[[i]] <- .get_classification(y$modifications[[i]])

  }

  unimods <- data.frame(title = titles,
                        mono_mass = mono_masses,
                        classification = classification)

  unimods <- unimods %>%
    dplyr::filter(!(classification %in% c("Chemical derivative",
                                          "Isotopic label",
                                          "Other"))) %>%
    dplyr::filter(
      !(title %in% c("Label:13C(1)2H(3)+Oxidation","AHA-SS_CAM","AHA-SS",
                     "BDMAPP","Glu->pyro-Glu+Methyl:2H(2)13C(1)",
                     "Iodoacetanilide:13C(6)","Iodoacetanilide",
                     "Label:2H(4)+GG",
                     "DMPO",
                     "Delta:H(10)C(8)O(1)",
                     "trifluoro",
                     "HN2_mustard",
                     "HN3_mustard",
                     "AEBS",
                     "CarbamidomethylDTT",
                     "Hydroxamic_acid"))
    ) %>%
    dplyr::select(mono_mass, title) %>%
    `colnames<-`(c("mass", "name"))


  return (unimods)

}

.get_classification <- function(x){
  if(class(x$specificity) == "list"){
    mod_class <- x$specificity$.attr["classification"]
  }else{
    mod_class <- x$specificity["classification"]
  }
  return(as.character(mod_class))
}
