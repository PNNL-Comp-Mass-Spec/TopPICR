#' Merge clusters using the overlap coefficient
#'
#' ...
#'
#' @param x ...
#'
#' @param oc_cutoff ...
#'
#' @return ...
#'
#' @importFrom magrittr %>%
#'
#' @export
#'
merge_clusters_oc <- function (x, oc_cutoff) {

  output <- x[[1]] %>%
    dplyr::group_by(Gene) %>%
    dplyr::do(out = merge_clusters_wrapper(., oc_cutoff = oc_cutoff)) %>%
    tidyr::unnest(cols=c(out))

  return (output)

}

oc <- function (xi, xj) {

  ci <- xi[, "Proteoform", drop = TRUE]
  cj <- xj[, "Proteoform", drop = TRUE]

  min(sum(ci %in% cj), sum(cj %in% ci))/min(length(ci), length(cj))

}

get_mergeable_clusters <- function(x, oc_cutoff){
  clusters <- unique(x$cluster)
  for(i in clusters){
    xi <- x[x$cluster == i,]
    for(j in clusters){
      if(i != j){
        xj <- x[x$cluster == j,]
        over_coef_val <- oc(xi, xj)
        # condition of a mergeable cluster
        if(over_coef_val >= oc_cutoff)
          return(c(i,j))
      }
    }
  }
  return(NULL)
}

# Recursion
# not an efficient algorithm, especially given the updates of cluster membership
# after every iteration. But it works.
merge_clusters2 <- function(clusters_df, oc_cutoff){
  ij <- get_mergeable_clusters(clusters_df, oc_cutoff)
  if(is.null(ij)){
    return(clusters_df)
  }else{
    # assumes that "i" (ij[1]) is larger cluster than "j" (ij[2]). This can be
    # achieved by pre-ordering.
    clusters_df[clusters_df$cluster == ij[2], "cluster"] <- ij[1]
    merge_clusters2(clusters_df, oc_cutoff)
  }
}

merge_clusters_wrapper <- function (x, oc_cutoff) {
  idx0 <- x$cluster == 0
  x$cluster_new <- 0
  y <- merge_clusters2(x[!idx0, c("Proteoform","cluster")], oc_cutoff)
  x$cluster_new[!idx0] <- y$cluster
  x$Gene <- NULL # removes duplicate gene due to the group/do functions.
  return (x)
}
