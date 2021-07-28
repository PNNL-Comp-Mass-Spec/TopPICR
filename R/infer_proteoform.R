#' ...
#'
#' ...
#'
#' @param x ...
#'
#' @return ...
#'
#' @import data.table
#'
#' @export
#'
infer_pf <- function (x) {

  # Keep three columns of the unique rows of x.
  y <- unique(x[,c("cleanSeq","UniProtAcc","AnnType")])
  y$AnnType <- ordered(y$AnnType, levels = c("SwissProt","VarSplic","TrEMBL"))

  # Convert y to a data table.
  setDT(y)

  # Order the data table by annotation, accession, and then cleanSeq. This order
  # will be carried through the call of the inference function.
  setorder(y, AnnType, UniProtAcc, cleanSeq)

  # Selects a set of accessions which maps to the most sequences. This is a
  # different method from accessions_seq which first selects sequences that only
  # map to one accession and then finds all sequences that those accessions map
  # to.
  res <- inference(y) # Occam's razor?

  # Compare x and res row-wise. Keep any row (all columns) in x that matches any
  # row in res. The row-wise comparison will be done with all columns that are
  # found in both x and res (with res being the limiting factor).
  x <- dplyr::semi_join(x, res)

  return (x)

}

# x_30 auxiliary functions -----------------------------------------------------

inference <- function (x) {

  res <- list()

  # Continue to run until there are no rows in x. At the end of the while loop
  # the number of rows in x will be reduced.
  while (nrow(x) > 0) {

    # First []: Counts the number of times each accession appears in x and
    # creates a column that contains the count for each UniProt accession.
    # Second []: Finds the row with the accession that has the highest count
    # (occurs the most number of times in x). $UniProtAcc: Extracts the UniProt
    # accession name with the highest count.
    top_prot <- x[, .N, by = UniProtAcc][which.max(N), , ]$UniProtAcc

    # Extracts the rows in x whose UniProtAcc value matches top_prot.
    top_peps <- subset(x, UniProtAcc == top_prot)

    # Add the data frame (top_peps which is a subset of x) to the res list.
    res <- c(res, list(top_peps))

    # Remove any rows in x whose cleanSeq value matches any of the sequences in
    # the top_peps data frame. In other words, any rows in x that contain the
    # sequences mapping to the accession with the highest count are removed.
    x <- subset(x, !(cleanSeq %in% top_peps[[1]]))

  }

  return(rbindlist(res,
                   use.names = FALSE,
                   fill = FALSE,
                   idcol = NULL))

}
