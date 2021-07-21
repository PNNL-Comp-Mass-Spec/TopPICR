# count_by: A character string indicating if observations should be counted by
# "Gene" or "ProteoForm".
# (Maybe ProteoForm can be replaced by a different variable such as cleanSeq?)
# threshold: An integer indicating the number of subjects (or project IDS) a
# Gene or ProteoRorm and UniProt accession combination must occur in to be kept.
# the default value is 10.
x_15 <- function (x, count_by, threshold = 10) {
  
  # Filter by the number of samples each ProteoForm occurs in ------------------
  
  # Remove from x any rows that contain values in the ProteoForm column that do
  # not occur in at least 10 subjects (ProjID).
  x <- x %>%
    distinct(ProjID, !!rlang::sym(count_by), UniProtAcc) %>%
    group_by(!!rlang::sym(count_by), UniProtAcc) %>%
    tally() %>%
    filter(n >= threshold) %>%
    select(-n) %>%
    semi_join(x, .) %>%
    ungroup()
  
  return (x)
  
}
