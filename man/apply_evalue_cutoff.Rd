% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/control_fdr.R
\name{apply_evalue_cutoff}
\alias{apply_evalue_cutoff}
\title{Filter the data by E-value}
\usage{
apply_evalue_cutoff(x, e_vals)
}
\arguments{
\item{x}{A \code{data.table} output from the \code{augment_annotation}
function.}

\item{e_vals}{A vector of three numeric values to be used as E-value cutoffs.
The E-value cutoffs must be in the order of unique(x$AnnType). For example,
if \code{unique(x$AnnType) = "SwissProt", "TrEMBL", "VarSplic"} then the
E-value cutoff for SwissProt must be the first element of e_vals, the
E-value cutoff for TrEMBL must be the second element of e_vals, and the
E-value cutoff for VarSplic must be the third element of e_vals.}
}
\value{
A \code{data.table} with the rows removed that fall below the E-value
cutoff.
}
\description{
The \code{control_fdr()} function filters the data by E-value for each
annotation type. If the E-values from the \code{eval_cutoff()} function are
used the overall FDR will fall below the given threshold.
}
\author{
Evan A Martin
}
