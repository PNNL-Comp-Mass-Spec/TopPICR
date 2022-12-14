% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster.R
\name{cluster}
\alias{cluster}
\title{Cluster data based on mass and retention time}
\usage{
cluster(x, errors, method, height, min_size)
}
\arguments{
\item{x}{A \code{data.table} output from the \code{recalibrate_mass}
function.}

\item{errors}{A \code{list} output from the \code{calc_error} function. The
first element of the list contains the standard deviation and median of the
mass measurement error for each data set. The second element is the
standard deviation of the mass measurement error across all data sets. The
third element is the standard deviation of the retention time in seconds
across all data sets.}

\item{method}{A character string indicating what agglomeration method should
be used in the hclust function. See \code{\link[stats]{hclust}} for more
details.}

\item{height}{An number specifying the height at which the tree created by
the hclust function should be cut. See \code{\link[stats]{hclust}} for more
details.}

\item{min_size}{An integer value indicating the minimum number of points a
cluster must have. All clusters with fewer members than min_size will be
reclassified as "noise" points. The cluster assignment for noise points is
0.}
}
\value{
A \code{data.table} with the cluster assignment for each observation.
}
\description{
Cluster the data using hierarchical clustering with the
\code{\link[stats]{hclust}} function. The data are clustered with the
normalized recalibrated mass and the normalized aligned retention times.
}
\author{
Evan A Martin
}
