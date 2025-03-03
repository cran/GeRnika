% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_instances.R
\name{calc_clone_proportions}
\alias{calc_clone_proportions}
\title{Calculate clone proportions for a tumor}
\usage{
calc_clone_proportions(B, selection)
}
\arguments{
\item{B}{A matrix representing the mutation relationships between the nodes in the phylogenetic tree (B matrix).}

\item{selection}{A character string representing the evolutionary mode the tumor follows. This should be either "positive" or "neutral".}
}
\value{
A data frame with two columns: 'clone_idx', which contains the clone identifiers, and 'proportion', which contains the calculated proportions of each clone.
}
\description{
This function calculates the proportions of each clone in a phylogenetic tree, following a given evolutionary mode (positive selection or neutral evolution).
}
\examples{
# Calculate clone proportions for a tumor phylogenetic tree represented by a B matrix 
# and following a positive selection model

# Create a mutation matrix for a phylogenetic tree with 10 nodes and k = 2
B_mat <- create_B(10, 2)

# Calculate the clone proportions following a positive selection model
clone_proportions <- calc_clone_proportions(B_mat, "positive")

}
