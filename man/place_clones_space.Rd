% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_instances.R
\name{place_clones_space}
\alias{place_clones_space}
\title{Create a model for the spatial distribution of the clones in a tumor}
\usage{
place_clones_space(B)
}
\arguments{
\item{B}{A matrix representing the mutation relationships between the nodes in the phylogenetic tree (B matrix).}
}
\value{
A list containing two elements: 'spatial_coords', a data frame where each column represents the density of a clone at different spatial coordinates, and 'x', a numeric vector representing the spatial coordinates.
}
\description{
This function creates a Gaussian mixture model for the spatial distribution of the clones in a tumor in a 1D space. In the model, each component represents a clone, and the mean of the component represents the position of the clone in the space. The standard deviation of the components is fixed to 1, while the mean values are random variables.
}
\examples{
# Create a model for the spatial distribution of the clones in a tumor represented by the 
# B matrix B_mat

# Create a mutation matrix for a phylogenetic tree with 10 nodes and k = 2
B_mat <- create_B(10, 2)

clone_placement <- place_clones_space(B_mat)

}
