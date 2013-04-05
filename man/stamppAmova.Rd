\name{stamppAmova}
\alias{stamppAmova}
\title{
Analysis of Molecular Variance
}
\description{
Calculates an AMOVA based on the genetic distance matrix from stamppNeisD() using the amova() function from the package PEGAS for exploring within and between population variation.
}
\usage{
stamppAmova(dist.mat, geno, perm = 100)
}
\arguments{
  \item{dist.mat}{
the matrix of genetic distances between individuals generated from stamppNeisD()
}
  \item{geno}{
a data frame containing allele frequency data generated from stamppConvert, or a genlight object containing genotype data, individual IDs, population IDs and ploidy levels.
}
  \item{perm}{
the number of permutations for the tests of hypotheses. 
}
}
\details{
Uses the formula distance ~ populations, to calculate an AMOVA for population differentiation and within & between population variation.

This function uses the amova function from the PEGAS package.
}
\value{
An object of class "amova" which is a list containing a table of sum of square deviations (SSD), mean square deviations (MSD) and the number of degrees of freedom as well as the variance components
}
\references{
Paradis E (2010) pegas: an R package for population genetics with an integrated-modular approach. Bioinformatics 26, 419-420.
}
\author{
LW Pembleton
}
\examples{
# import genotype data and convert to allele frequecies 
data(potato.mini, package="StAMPP")
potato.freq <- stamppConvert(potato.mini, 1, "r")

# Calculate genetic distance between individuals
potato.D.ind <- stamppNeisD(potato.freq, FALSE)

# Calculate AMOVA
stamppAmova(potato.D.ind, potato.freq, 100)
}
