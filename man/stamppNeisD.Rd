\name{stamppNeisD}
\alias{stamppNeisD}
\title{
Genetic Distance Calculation
}
\description{
This function calculates Nei's genetic distance (Nei 1972) between populations or individuals.
}
\usage{
stamppNeisD(geno, pop = TRUE)
}
\arguments{
  \item{geno}{
a data frame containing allele frequency data generated from stamppConvert, or a genlight object containing genotype data, individual IDs, population IDs and ploidy levels.
}
  \item{pop}{
logical. True if genetic distance should be calculated between populations, false if it should be calculated between individual.
}
}
\value{
A object of class matrix which contains the genetic distance between each population or individual.
}
\references{
Nei M (1972) Genetic Distance between Populations. The American Naturalist 106, 283-292.
}
\author{
LW Pembleton
}
\examples{
# import genotype data and convert to allele frequecies 
data(potato.mini, package="StAMPP")
potato.freq <- stamppConvert(potato.mini, "r")

# Calculate genetic distance between individuals
potato.D.ind <- stamppNeisD(potato.freq, FALSE)

# Calculate genetic distance between populations
potato.D.pop <- stamppNeisD(potato.freq, TRUE)
}
