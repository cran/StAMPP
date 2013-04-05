\name{stamppFst}
\alias{stamppFst}
\title{
Fst Computation
}
\description{
This function calculates pairwise Fst values along with confidence intervals and p-values between populations according to the method proposed by Wright(1949) and updated by Weir and Cockerham (1984).
}
\usage{
stamppFst(geno, nboots = 100, percent = 95, nclusters = 1)
}
\arguments{
  \item{geno}{
a data frame containing allele frequency data generated from stamppConvert, or a genlight object containing genotype data, individual IDs, population IDs and ploidy levels.
}
  \item{nboots}{
number of bootstraps to perform across loci to generate confidence intervals and p-values
}
  \item{percent}{
the percentile to calculate the confidence interval around
}
  \item{nclusters}{
the number of proccesor treads or cores to use during calculations.
}
}
\details{
If possible, using multiple processing threads or cores is recommended to assist in calculating Fst values over a large number of bootstraps. 
}
\value{
An object list with the components:
  \item{Fsts }{a matrix of pairwise Fst values between populations}
  \item{Pvalues }{a matrix of p-values for each of the pairwise Fst values containined in the 'Fsts' matrix}
  \item{Bootstraps }{a dataframe of each Fst value generated during Bootstrapping and the associated confidence intervals}

If nboots<2, no bootstrapping is performed and therefore only a matrix of Fst values is returned.  
}
\references{
Wright S (1949) The Genetical Structure of Populations. Annals of Human Genetics 15, 323-354.

Weir BS, Cockerham CC (1984) Estimating F Statistics for the ANalysis of Population Structure. Evolution 38, 1358-1370.
}
\author{
LW Pembleton
}
\examples{
# import genotype data and convert to allele frequecies 
data(potato.mini, package="StAMPP")
potato.freq <- stamppConvert(potato.mini, 1, "r")

# Calculate pairwise Fst values between each population
potato.fst <- stamppFst(potato.freq, 100, 95, 1)
}
