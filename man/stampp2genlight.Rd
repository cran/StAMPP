\name{stampp2genlight}
\alias{stampp2genlight}
\title{
Convert StAMPP genotype data to genlight object
}
\description{
Converts a StAMPP formated allele frequency data frame generated from the stamppConvert function to a genlight object for use in other packages
}
\usage{
stampp2genlight(geno, pop = TRUE)
}
\arguments{
  \item{geno}{
a data frame containing allele frequency data generated from stamppConvert
}
  \item{pop}{
logical. True if population IDs are present in the StAMPP genotype data, False if population IDs are absent.
}
}
\details{
StAMPP only exports to genlight objects as they are able to handle mixed ploidy datasets unlike genpop and genloci objects.

The genlight object allows the intergration between StAMPP and other common R packages such as ADEGENET
}
\value{
A object of class genlight which contains genotype data, individual IDs, population IDs (if present) and ploidy levels.
}
\author{
LW Pembleton
}
\examples{
# import genotype data and convert to allele frequecies 
data(potato.mini, package="StAMPP")
potato.freq <- stamppConvert(potato.mini, 1, "r")

# Convert the StAMPP formatted allele frequency data frame to a genlight object
potato.genlight <- stampp2genlight(potato.freq, TRUE)
}
