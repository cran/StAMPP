\name{stamppConvert}
\alias{stamppConvert}
\title{
Import and Convert
}
\description{
Imports biallelic AB formated or allele A frequency genotype data. If the data is in imported in biallelic AB format this function also converts it to allele frequencies.
}
\usage{
stamppConvert(genotype.file, nclusters = 1, type = "csv")
}
\arguments{
  \item{genotype.file}{
the genotype input file. This should be a R matrix object or a file path for a csv file containing the genotype data in either bialleleic AB format or allele 'A' frequency format, or a genlight object containing genotype data. 
}
  \item{nclusters}{
the number of proccesor treads or cores to use during calculations. 
}
  \item{type}{
the type of file the genotype data is being imported from; "csv" = comma seperated file, "r" = matrix in the R workspace, "genlight" = genlight object.
}
}
\details{
If importing from a csv file or the R workspace;
The first column of the genotype input file or matrix should contain genotype IDs, the second column should contain population IDs, 
the third column should contain a number indicating the ploidy level of each individual (ie 2=diploid, 4=tetraploid etc.), the fourth column should contain a code indicating the format of the genotype data ('BiA' for biallelic AB format, 'freq' for allele A frequency format) and each consecutive column should contain the genotype at a give locus coded in the coresponding format (for example Bialleleic format; AA or AB or BB or AAAB or AAAA etc. or allele A freq format; 1, 0.5, 0, 0.75 etc). Each column should have a header (eg. Inds, Pop, Format, Locus1, Locus2 etc.). 
}
\value{
An object of class data.frame which contains allele frequency data for use in other StAMPP functions
}
\author{
LW Pembleton
}
\examples{
# Import example data into the R workspace
data(potato.mini, package="StAMPP")

# Convert to allele frequencies
potato.freq <- stamppConvert(potato.mini, 1, "r")
}
