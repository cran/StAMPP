\name{potato.mini}
\alias{potato.mini}
\docType{data}
\title{
Smaller example genotype input format
}
\description{
A data frame containing Solcap potato genotype data in tetraploid and diploid format as an small example of the input format required by StAMPP
}
\usage{data(potato.mini)}
\format{
  A data frame with 6 observations on the following 48 variables.
  \describe{
    \item{\code{Sample}}{a factor with levels \code{ind1} \code{ind106} \code{ind140} \code{ind141} \code{ind142} \code{ind143} \code{ind144} \code{ind15} \code{ind16} \code{ind161} \code{ind162} \code{ind17} \code{ind18} \code{ind2} \code{ind257} \code{ind258} \code{ind259} \code{ind290} \code{ind291} \code{ind292} \code{ind293} \code{ind294} \code{ind295} \code{ind296} \code{ind3} \code{ind4} \code{ind5} \code{ind6} \code{ind94} \code{ind95}}
    \item{\code{Pop}}{a factor with levels \code{2x_A_Population} \code{4x_B_Population} \code{4x_C_Population}}
    \item{\code{Ploidy}}{a numeric vector}
    \item{\code{Format}}{a factor with levels \code{BiA}}
    \item{\code{solcap_snp_c1_1}}{a factor with levels \code{AA} \code{AAAA} \code{AAAB} \code{AABB} \code{AB}}
    \item{\code{solcap_snp_c1_1000}}{a factor with levels \code{AAAB} \code{AABB} \code{AB} \code{ABBB} \code{BB} \code{BBBB}}
    \item{\code{solcap_snp_c1_10000}}{a factor with levels \code{AA} \code{AAAB} \code{AABB} \code{AB} \code{ABBB} \code{BBBB}}
    \item{\code{solcap_snp_c1_10001}}{a factor with levels \code{AAAA} \code{AAAB} \code{AABB} \code{AB} \code{ABBB} \code{BB} \code{BBBB}}
    \item{\code{solcap_snp_c1_10011}}{a factor with levels \code{AA} \code{AB} \code{ABBB} \code{BBBB}}
    \item{\code{solcap_snp_c1_10012}}{a factor with levels \code{AAAA} \code{AAAB} \code{AB} \code{BB}}
    \item{\code{solcap_snp_c1_10031}}{a factor with levels \code{AA} \code{AAAA} \code{AAAB} \code{AABB} \code{AB} \code{ABBB} \code{BBBB}}
    \item{\code{solcap_snp_c1_10042}}{a factor with levels \code{AB} \code{ABBB} \code{BB} \code{BBBB}}
    \item{\code{solcap_snp_c1_10050}}{a factor with levels \code{AB} \code{BB} \code{BBBB}}
    \item{\code{solcap_snp_c1_10054}}{a factor with levels \code{AA} \code{AAAA} \code{AAAB} \code{AABB} \code{ABBB}}
    \item{\code{solcap_snp_c1_10109}}{a factor with levels \code{AA} \code{AAAA} \code{AAAB} \code{AABB} \code{AB} \code{ABBB} \code{BBBB}}
    \item{\code{solcap_snp_c1_10130}}{a factor with levels \code{-9} \code{AA} \code{AAAA} \code{AAAB} \code{AABB} \code{AB} \code{ABBB}}
    \item{\code{solcap_snp_c1_10157}}{a factor with levels \code{AABB} \code{AB} \code{ABBB} \code{BB} \code{BBBB}}
    \item{\code{solcap_snp_c1_10202}}{a factor with levels \code{AABB} \code{AB} \code{ABBB} \code{BB} \code{BBBB}}
    \item{\code{solcap_snp_c1_10252}}{a factor with levels \code{AA} \code{AAAA} \code{AAAB} \code{AABB} \code{ABBB}}
    \item{\code{solcap_snp_c1_10253}}{a factor with levels \code{AA} \code{AAAA} \code{AAAB} \code{AABB} \code{ABBB}}
    \item{\code{solcap_snp_c1_10255}}{a factor with levels \code{AAAB} \code{AABB} \code{AB} \code{ABBB} \code{BB} \code{BBBB}}
    \item{\code{solcap_snp_c1_1029}}{a factor with levels \code{AA} \code{AAAB} \code{AABB} \code{ABBB}}
    \item{\code{solcap_snp_c1_10295}}{a factor with levels \code{AA} \code{AAAA} \code{AAAB} \code{AABB}}
    \item{\code{solcap_snp_c1_10297}}{a factor with levels \code{-9} \code{AA} \code{AAAB} \code{AABB} \code{AB} \code{ABBB} \code{BBBB}}
    \item{\code{solcap_snp_c1_10351}}{a factor with levels \code{AAAB} \code{AABB} \code{ABBB} \code{BB} \code{BBBB}}
    \item{\code{solcap_snp_c1_10384}}{a factor with levels \code{AA} \code{AAAA} \code{AAAB} \code{AABB} \code{AB} \code{ABBB}}
    \item{\code{solcap_snp_c1_10397}}{a factor with levels \code{AA} \code{AAAA} \code{AAAB} \code{AABB} \code{ABBB}}
    \item{\code{solcap_snp_c1_10457}}{a factor with levels \code{AAAB} \code{AABB} \code{ABBB} \code{BB} \code{BBBB}}
    \item{\code{solcap_snp_c1_10491}}{a factor with levels \code{AABB} \code{AB} \code{ABBB} \code{BB} \code{BBBB}}
    \item{\code{solcap_snp_c1_10492}}{a factor with levels \code{AABB} \code{AB} \code{ABBB} \code{BB} \code{BBBB}}
    \item{\code{solcap_snp_c1_10494}}{a factor with levels \code{AAAA} \code{AAAB} \code{AB} \code{BB} \code{BBBB}}
    \item{\code{solcap_snp_c1_10579}}{a factor with levels \code{AA} \code{AAAB} \code{AABB} \code{AB} \code{ABBB} \code{BBBB}}
    \item{\code{solcap_snp_c1_10646}}{a factor with levels \code{AA} \code{AAAA} \code{AAAB} \code{AABB}}
    \item{\code{solcap_snp_c1_10669}}{a factor with levels \code{AA} \code{AAAA} \code{AAAB} \code{AABB} \code{AB}}
    \item{\code{solcap_snp_c1_10715}}{a factor with levels \code{AA} \code{AAAA} \code{AAAB} \code{AABB} \code{ABBB}}
    \item{\code{solcap_snp_c1_10737}}{a factor with levels \code{AABB} \code{AB} \code{ABBB} \code{BB} \code{BBBB}}
    \item{\code{solcap_snp_c1_10743}}{a factor with levels \code{-9} \code{AABB} \code{AB} \code{ABBB} \code{BB} \code{BBBB}}
    \item{\code{solcap_snp_c1_10762}}{a factor with levels \code{AA} \code{AAAA} \code{AAAB} \code{AABB} \code{AB} \code{ABBB}}
    \item{\code{solcap_snp_c1_10855}}{a factor with levels \code{AA} \code{AAAB} \code{AABB} \code{AB} \code{ABBB} \code{BBBB}}
    \item{\code{solcap_snp_c1_10873}}{a factor with levels \code{AA} \code{AAAA} \code{AAAB} \code{AABB} \code{AB} \code{ABBB} \code{BBBB}}
    \item{\code{solcap_snp_c1_10879}}{a factor with levels \code{-9} \code{AAAB} \code{AABB} \code{AB} \code{ABBB} \code{BB}}
    \item{\code{solcap_snp_c1_10900}}{a factor with levels \code{AA} \code{AAAA} \code{AAAB} \code{AABB} \code{AB} \code{ABBB}}
    \item{\code{solcap_snp_c1_10932}}{a factor with levels \code{AA} \code{AAAA} \code{AAAB} \code{AABB} \code{ABBB}}
    \item{\code{solcap_snp_c1_1094}}{a factor with levels \code{AAAA} \code{AAAB} \code{AABB} \code{AB} \code{ABBB} \code{BB}}
    \item{\code{solcap_snp_c1_11137}}{a factor with levels \code{AA} \code{AAAA} \code{AAAB} \code{AABB} \code{ABBB}}
    \item{\code{solcap_snp_c1_11144}}{a factor with levels \code{AA} \code{AAAA} \code{AAAB} \code{AABB} \code{ABBB}}
    \item{\code{solcap_snp_c1_11196}}{a factor with levels \code{AA} \code{AAAA} \code{AAAB} \code{AABB} \code{AB}}
    \item{\code{solcap_snp_c1_11206}}{a factor with levels \code{AAAB} \code{AABB} \code{ABBB} \code{BB} \code{BBBB}}
  }
}
\source{
The example genotype data is a subset of data from the publically avaliable Solcap potato dataset (http://solcap.msu.edu/potato_infinium.shtml) which was re-scored in GenomeStudio in diploid and tetraploid formats
}
\examples{
data(potato.mini, package="StAMPP")
}
\keyword{datasets}
