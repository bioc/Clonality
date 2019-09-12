\name{SNVtest2}
\alias{SNVtest2}

\title{Test for tumors from 2 different sites}
\description{This functions performs clonality testing of 2 tumors from the same patient that come from different types or organs. The test uses profiles of somatic mutations. The null hypothesis is that two tumors that come from different sites are independent. There are two alternative hypotheses: that they are clonal and from site 1, and clonal from site 2. This function calculates conditional maximum likelihood ratio relying only on loci where at least one of the tumors has a mutation, and p-value is calculated under the reference distribution under the hypothesis of independence.

}
\usage{
SNVtest2(tumor1, tumor2, pfreq, nrep = 999)
}
\arguments{
  \item{tumor1}{
Vector of the binary mutation calls from tumor 1, where 0 denotes no mutation, 1 denotes a mutation. Mutations should be in the same order as frequencies in  pfreq.
}
  \item{tumor2}{
Vector of the binary mutation calls from tumor 2, where 0 denotes no mutation, 1 denotes a mutation. Mutations should be in the same order as frequencies in  pfreq.
}
  \item{pfreq}{
Two column matrix of marginal frequencies of mutations in two sites known a priori. These can be obtained from TCGA or similar databases or calculated using function "get.mutation.frequencies" in this package. 
}
  \item{nrep}{
Number of simulations used for generating the reference distribution under the hypothesis of independence.
}
}
\details{
Test is related to the one described by "SNVtest" but it assumes that under null hypothesis two tumors come from 2 different sites with different marginal proabilities, thus the reference distribution is generated from 2 sets of frequencies. The p-value is significant when null hypothesis is rejected. It might be rejected sometimes in the absence of matches if the observed mutational profiles are unlikely to come from the declared tumor sites.

}
\value{
The output is a vector with 6 values:
"n.match", "n.site1only", "n.site2only", "xi.site1", "xi.site2", "p.value"

 \item{n.match }{Number of matches between two tumors}
 \item{n.site1only }{Number of mutations in the first tumor only}
\item{n.site2only }{Number of mutations in the second tumor only}
 \item{xi.site1 }{Maximum likelhood estimate of Ksi, parameter of the likelihood representing clonality strength, for the alternative hypothesis that two clonal tumors come from site 1. Value close to 0 indicates independence, value close to 1 indicates perfect concordance in mutational profiles.} 
  \item{xi.site2 }{Maximum likelhood estimate of Ksi, parameter of the likelihood representing clonality strength, for the alternative hypothesis that two clonal tumors come from site 2. Value close to 0 indicates independence, value close to 1 indicates perfect concordance in mutational profiles.} 
 \item{p.value }{p-value calculated using the null distribution generated using two independet tumors generated from two different sites}

}
\references{

Ostrovnaya I, Mauguen A, Seshan V, Begg CB "Testing Tumors from Different Anatomic Sites for Clonal Relatedness Using Somatic Mutation Data", submitted

}
\author{
Irina Ostrovnaya \email{ostrovni@mskcc.org}
}

\examples{


#Example of artificially generated independent tumors  from marginal mutation frequencies p1 and p2
n<-100
p1<-runif(n)/10
p2<-runif(n)/10
x1<-as.numeric(runif(n)<=p1)
x2<-as.numeric(runif(n)<=p2)
SNVtest2(x1,x2,cbind(p1,p2))


}