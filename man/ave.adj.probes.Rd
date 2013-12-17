\name{ave.adj.probes}
\alias{ave.adj.probes}

\title{
Averaging of adjacent probes in copy number arrays
}
\description{
For each sample the log-ratios at each consecutive K number of probes are averaged.
}
\usage{
ave.adj.probes(data, K)
}

\arguments{
  \item{data}{
Copy Number Array object (output of function CNA() from the package DNAcopy). First column contains chromosomes, second column contains genomic locations. Each remaining column contains log-ratios from a particular tumor or sample.
}
  \item{K}{
Number of markers to be averaged. Should be selected so that the final resolution of the averaged data would be 5,000-10,000 markers.
}
}
\details{
Averages log-ratios in every K consecutive markers. The purpose of this step is to reduce the noise in the data, eliminate possible very small germline copy number variations, and get rid of a possible wave effect.
}
\value{
Returns CNA object of reduced resolution
}

          
\examples{
# Same example as in clonality.analysis()
#Analysis of paired breast samples from study
#Hwang ES, Nyante SF, Chen YY, Moore D, DeVries S, Korkola JE, Esserman LJ, and Waldman FM. 
#Clonality of lobular carcinoma in situ and synchronous invasive lobular cancer. Cancer 100(12):2562-72, 2004.
#library(gdata)          #needed to read .xls files
#library(DNAcopy)
#arrayinfo<-read.xls("http://waldmanlab.org/Colon/nakao.data.xls")  #needed to extract  genomic locations
#data<-read.xls("http://waldmanlab.org/Breast/Hwang.data.xls")
#data<-data[!is.na(data[,2]),]
#data<-data[apply(is.na(data),1,sum)<=50,]
#data<-data[,apply(is.na(data),2,sum)<=1000]
#data$Position<-arrayinfo$Mb[match(toupper(as.character(data[,1])),toupper(as.character(arrayinfo[,1])))]
#data<-data[!is.na(data$Position),]
#data<-data[c(TRUE,data$Position[-1]!=data$Position[-1864]),] #discard probes with repeated genomic locations
#data<-data[data$Chromosome<=22,] #getting rid of X and Y chromosomes
#dataCNA<-CNA(as.matrix(data[,c(4:51)]),maploc=data$Position,chrom=data$Chromosome,sampleid=names(data)[c(4:51)])
#The averaging is not needed. Used only as example.
#dataAve<- ave.adj.probes(dataCNA,2)

}
