\name{splitChromosomes}
\alias{splitChromosomes}
\title{
Chromosome splitting
}
\description{
Divides the chromosomes into p and q arms.   }
\usage{
splitChromosomes(chrom,maploc)
}

\arguments{
  \item{chrom}{
Vector of chromosomes. They should be numeric 1 to 22.
}
  \item{maploc}{
Vector of genomic locations. They should be in Kilobases. 
}
}
\details{
The function returns the vector of chromosome arms labeled "chr01p", "chr01q", etc. The split into arms is accomplished using the following  centers (in Kb) for  chromosomes 1 through 22: (122356.96,  93189.90,  92037.54 , 50854.87  ,47941.40,  60438.12 , 59558.27,  45458.05 , 48607.50,  40434.94 , 52950.78,  35445.46 ,      16934.00,  16570.00,  16760.00 , 36043.30 , 22237.13,  16082.90  ,      28423.62 , 27150.40,  11760.00, 12830.00 ). 

}

\examples{
  # Same example as in clonality.analysis()
   #Analysis of paired breast samples from study
#Hwang ES, Nyante SF, Chen YY, Moore D, DeVries S, Korkola JE, Esserman LJ, and Waldman FM. 
#Clonality of lobular carcinoma in situ and synchronous invasive lobular cancer. Cancer 100(12):2562-72, 2004.
#library(gdata)    #needed to read .xls files
#library(DNAcopy) 
#arrayinfo<-read.xls("http://waldman.ucsf.edu/Colon/nakao.data.xls")  #needed to extract  genomic locations
#data<-read.xls("http://waldman.ucsf.edu/Breast/Hwang.data.xls")
#data<-data[!is.na(data[,2]),]
#data<-data[apply(is.na(data),1,sum)<=50,]
#data<-data[,apply(is.na(data),2,sum)<=1000]
#data$Position<-arrayinfo$Mb[match(toupper(as.character(data[,1])),toupper(as.character(arrayinfo[,1])))]
#data<-data[!is.na(data$Position),]
#dim(data)
#length(unique(paste(data$Chromosome, data$Position))) #there are repeated genomic locations
#data<-data[c(TRUE,data$Position[-1]!=data$Position[-1864]),] #discard probes with repeated genomic locations
#data<-data[data$Chromosome<=22,] #getting rid of X and Y chromosomes
#dataCNA<-CNA(as.matrix(data[,c(4:6,28:30)]),maploc=data$Position,chrom=data$Chromosome,sampleid=names(data)[c(4:6,28:30)]) #taking the first 3 patients only to shorten the computation time; use c(4:51) for the full dataset

#dataCNA$maploc<-dataCNA$maploc*1000 #transforming maploc to Kb scale
#dataCNA$chrom<- splitChromosomes( dataCNA$chrom,dataCNA$maploc)  #splits the chromosomes into arms


#ptlist<-substr(names(dataCNA)[-c(1,2)],1,4)
#samnms<-names(dataCNA)[-c(1,2)]

#results<-clonality.analysis(dataCNA, ptlist,  pfreq = NULL, refdata = NULL, nmad = 1.25, 
# reference = TRUE, allpairs = FALSE)
              
#genomewide plots of pairs of tumors from the same patient
#pdf(".pdf",height=7,width=11)
#for (i in unique(ptlist))
#{
#w<-which(ptlist==i) 
#ns<- length(w)
#if (ns>1)
#{
#for (p1 in c(1:(ns-1)))
#for (p2 in c((p1+1):ns))
#genomewidePlots(results$OneStepSeg, results$ChromClass,ptlist , ptpair=samnms[c(w[p1],w[p2])],results$LR,  plot.as.in.analysis = TRUE) 
#}
#}
#dev.off()



#pdf("hist.pdf",height=7,width=11)
#histogramPlot(results$LR[,4], results$refLR[,4])
#dev.off()


#for (i in unique(ptlist))
#{
#pdf(paste("pt",i,".pdf",sep=""),height=7,width=11)
#chromosomePlots(results$OneStepSeg, ptlist,ptname=i,nmad=1.25)
#dev.off()
#}
}
