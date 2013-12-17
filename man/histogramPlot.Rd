\name{histogramPlot}
\alias{histogramPlot}

\title{
Histrograms of Log-Likelihood Ratios
}
\description{
Function produces the histograms of the within-patient and between-patient log-Likelihood Ratios.
}
\usage{
histogramPlot(ptLRvec, refLRvec)
}

\arguments{
  \item{ptLRvec}{
Vector with the within-patient likelihood ratios -  output LR of clonality.analysis()
}
  \item{refLRvec}{
Vector with the between-patient likelihood ratios -  output refLR of clonality.analysis()
}
}
\details{
Functions plots two overlapping histograms:  within-patient log-likelihood ratios are in red and between-patient log-likelihood ratios (reference distribution under the hypothesis of independence) are in black.
}

\examples{
 # Same example as in clonality.analysis()
   #Analysis of paired breast samples from study
#Hwang ES, Nyante SF, Chen YY, Moore D, DeVries S, Korkola JE, Esserman LJ, and Waldman FM. 
#Clonality of lobular carcinoma in situ and synchronous invasive lobular cancer. Cancer 100(12):2562-72, 2004.
#library(gdata)    #needed to read .xls files
#library(DNAcopy) 
#arrayinfo<-read.xls("http://waldmanlab.org/Colon/nakao.data.xls")  #needed to extract  genomic locations
#data<-read.xls("http://waldmanlab.org/Breast/Hwang.data.xls")
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
#pdf("genomewideplots.pdf",height=7,width=11)
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
