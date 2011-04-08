\name{clonality.analysis}
\alias{clonality.analysis}

\title{  Clonality testing using copy number data

}
\description{
Function to test clonality of two tumors from the same patient based on their genomewide copy number profiles. This function calculates likelihood ratios and the reference distribution  under the hypothesis of independence.
}
\usage{
clonality.analysis(data, ptlist, pfreq = NULL, refdata = NULL, nmad = 1.25, reference = TRUE, allpairs = TRUE,segmethod="oneseg", segpar=NULL)
}

\arguments{
  \item{data}{
Copy Number Array object (output of function CNA() from package DNAcopy). First column contains chromosomes, second column contains genomic locations. Each remaining column contains log-ratios from a particular tumor or sample. Chromosomes X and Y should be removed prior to analysis, and chromosomes should be split into p and q arms to improve the power (use function splitChromosomes()).
}
  \item{ptlist}{
Vector of the patient IDs in the order of the samples appearing in the data.  For example, if the first three tumors (columns 3, 4, 5 of data) belong to patient A, and the following two (columns 6, 7 of data) belong to patient B, then ptlist=c('ptA', 'ptA', 'ptA', 'ptB', 'ptB'). Note that while sample names in data should be unique the ptlist should have repeated labels.
}

  \item{pfreq}{
Marginal frequencies of Gains, Losses and Normals for all the chromosomes. If it is not known, pfreq should be set to NULL and frequencies will be estimated from all the samples in the dataset.  If frequencies are known, pfreq should be a data frame with 4 columns: 1) chromosome arm in the format 'chr01p', probability of 2) gain, 3) loss and 4) normal.
}
  \item{refdata}{
If available,  additional cohort of patients with the same disease that should be used to estimate the marginal gain/loss frequencies. If NULL, the original set of tumors is used, otherwise, refdata should be a CNA object. It will be segmented with 1 step CBS and each chromosome will be classified as gain/loss as described in the manuscript, leading to frequency estimates. No averaging or chromosome splitting is done for this dataset, so users should make sure refdata has chromosomes in the format 'chr01p' and that its resolution is similar to the one of the original data.
}
  \item{nmad}{
Number of MADs (median absolute deviations) that is used for Gain/Loss calls.  For each array MAD of its residuals (that is, data minus segmentation means) is calculated. Residuals  represent the array's noise revel. Any segment of this array that has a mean at least nmad MADs above or below array's median is called a gain or a loss. We use value of 1.25, while values in the range of 0.5 to 2 can also be admissible depending on the resolution and presence of artifacts.
}

  \item{reference}{
If TRUE the reference distribution of likelihood ratios is created under  hypothesis of independence by pairing (independent) tumors from different patients.
}
  \item{allpairs}{
If TRUE  all possible pairs of tumors from different patients will be used for reference distribution. If two tumors in a pair are not  exchangeable, for example primary tumor  vs recurrence, or pre-cancerous lesion vs tumor, then allpairs should be set to FALSE and the 'first' tumor should always come earlier in the data before the 'second' tumor for all the patients. Then 'first' tumors of patients will only be paired with 'second' tumors of other patients for the reference distribution.
}
\item{segmethod}{
  The segmentation algorithm to be used. The default is "oneseg" which
  uses the built in function of the same name based on the CBS
  algorithm. An alternative segmentation algorithm can be used. A function
  should be created and the name passed as described in the vignette.
}
\item{segpar}{
  The parameters necessary for the segmentation algorithm as a list. For
  "oneseg" you can specify alpha (default = 0.01) and nperm (default =
  2000) necessary for the CBS algorithm.
}
}
\details{
The function implements the statistical procedure designed to distinguish whether the two tumors from the same patient are clonal (have the same progenitor cancer cell) or independent (developed from normal cells independently). At first  data are segmented with one step CBS (Olshen, A. B., Venkatraman, E. S., Lucito, R., Wigler, M. (2004). Circular binary segmentation for the analysis of array-based DNA copy number data. Biostatistics 5: 557-572) that picks at most one copy number change per chromosome arm. Then each chromosome arm is classified as Gain/Loss/Normal based  on a middle segment if there are 3 segments, or  based  on the most outstanding segment if there are 2 segments. The multinomial likelihood ratio comparing these classifications is computed (LR1). For each concordant partial arm gain or loss  we also calculate likelihood ratio that this change is exactly the same in both tumors. These likelihood ratios are multiplied by LR1 to obtain our final statistic, LR2. If LR2 is much greater than 1, that indicates clonality. If LR2 is much smaller than 1, it indicates independence. The reference distribution of LR2 under the hypothesis of independence is obtained  by pairing up tumors from different patients, which are independent by default.

Since only one gain/loss is admissible per chromosome arm it is highly recommended to apply this methodology to arrays with at most 10,000-15,000 markers.  We suggest averaging  blocks of consecutive probes for arrays with larger resolution, see function ave.adj.probes.
}
\value{
If the reference is TRUE, function returns the list with 4 elements: LR, OneStepSeg, ChromClass, refLR.

LR - matrix with the within patient comparisons. Each row corresponds to a pair of samples being compared. Columns are:
Sample1 - name of sample 1;
Sample2 - name of sample 2;
LR1 - likelihood ratio without comparisons of specific concordant gains/losses;
LR2  - final likelihood ratio with individual comparisons;
GGorLL - number of chromosome arms that are classified as Gains in both tumors or Losses in both tumors;
NN    - number of chromosome arms  that are classified as Normal in both tumors;
GL     - number of chromosome arms  that are classified as Gain in one tumors and Loss in another;
GNorLN    -  number of chromosome arms  that are classified as Gain(Loss) in one tumors and Normal in another;
IndividualComparisons     - list of chromosome arms that had comparisons of specific concordant gains/losses in both tumors and the corresponding likelihood ratio for them being exactly the same.
p-value - quantile of the reference distribution under the null hypothesis (refLR$LR2) that the value of LR2 match.

OneStepSeg  - is the output of one step  segmentation of the data. It has the same structure as the output of 'segment' from DNAcopy, but only one most prominent change per arm is allowed. 

ChromClass  - is the matrix of chromosome classifications based on the one step segmentation. Rows correspond to chromosome arms, columns correspond to samples. Chromosome arms are classified by the middle segment if there are 3 segments, and by the most outstanding segment if there are 2 segments.
 
refLR - matrix with the between patient comparisons (reference distribution under the hypothesis of independence). Has the same structure as LR but the pairs of tumors are selected from different patients.

Note that calculating the reference distribution might take along time.

If the  reference is FALSE,  there is no p-value column in LR and no refLR output.

}
\references{
Ostrovnaya, I., Olshen, A. B., Seshan, V.E., Orlow, I., Albertson, D. G. and Begg, C. B. (2010), A metastasis or a second independent 
cancer? Evaluating the clonal origin of tumors using array copy number data. Statistics in Medicine, 29: 1608-1621

Ostrovnaya, I. and Begg, C.  Testing Clonal Relatedness of Tumors Using Array Comparative Genomic Hybridization: A Statistical Challenge Clin Cancer Res March 1, 2010 16:1358-1367

Venkatraman, E. S. and Olshen, A. B. (2007). A faster circular binary segmentation algorithm for the analysis of array CGH data. 
Bioinformatics, 23:657 63.

Olshen, A. B., Venkatraman, E. S., Lucito, R., Wigler, M. (2004). Circular binary segmentation for the analysis of array-based DNA copy
number data. Biostatistics 5: 557-572.

}
\author{
Irina Ostrovnaya \email{ostrovni@mskcc.org}
}

\examples{
   #Analysis of paired breast samples from study
#Hwang ES, Nyante SF, Chen YY, Moore D, DeVries S, Korkola JE, Esserman LJ, and Waldman FM. 
#Clonality of lobular carcinoma in situ and synchronous invasive lobular cancer. Cancer 100(12):2562-72, 2004.
library(gdata)    #needed to read .xls files
library(DNAcopy) 
arrayinfo<-read.xls("http://waldman.ucsf.edu/Colon/nakao.data.xls")  #needed to extract  genomic locations
data<-read.xls("http://waldman.ucsf.edu/Breast/Hwang.data.xls")
data<-data[!is.na(data[,2]),]
data<-data[apply(is.na(data),1,sum)<=50,]
data<-data[,apply(is.na(data),2,sum)<=1000]
data$Position<-arrayinfo$Mb[match(toupper(as.character(data[,1])),toupper(as.character(arrayinfo[,1])))]
data<-data[!is.na(data$Position),]
dim(data)
length(unique(paste(data$Chromosome, data$Position))) #there are repeated genomic locations
data<-data[c(TRUE,data$Position[-1]!=data$Position[-1864]),] #discard probes with repeated genomic locations
data<-data[data$Chromosome<=22,] #getting rid of X and Y chromosomes
dataCNA<-CNA(as.matrix(data[,c(4:6,28:30)]),maploc=data$Position,chrom=data$Chromosome,sampleid=names(data)[c(4:6,28:30)]) #taking the first 3 patients only to shorten the computation time; use c(4:51) for the full dataset

dataCNA$maploc<-dataCNA$maploc*1000 #transforming maploc to Kb scale
dataCNA$chrom<- splitChromosomes( dataCNA$chrom,dataCNA$maploc)  #splits the chromosomes into arms


ptlist<-substr(names(dataCNA)[-c(1,2)],1,4)
samnms<-names(dataCNA)[-c(1,2)]

results<-clonality.analysis(dataCNA, ptlist,  pfreq = NULL, refdata = NULL, nmad = 1.25, 
 reference = TRUE, allpairs = FALSE)
              
#genomewide plots of pairs of tumors from the same patient
pdf("genomewideplots.pdf",height=7,width=11)
for (i in unique(ptlist))
{
w<-which(ptlist==i) 
ns<- length(w)
if (ns>1)
{
for (p1 in c(1:(ns-1)))
for (p2 in c((p1+1):ns))
genomewidePlots(results$OneStepSeg, results$ChromClass,ptlist , ptpair=samnms[c(w[p1],w[p2])],results$LR,  plot.as.in.analysis = TRUE) 
}
}
dev.off()



pdf("hist.pdf",height=7,width=11)
histogramPlot(results$LR[,4], results$refLR[,4])
dev.off()


for (i in unique(ptlist))
{
pdf(paste("pt",i,".pdf",sep=""),height=7,width=11)
chromosomePlots(results$OneStepSeg, ptlist,ptname=i,nmad=1.25)
dev.off()
}

}