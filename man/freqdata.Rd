\name{freqdata}
\alias{freqdata}

\title{
 TCGA pancancer mutation frequencies
}
\description{
Mutational data from TCGA
}
\usage{
  data(freqdata)
}

\details{
  Object 'freqdata' contains the frequencies of the mutations observed in the exome sequencing data in 3 cancer types, COAD, LUAD, and BRCA, profiled by TCGA. 
  
  
  Full set of 33 cancer types is available by loading a full object in GitHub 'load(url("https://github.com/IOstrovnaya/MutFreq/blob/master/freqdata.RData?raw=true"))'. After that proceed with the analysis in the same way.
  
  There are 3 columns for 3 cancer types abbreviated in TCGA as COAD, LUAD, and BRCA  (for the full list of abbreviations see https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations). The first column is the number of patients profiled in each cancer type. Each subsequent row is a mutation, where number of patients with this particular mutation in each cancer type is given. The mutation ID, contained in the row names of 'freqdata', is of the following format: {"Chromosome Location RefAllele AltAllele"}, each entry separated by space, where chromosome is a number 1-22 or "X" or "Y"; location is genomic location in  GRCh37 build; RefAllele is a reference allele and AltAllele is Alternative allele/Tumor_Seq_Allele2. For example "10 100003849 G A", is the mutation at chromosome 10, genomic location 100003849, where reference allele G is substituted with A, or "10 100011448 - CCGCTGCAAT" is the insertion of "CCGCTGCAAT" at chromosome 10, location 100011448. The ref and alt alleles follow standard TCGA maf file notations.
  
The code that was used to obtain freqdata is available in the vignette.
}
\references{
Ostrovnaya, Irina, Venkatraman E. Seshan, and Colin B. Begg. 2015. USING SOMATIC MUTATION DATA TO TEST TUMORS FOR CLONAL RELATEDNESS. The Annals of Applied Statistics 9 (3): 1533-48. 

Ellrott K, Bailey MH, Saksena G, Covington KR, Kandoth C, Stewart C, Hess J, Ma S, Chiotti KE, McLellan M, Sofia HJ, Hutter C, Getz G, Wheeler D, Ding L; MC3 Working Group; Cancer Genome Atlas Research Network,  Scalable Open Science Approach for Mutation Calling of Tumor Exomes Using Multiple Genomic Pipelines. Cell Syst. 2018 Mar 28;6(3):271-281.e7. doi: 10.1016/j.cels.2018.03.002. PubMed PMID: 29596782.
}

