\name{Clonality-internal}
\alias{calc.freq}
\alias{calculateLR}
\alias{class.all}
\alias{CM.pvalue}
\alias{complexts}
\alias{func}
\alias{getbdry}
\alias{GL}
\alias{grantLR}
\alias{indiv.test}
\alias{LRpv}
\alias{LRts}
\alias{oneseg}
\alias{prettyplot}
\alias{refvector}
\alias{run.oneseg.param}
\alias{segment1}
 \alias{tumors}
 
\title{Internal Clonality functions}
\description{
  Internal functions of package Clonality.
}
\usage{
LRpv(tum, pp,nref)
LRts (tum,p, MX)
tumors(p,cc,m)
refvector(p,cc,m,nref=1000)
CM.pvalue (rw)
calc.freq(pfreq,refdata, classall,nmad,sbdry=NULL) 
calculateLR(data.seg1,classall,ptlist,pfreq,reference,allpairs=TRUE,gainthres,lossthres,sbdry)
getbdry(eta, nperm, max.ones, tol = 0.01)                                                              
run.oneseg.param(alpha)
segment1(xcna,alpha=0.01,sbdry=NULL)
oneseg(x, alpha=0.01, nperm, sbdry=NULL, data.type="logratio") 
GL(seg,nmad)
func(ss1,ss2)
complexts(s1,s2,maplocs)
indiv.test(ss1,ss2,func,gainthres,lossthres,Nsim=100,sbdry=sbdry)
grantLR(tum1,tum2,p,pi,pc,cvalue=0.5,rescale=TRUE,prnfile=NULL,nm=NULL,chrlist)
class.all(seg1)
prettyplot(datseg,path,nm,lab.general="",t1lab="Tumor 1",t2lab="Tumor 2",
many=TRUE,gains=NULL,losses=NULL)




}
\details{These are not to be called directly by the user}
\keyword{internal}
