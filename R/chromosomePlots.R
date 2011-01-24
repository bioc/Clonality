chromosomePlots <-
function(data.seg1,ptlist,ptname,nmad)
{ 
samnms<-names(data.seg1$data)[-c(1,2)]
chrlist<-unique(data.seg1$data$chrom)
data.seg1GL<-GL(data.seg1,nmad)
gainthres<-data.seg1GL[[2]]
lossthres<-data.seg1GL[[3]]

datseg2<-subset(data.seg1,sample=samnms[ptlist==ptname])

ns2<-dim(datseg2$data)[2]-2
samnms2<-names(datseg2$data)[-c(1,2)]


upper<-gainthres[ptlist==ptname]
lower<-lossthres[ptlist==ptname]
   int.dev <- dev.interactive()
     parask <- par("ask")
    if (int.dev & !parask )
        par(ask = TRUE)

par(mfrow=c(ns2,1))
for (chr in chrlist)
{
for (p1 in c(1:ns2))
{
plot(subset(datseg2,sample=samnms2[p1],chrom=chr))
abline(h=c(upper[p1],lower[p1]),lty=2)
}
mtext(paste("Patient",ptname,"-",chr), side = 1, line = 2)

}
}

