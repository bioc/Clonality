run.oneseg.param <-
function(alpha)
 {
######parameters for one step seegmentation. Took them out of function because they take too long.
alpha.cbs=alpha<-0.01; nperm.cbs.ostat=2000;
max.ones.ostat <- floor(nperm.cbs.ostat * alpha.cbs) + 1
sbdry <- getbdry(eta=0.05, nperm.cbs.ostat, max.ones.ostat)

alpha.cbs=alpha<-0.001; nperm.cbs.ostat=2000;
max.ones.ostat <- floor(nperm.cbs.ostat * alpha.cbs) + 1
sbdry <- getbdry(eta=0.05, nperm.cbs.ostat, max.ones.ostat)

alpha.cbs=alpha<-0.05; nperm.cbs.ostat=2000;
max.ones.ostat <- floor(nperm.cbs.ostat * alpha.cbs) + 1
sbdry <- getbdry(eta=0.05, nperm.cbs.ostat, max.ones.ostat)
sbdry
 }

