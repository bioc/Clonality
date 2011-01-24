oneseg <-
function(x, alpha=0.01, nperm, sbdry=NULL, data.type="logratio") 
{
# setup data and parameters for fortran code
  current.genomdat <- x
  n <- current.n <- window.size <- wsize <- length(x)

if (is.null(sbdry))
{
 alpha.cbs=alpha; nperm.cbs.ostat=nperm;
max.ones.ostat <- floor(nperm.cbs.ostat * alpha.cbs) + 1
           sbdry <- getbdry(eta=0.05, nperm.cbs.ostat, max.ones.ostat)
}

  sbn <- length(sbdry)
  winnum <- 1
  winloc <- 0

# there are used only by hybrid code but need to set it
  kmax <- 25
  nmin <- 200
  ngrid <- 100
  tol <- 1e-6


if (n>200)  hybrid <- TRUE else  hybrid <- FALSE
 if (n>200) delta<-(kmax+1)/current.n else delta <- 0
            current.genomdat <- current.genomdat - mean(current.genomdat)
 current.tss <- sum(current.genomdat^2)
zzz<-NULL
#
# call the fortran segmentation cose
    zzz <- .Fortran("fndcpt",
n = as.integer(current.n),
             #   twon = as.integer(2 * n),
x = as.double(current.genomdat),
                tss = as.double(current.tss),
px = double(current.n),
                sx = double(n), #tx = double(2 * n),
 nperm = as.integer(nperm),
                cpval = as.double(alpha),
ncpt = integer(1),
                icpt = integer(2),
ibin = as.logical(data.type ==     "binary"),
hybrid = as.logical(hybrid),
al0 = as.integer(2),##
hk = as.integer(kmax),
                delta = as.double(delta),
ngrid = as.integer(ngrid),
                sbn = as.integer(sbn),
sbdry = as.integer(sbdry),
                tol = as.double(tol),
PACKAGE = "DNAcopy")


# number of changes detected at alpha level (could be 0, 1 or 2)
  switch(1+zzz$ncpt, c(zzz$ncpt, rep(0,4)),
         c(zzz$ncpt, 1, zzz$icpt[1], mean(x[1:zzz$icpt[1]]), mean(x[(zzz$icpt[1]+1):n])),
         c(zzz$ncpt, c(1,0)+zzz$icpt, mean(x[(zzz$icpt[1]+1):zzz$icpt[2]]), mean(x[-((zzz$icpt[1]+1):zzz$icpt[2])])))
}

