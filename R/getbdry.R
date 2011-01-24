getbdry <-
function (eta, nperm, max.ones, tol = 0.01)
{
    bdry <- rep(0, max.ones * (max.ones + 1)/2)
    zz <- .Fortran("getbdry", as.double(eta), as.integer(max.ones),
        as.integer(nperm), as.integer(max.ones * (max.ones +
            1)/2), bdry = as.integer(bdry), etastr = double(max.ones),
        as.double(tol), PACKAGE = "DNAcopy")
    zz$bdry
}

