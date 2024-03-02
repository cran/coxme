#
# Look at the use of Matrix for a sparse QR
#
tmat <- bdsmatrix(c(3,2,2,4), 
              c(22,1,2,21,3,20,19,4,18,17,5,16,15,6,7, 8,14,9,10,13,11,12),
              matrix(c(1,0,1,1,0,0,1,1,0,1,0,10,0,
                       0,1,1,0,1,1,0,1,1,0,1,0,10), ncol=2))
dimnames(tmat) <- list(NULL, letters[1:13])

bdsToM <- function(x) {
    nblock <- length(x@blocksize)
    blockend <- cumsum(x@blocksize)
    blockstart <- c(1, blockend[-nblock] +1)
    tlist <- vector("list", length=nblock)
    for (i in 1:nblock) {
        indx <- blockstart[i]:blockend[i]
        tlist[[i]] <- as.matrix(x[indx, indx])
    }
    new <- bdiag(tlist)

    if (length(tmat@rmat >0)) {
        temp1 <- Matrix(tmat@rmat)
        temp2 <- t(temp1[1:nrow(new),])

        new <- cBind(rBind(new, as(temp2, "dgCMatrix")),
                     as(temp1, "dgCMatrix"))
    }
    new
}
tmat2 <- bdsToM(tmat)
tmat3 <- as('dsCMatrix', tmat2)


temp <- Matrix(c(1,2,3,4,
                 2,5,7,0,
                 4,7,8,0,
                 0,0,0,12), 4)
temp2 <- as(temp, "dgCMatrix")
                           
cmat <- new("dgCMatrix", i=as.integer(c(0,1,2,3, 0,1,2, 0,1,2, 3)),
                         p=as.integer(c(0, 4, 7, 10,11)),
                         Dim=as.integer(c(4,4)),
                         Dimnames=list(NULL, NULL),
                         x=c(1,2,0,4, 2,5,7, 4,7,8,12),
            factors=list())
