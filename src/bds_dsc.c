/* Automatically generated from the noweb directory */
#include "coxmeS.h"
SEXP bds_dsc(SEXP blocksize2,  SEXP blocks2,  SEXP rmat2,
             SEXP dim2) {
    int i,j, iblock, bstart;
    int n, k, k2, rsize;
    int bs, nblock, rcol;
    
    /* pointers to input arguments */
    int *blocksize, dim;
    double *blocks, *rmat;
    
    /* output arguments */
    SEXP retlist, reti2, retp2, retx2;
    int *reti, *retp;
    double *retx;
    static const char *outnames[] = {"i", "p", "x", ""};
    
    /* Get sizes */
    blocksize = INTEGER(blocksize2);
    blocks = REAL(blocks2);
    rmat =  REAL(rmat2);
    dim = (INTEGER(dim2))[0];
    rcol = ncols(rmat2);
    
    nblock = LENGTH(blocksize2);  /* number of blocks */
    n = LENGTH(blocks2);          /* total number of non-zero elements */
    rsize = rcol*dim - (rcol*(rcol-1))/2;  /* the dense part */
 
    /* create output objects */
    PROTECT(reti2 = allocVector(INTSXP, n + rsize));
    reti = INTEGER(reti2);
    PROTECT(retp2 = allocVector(INTSXP, dim+1));
    retp = INTEGER(retp2);
    PROTECT(retx2 = allocVector(REALSXP, n + rsize));
    retx = REAL(retx2);
    
    k=0;  /* total elements processed */
    bstart =0;  /* row number for start of block */
    *retp =0;
    for (iblock=0; iblock<nblock; iblock++) {
        bs = blocksize[iblock];
        for (i=0; i<bs; i++) {   /* column in the output */
            retp[1] = *retp + i + 1; retp++;
            k2 = i+k;                        
            for (j=0; j<=i; j++) { /* row in the output*/
                 *retx++ = blocks[k2];
                *reti++ = bstart + j;
                k2 += bs - (j+1);
                }
            }
        bstart += bs;
        k += bs * (bs+1)/2;
        }
    
    /* Now do the rmat portion, if present 
       But not the lower right corner of rmat
    */
    k = 1+ dim - rcol;
    for (i=0; i<rcol; i++) {
        retp[1] = *retp + k; retp++;
        for (j=0; j<k; j++) {
            *retx++ = rmat[j];
            *reti++ = j;
            }
        rmat += dim;
        k++;
        }
                    
    retlist = PROTECT(mkNamed(VECSXP, outnames));
    SET_VECTOR_ELT(retlist, 0, reti2);
    SET_VECTOR_ELT(retlist, 1, retp2);
    SET_VECTOR_ELT(retlist, 2, retx2);
    UNPROTECT(4);
    return(retlist);
    }
