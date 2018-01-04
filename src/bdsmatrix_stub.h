/*
** Define all the bdsmatrix routines
*/
#include "bdsmatrix.h"
#include <R.h>
#include <R_ext/Rdynload.h>

/*
** Rdynload defines "typedef void * (*DL_FUNC) (), which is just right
**  for almost all the routines.  Add two more.
*/
typedef int * (*INT_FUNC) ();
typedef double* (*DBL_FUNC) ();

void bdsmatrix_prod2(int nblock,     int *bsize,     int nrow,
                     double *bmat,   double *rmat,  
                     double *y,      double *result, int *itemp) {
    static DL_FUNC fun = NULL;
    if (fun==NULL) {
	fun = R_GetCCallable("bdsmatrix", "bdsmatrix_prod2");
	if (fun==NULL) Rf_error("cannot find function 'bsdmatrix_prod2'");
	}
    fun(nblock, bsize, nrow, bmat, rmat, y, result, itemp);
    }

void bdsmatrix_prod4(int nrow,    int nblock,   int *bsize, 
                    double *bmat, double *rmat,    
                    int nfrail,   double *y) {
    static DL_FUNC fun = NULL;
    if (fun==NULL) {
	fun = R_GetCCallable("bdsmatrix", "bdsmatrix_prod4");
	if (fun==NULL) Rf_error("cannot find function 'bsdmatrix_prod4'");
	}
    fun(nrow, nblock, bsize, bmat, rmat, nfrail, y);
    }

int cholesky4(double **matrix,  int n,          int nblock,     int *bsize,
              double *bd,       double toler) {
    static INT_FUNC fun = NULL;
    if (fun==NULL) {
	fun= R_GetCCallable("bdsmatrix", "cholesky4");
 	if (fun==NULL) Rf_error("cannot find function 'cholesky4'");
	}
    return(fun(matrix, n, nblock, bsize, bd, toler));
    }

int cholesky5(double **matrix,  int n,          double toler){
    static INT_FUNC fun =NULL;
    if (fun==NULL) {
	fun=  R_GetCCallable("bdsmatrix", "cholesky5");
  	if (fun==NULL) Rf_error("cannot find function 'cholesky5'");
	}
    return(fun(matrix, n, toler));
    }
    
void chinv4(double **matrix,    int n,          int nblock,     int *bsize, 
            double *bd,         int flag) {
    static DL_FUNC fun = NULL;
    if (fun==NULL){
	fun=  R_GetCCallable("bdsmatrix", "chinv4");
  	if (fun==NULL) Rf_error("cannot find function 'chinv4'");
	}
    fun(matrix, n, nblock, bsize, bd, flag);
    }

void chinv5(double **matrix ,   int n, int flag) {
    static DL_FUNC fun = NULL;
    if (fun==NULL){
	fun= R_GetCCallable("bdsmatrix", "chinv5");
  	if (fun==NULL) Rf_error("cannot find function 'chinv5'");
	}
    fun(matrix, n, flag);
    }
 
void chsolve4(double **matrix,  int n,          int nblock,     int *bsize,
              double *bd,       double *y,      int flag){
    static DL_FUNC fun = NULL;
    if (fun==NULL){
	fun= R_GetCCallable("bdsmatrix", "chsolve4");
  	if (fun==NULL) Rf_error("cannot find function 'chsolve4'");
	}
     fun(matrix, n, nblock, bsize, bd, y, flag);
    }

void chsolve5(double **matrix,  int n, double *y,int flag){
    static DL_FUNC fun = NULL;
    if (fun==NULL) {
	fun= R_GetCCallable("bdsmatrix", "chsolve5");
  	if (fun==NULL) Rf_error("cannot find function 'chsolve5'");
	}
    fun(matrix, n, y, flag);
    }

double **dmatrix(double *array, int ncol, int nrow){
    static DBL_FUNC fun = NULL;
    if (fun==NULL) {
	fun= R_GetCCallable("bdsmatrix", "dmatrix");
  	if (fun==NULL) Rf_error("cannot find function 'dmatrix'");
	}
    return(fun(array, ncol, nrow));
    }
