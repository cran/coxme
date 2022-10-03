/*
** Define all the bdsmatrix routines
*/
#include "bdsmatrix.h"
#include <R.h>
#include <R_ext/Rdynload.h>

typedef void (*bdsmatrix_prod2_func)
                    (int nblock,     int *bsize,     int nrow,
                     double *bmat,   double *rmat,  
                     double *y,      double *result, int *itemp);

void bdsmatrix_prod2(int nblock,     int *bsize,     int nrow,
                     double *bmat,   double *rmat,  
                     double *y,      double *result, int *itemp) {
    static bdsmatrix_prod2_func fun = NULL;
    if (fun==NULL) {
	fun = (bdsmatrix_prod2_func) R_GetCCallable("bdsmatrix", "bdsmatrix_prod2");
	if (fun==NULL) Rf_error("cannot find function 'bsdmatrix_prod2'");
	}
    fun(nblock, bsize, nrow, bmat, rmat, y, result, itemp);
    }


typedef void (*bdsmatrix_prod4_func)
                   (int nrow,    int nblock,   int *bsize, 
                    double *bmat, double *rmat,    
                    int nfrail,   double *y);

void bdsmatrix_prod4(int nrow,    int nblock,   int *bsize, 
                    double *bmat, double *rmat,    
                    int nfrail,   double *y) {
    static bdsmatrix_prod4_func fun = NULL;
    if (fun==NULL) {
	fun = (bdsmatrix_prod4_func) R_GetCCallable("bdsmatrix", "bdsmatrix_prod4");
	if (fun==NULL) Rf_error("cannot find function 'bsdmatrix_prod4'");
	}
    fun(nrow, nblock, bsize, bmat, rmat, nfrail, y);
    }

typedef int (*cholesky4_func)
             (double **matrix,  int n,          int nblock,     int *bsize,
              double *bd,       double toler);

int cholesky4(double **matrix,  int n,          int nblock,     int *bsize,
              double *bd,       double toler) {
    static cholesky4_func fun = NULL;
    if (fun==NULL) {
	fun = (cholesky4_func) R_GetCCallable("bdsmatrix", "cholesky4");
 	if (fun==NULL) Rf_error("cannot find function 'cholesky4'");
	}
    return(fun(matrix, n, nblock, bsize, bd, toler));
    }

typedef int (*cholesky5_func)
             (double **matrix,  int n,          double toler);

int cholesky5(double **matrix,  int n,          double toler){
    static cholesky5_func fun =NULL;
    if (fun==NULL) {
	fun = (cholesky5_func) R_GetCCallable("bdsmatrix", "cholesky5");
  	if (fun==NULL) Rf_error("cannot find function 'cholesky5'");
	}
    return(fun(matrix, n, toler));
    }

typedef void (*chinv4_func)
           (double **matrix,    int n,          int nblock,     int *bsize, 
            double *bd,         int flag);

void chinv4(double **matrix,    int n,          int nblock,     int *bsize, 
            double *bd,         int flag) {
    static chinv4_func fun = NULL;
    if (fun==NULL){
	fun = (chinv4_func) R_GetCCallable("bdsmatrix", "chinv4");
  	if (fun==NULL) Rf_error("cannot find function 'chinv4'");
	}
    fun(matrix, n, nblock, bsize, bd, flag);
    }

typedef void (*chinv5_func)
           (double **matrix ,   int n, int flag);

void chinv5(double **matrix ,   int n, int flag) {
    static chinv5_func fun = NULL;
    if (fun==NULL){
	fun = (chinv5_func) R_GetCCallable("bdsmatrix", "chinv5");
  	if (fun==NULL) Rf_error("cannot find function 'chinv5'");
	}
    fun(matrix, n, flag);
    }
 
typedef void (*chsolve4_func)
             (double **matrix,  int n,          int nblock,     int *bsize,
              double *bd,       double *y,      int flag);

void chsolve4(double **matrix,  int n,          int nblock,     int *bsize,
              double *bd,       double *y,      int flag){
    static chsolve4_func fun = NULL;
    if (fun==NULL){
	fun = (chsolve4_func) R_GetCCallable("bdsmatrix", "chsolve4");
  	if (fun==NULL) Rf_error("cannot find function 'chsolve4'");
	}
    fun(matrix, n, nblock, bsize, bd, y, flag);
    }

typedef void (*chsolve5_func)
             (double **matrix,  int n, double *y,int flag);

void chsolve5(double **matrix,  int n, double *y,int flag){
    static chsolve5_func fun = NULL;
    if (fun==NULL) {
	fun= (chsolve5_func) R_GetCCallable("bdsmatrix", "chsolve5");
  	if (fun==NULL) Rf_error("cannot find function 'chsolve5'");
	}
    fun(matrix, n, y, flag);
    }

typedef double **(*dmatrix_func)
                (double *array, int ncol, int nrow);

double **dmatrix(double *array, int ncol, int nrow){
    static dmatrix_func fun = NULL;
    if (fun==NULL) {
	fun = (dmatrix_func) R_GetCCallable("bdsmatrix", "dmatrix");
  	if (fun==NULL) Rf_error("cannot find function 'dmatrix'");
	}
    return(fun(array, ncol, nrow));
    }
