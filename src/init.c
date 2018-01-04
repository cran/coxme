/*
** This file causes the entry points of my .C routines to be preloaded
** It adds one more layer of protection by declaring the number of arguments,
**  and perhaps a tiny bit of speed.  
*/
#include "coxmeS.h"
#include "R_ext/Rdynload.h"
#include "Rversion.h"

SEXP agfit6b (SEXP maxiter2,   SEXP beta2,   SEXP pmatb2,  SEXP pmatr2);
void agfit6d (Sint *nrefine,   double *beta, double *bhat, double *loglik);
SEXP bds_dsc (SEXP blocksize2, SEXP blocks2, SEXP rmat2,   SEXP dim2);
void coxfit6a(Sint *nused,      Sint *nvar,      Sint *ny,     
	      double *yy,       double *covar2,  double *offset2,
	      double *weights2, Sint  *nstrat,   Sint *strata,    
	      Sint *sorts,      Sint *fcol,      Sint *fx,
	      Sint *findex,
	      Sint *nblock,     Sint *bsize,     Sint *rsize,
	      double *means,    double *xscale,  Sint *method,
	      double *tolerch,  double *eps,     Sint *standard);
SEXP coxfit6b(SEXP maxiter2,    SEXP beta2,      SEXP pmatb2,  SEXP pmatr2);
void coxfit6c(double *u,        double *imatb,   double *imatr,
	      double *hinvb,    double *hinvr,   Sint  *rank,   Sint *ny);
void coxfit6d(Sint *nrefine,    double *beta,    double *bhat,  double *loglik); 
void coxfit6e(Sint *ny);

static const R_CMethodDef Centries[] = {
    {"Cagfit6d",  (DL_FUNC) &agfit6d,   4},
    {"Ccoxfit6a", (DL_FUNC) &coxfit6a, 22},
    {"Ccoxfit6c", (DL_FUNC) &coxfit6c,  7},
    {"Ccoxfit6d", (DL_FUNC) &coxfit6d,  4},
    {"Ccoxfit6e", (DL_FUNC) &coxfit6e,  1},
    {NULL, NULL, 0} };

static const R_CallMethodDef Callentries[] = {
    {"Cagfit6b",  (DL_FUNC) &agfit6b,  4},
    {"Ccoxfit6b", (DL_FUNC) &coxfit6b, 4},
    {"Cbds_dsc" , (DL_FUNC) &bds_dsc,  4},
    {NULL, NULL, 0} };

void R_init_coxme(DllInfo *dll) {
    R_registerRoutines(dll, Centries, Callentries, NULL, NULL);
 
   /* The following line makes only those routines defined above
       available to outside packages, i.e., internal things like
       dmatrix() are now invisible.
    */
    R_useDynamicSymbols(dll, FALSE); 

    /*
    ** This line makes them only be available via the symbols above
    **  i.e., .Call("agfit6b", ) won't work but .Call(Cagfit6b, )  will
    ** This feature was added in version 2.16
    */
#if defined(R_VERSION) && R_VERSION >= R_Version(2, 16, 0)
    R_forceSymbols(dll, TRUE);
#endif
}
