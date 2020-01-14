/*
** This call is for the monte carlo computation of the error in the
**   Laplace transform estimate.  The basic input is an estimate 
**   beta-hat of the coefficients, along with a matrix of trial values
**   for the random effects b-hat.  For each of the trial estimates, it
**   computes the partial likelihood.  No need for a first derivative, a
**   variance, or iteration --- so we can be very fast.
**   
**    Input
** beta         : vector of coefficients, random first, then others
** bhat         : matrix of trial values for the random coefficients
**
**    Output
** loglik       : vector of partial likelihoods
*/
#include "coxmeS.h"
#include "coxfit6.h"
#include <math.h>

/* the next line is just so that I can use "c6.n" instead of "coxfit6.n", etc*/
#define c6 coxfit6  

void agfit6d(Sint *nrefine,  double *beta,  double *bhat,
	      double *loglik) {
    int i,j,p;
    int person;
    int ii, istrat;
    int    nvar, nvar2;
    int    nfrail, nfac;
    int    nfns;    /* number of factors that are not sparse */
    int    indx2;
    double time;

    double  newlik;
    double  denom, zbeta, risk;
    double  temp, temp2;
    double  d2, efron_wt;
    double  ndead;

    nfrail = c6.nfrail;    /* number of penalized coefficients */
    nvar   = c6.nvar;      /* number of unpenalized coefficients */
    nfac   = c6.nfactor;   /* number of factor levels (penalized) */
    nvar2  = nvar + (nfrail - nfac);  /* number of cols of X */
    nfns   = nfrail - nfac;  /* penalized terms that are not factors */

    for (ii=0; ii< *nrefine; ii++) {
	/*
	** Loop through the data, and compute the loglik
	*/ 
	istrat=0;
	denom =0;
	efron_wt =0;
	newlik =0;
	indx2 =0;

	for (person=0; person<c6.n; person++) {  /* p = person */
	    p = c6.sort1[person];
	    if (person==c6.strata[istrat]) {
		indx2 =person;
		istrat++;
		efron_wt =0;
		denom = 0;
		}
	    /*
	    ** Form the linear predictor zbeta, and the risk score
	    **   For the sparse coefs use bhat, and beta for the others
	    */
	    zbeta = c6.offset[p];
	    for (i=0; i<c6.nfx; i++) {
		j = c6.fx[p + i*c6.n];  /* level of covariate i */
		zbeta = zbeta + bhat[j];
		}
	    for (i=0; i<nfns; i++)
		zbeta += bhat[i+nfac] * c6.x[i][p];
	    for (i=nfns; i<nvar2; i++)
		zbeta += beta[i+nfac]* c6.x[i][p];
	    risk = exp(zbeta) * c6.weights[p];
	    denom += risk;

	    /*
	    ** Extra terms for the deaths
	    */
	    if (c6.status[p]==1) {
		newlik += c6.weights[p] *zbeta;
		efron_wt += risk;

		/*
		** Take people out of the sum who are no longer at risk
		*/
		time = c6.stop[p];
		for (; indx2<c6.strata[istrat]; indx2++) {
		    p = c6.sort2[indx2];
		    if (c6.start[p] < time) break;
		    zbeta = c6.offset[p];
		    for (i=0; i<c6.nfx; i++) {
			j = c6.fx[p + i*c6.n];  /* level of covariate i */
			zbeta = zbeta + bhat[j];
			}
		    for (i=0; i<nfns; i++)
			zbeta += bhat[i+nfac]* c6.x[i][p];
		    for (i=nfns; i<nvar2; i++)
			zbeta += beta[i+nfac]* c6.x[i][p];
		    risk = exp(zbeta) * c6.weights[p];
		    denom -= risk;
		    }
		p = c6.sort1[person];  /* restore the pointer */
		}

	    if (c6.mark[p] >0) {  /* once per unique death time */
		ndead = c6.mark[p];
		if (c6.method==0 || ndead==1)  {
		    /*
		    ** Breslow approx 
		    */
		    temp = c6.wtave[p] * ndead;
		    newlik -= temp *log(denom);
		    }
		
		else {
		    /* 
		    ** Do the Efron approx 
		    */
		    for (temp2=0; temp2<ndead; temp2++) {
			temp = temp2/ ndead;
			d2= denom - temp*efron_wt;
			newlik -= c6.wtave[p] *log(d2);
			}
		    } /* end of Efron loop */
		
		efron_wt =0;
	        }   /* matches "if (mark[p] >0)"  */
	    } /* end  of accumulation loop  */

	loglik[ii] = newlik;
	bhat += nfrail;
	}   /* return for another iteration */

    return;
    }
