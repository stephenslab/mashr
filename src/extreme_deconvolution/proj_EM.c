/*
  NAME:
     proj_EM
  PURPOSE:
     goes through proj_EM
  CALLING SEQUENCE:
     proj_EM(struct datapoint * data, int N, struct gaussian * gaussians, 
     int K, bool * fixamp, bool * fixmean, bool * fixcovar, 
     double * avgloglikedata, double tol,long long int maxiter, 
     bool likeonly, double w,int partial_indx[3],double * qstarij,
     bool keeplog, FILE *logfile, FILE *tmplogfile, bool noproj, 
     bool diagerrs, bool noweight)
  INPUT:
     data         - the data
     N            - number of data points
     gaussians    - model gaussians
     K            - number of gaussians
     fixamp       - fix the amplitude?
     fixmean      - fix the mean?
     fixcovar     - fix the covariance?
     tol          - convergence limit
     maxiter      - maximum number of iterations
     likeonly     - only compute the likelihood?
     w            - regularization parameter
     keeplog      - keep a log in a logfile?
     logfile      - pointer to the logfile
     tmplogfile   - pointer to a tmplogfile to which the log likelihoods 
                    are written
     noproj       - don't perform any projections
     diagerrs     - the data->SS errors-squared are diagonal
     noweight     - don't use data-weights
  OUTPUT:
     avgloglikedata - average log likelihood of the data
  REVISION HISTORY:
     2008-09-21 - Written Bovy
     2010-03-01 Added noproj option - Bovy
     2010-04-01 Added noweight option - Bovy
*/
#include <stdio.h>
#include <math.h>
#include <proj_gauss_mixtures.h>

void proj_EM(struct datapoint * data, int N, struct gaussian * gaussians, 
	     int K,bool * fixamp, bool * fixmean, bool * fixcovar, 
	     double * avgloglikedata, double tol,long long int maxiter, 
	     bool likeonly, double w, bool keeplog, FILE *logfile,
	     FILE *tmplogfile, bool noproj, bool diagerrs, bool noweight){
  double diff = 2. * tol, oldavgloglikedata;
  int niter = 0;
  int d = (gaussians->mm)->size;
  halflogtwopi  = 0.5 * log(8. * atan(1.0));
  while ( diff > tol && niter < maxiter){
    proj_EM_step(data,N,gaussians,K,fixamp,fixmean,fixcovar,avgloglikedata,
		 likeonly,w,noproj,diagerrs,noweight);
    if (keeplog){
      fprintf(logfile,"%f\n",*avgloglikedata);
      fprintf(tmplogfile,"%f\n",*avgloglikedata);
      fflush(logfile);
      fflush(tmplogfile);
      //printf("%f\n",*avgloglikedata);
    }
    if (niter > 0){
      diff = *avgloglikedata - oldavgloglikedata;
      if (diff < 0){
	printf("Warning: log likelihood decreased by %g\n",diff);
	//fprintf(logfile,"oldavgloglike was %g\navgloglike is %g\n",oldavgloglikedata,*avgloglikedata);
      }
    }
    oldavgloglikedata = *avgloglikedata;
    if (likeonly) break;
    ++niter;
    //write_model("result.dat");
  }
  
 //post-processing: only the upper right of VV was computed, copy this to the lower left of VV
  int dd1,dd2,kk;
  for (kk = 0; kk != K; ++kk){
    for (dd1 = 0; dd1 != d; ++dd1)
      for (dd2 = dd1+1; dd2 != d ; ++dd2)
	gsl_matrix_set(gaussians->VV,dd2,dd1,
		       gsl_matrix_get(gaussians->VV,dd1,dd2));
    ++gaussians;
  }
  gaussians -= K;



  return;
}
