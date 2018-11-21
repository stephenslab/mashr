/*
  NAME:
     proj_gauss_mixtures_IDL
  PURPOSE:
     run the projected gaussian mixtures algorithm from R
	[Patched from the original IDL (or python) interface by Gao Wang]
  CALLING SEQUENCE:
     see R wrapper
  INPUT:
     from R wrapper
  OUTPUT:
     updated model gaussians and average loglikelihood, see R WRAPPER
  REVISION HISTORY:
     2008-09-21 - Written Bovy
     2010-03-01 Added noproj option - Bovy
     2010-04-01 Added noweight option and logweights - Bovy
     2015-08-08 Patched by Gao Wang to interface with R instead
*/
#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "proj_gauss_mixtures.h"

bool * int2bool(int * a, int K)
{
  bool * x = (bool *) malloc(K * sizeof (bool));
  int kk;
  for (kk = 0; kk != K; ++kk) *(x++) = (bool) *(a++);
  x -=K;
  return x;
}

int proj_gauss_mixtures_IDL(double * ydata, double * ycovar, 
			    double * projection, double * logweights,
			    int *p_N, int *p_dy,
			    double * amp, double * xmean,
			    double * xcovar, int *p_d, int *p_K,
			    int * fixamp_int, int * fixmean_int, int * fixcovar_int,
			    double * avgloglikedata, double *p_tol,
			    int *p_maxiter, int *p_likeonly, double *p_w,
			    int * logfilename, int *p_slen, int *p_splitnmerge,
			    int * convlogfilename, int *p_convloglen,
			    int *p_noprojection, int *p_diagerrors,
			    int *p_noweights){
  // convert variables from R interface
  int N = *p_N, dy = *p_dy, d = *p_d, K = *p_K,
    maxiter = *p_maxiter, slen = *p_slen,
    splitnmerge = *p_splitnmerge, convloglen = *p_convloglen,
    likeonly = *p_likeonly, noprojection = *p_noprojection,
    diagerrors = *p_diagerrors, noweights = *p_noweights;
  double tol = *p_tol, w = *p_w;
  bool * fixamp = int2bool(fixamp_int, K);
  bool * fixmean = int2bool(fixmean_int, K);
  bool * fixcovar = int2bool(fixcovar_int, K);
  //Set up logfiles  
  bool keeplog = true;
  char logname[slen+1];
  char convlogname[convloglen+1];
  int ss;
  if (*logfilename == 0 || likeonly != 0 || slen == 0)
    keeplog = false;
  else {
    for (ss = 0; ss != slen; ++ss)
      logname[ss] = (char) *(logfilename++);
    for (ss = 0; ss != convloglen; ++ss)
      convlogname[ss] = (char) *(convlogfilename++);
    logfilename -= slen;
    convlogfilename -= convloglen;
    logname[slen] = '\0';
    convlogname[convloglen] = '\0';
  }

  if (keeplog) {
    logfile = fopen(logname,"a");
    if (logfile == NULL) return -1;
    convlogfile = fopen(convlogname,"w");
    if (convlogfile == NULL) return -1;
  }


  if (keeplog){
    time_t now;
    time(&now);
    fprintf(logfile,"#----------------------------------\n");
    fprintf(logfile,"#\n#%s\n",asctime(localtime(&now)));
    fprintf(logfile,"#----------------------------------\n");
    fflush(logfile);
  }
  
  //Copy everything into the right formats
  struct datapoint * data = (struct datapoint *) malloc( N * sizeof (struct datapoint) );
  struct gaussian * gaussians = (struct gaussian *) malloc (K * sizeof (struct gaussian) );

  bool noproj= (bool) noprojection;
  bool noweight= (bool) noweights;
  bool diagerrs= (bool) diagerrors;
  int ii, jj,dd1,dd2;
  for (ii = 0; ii != N; ++ii){
    data->ww = gsl_vector_alloc(dy);
    if ( ! noweight ) data->logweight = *(logweights++);
    if ( diagerrs ) data->SS = gsl_matrix_alloc(dy,1);
    else data->SS = gsl_matrix_alloc(dy,dy);
    if ( ! noproj ) data->RR = gsl_matrix_alloc(dy,d);
    for (dd1 = 0; dd1 != dy;++dd1)
      gsl_vector_set(data->ww,dd1,*(ydata++));
    if ( diagerrs)
      for (dd1 = 0; dd1 != dy; ++dd1)
	  gsl_matrix_set(data->SS,dd1,0,*(ycovar++));
    else
      for (dd1 = 0; dd1 != dy; ++dd1)
	for (dd2 = 0; dd2 != dy; ++dd2)
	  gsl_matrix_set(data->SS,dd1,dd2,*(ycovar++));
    if ( ! noproj )
      for (dd1 = 0; dd1 != dy; ++dd1)
	for (dd2 = 0; dd2 != d; ++dd2)
	  gsl_matrix_set(data->RR,dd1,dd2,*(projection++));
    else data->RR= NULL;
    ++data;
  }
  data -= N;
  ydata -= N*dy;
  if ( diagerrs ) ycovar -= N*dy;
  else ycovar -= N*dy*dy;
  if ( ! noproj ) projection -= N*dy*d;

  for (jj = 0; jj != K; ++jj){
    gaussians->mm = gsl_vector_alloc(d);
    gaussians->VV = gsl_matrix_alloc(d,d);
    gaussians->alpha = *(amp++);
    for (dd1 = 0; dd1 != d; ++dd1)
      gsl_vector_set(gaussians->mm,dd1,*(xmean++));
    for (dd1 = 0; dd1 != d; ++dd1)
      for (dd2 = 0; dd2 != d; ++dd2)
	gsl_matrix_set(gaussians->VV,dd1,dd2,*(xcovar++));
    ++gaussians;
  }
  gaussians -= K;
  amp -= K;
  xmean -= K*d;
  xcovar -= K*d*d;


  //Print the initial model parameters to the logfile
  int kk;
  if (keeplog){
    fprintf(logfile,"#\n#Using %i Gaussians and w = %f\n\n",K,w);
    fprintf(logfile,"#\n#Initial model parameters used:\n\n");
    for (kk=0; kk != K; ++kk){
      fprintf(logfile,"#Gaussian ");
      fprintf(logfile,"%i",kk);
      fprintf(logfile,"\n");
      fprintf(logfile,"#amp\t=\t");
      fprintf(logfile,"%f",(*gaussians).alpha);
      fprintf(logfile,"\n");
      fprintf(logfile,"#mean\t=\t");
      for (dd1=0; dd1 != d; ++dd1){
	fprintf(logfile,"%f",gsl_vector_get(gaussians->mm,dd1));
	if (dd1 < d-1) fprintf(logfile,"\t");
      }
      fprintf(logfile,"\n");
      fprintf(logfile,"#covar\t=\t");
      for (dd1=0; dd1 != d; ++dd1)
	fprintf(logfile,"%f\t",gsl_matrix_get(gaussians->VV,dd1,dd1));
      for (dd1=0; dd1 != d-1; ++dd1)
	for (dd2=dd1+1; dd2 != d; ++dd2){
	  fprintf(logfile,"%f\t",gsl_matrix_get(gaussians->VV,dd1,dd2));
	}
      ++gaussians;
      fprintf(logfile,"\n#\n");
    }
    gaussians -= K;
    fflush(logfile);
  }



  //Then run projected_gauss_mixtures
  proj_gauss_mixtures(data,N,gaussians,K,(bool *) fixamp,
		      (bool *) fixmean, (bool *) fixcovar,avgloglikedata,
		      tol,(long long int) maxiter, (bool) likeonly, w,
		      splitnmerge,keeplog,logfile,convlogfile,noproj,diagerrs,
		      noweight);


  //Print the final model parameters to the logfile
  if (keeplog){
    fprintf(logfile,"\n#Final model parameters obtained:\n\n");
    for (kk=0; kk != K; ++kk){
      fprintf(logfile,"#Gaussian ");
      fprintf(logfile,"%i",kk);
      fprintf(logfile,"\n");
      fprintf(logfile,"#amp\t=\t");
      fprintf(logfile,"%f",(*gaussians).alpha);
      fprintf(logfile,"\n");
      fprintf(logfile,"#mean\t=\t");
      for (dd1=0; dd1 != d; ++dd1){
	fprintf(logfile,"%f",gsl_vector_get(gaussians->mm,dd1));
	if (dd1 < d-1) fprintf(logfile,"\t");
      }
      fprintf(logfile,"\n");
      fprintf(logfile,"#covar\t=\t");
      for (dd1=0; dd1 != d; ++dd1)
	fprintf(logfile,"%f\t",gsl_matrix_get(gaussians->VV,dd1,dd1));
      for (dd1=0; dd1 != d-1; ++dd1)
	for (dd2=dd1+1; dd2 != d; ++dd2){
	  fprintf(logfile,"%f\t",gsl_matrix_get(gaussians->VV,dd1,dd2));
	}
      ++gaussians;
      fprintf(logfile,"\n#\n");
    }
    gaussians -= K;
    fflush(logfile);
  }



  //Then update the arrays given to us by IDL
  for (jj = 0; jj != K; ++jj){
    *(amp++) = gaussians->alpha;
    for (dd1 = 0; dd1 != d; ++dd1)
      *(xmean++) = gsl_vector_get(gaussians->mm,dd1);
    for (dd1 = 0; dd1 != d; ++dd1)
      for (dd2 = 0; dd2 != d; ++dd2)
	*(xcovar++) = gsl_matrix_get(gaussians->VV,dd1,dd2);
    ++gaussians;
  }
  gaussians -= K;
  amp -= K;
  xmean -= K*d;
  xcovar -= K*d*d;
  
  //And free any memory we allocated
  for (ii = 0; ii != N; ++ii){
    gsl_vector_free(data->ww);
    gsl_matrix_free(data->SS);
    if ( ! noproj )  gsl_matrix_free(data->RR);
    ++data;
  }
  data -= N;
  free(data);
  
  for (jj = 0; jj != K; ++jj){
    gsl_vector_free(gaussians->mm);
    gsl_matrix_free(gaussians->VV);
    ++gaussians;
  }
  gaussians -= K;
  free(gaussians);

  if (keeplog){
    fclose(logfile);
    fclose(convlogfile);
  }

  return 0;
}
