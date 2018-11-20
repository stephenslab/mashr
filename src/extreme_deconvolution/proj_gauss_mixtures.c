/*
  NAME:
     proj_gauss_mixtures
  PURPOSE:
     runs the full projected gaussian mixtures algorithms, with 
     regularization and split and merge
  CALLING SEQUENCE:
     proj_gauss_mixtures(struct datapoint * data, int N, 
     struct gaussian * gaussians, int K,bool * fixamp, bool * fixmean, 
     bool * fixcovar, double * avgloglikedata, double tol,
     long long int maxiter, bool likeonly, double w, int splitnmerge, 
     bool keeplog, FILE *logfile, FILE *convlogfile, bool noproj, 
     bool diagerrs,noweight)
  INPUT:
     data        - the data
     N           - number of datapoints
     gaussians   - model gaussians (initial conditions)
     K           - number of model gaussians
     fixamp      - fix the amplitude?
     fixmean     - fix the mean?
     fixcovar    - fix the covar?
     tol         - proj_EM convergence limit
     maxiter     - maximum number of iterations in each proj_EM
     likeonly    - only compute the likelihood?
     w           - regularization paramter
     splitnmerge - split 'n' merge depth (how far down the list to go)
     keeplog     - keep a log in a logfile?
     logfile     - pointer to the logfile
     convlogfile - pointer to the convlogfile
     noproj      - don't perform any projections
     diagerrs    - the data->SS errors-squared are diagonal
     noweight    - don't use data-weights
  OUTPUT:
     updated model gaussians
     avgloglikedata - average log likelihood of the data
  REVISION HISTORY:
     2008-08-21 Written Bovy
     2010-03-01 Added noproj option - Bovy
     2010-04-01 Added noweight option - Bovy
*/
#ifdef _OPENMP
#include <omp.h>
#endif
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <proj_gauss_mixtures.h>

void proj_gauss_mixtures(struct datapoint * data, int N, 
			 struct gaussian * gaussians, int K,
			 bool * fixamp, bool * fixmean, bool * fixcovar, 
			 double * avgloglikedata, double tol,
			 long long int maxiter, bool likeonly, double w, 
			 int splitnmerge, bool keeplog, FILE *logfile, 
			 FILE *convlogfile, bool noproj, bool diagerrs,
			 bool noweight){
  //Allocate some memory
  struct gaussian * startgaussians;
  startgaussians = gaussians;
  int d = (gaussians->VV)->size1;//dim of mm
  //Only give copies of the fix* vectors to the EM algorithm
  bool *fixamp_tmp, *fixmean_tmp,*fixcovar_tmp;
  fixamp_tmp = (bool *) malloc(K * sizeof (bool) );
  fixmean_tmp = (bool *) malloc(K * sizeof (bool) );
  fixcovar_tmp = (bool *) malloc(K * sizeof (bool) );
  int kk;
  for (kk = 0; kk != K; ++kk){
    *(fixamp_tmp++) = *(fixamp++);
    *(fixmean_tmp++) = *(fixmean++);
    *(fixcovar_tmp++) = *(fixcovar++);
  }
  fixamp -= K;
  fixamp_tmp -= K;
  fixmean -= K;
  fixmean_tmp -= K;
  fixcovar -= K;
  fixcovar_tmp -= K;
  //allocate the newalpha, newmm and newVV matrices
#ifdef _OPENMP
  nthreads = omp_get_max_threads();
#else
  nthreads = 1;
#endif
  newgaussians = (struct gaussian *) malloc(K * nthreads * sizeof (struct gaussian) );
  startnewgaussians = newgaussians;
  int ll;
  for (kk=0; kk != K*nthreads; ++kk){
    newgaussians->alpha = 0.0;
    newgaussians->mm = gsl_vector_calloc (d);
    newgaussians->VV = gsl_matrix_calloc (d,d);
    ++newgaussians;
  }
  newgaussians= startnewgaussians;
  double oldavgloglikedata;
  //allocate the q_ij matrix
  qij = gsl_matrix_alloc(N,K);
  double * qstarij = (double *) malloc(N * sizeof (double) );
  I = gsl_matrix_alloc(d,d);
  gsl_matrix_set_identity(I);//Unit matrix
  gsl_matrix_scale(I,w);//scaled to w
  //Also take care of the bbij's and the BBij's
  bs = (struct modelbs *) malloc(nthreads * K * sizeof (struct modelbs) );
  for (kk = 0; kk != nthreads*K; ++kk){
    bs->bbij = gsl_vector_alloc (d);
    bs->BBij = gsl_matrix_alloc (d,d);
    ++bs;
  }
  bs -= nthreads*K;
  //splitnmerge
  int maxsnm = K*(K-1)*(K-2)/2;
  int * snmhierarchy = (int *) malloc(maxsnm*3* sizeof (int) );
  int j,k,l;
  struct gaussian * oldgaussians = (struct gaussian *) malloc(K * sizeof (struct gaussian) );
  gsl_matrix * oldqij = gsl_matrix_alloc(N,K);
  for (kk=0; kk != K; ++kk){
    oldgaussians->mm = gsl_vector_calloc (d);
    oldgaussians->VV = gsl_matrix_calloc (d,d);
    ++oldgaussians;
  }
  oldgaussians -= K;

  //create temporary file to hold convergence info
  FILE *tmpconvfile;
  char * tmpfilename= tmpnam(NULL);;
  if (keeplog)
    tmpconvfile= fopen(tmpfilename,"w+");


  //proj_EM
  if (keeplog)
    fprintf(logfile,"#Initial proj_EM\n");
  //printf("Where's the segmentation fault?\n");
  proj_EM(data,N,gaussians,K,fixamp_tmp,fixmean_tmp,fixcovar_tmp,
	  avgloglikedata,tol,maxiter,likeonly,w,
	  keeplog,logfile,convlogfile,noproj,diagerrs,noweight);
  if (keeplog){
    fprintf(logfile,"\n");
    fprintf(convlogfile,"\n");
  }
  //reset fix* vectors
  for (kk = 0; kk != K; ++kk){
    *(fixamp_tmp++) = *(fixamp++);
    *(fixmean_tmp++) = *(fixmean++);
    *(fixcovar_tmp++) = *(fixcovar++);
  }
  fixamp -= K;
  fixamp_tmp -= K;
  fixmean -= K;
  fixmean_tmp -= K;
  fixcovar -= K;
  fixcovar_tmp -= K;

  //Run splitnmerge
  randgen = gsl_rng_alloc(gsl_rng_mt19937);
  bool weretrying = true;
  if (likeonly || splitnmerge == 0 || K < 3)
    ;
  else {
    while (weretrying){
      weretrying = false; /* this is set back to true if an improvement is found */
      //store avgloglike from normal EM and model parameters
      oldavgloglikedata = *avgloglikedata;
      gsl_matrix_memcpy(oldqij,qij);
      for (kk=0; kk != K; ++kk){
	oldgaussians->alpha = gaussians->alpha;
	gsl_vector_memcpy(oldgaussians->mm,gaussians->mm);
	gsl_matrix_memcpy(oldgaussians->VV,gaussians->VV);
	++oldgaussians;
	++gaussians;
      }
      gaussians -= K;
      oldgaussians -= K;
      //Then calculate the splitnmerge hierarchy
      calc_splitnmerge(data,N,gaussians,K,qij,snmhierarchy);
      //Then go through this hierarchy
      kk=0;
      while (kk != splitnmerge && kk != maxsnm){
	//splitnmerge
	j = *(snmhierarchy++);
	k = *(snmhierarchy++);
	l = *(snmhierarchy++);
	splitnmergegauss(gaussians,K,oldqij,j,k,l);
	//partial EM
	//Prepare fixed vectors for partial EM
	for (ll = 0; ll != K; ++ll){
	  *(fixamp_tmp++) = (ll != j && ll != k && ll != l);
	  *(fixmean_tmp++) = (ll != j && ll != k && ll != l);
	  *(fixcovar_tmp++) = (ll != j && ll != k && ll != l);
	}
	fixamp_tmp -= K;
	fixmean_tmp -= K;
	fixcovar_tmp -= K;
	if (keeplog)
	  fprintf(logfile,"#Merging %i and %i, splitting %i\n",j,k,l);
	proj_EM(data,N,gaussians,K,fixamp_tmp,fixmean_tmp,fixcovar_tmp,
		avgloglikedata,tol,maxiter,likeonly,w,keeplog,logfile,
		tmpconvfile,noproj,diagerrs,noweight);
	//reset fix* vectors
	for (ll = 0; ll != K; ++ll){
	    *(fixamp_tmp++) = *(fixamp++);
	    *(fixmean_tmp++) = *(fixmean++);
	    *(fixcovar_tmp++) = *(fixcovar++);
	}
	fixamp -= K;
	fixamp_tmp -= K;
	fixmean -= K;
	fixmean_tmp -= K;
	fixcovar -= K;
	fixcovar_tmp -= K;
	//Full EM
	if (keeplog){
	  fprintf(logfile,"#full EM:\n");
	  fprintf(tmpconvfile,"\n");
	}
	proj_EM(data,N,gaussians,K,fixamp_tmp,fixmean_tmp,fixcovar_tmp,
		avgloglikedata,tol,maxiter,likeonly,w,keeplog,logfile,
		tmpconvfile,noproj,diagerrs,noweight);
	if (keeplog){
	  fprintf(logfile,"\n");
	  fprintf(tmpconvfile,"\n");
	}
	//reset fix* vectors
	for (ll = 0; ll != K; ++ll){
	    *(fixamp_tmp++) = *(fixamp++);
	    *(fixmean_tmp++) = *(fixmean++);
	    *(fixcovar_tmp++) = *(fixcovar++);
	  }
	  fixamp -= K;
	  fixamp_tmp -= K;
	  fixmean -= K;
	  fixmean_tmp -= K;
	  fixcovar -= K;
	  fixcovar_tmp -= K;
	//Better?
	if (*avgloglikedata > oldavgloglikedata){
	  if (keeplog){
	    fprintf(logfile,"#accepted\n");
	    //Copy tmpfile into convfile
	    fseek(tmpconvfile,0,SEEK_SET);
	    while (feof(tmpconvfile) == 0)
	      fputc(fgetc(tmpconvfile),convlogfile);
	    fseek(convlogfile,-1,SEEK_CUR);
	    fclose(tmpconvfile);
	    tmpconvfile= fopen(tmpfilename,"w+");
	  }
	  weretrying = true;
	  ++kk;
	  break;
	}
	else {
	  if (keeplog){
	    fprintf(logfile,"#didn't improve likelihood\n");
	    fclose(tmpconvfile);
	    tmpconvfile= fopen(tmpfilename,"w+");
	  }
	  //revert back to the older solution
	  *avgloglikedata = oldavgloglikedata;
	  /*  gsl_matrix_memcpy(qij,oldqij);*/
	  for (ll=0; ll != K; ++ll){
	    gaussians->alpha = oldgaussians->alpha;
	    gsl_vector_memcpy(gaussians->mm,oldgaussians->mm);
	    gsl_matrix_memcpy(gaussians->VV,oldgaussians->VV);
	    ++oldgaussians;
	    ++gaussians;
	  }
	  gaussians -= K;
	  oldgaussians -= K;
	}
	++kk;
      }
      snmhierarchy -= 3*kk;
    }
  }

  if (keeplog){
    fclose(tmpconvfile);
    fprintf(convlogfile,"\n");
  }


  //Compute some criteria to set the number of Gaussians and print these to the logfile
  int ii;
  int npc,np;
  double pc,aic,mdl;
  if (keeplog){
    //Partition coefficient
    pc=0.;
    for (ii=0; ii != N; ++ii)
      for (kk=0; kk != K; ++kk)
	pc += pow(exp(gsl_matrix_get(qij,ii,kk)),2);
    pc /= N;
    fprintf(logfile,"Partition coefficient \t=\t%f\n",pc);
    //Akaike's information criterion
    npc = 1 + d + d * (d - 1) / 2;
    np = K * npc;
    aic= -2 * ((double) N - 1. - npc - 100)* *avgloglikedata + 3 * np;
    fprintf(logfile,"AIC \t\t\t=\t%f\n",aic);
    //MDL
    mdl = - *avgloglikedata * N + 0.5 * np * log(N);
    fprintf(logfile,"MDL \t\t\t=\t%f\n",mdl);
  }



  
  //Free memory
  gsl_matrix_free(I);
  gsl_matrix_free(qij);
  for (kk = 0; kk != nthreads*K; ++kk){
    gsl_vector_free(bs->bbij);
    gsl_matrix_free(bs->BBij);
    ++bs;
  }
  bs -= nthreads*K;
  free(bs);
  for (kk=0; kk != K*nthreads; ++kk){
    gsl_vector_free(newgaussians->mm);
    gsl_matrix_free(newgaussians->VV);
    ++newgaussians;
  }
  newgaussians= startnewgaussians;
  free(newgaussians);
  gsl_rng_free(randgen);

  gsl_matrix_free(oldqij);
  for (kk=0; kk != K; ++kk){
    gsl_vector_free(oldgaussians->mm);
    gsl_matrix_free(oldgaussians->VV);
    ++oldgaussians;
  }
  oldgaussians -= K;
  free(oldgaussians);
  free(qstarij);
  free(snmhierarchy);
  free(fixamp_tmp);
  free(fixmean_tmp);
  free(fixcovar_tmp);

  return;
}

