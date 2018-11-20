/*
  NAME:
     calc_splitnmerge
  PURPOSE:
     calculates the split and merge hierarchy after an proj_EM convergence
  CALLING SEQUENCE:
     calc_splitnmerge(struct datapoint * data,int N,
     struct gaussian * gaussians, int K, gsl_matrix * qij, 
     int * snmhierarchy){
  INPUT:
     data      - the data
     N         - number of data points
     gaussians - model gaussians
     K         - number of gaussians
     qij       - matrix of log(posterior likelihoods)
     bs        - bs from previous EM estimate
  OUTPUT:
     snmhierarchy - the hierarchy, first row has the highest prioriry, 
                    goes down from there
  REVISION HISTORY:
     2008-09-21 - Written Bovy
*/
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <proj_gauss_mixtures.h>

void calc_splitnmerge(struct datapoint * data,int N,
		      struct gaussian * gaussians, int K, 
		      gsl_matrix * qij, int * snmhierarchy){
  gsl_matrix * tempqij = gsl_matrix_alloc(N,K);
  gsl_matrix_memcpy(tempqij,qij);
  gsl_matrix * Jmerge = gsl_matrix_alloc(K,K);
  gsl_matrix_set_all(Jmerge,-1.);
  int kk1, kk2, kk, ii,maxsnm= K*(K-1)*(K-2)/2;
  int d = (gaussians->VV)->size1;//dim of mm
  double temp1,temp2,temp;
  for (kk1 = 0; kk1 != K; ++kk1)
    for (kk2 = kk1+1; kk2 != K; ++kk2){
      temp = 0.;
      for (ii=0; ii != N; ++ii){
	//make them all exps
	temp1 = exp(gsl_matrix_get(qij,ii,kk1));
	temp2 = exp(gsl_matrix_get(qij,ii,kk2));
	temp += temp1*temp2;
      }
      gsl_matrix_set(Jmerge,kk1,kk2,temp);
    }

  //Then calculate Jsplit
  gsl_vector * Jsplit = gsl_vector_alloc(K);
  gsl_vector * Jsplit_temp = gsl_vector_alloc(K);
  gsl_vector_set_all(Jsplit,-1.);
  //if there is missing data, fill in the missing data
  struct missingdatapoint{
    gsl_vector *ww;
  };
  struct missingdatapoint * missingdata;
  missingdata = (struct missingdatapoint *) malloc(N * sizeof (struct missingdatapoint) );
  for (ii=0; ii != N; ++ii){
    missingdata->ww = gsl_vector_alloc(d);
    ++missingdata;
  }
  missingdata -= N;

  gsl_matrix * tempRR,* tempVV;
  gsl_vector * tempSS,* tempwork;
  gsl_vector * expectedww = gsl_vector_alloc(d);
  //gsl_vector_view tempUcol;
  double lambda;
  int di,signum;
  for (ii=0; ii != N; ++ii){
    //First check whether there is any missing data
    if ((data->ww)->size == d){
      gsl_vector_memcpy(missingdata->ww,data->ww);
      ++missingdata;
      ++data;
      continue;
    }

    /*AS IT STANDS THE MISSING DATA PART IS *NOT* IMPLEMENTED CORRECTLY: 
      A CORRECT IMPLEMENTATION NEEDS THE NULL SPACE OF THE PROJECTION 
      MATRIX WHICH CAN BE FOUND FROM THE FULL SINGULAR VALUE DECOMPOSITIIN, 
      UNFORTUNATELY GSL DOES NOT COMPUTE THE FULL SVD, BUT ONLY THE THIN SVD. 
      LAPACK MIGHT DO, BUT MIGHT NOT BE INSTALLED (?) AND THIS MIGHT BE HARD TO IMPLEMENT.

      INDICATED BELOW ARE THE SECTION THAT WOULD HAVE TO BE FIXED TO MAKE THIS WORK
     */

    //calculate expectation, for this we need to calculate the bbijs (EXACTLY THE SAME AS IN PROJ_EM, SHOULD WRITE GENERAL FUNCTION TO DO THIS)
    gsl_vector_set_zero(expectedww);
    for (kk = 0; kk != K; ++kk){
      //prepare...
      di = (data->SS)->size1;
      p = gsl_permutation_alloc (di);
      wminusRm = gsl_vector_alloc (di);
      gsl_vector_memcpy(wminusRm,data->ww);
      TinvwminusRm = gsl_vector_alloc (di);
      Tij = gsl_matrix_alloc(di,di);
      gsl_matrix_memcpy(Tij,data->SS);
      Tij_inv = gsl_matrix_alloc(di,di);
      VRT = gsl_matrix_alloc(d,di);
      VRTTinv = gsl_matrix_alloc(d,di);
      Rtrans = gsl_matrix_alloc(d,di);
      //Calculate Tij
      gsl_matrix_transpose_memcpy(Rtrans,data->RR);
      gsl_blas_dsymm(CblasLeft,CblasUpper,1.0,gaussians->VV,Rtrans,0.0,VRT);//Only the upper right part of VV is calculated
      gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,data->RR,VRT,1.0,Tij);//This is Tij
      //Calculate LU decomp of Tij and Tij inverse
      gsl_linalg_LU_decomp(Tij,p,&signum);
      gsl_linalg_LU_invert(Tij,p,Tij_inv);
      //Calculate Tijinv*(w-Rm)
      gsl_blas_dgemv(CblasNoTrans,-1.0,data->RR,gaussians->mm,1.0,wminusRm);
      gsl_blas_dsymv(CblasUpper,1.0,Tij_inv,wminusRm,0.0,TinvwminusRm);
      //Now calculate bij and Bij
      gsl_vector_memcpy(bs->bbij,gaussians->mm);
      gsl_blas_dgemv(CblasNoTrans,1.0,VRT,TinvwminusRm,1.0,bs->bbij);
      //..and add the result to expectedww
      gsl_vector_scale(bs->bbij,exp(gsl_matrix_get(qij,ii,kk)));
      gsl_vector_add(expectedww,bs->bbij);
      //Clean up
      gsl_permutation_free (p);
      gsl_vector_free(wminusRm);
      gsl_vector_free(TinvwminusRm);
      gsl_matrix_free(Tij);
      gsl_matrix_free(Tij_inv);
      gsl_matrix_free(VRT);
      gsl_matrix_free(VRTTinv);
      gsl_matrix_free(Rtrans);
      ++gaussians;
    }
    gaussians -= K;
    //if missing, fill in the missing data
    tempRR = gsl_matrix_alloc((data->RR)->size2,(data->RR)->size1);//will hold the transpose of RR
    //tempVV = gsl_matrix_alloc((data->RR)->size1,(data->RR)->size1);
    //tempSS = gsl_vector_alloc((data->RR)->size1);
    //tempwork = gsl_vector_alloc((data->RR)->size1);
    gsl_matrix_transpose_memcpy(tempRR,data->RR);
    gsl_blas_dgemv(CblasNoTrans,1.,tempRR,data->ww,0.,missingdata->ww);
    //gsl_linalg_SV_decomp(tempRR,tempVV,tempSS,tempwork);
    //compute the missing data THIS PART IS NOT IMPLEMENTED CORRECTLY
    //for (kk = 0; kk != d-(data->ww)->size; ++kk){
      //tempUcol = gsl_matrix_column(tempRR,d-1-kk);
    //gsl_blas_ddot(&(tempUcol.vector),expectedww,&lambda);
    //gsl_vector_scale(&(tempUcol.vector),lambda);
    //gsl_vector_add(missingdata->ww,&(tempUcol.vector));
    //}
    ++missingdata;
    ++data;
    //free
    gsl_matrix_free(tempRR);
    //gsl_matrix_free(tempVV);
    //gsl_vector_free(tempSS);
    //gsl_vector_free(tempwork);
  }
  data -= N;
  missingdata -= N;

  //then for every gaussian, calculate the KL divergence between the local data density and the l-th gaussian
  double tempsplit;
  p = gsl_permutation_alloc (d);
  tempVV= gsl_matrix_alloc(d,d);
  tempRR= gsl_matrix_alloc(d,d);
  tempSS= gsl_vector_alloc(d);
  tempwork= gsl_vector_alloc(d);
  for (kk = 0; kk != K; ++kk){
    //calculate qil/ql factors
    normalize_row(tempqij,kk,false,true,0.);
    //calculate inverse of V and det(V)
    gsl_matrix_memcpy(tempVV,gaussians->VV);
    gsl_linalg_LU_decomp(tempVV,p,&signum);
    gsl_linalg_LU_invert(tempVV,p,tempRR);//tempRR now has the inverse of VV
    tempsplit = d * halflogtwopi + 0.5 * gsl_linalg_LU_lndet(tempVV);
    for (ii=0; ii != N; ++ii){
      if (exp(gsl_matrix_get(tempqij,ii,kk)) == 0.){
	++missingdata;
	continue;
      }
      tempsplit += gsl_matrix_get(tempqij,ii,kk) * exp(gsl_matrix_get(tempqij,ii,kk));
      gsl_vector_memcpy(tempSS,gaussians->mm);
      gsl_vector_scale(tempSS,-1.);
      gsl_vector_add(tempSS,missingdata->ww);
      gsl_blas_dgemv(CblasNoTrans,1.0,tempRR,tempSS,0.,tempwork);
      gsl_blas_ddot(tempSS,tempwork,&lambda);
      tempsplit += 0.5 * exp(gsl_matrix_get(tempqij,ii,kk)) * lambda;
      ++missingdata;
    }
    gsl_vector_set(Jsplit,kk,tempsplit);
    //printf("Jsplit for gaussian %i = %f\n",kk,tempsplit);
    missingdata -= N;
    ++gaussians;
  }
  gaussians -= K;

  //free
  gsl_permutation_free(p);
  gsl_matrix_free(tempRR);
  gsl_matrix_free(tempVV);
  gsl_vector_free(tempSS);
  gsl_vector_free(tempwork);
  

  //and put everything in the hierarchy
  size_t maxj, maxk, maxl;
  for (kk1 = 0; kk1 != maxsnm; kk1 += (K-2)){
    gsl_matrix_max_index(Jmerge,&maxj,&maxk);
    gsl_vector_memcpy(Jsplit_temp,Jsplit);
    gsl_vector_set(Jsplit_temp,maxj,-1.);
    gsl_vector_set(Jsplit_temp,maxk,-1.);
    for (kk2=0; kk2 != K-2; ++kk2){
      maxl = gsl_vector_max_index(Jsplit_temp);
      gsl_vector_set(Jsplit_temp,maxl,-1.);
      *(snmhierarchy++)= maxj;
      *(snmhierarchy++)= maxk;
      *(snmhierarchy++)= maxl;
      //printf("j = %i, k = %i, l = %i\n",(int)maxj,(int)maxk,(int)maxl);
    }
    //then set it to zero and find the next
    gsl_matrix_set(Jmerge,maxj,maxk,-1.);
  }
  snmhierarchy -= 3*maxsnm;
    


  //clean up
  gsl_matrix_free(Jmerge);
  gsl_vector_free(Jsplit);
  gsl_vector_free(Jsplit_temp);

  return ;
}
