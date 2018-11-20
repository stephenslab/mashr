/*
  NAME:
     splitnmergegauss
  PURPOSE:
     split one gaussian and merge two other gaussians
  CALLING SEQUENCE:
     splitnmergegauss(struct gaussian * gaussians,int K, gsl_matrix * qij, 
     int j, int k, int l)
  INPUT:
     gaussians   - model gaussians
     K           - number of gaussians
     qij         - matrix of log(posterior likelihoods)
     j,k         - gaussians that need to be merged
     l           - gaussian that needs to be split
  OUTPUT:
     updated gaussians
  REVISION HISTORY:
     2008-09-21 - Written Bovy
*/
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <proj_gauss_mixtures.h>

void splitnmergegauss(struct gaussian * gaussians,int K, 
		      gsl_matrix * qij, int j, int k, int l){
  //get the gaussians to be split 'n' merged
  int d = (gaussians->VV)->size1;//dim of mm
  //int partial_indx[]= {-1,-1,-1};/* dummy argument for logsum */
  //bool * dummy_allfixed = (bool *) calloc(K,sizeof(bool));
  //j,k,l gaussians
  struct gaussian gaussianj, gaussiank, gaussianl;
  gaussianj.mm = gsl_vector_alloc(d);
  gaussianj.VV = gsl_matrix_alloc(d,d);
  gaussiank.mm = gsl_vector_alloc(d);
  gaussiank.VV = gsl_matrix_alloc(d,d);
  gaussianl.mm = gsl_vector_alloc(d);
  gaussianl.VV = gsl_matrix_alloc(d,d);
  
  gsl_matrix * unitm = gsl_matrix_alloc(d,d);
  gsl_matrix_set_identity(unitm);
  gsl_vector * eps = gsl_vector_alloc(d);
  double qjj,qjk,detVVjl;
  int kk;
  for (kk = 0; kk != K; ++kk){
    if (kk == j){
      gaussianj.alpha = gaussians->alpha;
      gsl_vector_memcpy(gaussianj.mm,gaussians->mm);
      gsl_matrix_memcpy(gaussianj.VV,gaussians->VV);
      qjj = exp(logsum(qij,j,false));//,dummy_allfixed));
    }
    if (kk == k){
      gaussiank.alpha = gaussians->alpha;
      gsl_vector_memcpy(gaussiank.mm,gaussians->mm);
      gsl_matrix_memcpy(gaussiank.VV,gaussians->VV);
      qjk = exp(logsum(qij,k,false));//,dummy_allfixed));
    }
    if (kk == l){
      gaussianl.alpha = gaussians->alpha;
      gsl_vector_memcpy(gaussianl.mm,gaussians->mm);
      gsl_matrix_memcpy(gaussianl.VV,gaussians->VV);
    }
    ++gaussians;
  }
  gaussians -= K;

  //merge j & k
  gaussianj.alpha += gaussiank.alpha;
  if (qjk == 0. && qjj == 0){
    gsl_vector_add(gaussianj.mm,gaussiank.mm);
    gsl_vector_scale(gaussianj.mm,0.5);
    gsl_matrix_add(gaussianj.VV,gaussiank.VV);
    gsl_matrix_scale(gaussianj.VV,0.5);
  }
  else{
    gsl_vector_scale(gaussianj.mm,qjj/(qjj+qjk));
    gsl_vector_scale(gaussiank.mm,qjk/(qjj+qjk));
    gsl_vector_add(gaussianj.mm,gaussiank.mm);
    gsl_matrix_scale(gaussianj.VV,qjj/(qjj+qjk));
    gsl_matrix_scale(gaussiank.VV,qjk/(qjj+qjk));
    gsl_matrix_add(gaussianj.VV,gaussiank.VV);
  }

  //split l
  gaussianl.alpha /= 2.;
  gaussiank.alpha = gaussianl.alpha;
  detVVjl = bovy_det(gaussianl.VV);
  detVVjl= pow(detVVjl,1./d);
  gsl_matrix_scale(unitm,detVVjl);
  gsl_matrix_memcpy(gaussiank.VV,unitm);
  gsl_matrix_memcpy(gaussianl.VV,unitm);
  gsl_vector_memcpy(gaussiank.mm,gaussianl.mm);
  bovy_randvec(eps,d,sqrt(detVVjl));
  gsl_vector_add(gaussiank.mm,eps);
  bovy_randvec(eps,d,sqrt(detVVjl));
  gsl_vector_add(gaussianl.mm,eps);
  
  //copy everything back into the right gaussians
  for (kk = 0; kk != K; ++kk){
    if (kk == j){
      gaussians->alpha = gaussianj.alpha;
      gsl_vector_memcpy(gaussians->mm,gaussianj.mm);
      gsl_matrix_memcpy(gaussians->VV,gaussianj.VV);
    }
    if (kk == k){
      gaussians->alpha = gaussiank.alpha;
      gsl_vector_memcpy(gaussians->mm,gaussiank.mm);
      gsl_matrix_memcpy(gaussians->VV,gaussiank.VV);
    }
    if (kk == l){
      gaussians->alpha = gaussianl.alpha;
      gsl_vector_memcpy(gaussians->mm,gaussianl.mm);
      gsl_matrix_memcpy(gaussians->VV,gaussianl.VV);
    }
    ++gaussians;
  }
  gaussians -= K;

  //cleanup
  gsl_matrix_free(unitm);
  gsl_vector_free(eps);
  //free(dummy_allfixed);

  return ;
}
