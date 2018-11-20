/*
  NAME:
     write_model
  PURPOSE:
     writes the final model parameters to a file
  CALLING SEQUENCE:
     write_model(char outputfilename[])
  INPUT:
     outputfilename
  OUTPUT:
     -
  REVISION HISTORY:
     2008-09-21 - Written Bovy
*/
#include <stdbool.h>
#include <stdio.h>
#include <proj_gauss_mixtures.h>
#include <proj_gauss_main.h>

bool write_model(char outputfilename[]){

  FILE *outputfile;
  if ( (outputfile= fopen(outputfilename,"w+")) == NULL){
    printf ("Opening the data file failed...\n");
    return false;
  }  
  
  //loop over the gaussians and write the results to a file
  int kk,dd,dd1,dd2;
  for (kk=0; kk != K; ++kk){
    printf("#K=%i\n",kk+1);
    printf("%g\n",gaussians->alpha);
    for (dd=0; dd != d; ++dd)
      printf("%g\n",gsl_vector_get(gaussians->mm,dd));
    for (dd1=0; dd1 != d; ++dd1)
      printf("%g\n",gsl_matrix_get(gaussians->VV,dd1,dd1));
    for (dd1=0; dd1 != d-1; ++dd1)
      for (dd2=dd1+1; dd2 != d; ++dd2){
	printf("%g\n",gsl_matrix_get(gaussians->VV,dd1,dd2));
      }
    printf("\n");
    ++gaussians;
  }
  
  gaussians -= K;

  return true;
}
