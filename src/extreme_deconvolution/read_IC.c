/*
  NAME:
     read_IC
  PURPOSE:
     read the initial conditions file
  CALLING SEQUENCE:
     read_IC(char ICfilename[])
  INPUT:
     ICfilename   - initial conditions filename
  OUTPUT:
     sets the options and the initial conditions
  REVISION HISTORY:
     2008-09-21 - Written Bovy
*/
#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <proj_gauss_mixtures.h>
#include <proj_gauss_main.h>

bool read_IC(char ICfilename[]){
  
  FILE *ICfile;
  if ( (ICfile= fopen(ICfilename,"r")) == NULL){
    printf ("Opening the initial conditions file failed...\n");
    return false;
  }
  char line[100];
  printf("Reading the options section of the initial conditions file...\n");
  while (fgets(line,100,ICfile) != NULL){
    if (line[0] != '#'){
      if (line[0] == '\n') break;
      if (parse_option(line) == false){
	printf("One of the lines in the options section of the initial conditions file is corrupted\n");
	printf("Please check the initial conditions file and try again\n");
	return false;
      }
    }
  }
  
  printf("Successfully read the options\n");

  printf("Reading the initial model parameters...\n");
  

  //Read first block, establish the dimension d of the modeled quantities
  int countd=0;
  double *mmtemp = (double *) malloc (1000000 * sizeof (double));
  while (fgets(line,100,ICfile) != NULL){
    if (line[0] != '#'){
      if (line[0] == '\n') break;
      *(mmtemp++) = atof(line);
      ++countd;
    }
  }
  //now determine d
  d = (int) (-3 + sqrt(9 + 8 * (countd-1)))/2 ; 
  dV = (int) (d*(d+1)/2);

  //allocate the alpha, mm and VV matrices
  int kk;
  for (kk=0; kk != K; ++kk){
    gaussians->mm = gsl_vector_alloc (d);
    gaussians->VV = gsl_matrix_alloc (d,d);;
    ++gaussians;
  }
  gaussians -= K;
  
  //first map the mmtemp values on the right alpha, mm, VV
  mmtemp -= countd;
  (*gaussians).alpha = *(mmtemp++);
  int dd;
  for (dd=0; dd != d; ++dd)
    gsl_vector_set(gaussians->mm,dd,*(mmtemp++));
  int dd1,dd2;
  for (dd1=0; dd1 != d; ++dd1)
    gsl_matrix_set(gaussians->VV,dd1,dd1,*(mmtemp++));
  for (dd1=0; dd1 != d-1; ++dd1)
    for (dd2=dd1+1; dd2 != d; ++dd2){
      gsl_matrix_set(gaussians->VV,dd1,dd2,*mmtemp);
      gsl_matrix_set(gaussians->VV,dd2,dd1,*mmtemp);
      mmtemp++;
    }

  ++gaussians;
  
  //reallocate mmtemp
  mmtemp -= countd;
  mmtemp = (double *) realloc (mmtemp,countd * sizeof (double) );
  if (mmtemp == NULL){
    printf("Error reallocating memory\n");
    printf("Returning\n");
    return false;
  }

  //Then read the rest of the Gaussians.
  for (kk=1; kk != K; ++kk){
    while (fgets(line,100,ICfile) != NULL){
      if (line[0] != '#'){
	if (line[0] == '\n') break;
	*(mmtemp++) = atof(line);
      }
    }
    mmtemp -=countd;
    (*gaussians).alpha = *(mmtemp++);
    for (dd=0; dd != d; ++dd)
      gsl_vector_set(gaussians->mm,dd,*(mmtemp++));
    for (dd1=0; dd1 != d; ++dd1)
      gsl_matrix_set(gaussians->VV,dd1,dd1,*(mmtemp++));
    for (dd1=0; dd1 != d-1; ++dd1)
      for (dd2=dd1+1; dd2 != d; ++dd2){
	gsl_matrix_set(gaussians->VV,dd1,dd2,*mmtemp);
	gsl_matrix_set(gaussians->VV,dd2,dd1,*mmtemp);
	mmtemp++;
      }
    ++gaussians;
    mmtemp -= countd;
  }
  gaussians -= K;
    
  free(mmtemp);

  fclose(ICfile);

  printf("Successfully read initial model parameters from the initial conditions file\n");


  //Print options
  printf("\nThe options are set to:\n");
  printf("K\t\t=\t");
  printf("%i",K);
  printf("\n");
  printf("maxiter\t\t=\t");
  printf("%lli",maxiter);
  printf("\n");
  printf("tol\t\t=\t");
  printf("%f",tol);
  printf("\n");
  printf("splitnmerge\t=\t");
  printf("%i",splitnmerge);
  printf("\n");
  printf("likeonly\t=\t");
  printf("%i",likeonly);
  printf("\n");
  printf("w\t\t=\t");
  printf("%f",w);
  printf("\n");
  printf("fixamp\t\t=\t");
  int ii;
  for (ii=0; ii != K; ++ii){
    printf("%i",*fixampP);
    if (ii < K-1) printf("\t");
    fixampP++;
  }
  fixampP -= K;
  printf("\n");
  printf("fixmean\t\t=\t");
  for (ii=0; ii != K; ++ii){
    printf("%i",*fixmeanP);
    if (ii < K-1) printf("\t");
    fixmeanP++;
  }
  fixmeanP -= K;
  printf("\n");
  printf("fixcovar\t=\t");
  for (ii=0; ii != K; ++ii){
    printf("%i",*fixcovarP);
    if (ii < K-1) printf("\t");
    fixcovarP++;
  }
  fixcovarP -= K;
  printf("\n");


  //Print the initial model parameters
  printf("\nInitial model parameters used:\n\n");
  for (kk=0; kk != K; ++kk){
    printf("Gaussian ");
    printf("%i",kk);
    printf("\n");
    printf("amp\t=\t");
    printf("%f",(*gaussians).alpha);
    printf("\n");
    printf("mean\t=\t");
    for (dd=0; dd != d; ++dd){
      printf("%f",gsl_vector_get(gaussians->mm,dd));
      if (dd < d-1) printf("\t");
    }
    printf("\n");
    printf("covar\t=\t");
    for (dd1=0; dd1 != d; ++dd1)
      printf("%f\t",gsl_matrix_get(gaussians->VV,dd1,dd1));
    for (dd1=0; dd1 != d-1; ++dd1)
      for (dd2=dd1+1; dd2 != d; ++dd2){
	printf("%f\t",gsl_matrix_get(gaussians->VV,dd1,dd2));
      }
    ++gaussians;
    printf("\n\n");
  }
  gaussians -= K;

  return true;
  
}
