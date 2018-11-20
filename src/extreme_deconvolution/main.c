/*
  NAME:
     main
  PURPOSE:
     runs the proj_gauss_mixtures using input from files
  CALLING SEQUENCE:
     ./proj_gauss_mixtures
  INPUT:
  OUTPUT:
  REVISION HISTORY:
     2008-09-21 - Written Bovy
*/
#include <stdio.h>
#include <stdbool.h>
#include <proj_gauss_main.h>
#include <proj_gauss_mixtures.h>

int main (void){

  //Initialize some of the arrays (will be reallocated once more info is available)


  //Loop over commandline arguments to extract the different filenames


  //Read initial conditions file
  printf ("Reading initial conditions file...\n");
  if ( read_IC("IC.dat") )
    printf ("Successfully read initial conditions file\n");
  else {
    printf ("Reading initial conditions failed\n");
    printf("Returning");
    return -1;
  }

  //Successfully implemented reading the initial conditions file

  //Next up: reading the datafile
  //put the data in a structure that has [d_i,ww_i,SS_i,PP_i]
  //[d_i is a byte array which contains which values are missing


  //Read datafile  
  printf("\nReading data...\n");
  if (read_data("inputdata.dat"))
    printf("Succesfully read the data\n");
  else {
    printf("Couldn't read the data...\n");
    printf("Returning\n");
    return -1;
  }


  //Try showing the data
  int ii;
  for (ii=0; ii != N-1; ++ii){
    printf("datapoint %i:\n",ii+1);
    printf("Mean\t\t\t=\t%f\t%f\n",gsl_vector_get(data->ww,0),gsl_vector_get(data->ww,1));
    fflush(stdout);
    printf("Variance\t\t=\t%f\t%f\t%f\n",gsl_matrix_get(data->SS,0,0),gsl_matrix_get(data->SS,1,1),gsl_matrix_get(data->SS,0,1));
    printf("Projection matrix\t=\t%f\t%f\t%f\t%f\n",gsl_matrix_get(data->RR,0,0),gsl_matrix_get(data->RR,0,1),gsl_matrix_get(data->RR,1,0),gsl_matrix_get(data->RR,1,1));
    ++data;
  }   
  for (ii=N-1; ii != N; ++ii){
    printf("datapoint %i:\n",ii+1);
    printf("Mean\t\t\t=\t%f\n",gsl_vector_get(data->ww,0));
    fflush(stdout);
    printf("Variance\t\t=\t%f\n",gsl_matrix_get(data->SS,0,0));
    printf("Projection matrix\t=\t%f\t%f\n",gsl_matrix_get(data->RR,0,0),gsl_matrix_get(data->RR,0,1));
    ++data;
  } 
  data -= N;
  //data = startdata;



  //Successfully implemented reading the data into a structure!!

  //Next up: ouputting the result to a file (should basically be the same as the code at the end of read_IC
  //Done!


  //Then decide on the function EM
  //Run EM
  double avgloglikedata;
  proj_gauss_mixtures(data, N, gaussians, K,fixampP, fixmeanP, fixcovarP, &avgloglikedata, tol, maxiter, likeonly, w,splitnmerge,false,NULL,NULL);





  //Write the result
  printf("Writing the results...\n");
  if (write_model("result.dat") )
    printf("Succesfully wrote the results...\n");
  else {
    printf("Couldn't write the results...\n");
    printf("Returning\n");
    return -1;
  }


  //clean up
  if (cleanup());
  else {
    printf("Couldn't quite clean everything up...\n");
    printf("Sorry!\n");
    return -1;
  }

  return 0;
}
