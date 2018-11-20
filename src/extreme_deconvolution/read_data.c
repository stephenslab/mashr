/*
  NAME:
      read_data
  PURPOSE:
      read the data from the data input file
  CALLING SEQUENCE:
      read_data(inputfilename)
  INPUT:
      inputfilename - name of the data file
  OUPUT:
      reads everything into the right structures
  REVISION HISTORY:
     2008-09-21 - Written Bovy
*/
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <proj_gauss_main.h>
#include <proj_gauss_mixtures.h>

bool read_data(char inputfilename[]){
  
  FILE *inputfile;
  if ( (inputfile= fopen(inputfilename,"r")) == NULL){
    printf ("Opening the data file failed...\n");
    return false;
  }
  
  //First count the number of datapoints
  N=0;
  while (feof(inputfile) ==0) if (getc(inputfile) == '\n') N++;
  
  printf("%i datapoints found in ",N);
  printf("%s", inputfilename);
  printf("\n");
  
  //Then read the actual data
  fseek(inputfile,0,SEEK_SET);

  //allocate data
  data = (struct datapoint *) malloc(N*sizeof (struct datapoint));
  if (data == NULL){
    printf("Allocation of arrays failed, not enough free memory?\n");
    exit(-1);
  } 
  startdata = data;
    
  //current line holds a maximum of d+dV+d*d elements
  //First read the whole line
  double *curr_line = (double *) calloc(d+dV+d*d, sizeof (double) );
  bool end_line;
  int ii=0,dd=0,di;//di is the dimension of the individual datapoints
  while (feof(inputfile) == 0){
    char curr_value[20]="";
    end_line = read_till_sep(curr_value,inputfile,'|');
    *(curr_line++) = atof(curr_value);
    ++dd;
    if (end_line){
      ++ii;
      curr_line -= dd;
      //Determine the dimension of the datapoint, di, from the equation di+di*(di+1)/2+d*di=dd, or, 
      di = (int) ((-(3 + 2 * d) + sqrt((3 + 2 * d)*(3 + 2 * d)+8 * dd))/2);
      //then write data values to memory
      //first allocate the space for the data
      data->ww = gsl_vector_alloc(di);
      data->SS = gsl_matrix_alloc(di,di);
      data->RR = gsl_matrix_alloc(di,d);
      int dd1,dd2;
      for (dd1=0; dd1 != di; ++dd1)
	gsl_vector_set(data->ww,dd1,*(curr_line++));
      for (dd1=0; dd1 != di; ++dd1)
	gsl_matrix_set(data->SS,dd1,dd1,*(curr_line++));
      for (dd1=0; dd1 != di-1; ++dd1)
	for (dd2=dd1+1; dd2 != di; ++dd2){
	  gsl_matrix_set(data->SS,dd1,dd2,*curr_line);
	  gsl_matrix_set(data->SS,dd2,dd1,*curr_line);
	  curr_line++;
	}
      for (dd1=0; dd1 != di; ++dd1)
	for (dd2=0; dd2 != d; ++dd2)
	  gsl_matrix_set(data->RR,dd1,dd2,*(curr_line++));
      curr_line -= dd;
      dd = 0;
      data++;
      if (ii == N) break;
    }
  }
  
  data = startdata;

  fclose(inputfile);
  free(curr_line);

  return true;
}
