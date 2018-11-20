/*
  NAME:
     logsum
  PURPOSE:
     calculate the log of the sum of numbers which are given as logs using as much of the dynamic range as possible
  CALLING SEQUENCE:
     logsum(gsl_matrix * q, int row, bool isrow,bool partial, int partial_indx[3])
  INPUT:
     q            - matrix of which we will sum a row/column
     row          - row/column number
     isrow        - are we summing a row or a column
  OUTPUT:
     log of the sum
  REVISION HISTORY:
     2008-09-21 - Written Bovy
*/
#include <stdbool.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_matrix.h>
#include <proj_gauss_mixtures.h>

double logsum(gsl_matrix * q, int row, bool isrow){
  double logxmin = log(DBL_MIN);
  double logxmax = log(DBL_MAX);
  int l = (isrow) ? q->size2 : q->size1;

  /* First find the maximum and mininum */
  double max, min;
  minmax(q,row,isrow,&min,&max);//,allfixed);


  min *= -1.;
  min += logxmin;
  max *= -1.;
  max += logxmax - log(l);
  (min >  max) ? (max=max) : (max=min);
  double loglike=0.0;
  int dd;
  if (isrow)
    for (dd = 0; dd != q->size2; ++dd)
	loglike += exp(gsl_matrix_get(q,row,dd)+max);
  else
    for (dd = 0; dd != q->size1; ++dd)
	loglike += exp(gsl_matrix_get(q,dd,row)+max);
  

  return log(loglike)-max;
}
