/*
  NAME:
     minmax
  PURPOSE:
     returns the (finite) minimum and the maximum of a vector which is the row/column of a matrix
  CALLING SEQUENCE:
     minmax(gsl_matrix * q, int row, bool isrow, double * min, double * max, bool partial, int partial_indx[3])
  INPUT:
     q            - matrix which holds the vector
     row          - row/column which holds the vector
     isrow        - is it the row or is it the column
  OUTPUT:
     min   - finite minimum
     max   - finite maximum
  REVISION HISTORY:
     2008-09-21 - Written Bovy
*/
#include <gsl/gsl_matrix.h>
#include <float.h>
#include <stdbool.h>
#include <proj_gauss_mixtures.h>

void minmax(gsl_matrix * q, int row, bool isrow, double * min, 
	    double * max){
  *max = -DBL_MAX;
  *min = DBL_MAX;
  int dd;
  double temp;
  if (isrow) {
    for (dd = 0; dd != q->size2; ++dd){
	temp  = gsl_matrix_get(q,row,dd);
	if (temp > *max && bovy_isfin(temp))
	  *max = temp;
	if (temp < *min && bovy_isfin(temp))
	  *min = temp;
    }
  }
  else {
    for (dd = 0; dd != q->size1; ++dd){
	temp  = gsl_matrix_get(q,dd,row);
	if (temp > *max && bovy_isfin(temp))
	  *max = temp;
	if (temp < *min && bovy_isfin(temp))
	  *min = temp;
    }
  }



  return ;
}
