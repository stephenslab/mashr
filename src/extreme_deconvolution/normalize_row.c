/*
  NAME:
     normalize_row
  PURPOSE:
     normalize a row (or column) of a matrix given as logs
  CALLING SEQUENCE:
     normalize_row(gsl_matrix * q, int row,bool isrow,bool noweight,
     double weight)
  INPUT:
     q            - matrix
     row          - row to be normalized
     isrow        - is it a row or a column
     noweight    - add a weight to all of the values?
     weight       - weight to be added to all of the values
  OUTPUT:
     normalization factor (i.e. logsum)
  REVISION HISTORY:
     2008-09-21 - Written Bovy
     2010-04-01 - Added noweight and weight inputs to allow the qij to have 
                  weights - Bovy
*/
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <proj_gauss_mixtures.h>

double normalize_row(gsl_matrix * q, int row, bool isrow,
		     bool noweight, double weight){
  double loglike;
  if (isrow)
    loglike = logsum(q,row,true);
  else
    loglike = logsum(q,row,false);

  int dd;
  if (isrow)
    for (dd = 0; dd != q->size2; ++dd) {
      if ( noweight ) 
	gsl_matrix_set(q,row,dd,gsl_matrix_get(q,row,dd)-loglike);
      else
	gsl_matrix_set(q,row,dd,gsl_matrix_get(q,row,dd)-loglike+weight);
    }
  else
    for (dd = 0; dd != q->size1; ++dd) {
      if ( noweight ) 
	gsl_matrix_set(q,dd,row,gsl_matrix_get(q,dd,row)-loglike);
      else
	gsl_matrix_set(q,dd,row,gsl_matrix_get(q,dd,row)-loglike+weight);
    }
  if ( ! noweight ) loglike*= exp(weight);

  return loglike;
}
