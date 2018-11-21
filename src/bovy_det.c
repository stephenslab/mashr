/*
  NAME:
     bovy_det
  PURPOSE:
     returns the determinant of a matrix
  CALLING SEQUENCE:
     bovy_det(gsl_matrix * A)
  INPUT:
     A - matrix
  OUTPUT:
     det(A)
  REVISION HISTORY:
     2008-09-21 - Written Bovy
 */
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

double bovy_det(gsl_matrix * A){
  gsl_permutation * p = gsl_permutation_alloc (A->size1);
  int signum;
  double det;
  gsl_linalg_LU_decomp(A,p,&signum);
  det= gsl_linalg_LU_det(A,signum);
  gsl_permutation_free(p);

  return det;
}
