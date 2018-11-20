/* 
   NAME:
      bovy_randvec
   PURPOSE:
      returns a (uniform) random vector with a maximum length
   CALLING SEQUENCE:
      bovy_randvec(gsl_vector * eps, int d, double length)
   INPUT:
      d      - dimension of the vector
      length - maximum length of the random vector
   OUTPUT:
      eps    - random vector
   REVISION HISTORY:
      2008-09-21
*/
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <proj_gauss_mixtures.h>

void bovy_randvec(gsl_vector * eps, int d, double length){
  length /= sqrt((double)d);
  int dd;
  for (dd = 0; dd != d; ++dd)
    gsl_vector_set(eps,dd,(2.*gsl_rng_uniform(randgen)-1.)*length);

  return;
}
