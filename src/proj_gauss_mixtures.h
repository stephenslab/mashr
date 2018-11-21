#ifndef __PROJ_GAUSS_MIXTURES_H__
#define __PROJ_GAUSS_MIXTURES_H__


/* include */
#include <stdbool.h>
#include <float.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_linalg.h>

/* variable declarations and definitions */

double halflogtwopi; /* constant used in calculation */

int nthreads;
int tid;

struct gaussian{
  double alpha;
  gsl_vector *mm;
  gsl_matrix *VV;
};

struct datapoint{
  gsl_vector *ww;
  gsl_matrix *SS;
  gsl_matrix *RR;
  double logweight;
};

struct gaussian * newgaussians, * startnewgaussians;
gsl_matrix * qij;
gsl_permutation * p;
gsl_vector * wminusRm, * TinvwminusRm;
gsl_matrix * Tij,* Tij_inv, * VRT, * VRTTinv, * Rtrans, * I;
struct modelbs{
  gsl_vector *bbij;
  gsl_matrix *BBij;
};

struct modelbs * bs;

FILE *logfile,*convlogfile;

gsl_rng * randgen; /* global random number generator */




/* function declarations */

static inline bool bovy_isfin(double x){ /* true if x is finite */
  return (x > DBL_MAX || x < -DBL_MAX) ? false : true;}
void minmax(gsl_matrix * q, int row, bool isrow, double * min, double * max);
double logsum(gsl_matrix * q, int row, bool isrow);
double normalize_row(gsl_matrix * q, int row,bool isrow,bool noweight, double weight);
void bovy_randvec(gsl_vector * eps, int d, double length); /* returns random vector */
double bovy_det(gsl_matrix * A);/* determinant of matrix A */
void calc_splitnmerge(struct datapoint * data,int N,struct gaussian * gaussians, int K, gsl_matrix * qij, int * snmhierarchy);
void splitnmergegauss(struct gaussian * gaussians,int K, gsl_matrix * qij, int j, int k, int l);
void proj_EM_step(struct datapoint * data, int N, struct gaussian * gaussians, int K,bool * fixamp, bool * fixmean, bool * fixcovar, double * avgloglikedata, bool likeonly, double w,bool noproj, bool diagerrs, bool noweight);
void proj_EM(struct datapoint * data, int N, struct gaussian * gaussians, int K,bool * fixamp, bool * fixmean, bool * fixcovar, double * avgloglikedata, double tol,long long int maxiter, bool likeonly, double w,bool keeplog, FILE *logfile,FILE *tmplogfile, bool noproj, bool diagerrs, bool noweight);
void proj_gauss_mixtures(struct datapoint * data, int N, struct gaussian * gaussians, int K,bool * fixamp, bool * fixmean, bool * fixcovar, double * avgloglikedata, double tol,long long int maxiter, bool likeonly, double w, int splitnmerge, bool keeplog, FILE *logfile,FILE *convlogfile, bool noproj, bool diagerrs, bool noweight);
void calc_qstarij(double * qstarij, gsl_matrix * qij, int partial_indx[3]);

#endif /* proj_gauss_mixtures.h */
