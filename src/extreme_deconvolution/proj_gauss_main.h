#ifndef __PROJ_GAUSS_MAIN_H__
#define __PROJ_GAUSS_MAIN_H__

#include <stdbool.h>
#include <stdio.h>
#include <proj_gauss_mixtures.h>



/* variable declarations */

//Options:
int K;
bool *fixampP, *fixmeanP, *fixcovarP;
double avgloglikedata;
double tol;
long long int maxiter;
bool likeonly;
double w;
int splitnmerge;

//Model parameters:
int d;
int dV;//dim of VV, i.e. d(d+1)/2
struct gaussian * gaussians;

//Data parameters:
int N;//number of datapoints
struct datapoint *data, *startdata;



/* function declarations */
bool parse_option(char line[]);
bool read_IC(char line[]);
bool read_data(char line[]);
bool read_till_sep(char curr_value[],FILE *file,char sep);
bool write_model(char outputfilename[]);
bool cleanup();

#endif  /* proj_gauss_main.h */
