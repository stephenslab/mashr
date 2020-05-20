// Wrapper to the projected gaussian mixture algorithm (Bovy 2009)
// Gao Wang (c) 2018-2020 wang.gao@columbia.edu
#include "RcppGSL.h"
#include <iostream>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "extreme_deconvolution.h"

using Rcpp::List;
using Rcpp::Named;

void
int2bool(RcppGSL::vector<int> & a, int K, bool* x)
{
	int kk;

	for (kk = 0; kk != K; ++kk) *(x++) = (bool) a[kk];
	x -= K;
}

// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::export]]
List
extreme_deconvolution_rcpp(
	RcppGSL::matrix<double> & ydata,
	RcppGSL::vector<double> & ycovar,
	RcppGSL::vector<double> & projection,
	RcppGSL::vector<double> & logweights,
	RcppGSL::vector<double> & amp,
	RcppGSL::matrix<double> & xmean,
	RcppGSL::vector<double> & xcovar,
	RcppGSL::vector<int> & fixamp_int,
	RcppGSL::vector<int> & fixmean_int,
	RcppGSL::vector<int> & fixcovar_int,
	double tol,
	int maxiter, int likeonly, double w,
	RcppGSL::vector<int> & logfilename,
	int splitnmerge,
	RcppGSL::vector<int> & convlogfilename,
	bool noproj, bool diagerrs, bool noweight)
{
	// convert variables from R interface
	int N = ydata.nrow(), dy = ydata.ncol(), d = xmean.ncol(), K = xmean.nrow(),
	    slen = logfilename.size(), convloglen = convlogfilename.size();

        bool fixamp_array[K];
	bool fixmean_array[K];
	bool fixcovar_array[K];
	bool * fixamp   = &fixamp_array[0];
	bool * fixmean  = &fixmean_array[0];
	bool * fixcovar = &fixcovar_array[0];
	int2bool(fixamp_int,K,fixamp);
	int2bool(fixmean_int,K,fixmean);
	int2bool(fixcovar_int,K,fixcovar);
	
	// Set up logfiles
	bool keeplog = true;
	char logname[slen + 1];
	char convlogname[convloglen + 1];
	int ss;

	if (likeonly != 0 || slen <= 1) {
		keeplog = false;
	} else {
		for (ss = 0; ss != slen; ++ss)
			logname[ss] = (char) logfilename[ss];
		for (ss = 0; ss != convloglen; ++ss)
			convlogname[ss] = (char) convlogfilename[ss];
		logname[slen] = '\0';
		convlogname[convloglen] = '\0';
	}

	if (keeplog) {
		logfile = fopen(logname, "a");
		if (logfile == NULL) return -1;

		convlogfile = fopen(convlogname, "w");
		if (convlogfile == NULL) return -1;
	}

	if (keeplog) {
		time_t now;
		time(&now);
		fprintf(logfile, "#----------------------------------\n");
		fprintf(logfile, "#\n#%s\n", asctime(localtime(&now)));
		fprintf(logfile, "#----------------------------------\n");
		fflush(logfile);
	}

	// Copy everything into the right formats
	struct datapoint * data      = (struct datapoint *) malloc(N * sizeof(struct datapoint) );
	struct gaussian  * gaussians = (struct gaussian *) malloc(K * sizeof(struct gaussian) );

	int ii, jj, dd1, dd2;
	for (ii = 0; ii != N; ++ii) {
		data->ww = gsl_vector_alloc(dy);
		if (!noweight) data->logweight = logweights[ii];
		if (diagerrs) data->SS = gsl_matrix_alloc(dy, 1);
		else data->SS = gsl_matrix_alloc(dy, dy);
		if (!noproj) data->RR = gsl_matrix_alloc(dy, d);
		for (dd1 = 0; dd1 != dy; ++dd1)
			gsl_vector_set(data->ww, dd1, ydata(ii, dd1));
		if (diagerrs) {
			for (dd1 = 0; dd1 != dy; ++dd1)
				gsl_matrix_set(data->SS, dd1, 0, ycovar[ii * dy + dd1]);
		} else {
			for (dd1 = 0; dd1 != dy; ++dd1)
				for (dd2 = 0; dd2 != dy; ++dd2)
					gsl_matrix_set(data->SS, dd1, dd2, ycovar[ii * dy * dy + dd1 * dy + dd2]);
		}
		if (!noproj) {
			for (dd1 = 0; dd1 != dy; ++dd1)
				for (dd2 = 0; dd2 != d; ++dd2)
					gsl_matrix_set(data->RR, dd1, dd2, projection[ii * dy * d + dd1 * dy + dd2]);
		} else { data->RR = NULL; }
		++data;
	}
	data -= N;

	for (jj = 0; jj != K; ++jj) {
		gaussians->mm    = gsl_vector_alloc(d);
		gaussians->VV    = gsl_matrix_alloc(d, d);
		gaussians->alpha = amp[jj];
		for (dd1 = 0; dd1 != d; ++dd1)
			gsl_vector_set(gaussians->mm, dd1, xmean(jj, dd1));
		for (dd1 = 0; dd1 != d; ++dd1)
			for (dd2 = 0; dd2 != d; ++dd2)
				gsl_matrix_set(gaussians->VV, dd1, dd2, xcovar[jj * d * d + dd1 * d + dd2]);
		++gaussians;
	}
	gaussians -= K;

	// Print the initial model parameters to the logfile
	int kk;
	if (keeplog) {
		fprintf(logfile, "#\n#Using %i Gaussians and w = %f\n\n", K, w);
		fprintf(logfile, "#\n#Initial model parameters used:\n\n");
		for (kk = 0; kk != K; ++kk) {
			fprintf(logfile, "#Gaussian ");
			fprintf(logfile, "%i", kk);
			fprintf(logfile, "\n");
			fprintf(logfile, "#amp\t=\t");
			fprintf(logfile, "%f", (*gaussians).alpha);
			fprintf(logfile, "\n");
			fprintf(logfile, "#mean\t=\t");
			for (dd1 = 0; dd1 != d; ++dd1) {
				fprintf(logfile, "%f", gsl_vector_get(gaussians->mm, dd1));
				if (dd1 < d - 1) fprintf(logfile, "\t");
			}
			fprintf(logfile, "\n");
			fprintf(logfile, "#covar\t=\t");
			for (dd1 = 0; dd1 != d; ++dd1)
				fprintf(logfile, "%f\t", gsl_matrix_get(gaussians->VV, dd1, dd1));
			for (dd1 = 0; dd1 != d - 1; ++dd1)
				for (dd2 = dd1 + 1; dd2 != d; ++dd2) {
					fprintf(logfile, "%f\t", gsl_matrix_get(gaussians->VV, dd1, dd2));
				}
			++gaussians;
			fprintf(logfile, "\n#\n");
		}
		gaussians -= K;
		fflush(logfile);
	}

	double avgloglikedata_np = 0.0;
	double * avgloglikedata;
	avgloglikedata = &avgloglikedata_np;

	// Then run projected_gauss_mixtures
	proj_gauss_mixtures(data, N, gaussians, K, fixamp, fixmean, fixcovar,
	                    avgloglikedata, tol, (long long int) maxiter, (bool) likeonly, w,
	                    splitnmerge, keeplog, logfile, convlogfile, noproj, diagerrs, noweight);

	// Print the final model parameters to the logfile
	if (keeplog) {
		fprintf(logfile, "\n#Final model parameters obtained:\n\n");
		for (kk = 0; kk != K; ++kk) {
			fprintf(logfile, "#Gaussian ");
			fprintf(logfile, "%i", kk);
			fprintf(logfile, "\n");
			fprintf(logfile, "#amp\t=\t");
			fprintf(logfile, "%f", (*gaussians).alpha);
			fprintf(logfile, "\n");
			fprintf(logfile, "#mean\t=\t");
			for (dd1 = 0; dd1 != d; ++dd1) {
				fprintf(logfile, "%f", gsl_vector_get(gaussians->mm, dd1));
				if (dd1 < d - 1) fprintf(logfile, "\t");
			}
			fprintf(logfile, "\n");
			fprintf(logfile, "#covar\t=\t");
			for (dd1 = 0; dd1 != d; ++dd1)
				fprintf(logfile, "%f\t", gsl_matrix_get(gaussians->VV, dd1, dd1));
			for (dd1 = 0; dd1 != d - 1; ++dd1)
				for (dd2 = dd1 + 1; dd2 != d; ++dd2) {
					fprintf(logfile, "%f\t", gsl_matrix_get(gaussians->VV, dd1, dd2));
				}
			++gaussians;
			fprintf(logfile, "\n#\n");
		}
		gaussians -= K;
		fflush(logfile);
	}

	// Then update the arrays given to us by IDL
	for (jj = 0; jj != K; ++jj) {
		amp[jj] = gaussians->alpha;
		for (dd1 = 0; dd1 != d; ++dd1)
			xmean(jj, dd1) = gsl_vector_get(gaussians->mm, dd1);
		for (dd1 = 0; dd1 != d; ++dd1)
			for (dd2 = 0; dd2 != d; ++dd2)
				xcovar[jj * d * d + dd1 * d + dd2] = gsl_matrix_get(gaussians->VV, dd1, dd2);
		++gaussians;
	}
	gaussians -= K;

	// And free any memory we allocated
	for (ii = 0; ii != N; ++ii) {
		gsl_vector_free(data->ww);
		gsl_matrix_free(data->SS);
		if (!noproj) gsl_matrix_free(data->RR);
		++data;
	}
	data -= N;
	free(data);

	for (jj = 0; jj != K; ++jj) {
		gsl_vector_free(gaussians->mm);
		gsl_matrix_free(gaussians->VV);
		++gaussians;
	}
	gaussians -= K;
	free(gaussians);

	if (keeplog) {
		fclose(logfile);
		fclose(convlogfile);
	}

	return List::create(Named("xmean") = xmean,
	                    Named("xcovar")         = xcovar,
	                    Named("xamp")           = amp,
	                    Named("avgloglikedata") = avgloglikedata_np);
} // extreme_deconvolution_rcpp
