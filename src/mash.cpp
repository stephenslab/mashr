// Wrapper to various C++ functions/objects for inference in MASH
// Gao Wang (c) 2017-2020 wang.gao@columbia.edu
#include <iostream>
#include <stdexcept>
#ifdef _OPENMP
# include <omp.h>
#endif
#include "RcppArmadillo.h"
#include "mash.h"

using Rcpp::List;
using Rcpp::Named;
using Rcpp::NumericVector;
using Rcpp::IntegerVector;

using arma::vectorise;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List
inv_chol_tri_rcpp(const arma::mat & x_mat)
{
	mat res = trans(inv(trimatu(chol(x_mat))));

	return List::create(Named("data") = res,
	                          Named("status") = 0);
}

// [[Rcpp::export]]
List
calc_lik_rcpp(const arma::mat &   b_mat,
              const arma::mat &   s_mat,
              const arma::mat &   v_mat,
              const arma::mat &   l_mat,
              NumericVector U_3d,
              NumericVector sigma_3d,
              bool                logd,
              bool                common_cov,
              int                 n_thread = 1)
{
	// hide armadillo warning / error messages
	mat res;
	if (!Rf_isNull(U_3d.attr("dim"))) {
		// matrix version
		// set cube data from R 3D array
		IntegerVector dimU = U_3d.attr("dim");
		cube U_cube(U_3d.begin(), dimU[0], dimU[1], dimU[2]);
		cube sigma_cube;
		if (!Rf_isNull(sigma_3d.attr("dim"))) {
			IntegerVector dimSigma = sigma_3d.attr("dim");
			cube tmp_cube(sigma_3d.begin(), dimSigma[0], dimSigma[1], dimSigma[2]);
			sigma_cube = tmp_cube;
		}
		res = calc_lik(b_mat, s_mat, v_mat, l_mat, U_cube, sigma_cube, logd, common_cov, n_thread);
	} else {
		// vector version
		res = calc_lik(vectorise(b_mat), vectorise(s_mat), v_mat(0, 0), Rcpp::as<arma::vec>(U_3d), logd);
	}
	return List::create(Named("data") = res,
	                          Named("status") = 0);
} // calc_lik_rcpp

// [[Rcpp::export]]
List
calc_lik_precomputed_rcpp(const arma::mat &   b_mat,
                          NumericVector rooti_3d,
                          bool                logd,
                          bool                common_cov,
                          int                 n_thread = 1)
{
	// hide armadillo warning / error messages
	mat res;
	// set cube data from R 3D array
	IntegerVector dimR = rooti_3d.attr("dim");
	cube rooti_cube(rooti_3d.begin(), dimR[0], dimR[1], dimR[2]);

	res = calc_lik(b_mat, rooti_cube, logd, common_cov, n_thread);
	return List::create(Named("data") = res,
	                          Named("status") = 0);
}

// [[Rcpp::export]]
List
calc_post_rcpp(const arma::mat &   b_mat,
               const arma::mat &   s_mat,
               const arma::mat &   s_alpha_mat,
               const arma::mat &   s_orig_mat,
               const arma::mat &   v_mat,
               const arma::mat &   l_mat,
               const arma::mat &   a_mat,
               NumericVector U_3d,
               const arma::mat &   posterior_weights,
               bool                common_cov,
               int                 report_type,
               int                 n_thread = 1)
{
	// hide armadillo warning / error messages

	if (!Rf_isNull(U_3d.attr("dim"))) {
		// set cube data from R 3D array
		IntegerVector dimU = U_3d.attr("dim");
		cube U_cube(U_3d.begin(), dimU[0], dimU[1], dimU[2]);
		PosteriorMASH pc(b_mat, s_mat, s_alpha_mat, s_orig_mat, v_mat, l_mat, a_mat, U_cube);
		pc.set_thread(n_thread);
		if (!common_cov) pc.compute_posterior(posterior_weights, report_type);
		else pc.compute_posterior_comcov(posterior_weights, report_type);
		return List::create(
			Named("post_mean") = pc.PosteriorMean(),
			Named("post_sd")   = pc.PosteriorSD(),
			Named("post_cov")  = pc.PosteriorCov(),
			Named("post_zero") = pc.ZeroProb(),
			Named("post_neg")  = pc.NegativeProb());
	} else {
		// U_3d is in fact a vector
		PosteriorASH pc(vectorise(b_mat),
		                vectorise(s_mat),
		                vectorise(s_alpha_mat),
		                v_mat(0, 0),
		                Rcpp::as<arma::vec>(U_3d));

		pc.compute_posterior(posterior_weights);
		return List::create(
			Named("post_mean") = pc.PosteriorMean(),
			Named("post_cov")  = pc.PosteriorCov(),
			Named("post_sd")   = pc.PosteriorSD(),
			Named("post_zero") = pc.ZeroProb(),
			Named("post_neg")  = pc.NegativeProb());
	}
} // calc_post_rcpp

// [[Rcpp::export]]
List
calc_sermix_rcpp(const arma::mat &   b_mat,
                 const arma::mat &   s_mat,
                 const arma::mat &   v_mat,
                 NumericVector vinv_3d,
                 NumericVector U_3d,
                 NumericVector Uinv_3d_drank,
                 NumericVector U0_3d,
                 const arma::mat &   posterior_mixture_weights,
                 const arma::mat &   posterior_variable_weights,
                 bool                common_cov,
                 int                 n_thread = 1)
{
	// hide armadillo warning / error messages
	if (Rf_isNull(U_3d.attr("dim")) && Rf_isNull(U0_3d.attr("dim"))) {
		throw std::invalid_argument(
			      "Either U_3d (prior matrices) or U0_3d (precomputed prior quantaties) has to be specified");
	}
	// set cube data from R 3D array
	cube U_cube;
	IntegerVector dimU = (Rf_isNull(U_3d.attr("dim"))) ? U0_3d.attr("dim") : U_3d.attr("dim");
	if (Rf_isNull(U_3d.attr("dim"))) {
		if (!common_cov) U_cube.resize(dimU[0], dimU[1], dimU[2] / b_mat.n_cols);
		else U_cube.resize(dimU[0], dimU[1], dimU[2]);
	} else {
		U_cube = cube(U_3d.begin(), dimU[0], dimU[1], dimU[2]);
	}
	MVSERMix pc(b_mat, s_mat, v_mat, U_cube);
	pc.set_thread(n_thread);
	if (!Rf_isNull(U0_3d.attr("dim"))) {
		IntegerVector dimU0 = U0_3d.attr("dim");
		cube U0_cube(U0_3d.begin(), dimU0[0], dimU0[1], dimU0[2]);
		pc.set_U0(U0_cube);
	}
	if (!Rf_isNull(vinv_3d.attr("dim"))) {
		// should be dimV = dimU
		IntegerVector dimV = vinv_3d.attr("dim");
		cube vinv_cube(vinv_3d.begin(), dimV[0], dimV[1], dimV[2]);
		pc.set_Vinv(vinv_cube);
	}
	if (!Rf_isNull(Uinv_3d_drank.attr("dim"))) {
		// inverse of prior matrices
		// relevent to computing updated scalar for prior variances
		// a feature used in mmbr package
		cube Uinv_cube_drank(Uinv_3d_drank.begin(), dimU[0], dimU[1], dimU[2]);
		pc.set_Uinv(Uinv_cube_drank);
	}
	if (!common_cov) pc.compute_posterior(posterior_mixture_weights, posterior_variable_weights);
	else pc.compute_posterior_comcov(posterior_mixture_weights,
		                         posterior_variable_weights);
	List res = List::create(
		Named("post_mean") = pc.PosteriorMean(),
		Named("post_sd")   = pc.PosteriorSD(),
		Named("post_cov")  = pc.PosteriorCov(),
		Named("post_zero") = pc.ZeroProb(),
		Named("post_neg")  = pc.NegativeProb());
	if (posterior_variable_weights.n_rows > 0) res.push_back(pc.PriorScalar(), "prior_scale_em_update");
	return res;
} // calc_sermix_rcpp

// [[Rcpp::export]]
List
fit_teem_rcpp(const arma::mat &   x_mat,
              const arma::vec &   w_vec,
              NumericVector U_3d,
              int                 maxiter,
              double              converge_tol,
              double              eigen_tol,
              bool                verbose)
{
	// Convert R 3d array to Rcpp cube
	if (Rf_isNull(U_3d.attr("dim"))) {
		throw std::invalid_argument(
			      "U_3d has to be a 3D array");
	}
	IntegerVector dimU = U_3d.attr("dim");
	cube U_cube(U_3d.begin(), dimU[0], dimU[1], dimU[2]);
	TEEM teem(x_mat, w_vec, U_cube);
	teem.fit(maxiter, converge_tol, eigen_tol, verbose);
	List res = List::create(
		Named("w")         = teem.get_w(),
		Named("U")         = teem.get_U(),
		Named("objective") = teem.get_objective(),
		Named("maxd")      = teem.get_maxd());
	return res;
}
