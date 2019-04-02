// Gao Wang (c) 2017 gaow@uchicago.edu
#include <iostream>
#include "RcppArmadillo.h"
#include "mash.h"
// [[Rcpp::depends(RcppArmadillo)]]

// Enable C++11.
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
Rcpp::List calc_lik_rcpp(Rcpp::NumericMatrix b_mat,
                         Rcpp::NumericMatrix s_mat,
                         Rcpp::NumericMatrix v_mat,
                         Rcpp::NumericMatrix l_mat,
                         Rcpp::NumericVector U_3d,
                         bool logd,
                         bool common_cov)
{

	// hide armadillo warning / error messages
	// std::ostream nullstream(0);
	// arma::set_stream_err2(nullstream);
	arma::mat res;
	if (!Rf_isNull(U_3d.attr("dim"))) {
		// set cube data from R 3D array
		Rcpp::IntegerVector dimU = U_3d.attr("dim");
		arma::cube U_cube(U_3d.begin(), dimU[0], dimU[1], dimU[2]);
		res = calc_lik(Rcpp::as<arma::mat>(b_mat),
			Rcpp::as<arma::mat>(s_mat),
			Rcpp::as<arma::mat>(v_mat),
			Rcpp::as<arma::mat>(l_mat),
			U_cube,
			logd,
			common_cov);
	} else {
		res = calc_lik(Rcpp::as<arma::vec>(b_mat),
			Rcpp::as<arma::vec>(s_mat),
			v_mat(0, 0),
			Rcpp::as<arma::vec>(U_3d),
			logd);
	}
	return Rcpp::List::create(Rcpp::Named("data") = res,
		Rcpp::Named("status") = 0);
}

// [[Rcpp::export]]
// This only works for multivariate case (U is 3D array)
Rcpp::List calc_lik_v_known_rcpp(Rcpp::NumericMatrix b_mat,
                         Rcpp::NumericVector V_3d,
                         Rcpp::NumericVector U_3d,
                         bool logd)
{

	// hide armadillo warning / error messages
	// std::ostream nullstream(0);
	// arma::set_stream_err2(nullstream);
	arma::mat res;
	// set cube data from R 3D array
	Rcpp::IntegerVector dimU = U_3d.attr("dim");
	arma::cube U_cube(U_3d.begin(), dimU[0], dimU[1], dimU[2]);
	Rcpp::IntegerVector dimV = V_3d.attr("dim");
	if (dimV.size() == 3) {
		arma::cube V_cube(V_3d.begin(), dimV[0], dimV[1], dimV[2]);
		res = calc_lik(Rcpp::as<arma::mat>(b_mat),
			V_cube,
			U_cube,
			logd);
	} else {
		// common V
		arma::mat V_mat(V_3d.begin(), dimV[0], dimV[1]);
		res = calc_lik(Rcpp::as<arma::vec>(b_mat),
			V_mat,
			U_cube,
			logd);
	}
	return Rcpp::List::create(Rcpp::Named("data") = res,
		Rcpp::Named("status") = 0);
}


// [[Rcpp::export]]
Rcpp::List calc_post_rcpp(Rcpp::NumericMatrix b_mat,
                          Rcpp::NumericMatrix s_mat,
                          Rcpp::NumericMatrix s_alpha_mat,
                          Rcpp::NumericMatrix s_orig_mat,
                          Rcpp::NumericMatrix v_mat,
                          Rcpp::NumericMatrix l_mat,
                          Rcpp::NumericMatrix a_mat,
                          Rcpp::NumericVector U_3d,
                          Rcpp::NumericMatrix posterior_weights,
                          bool common_cov,
                          int report_type)
{

	// hide armadillo warning / error messages
	// std::ostream nullstream(0);
	// arma::set_stream_err2(nullstream);

	if (!Rf_isNull(U_3d.attr("dim"))) {
		// set cube data from R 3D array
		Rcpp::IntegerVector dimU = U_3d.attr("dim");
		arma::cube U_cube(U_3d.begin(), dimU[0], dimU[1], dimU[2]);
		PosteriorMASH pc(Rcpp::as<arma::mat>(b_mat),
		                 Rcpp::as<arma::mat>(s_mat),
		                 Rcpp::as<arma::mat>(s_alpha_mat),
		                 Rcpp::as<arma::mat>(s_orig_mat),
		                 Rcpp::as<arma::mat>(v_mat),
		                 Rcpp::as<arma::mat>(l_mat),
		                 Rcpp::as<arma::mat>(a_mat),
		                 U_cube);
		if (!common_cov) pc.compute_posterior(Rcpp::as<arma::mat>(posterior_weights), report_type);
		else pc.compute_posterior_comcov(Rcpp::as<arma::mat>(posterior_weights), report_type);
		return Rcpp::List::create(
			Rcpp::Named("post_mean") = pc.PosteriorMean(),
			Rcpp::Named("post_sd") = pc.PosteriorSD(),
			Rcpp::Named("post_cov") = pc.PosteriorCov(),
			Rcpp::Named("post_zero") = pc.ZeroProb(),
			Rcpp::Named("post_neg") = pc.NegativeProb());
	} else {
		// U_3d is in fact a vector
		PosteriorASH pc(Rcpp::as<arma::vec>(b_mat),
		                Rcpp::as<arma::vec>(s_mat),
		                Rcpp::as<arma::vec>(s_alpha_mat),
		                v_mat(0, 0),
		                Rcpp::as<arma::vec>(U_3d));

		pc.compute_posterior(Rcpp::as<arma::mat>(posterior_weights));
		return Rcpp::List::create(
			Rcpp::Named("post_mean") = pc.PosteriorMean(),
			Rcpp::Named("post_cov") = pc.PosteriorCov(),
			Rcpp::Named("post_sd") = pc.PosteriorSD(),
			Rcpp::Named("post_zero") = pc.ZeroProb(),
			Rcpp::Named("post_neg") = pc.NegativeProb());
	}
}


// [[Rcpp::export]]
// This only works for multivariate case (U is 3D array)
Rcpp::List calc_post_v_known_rcpp(Rcpp::NumericMatrix b_mat,
                          Rcpp::NumericMatrix a_mat,
                          Rcpp::NumericVector V_3d,
                          Rcpp::NumericVector Vinv_3d,
                          Rcpp::NumericVector U_3d,
                          Rcpp::NumericMatrix posterior_weights,
                          int report_type)
{

	// hide armadillo warning / error messages
	// std::ostream nullstream(0);
	// arma::set_stream_err2(nullstream);
	// set cube data from R 3D array
	Rcpp::IntegerVector dimU = U_3d.attr("dim");
	arma::cube U_cube(U_3d.begin(), dimU[0], dimU[1], dimU[2]);
	Rcpp::IntegerVector dimV = V_3d.attr("dim");
	if (dimV.size() == 3) {
		arma::cube V_cube(V_3d.begin(), dimV[0], dimV[1], dimV[2]);
		arma::cube Vinv_cube(V_3d.begin(), dimV[0], dimV[1], dimV[2]);
		PosteriorMASH_V_known pc(Rcpp::as<arma::mat>(b_mat),
		                 Rcpp::as<arma::mat>(a_mat),
		                 V_cube, Vinv_cube,
		                 U_cube);
		pc.compute_posterior(Rcpp::as<arma::mat>(posterior_weights), report_type);
		return Rcpp::List::create(
			Rcpp::Named("post_mean") = pc.PosteriorMean(),
			Rcpp::Named("post_sd") = pc.PosteriorSD(),
			Rcpp::Named("post_cov") = pc.PosteriorCov(),
			Rcpp::Named("post_zero") = pc.ZeroProb(),
			Rcpp::Named("post_neg") = pc.NegativeProb());
	} else {
		// common V
		arma::mat V_mat(V_3d.begin(), dimV[0], dimV[1]);
		arma::mat Vinv_mat(Vinv_3d.begin(), dimV[0], dimV[1]);
		PosteriorMASH_V_known pc(Rcpp::as<arma::mat>(b_mat),
		                 Rcpp::as<arma::mat>(a_mat),
		                 V_mat, Vinv_mat,
		                 U_cube);
		pc.compute_posterior_comcov(Rcpp::as<arma::mat>(posterior_weights), report_type);
		return Rcpp::List::create(
			Rcpp::Named("post_mean") = pc.PosteriorMean(),
			Rcpp::Named("post_sd") = pc.PosteriorSD(),
			Rcpp::Named("post_cov") = pc.PosteriorCov(),
			Rcpp::Named("post_zero") = pc.ZeroProb(),
			Rcpp::Named("post_neg") = pc.NegativeProb());
	}
}