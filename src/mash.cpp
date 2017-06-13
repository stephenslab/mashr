// Gao Wang (c) 2017 gaow@uchicago.edu
#include <iostream>
#include "RcppArmadillo.h"
#include "mashutils.hpp"

//' @export
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List calc_lik_rcpp (Rcpp::NumericMatrix b_mat,
			  Rcpp::NumericMatrix s_mat,
			  Rcpp::NumericMatrix v_mat,
			  Rcpp::NumericVector U_3d,
			  bool logd = false) {
  
	// hide armadillo warning / error messages
	std::ostream nullstream(0);
	arma::set_stream_err2(nullstream);
	arma::mat res;

	if (!Rf_isNull(U_3d.attr("dim"))) {
		// set cube data from R 3D array
		Rcpp::IntegerVector dimU = U_3d.attr("dim");
		arma::cube U_cube(U_3d.begin(), dimU[0], dimU[1], dimU[2]);
		res = calc_lik(Rcpp::as<arma::mat>(b_mat),
			Rcpp::as<arma::mat>(s_mat),
			Rcpp::as<arma::mat>(v_mat),
			U_cube,
			logd);
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

//' @export
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List calc_post_rcpp(Rcpp::NumericMatrix b_mat,
                          Rcpp::NumericMatrix s_mat,
                          Rcpp::NumericMatrix v_mat,
                          Rcpp::NumericVector U_3d,
                          Rcpp::NumericMatrix posterior_weights) {

	// hide armadillo warning / error messages
	std::ostream nullstream(0);
	arma::set_stream_err2(nullstream);

	if (!Rf_isNull(U_3d.attr("dim"))) {
		// set cube data from R 3D array
		Rcpp::IntegerVector dimU = U_3d.attr("dim");
		arma::cube U_cube(U_3d.begin(), dimU[0], dimU[1], dimU[2]);
		PosteriorMASH pc(Rcpp::as<arma::mat>(b_mat),
		                 Rcpp::as<arma::mat>(s_mat),
		                 Rcpp::as<arma::mat>(v_mat),
		                 U_cube);

		pc.compute_posterior(Rcpp::as<arma::mat>(posterior_weights));
		return Rcpp::List::create(Rcpp::Named("post_mean") = 
					    pc.PosteriorMean(),
			Rcpp::Named("post_sd") = pc.PosteriorSD(),
			Rcpp::Named("post_zero") = pc.ZeroProb(),
			Rcpp::Named("post_neg") = pc.NegativeProb());
	} else {
		// U_3d is in fact a vector
		PosteriorASH pc(Rcpp::as<arma::vec>(b_mat),
		                Rcpp::as<arma::vec>(s_mat),
		                v_mat(0, 0),
		                Rcpp::as<arma::vec>(U_3d));

		pc.compute_posterior(Rcpp::as<arma::mat>(posterior_weights));
		return Rcpp::List::create(Rcpp::Named("post_mean") = 
					    pc.PosteriorMean(),
			Rcpp::Named("post_sd") = pc.PosteriorSD(),
			Rcpp::Named("post_zero") = pc.ZeroProb(),
			Rcpp::Named("post_neg") = pc.NegativeProb());
	}
}


