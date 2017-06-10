// Gao Wang (c) 2017 gaow@uchicago.edu
#ifndef _UTILS_HPP
#define _UTILS_HPP
#include <armadillo>

const double log2pi = std::log(2.0 * M_PI);

arma::vec dmvnorm(arma::mat x,
                  arma::rowvec mean,
                  arma::mat sigma,
                  bool logd = false)
{
	int n = x.n_rows;
	int xdim = x.n_cols;
	arma::vec out(n);
	arma::mat rooti;

	try {
		rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
	} catch (const std::runtime_error & error) {
		if (logd) out.fill(-arma::datum::inf);
		else out.fill(0.0);
		return out;
	}
	double rootisum = arma::sum(arma::log(rooti.diag()));
	double constants = -(static_cast<double>(xdim) / 2.0) * log2pi;

	for (int i = 0; i < n; i++) {
		arma::vec z = rooti * arma::trans(x.row(i) - mean) ;
		out(i) = constants - 0.5 * arma::sum(z % z) + rootisum;
	}

	if (logd == false) {
		out = arma::exp(out);
	}
	return out;
}


double dmvnorm(arma::vec x,
               arma::vec mean,
               arma::mat sigma,
               bool logd = false)
{
	arma::mat rooti;

	try {
		rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
	} catch (const std::runtime_error & error) {
		if (logd) return -arma::datum::inf;
		else return 0.0;
	}
	double rootisum = arma::sum(arma::log(rooti.diag()));
	double constants = -(static_cast<double>(x.n_elem) / 2.0) * log2pi;

	arma::vec z = rooti * (x - mean) ;
	double out = constants - 0.5 * arma::sum(z % z) + rootisum;

	if (logd == false) {
		out = std::exp(out);
	}
	return out;
}


// @title calc_lik
// @description computes matrix of likelihoods for each of J cols of Bhat for each of P prior covariances
// @param b_mat R by J
// @param s_mat R by J
// @param v_mat R by R
// @param U_cube list of prior covariance matrices
// @param logd if true computes log-likelihood
// @return J x P matrix of multivariate normal likelihoods, p(bhat | U[p], V)
arma::mat calc_lik(arma::mat b_mat,
                   arma::mat s_mat,
                   arma::mat v_mat,
                   arma::cube U_cube,
                   bool logd = false)
{
	// In armadillo data are stored with column-major ordering
	// slicing columns are therefore faster than rows
	// lik is a J by P matrix
	arma::mat lik(b_mat.n_cols, U_cube.n_slices, arma::fill::zeros);
	arma::vec mean(b_mat.n_rows, arma::fill::zeros);

	for (unsigned j = 0; j < lik.n_rows; ++j) {
		arma::mat sigma = (v_mat.each_col() % s_mat.col(j)).each_row() % s_mat.col(j).t(); // quicker than diagmat(s) * v diagmat(s)
		for (unsigned p = 0; p < lik.n_cols; ++p) {
			lik.at(j, p) = dmvnorm(b_mat.col(j), mean, sigma + U_cube.slice(p), logd);
		}
	}
	return lik;
}


#endif
