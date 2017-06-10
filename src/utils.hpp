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
	arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
	double rootisum = arma::sum(log(rooti.diag()));
	double constants = -(static_cast<double>(xdim) / 2.0) * log2pi;

	for (int i = 0; i < n; i++) {
		arma::vec z = rooti * arma::trans(x.row(i) - mean) ;
		out(i) = constants - 0.5 * arma::sum(z % z) + rootisum;
	}

	if (logd == false) {
		out = exp(out);
	}
	return(out);
}

// @title calc_lik
// @description computes matrix of likelihoods for each of J rows of Bhat for each of P prior covariances
// @param b_mat
// @param s_mat
// @param v_mat
// @param U_cube list of prior covariance matrices
// @param logd if true computes log-likelihood
// @return J x P vector of multivariate normal likelihoods, p(bhat | U[p], V)
arma::mat calc_lik(arma::mat b_mat,
                   arma::mat s_mat,
                   arma::mat v_mat,
                   arma::cube U_cube,
                   bool logd = false) {
  b_mat.print("b_mat");
  s_mat.print("s_mat");
  v_mat.print("v_mat");
  U_cube.print("U_cube");
  return b_mat;
}

#endif
