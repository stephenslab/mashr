// Gao Wang (c) 2017 gaow@uchicago.edu
#ifndef _UTILS_HPP
#define _UTILS_HPP
#include <cmath>
#include <armadillo>

const double LOG_2PI = std::log(2.0 * M_PI);
static const double INV_SQRT_2PI = 1.0 / std::sqrt(2.0 * M_PI);
static const double LOG_INV_SQRT_2PI = std::log(INV_SQRT_2PI);

inline arma::vec dmvnorm(const arma::mat & x,
                         const arma::rowvec & mean,
                         const arma::mat & sigma,
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
	double constants = -(static_cast<double>(xdim) / 2.0) * LOG_2PI;

	for (int i = 0; i < n; i++) {
		arma::vec z = rooti * arma::trans(x.row(i) - mean) ;
		out(i) = constants - 0.5 * arma::sum(z % z) + rootisum;
	}

	if (logd == false) {
		out = arma::exp(out);
	}
	return out;
}


inline double dmvnorm(const arma::vec & x,
                      const arma::vec & mean,
                      const arma::mat & sigma,
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
	double constants = -(static_cast<double>(x.n_elem) / 2.0) * LOG_2PI;

	arma::vec z = rooti * (x - mean) ;
	double out = constants - 0.5 * arma::sum(z % z) + rootisum;

	if (logd == false) {
		out = std::exp(out);
	}
	return out;
}


inline arma::vec pnorm(const arma::vec & x, const arma::vec & m, const arma::vec & s,
                       bool logd = false, bool lower_tail = true)
{
	arma::vec res(x.n_elem);
	for (unsigned i = 0; i < x.n_elem; i++) {
		// FIXME: not sure if erfc approximation is accurate
		// see `normalCDF` function at:
		// http://en.cppreference.com/w/cpp/numeric/math/erfc
		res(i) = 0.5 * std::erfc(-(x(i) - m(i)) / s(i) * M_SQRT1_2);
	}
	if (!lower_tail & !logd) {
		return 1 - res;
	} else if (lower_tail & !logd) {
		return res;
	} else if (!lower_tail & logd) {
		return arma::log(1 - res);
	} else {  // (lower_tail & logd)
		return arma::log(res);
	}
}


// a quicker way to compute diag(s) %*% V %*% diag(s)
inline arma::mat get_cov(const arma::vec & s, const ::arma::mat & V)
{
	return (V.each_col() % s).each_row() % s.t();
}


// @title calc_lik
// @description computes matrix of likelihoods for each of J cols of Bhat for each of P prior covariances
// @param b_mat R by J
// @param s_mat R by J
// @param v_mat R by R
// @param U_cube list of prior covariance matrices
// @param logd if true computes log-likelihood
// @return J x P matrix of multivariate normal likelihoods, p(bhat | U[p], V)
arma::mat calc_lik(const arma::mat & b_mat,
                   const arma::mat & s_mat,
                   const arma::mat & v_mat,
                   const arma::cube & U_cube,
                   bool logd = false)
{
	// In armadillo data are stored with column-major ordering
	// slicing columns are therefore faster than rows
	// lik is a J by P matrix
	arma::mat lik(b_mat.n_cols, U_cube.n_slices, arma::fill::zeros);
	arma::vec mean(b_mat.n_rows, arma::fill::zeros);

	for (unsigned j = 0; j < lik.n_rows; ++j) {
		arma::mat sigma = get_cov(s_mat.col(j), v_mat);
		for (unsigned p = 0; p < lik.n_cols; ++p) {
			lik.at(j, p) = dmvnorm(b_mat.col(j), mean, sigma + U_cube.slice(p), logd);
		}
	}
	return lik;
}


// @param b_mat R by J
// @param s_mat R by J
// @param v_mat R by R
// @param U_cube list of prior covariance matrices
// @param logd if true computes log-likelihood
class PosteriorCalculator
{
public:
	PosteriorCalculator(const arma::mat & b_mat,
	                    const arma::mat & s_mat,
	                    const arma::mat & v_mat,
	                    const arma::cube & U_cube) :
		b_mat(b_mat), s_mat(s_mat), v_mat(v_mat), U_cube(U_cube) {}
	~PosteriorCalculator() {}

private:
	arma::mat b_mat;
	arma::mat s_mat;
	arma::mat v_mat;
	arma::mat U_cube;
};


#endif
