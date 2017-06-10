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
		// FIXME: not sure if erfc approximation is accurate enough compared to R's pnorm()
		// see `normalCDF` function at:
		// http://en.cppreference.com/w/cpp/numeric/math/erfc
		res(i) = 0.5 * std::erfc(-(x(i) - m(i)) / s(i) * M_SQRT1_2);
	}
	if (!lower_tail & !logd) {
		return 1.0 - res;
	} else if (lower_tail & !logd) {
		return res;
	} else if (!lower_tail & logd) {
		return arma::log(1.0 - res);
	} else {  // (lower_tail & logd)
		return arma::log(res);
	}
}


// a quicker way to compute diag(s) %*% V %*% diag(s)
inline arma::mat get_cov(const arma::vec & s, const ::arma::mat & V)
{
	return (V.each_col() % s).each_row() % s.t();
}


// @title posterior_cov
// @param Vinv R x R inverse covariance matrix for the likelihood
// @param U R x R prior covariance matrix
// @return R x R posterior covariance matrix
// @description If bhat is N(b,V) and b is N(0,U) then b|bhat N(mu1,U1). This function returns U1.
inline arma::mat get_posterior_cov(const arma::mat & Vinv, const arma::mat & U)
{
	// U %*% solve(Vinv %*% U + diag(nrow(U)))
	arma::mat S = Vinv * U;
	S.diag() += static_cast<double>(U.n_rows);
	return U * S.i();
}


// @title posterior_mean
// @param bhat R vector of observations
// @param Vinv R x R inverse covariance matrix for the likelihood
// @param U1 R x R posterior covariance matrix, computed using posterior_cov
// @return R vector of posterior mean
// @description If bhat is N(b,V) and b is N(0,U) then b|bhat N(mu1,U1). This function returns mu1.
inline arma::vec get_posterior_mean(const arma::vec & bhat, const arma::mat & Vinv,
                                    const arma::mat & U1)
{
	return U1 * Vinv * bhat;
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
// @param U_cube list of prior covariance matrices, for each mixture component P by R by R
// @param logd if true computes log-likelihood
class PosteriorCalculator
{
public:
	PosteriorCalculator(const arma::mat & b_mat,
	                    const arma::mat & s_mat,
	                    const arma::mat & v_mat,
	                    const arma::cube & U_cube) :
		b_mat(b_mat), s_mat(s_mat), v_mat(v_mat), U_cube(U_cube)
	{
		int J = b_mat.n_cols, R = b_mat.n_rows;

		post_mean.set_size(R, J);
		neg_prob.set_size(R, J);
		zero_prob.set_size(R, J);
	}


	~PosteriorCalculator() {}

	// @title Compute posterior matrices
	// @description More detailed description of function goes here.
	// @param posterior_weights P X J matrix, the posterior probabilities of each mixture component for each effect
	int compute_posterior(const arma::mat & posterior_weights)
	{
		for (unsigned j = 0; j < b_mat.n_cols; ++j) {
      // FIXME: improved math may help here
			arma::mat Vinv = get_cov(s_mat.col(j), v_mat).i();
      arma::vec mean(b_mat.n_rows, arma::fill::zeros);
      // R X P matrices
      arma::mat mu1_mat(b_mat.n_rows, U_cube.n_slices, arma::fill::zeros);
      arma::mat mu2_mat(b_mat.n_rows, U_cube.n_slices, arma::fill::zeros);
      arma::mat zero_mat(b_mat.n_rows, U_cube.n_slices, arma::fill::zeros);
      arma::mat neg_mat(b_mat.n_rows, U_cube.n_slices, arma::fill::zeros);
      for (unsigned p = 0; p < U_cube.n_slices; ++p) {
        arma::mat U1 = get_posterior_cov(Vinv, U_cube.slice(p));
        mu1_mat.col(p) = get_posterior_mean(b_mat.col(j), Vinv, U1);
        arma::vec sigma = U1.diag(); // U1.diag() is the posterior covariance
        mu2_mat.col(p) = arma::pow(mu1_mat.col(p), 2) + sigma;
        // FIXME: please double-check the implementation logic here
        // I'm not sure if it is correct
        // against https://github.com/stephenslab/mashr/blob/master/R/posterior.R#L83
        neg_mat.col(p) = pnorm(mean, mu1_mat.col(p), arma::sqrt(sigma));
        for (unsigned r = 0; r < sigma.n_elem; ++r) {
          if (sigma.at(r) == 0) {
            zero_mat.at(r, p) = 1.0;
            neg_mat.at(r, p) = 0.0;
          }
        }
      }
      // compute weighted means of posterior arrays
      post_mean.col(j) = mu1_mat * posterior_weights.col(j);
      post_mean2.col(j) = mu2_mat * posterior_weights.col(j);
      neg_prob.col(j) = neg_mat * posterior_weights.col(j);
      zero_prob.col(j) = zero_prob * posterior_weights.col(j);
		}
    return 0;
	}


	// @return PosteriorMean JxR matrix of posterior means
	// @return PosteriorSD JxR matrix of posterior (marginal) standard deviations
	// @return NegativeProb JxR matrix of posterior (marginal) probability of being negative
	// @return ZeroProb JxR matrix of posterior (marginal) probability of being zero
	arma::mat PosteriorMean() { return post_mean; }
	arma::mat PosteriorSD() { return arma::sqrt(post_mean2 - arma::pow(post_mean, 2)); }
	arma::mat NegativeProb() { return neg_prob; }
	arma::mat ZeroProb() { return zero_prob; }

private:
	// input
	arma::mat b_mat;
	arma::mat s_mat;
	arma::mat v_mat;
	arma::cube U_cube;
	// output
	// all R X J mat
	arma::mat post_mean;
	arma::mat post_mean2;
	arma::mat neg_prob;
	arma::mat zero_prob;
};


#endif
