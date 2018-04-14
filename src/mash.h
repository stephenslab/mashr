// Gao Wang (c) 2017 gaow@uchicago.edu
#ifndef _MASH_H
#define _MASH_H
#include <cmath>
#include <armadillo>

const double LOG_2PI = std::log(2.0 * M_PI);
static const double INV_SQRT_2PI = 1.0 / std::sqrt(2.0 * M_PI);
static const double LOG_INV_SQRT_2PI = std::log(INV_SQRT_2PI);

inline arma::vec dnorm(const arma::vec & x,
                       const arma::vec & mu,
                       const arma::vec & sigma2,
                       bool logd = false)
{
	arma::vec res = LOG_INV_SQRT_2PI -
	                arma::log(arma::sqrt(sigma2)) -
	                arma::pow(x - mu, 2.0) / (2.0 * sigma2);

	if (logd) return res;
	else return arma::exp(res);
}


inline arma::vec dmvnorm_mat(const arma::mat & x,
                             const arma::vec & mean,
                             const arma::mat & sigma,
                             bool logd = false)
{
	double xdim = static_cast<double>(x.n_rows);

	arma::vec out(x.n_cols);
	arma::mat rooti;

	try {
		rooti = arma::trans(arma::inv(arma::trimatu(arma::chol(sigma))));
	} catch (const std::runtime_error & error) {
		if (logd) out.fill(-arma::datum::inf);
		else out.fill(0.0);
		for (arma::uword i = 0; i < x.n_cols; ++i)
			if (arma::accu(arma::abs(x.col(i) - mean)) < 1e-6) out.at(i) = arma::datum::inf;
		return out;
	}
	double rootisum = arma::sum(arma::log(rooti.diag()));
	double constants = -(xdim / 2.0) * LOG_2PI;

	for (unsigned i = 0; i < x.n_cols; i++) {
		arma::vec z = rooti * (x.col(i) - mean) ;
		out.at(i) = constants - 0.5 * arma::sum(z % z) + rootisum;
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
		rooti = arma::trans(arma::inv(arma::trimatu(arma::chol(sigma))));
	} catch (const std::runtime_error & error) {
		double diff = arma::accu(arma::abs(x - mean));
		if (logd) return (diff < 1e-6) ? arma::datum::inf : -arma::datum::inf;
		else return (diff < 1e-6) ? arma::datum::inf : 0.0;
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


template <class T, class U>
inline T pnorm(const U & x, const T & m, const T & s,
               bool logd = false, bool lower_tail = true)
{
	// FIXME: not sure if erfc approximation is accurate enough compared to R's pnorm()
	// see `normalCDF` function at:
	// http://en.cppreference.com/w/cpp/numeric/math/erfc
	T res = 0.5 * arma::erfc((x - m) / s * M_SQRT1_2);

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
	/* return arma::diagmat(s) * V * arma::diagmat(s); */
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

	S.diag() += 1.0;
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


inline arma::mat get_posterior_mean_mat(const arma::mat & bhat, const arma::mat & Vinv,
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
// @param common_cov if true use version for common covariance
// @return J x P matrix of multivariate normal likelihoods, p(bhat | U[p], V)
arma::mat calc_lik(const arma::mat & b_mat,
                   const arma::mat & s_mat,
                   const arma::mat & v_mat,
                   const arma::cube & U_cube,
                   bool logd,
                   bool common_cov)
{
	// In armadillo data are stored with column-major ordering
	// slicing columns are therefore faster than rows
	// lik is a J by P matrix
	arma::mat lik(b_mat.n_cols, U_cube.n_slices, arma::fill::zeros);

	if (common_cov) {
		arma::vec mean(b_mat.n_rows, arma::fill::zeros);
		arma::mat sigma = get_cov(s_mat.col(0), v_mat);
		for (arma::uword p = 0; p < lik.n_cols; ++p) {
			lik.col(p) = dmvnorm_mat(b_mat, mean, sigma + U_cube.slice(p), logd);
		}
	} else {
		arma::vec mean(b_mat.n_rows, arma::fill::zeros);
		for (arma::uword j = 0; j < lik.n_rows; ++j) {
			arma::mat sigma = get_cov(s_mat.col(j), v_mat);
			for (arma::uword p = 0; p < lik.n_cols; ++p) {
				lik.at(j, p) = dmvnorm(b_mat.col(j), mean, sigma + U_cube.slice(p), logd);
			}
		}
	}
	return lik;
}


// @title calc_lik univariate version
// @description computes matrix of likelihoods for each of J cols of Bhat for each of P prior sigma
// @param b_vec of J
// @param s_vec of J
// @param v numeric
// @param U_vec P vector
// @param logd if true computes log-likelihood
// @return J x P matrix of multivariate normal likelihoods, p(bhat | U[p], V)
arma::mat calc_lik(const arma::vec & b_vec,
                   const arma::vec & s_vec,
                   double v,
                   const arma::vec & U_vec,
                   bool logd)
{
	arma::mat lik(b_vec.n_elem, U_vec.n_elem, arma::fill::zeros);
	arma::vec sigma = s_vec % s_vec * v;
	arma::vec mean(b_vec.n_elem, arma::fill::zeros);

	for (arma::uword p = 0; p < lik.n_cols; ++p) {
		lik.col(p) = dnorm(b_vec, mean, sigma + U_vec.at(p), logd);
	}
	return lik;
}


// @param b_mat R by J
// @param s_mat R by J
// @param v_mat R by R
// @param U_cube list of prior covariance matrices, for each mixture component P by R by R
class PosteriorMASH
{
public:
	PosteriorMASH(const arma::mat & b_mat,
	              const arma::mat & s_mat,
	              const arma::mat & v_mat,
	              const arma::cube & U_cube) :
		b_mat(b_mat), s_mat(s_mat), v_mat(v_mat), U_cube(U_cube)
	{
		int J = b_mat.n_cols, R = b_mat.n_rows;

		post_mean.set_size(R, J);
		post_var.set_size(R, J);
		post_cov.set_size(J, R, R);
		neg_prob.set_size(R, J);
		zero_prob.set_size(R, J);
		post_mean.zeros();
		post_var.zeros();
		post_cov.zeros();
		neg_prob.zeros();
		zero_prob.zeros();
	}


	~PosteriorMASH() {}

	// @title Compute posterior matrices
	// @description More detailed description of function goes here.
	// @param posterior_weights P X J matrix, the posterior probabilities of each mixture component for each effect
	// @param report_post_cov Boolean variable that decides whether or not posterior covariance should be computed
	int compute_posterior(const arma::mat & posterior_weights, const bool report_post_cov)
	{
		arma::vec mean(b_mat.n_rows, arma::fill::zeros);

		for (arma::uword j = 0; j < b_mat.n_cols; ++j) {
			// FIXME: improved math may help here
			arma::mat Vinv = arma::inv_sympd(get_cov(s_mat.col(j), v_mat));
			// R X P matrices
			arma::mat mu1_mat(b_mat.n_rows, U_cube.n_slices, arma::fill::zeros);
			arma::mat mu2_mat(b_mat.n_rows, U_cube.n_slices, arma::fill::zeros);
			arma::mat zero_mat(b_mat.n_rows, U_cube.n_slices, arma::fill::zeros);
			arma::mat neg_mat(b_mat.n_rows, U_cube.n_slices, arma::fill::zeros);
			for (arma::uword p = 0; p < U_cube.n_slices; ++p) {
				arma::mat U1 = get_posterior_cov(Vinv, U_cube.slice(p));
				mu1_mat.col(p) = get_posterior_mean(b_mat.col(j), Vinv, U1);
				if (report_post_cov) {
					post_cov.slice(j) +=
						posterior_weights.at(p, j) * (U1 + mu1_mat.col(p) * mu1_mat.col(p).t());
				}
				arma::vec sigma = arma::sqrt(U1.diag()); // U1.diag() is the posterior covariance
				mu2_mat.col(p) = arma::pow(mu1_mat.col(p), 2.0) + U1.diag();
				neg_mat.col(p) = pnorm(mu1_mat.col(p), mean, sigma);
				for (arma::uword r = 0; r < sigma.n_elem; ++r) {
					if (sigma.at(r) == 0) {
						zero_mat.at(r, p) = 1.0;
						neg_mat.at(r, p) = 0.0;
					}
				}
			}
			// compute weighted means of posterior arrays
			post_mean.col(j) = mu1_mat * posterior_weights.col(j);
			post_var.col(j) = mu2_mat * posterior_weights.col(j);
			neg_prob.col(j) = neg_mat * posterior_weights.col(j);
			zero_prob.col(j) = zero_mat * posterior_weights.col(j);
		}
		post_var -= arma::pow(post_mean, 2.0);
		return 0;
	}


	// @title Compute posterior matrices when covariance SVS is the same for all J conditions
	// @description More detailed description of function goes here.
	// @param posterior_weights P X J matrix, the posterior probabilities of each mixture component for each effect
	// @param report_post_cov Boolean variable that decides whether or not posterior covariance should be computed
	int compute_posterior_comcov(const arma::mat & posterior_weights, const bool report_post_cov)
	{
		arma::mat mean(b_mat.n_rows, b_mat.n_cols, arma::fill::zeros);
		// R X R
		arma::mat Vinv = arma::inv_sympd(get_cov(s_mat.col(0), v_mat));
		arma::rowvec ones(b_mat.n_cols, arma::fill::ones);
		arma::rowvec zeros(b_mat.n_cols, arma::fill::zeros);
		arma::mat sigma(b_mat.n_rows, b_mat.n_cols, arma::fill::zeros);
		for (arma::uword p = 0; p < U_cube.n_slices; ++p) {
			arma::mat zero_mat(b_mat.n_rows, b_mat.n_cols, arma::fill::zeros);
			// R X R
			arma::mat U1 = get_posterior_cov(Vinv, U_cube.slice(p));
			// R X J
			arma::mat mu1_mat = get_posterior_mean_mat(b_mat, Vinv, U1);
			// FIXME: better initialization?
			arma::vec Svec = arma::sqrt(U1.diag()); // U1.diag() is the posterior covariance
			for (arma::uword j = 0; j < sigma.n_cols; ++j) {
				sigma.col(j) = Svec;
				if (report_post_cov) {
					post_cov.slice(j) +=
						posterior_weights.at(p, j) * (U1 + mu1_mat.col(j) * mu1_mat.col(j).t());
				}
			}
			// R X J
			arma::mat mu2_mat = arma::pow(mu1_mat, 2.0);
			mu2_mat.each_col() += U1.diag();
			// R X J
			arma::mat neg_mat = pnorm(mu1_mat, mean, sigma);
			for (arma::uword r = 0; r < Svec.n_elem; ++r) {
				if (Svec.at(r) == 0) {
					zero_mat.row(r) = ones;
					neg_mat.row(r) = zeros;
				}
			}
			// compute weighted means of posterior arrays
			post_mean += mu1_mat.each_row() % posterior_weights.row(p);
			post_var += mu2_mat.each_row() % posterior_weights.row(p);
			neg_prob += neg_mat.each_row() % posterior_weights.row(p);
			zero_prob += zero_mat.each_row() % posterior_weights.row(p);
		}
		post_var -= arma::pow(post_mean, 2.0);
		return 0;
	}


	// @return PosteriorMean JxR matrix of posterior means
	// @return PosteriorSD JxR matrix of posterior (marginal) standard deviations
	// @return NegativeProb JxR matrix of posterior (marginal) probability of being negative
	// @return ZeroProb JxR matrix of posterior (marginal) probability of being zero
	arma::mat PosteriorMean() { return post_mean.t(); }
	arma::mat PosteriorSD() { return arma::sqrt(post_var).t(); }
	arma::cube PosteriorCov() { return post_cov; }
	arma::mat NegativeProb() { return neg_prob.t(); }
	arma::mat ZeroProb() { return zero_prob.t(); }

private:
	// input
	arma::mat b_mat;
	arma::mat s_mat;
	arma::mat v_mat;
	arma::cube U_cube;
	// output
	// all R X J mat
	arma::mat post_mean;
	arma::mat post_var;
	arma::mat neg_prob;
	arma::mat zero_prob;
	// J X R X R cube
	arma::cube post_cov;
};

// @param b_vec of J
// @param s_vec of J
// @param v double
// @param U_vec of P
class PosteriorASH
{
public:
	PosteriorASH(const arma::vec & b_vec,
	             const arma::vec & s_vec,
	             double v,
	             const arma::vec & U_vec) :
		b_vec(b_vec), s_vec(s_vec), v(v), U_vec(U_vec)
	{
		int J = b_vec.n_elem;

		post_mean.set_size(J);
		post_var.set_size(J);
		neg_prob.set_size(J);
		zero_prob.set_size(J);
	}


	~PosteriorASH() {}

	// @title Compute posterior matrices
	// @description univariate version of PosteriorMASH::compute_posterior(), same logic
	// @param posterior_weights P X J matrix, the posterior probabilities of each mixture component for each effect
	int compute_posterior(const arma::mat & posterior_weights)
	{
		arma::vec sv = s_vec % s_vec * v;
		unsigned J = b_vec.n_elem;
		unsigned P = U_vec.n_elem;
		arma::vec mean(J, arma::fill::zeros);
		// J X P matrices
		arma::mat mu1_mat(J, P, arma::fill::zeros);
		arma::mat mu2_mat(J, P, arma::fill::zeros);
		arma::mat zero_mat(J, P, arma::fill::zeros);
		arma::mat neg_mat(J, P, arma::fill::zeros);

		for (arma::uword p = 0; p < P; ++p) {
			arma::vec U1 = U_vec.at(p) / (sv * U_vec.at(p) + 1);
			mu1_mat.col(p) = U1 / sv % b_vec;
			mu2_mat.col(p) = arma::pow(mu1_mat.col(p), 2.0) + U1;
			arma::vec sigma = arma::sqrt(U1);
			neg_mat.col(p) = pnorm(mu1_mat.col(p), mean, sigma);
			for (arma::uword j = 0; j < J; ++j) {
				if (U1.at(j) == 0) {
					zero_mat.at(j, p) = 1.0;
					neg_mat.at(j, p) = 0.0;
				}
			}
		}
		// compute weighted means of posterior arrays
		for (arma::uword j = 0; j < J; ++j) {
			post_mean.at(j) = arma::dot(mu1_mat.row(j), posterior_weights.col(j));
			post_var.at(j) = arma::dot(mu2_mat.row(j), posterior_weights.col(j));
			neg_prob.at(j) = arma::dot(neg_mat.row(j), posterior_weights.col(j));
			zero_prob.at(j) = arma::dot(zero_mat.row(j), posterior_weights.col(j));
		}
		post_var -= arma::pow(post_mean, 2.0);
		return 0;
	}


	// @return PosteriorMean J vec of posterior means
	// @return PosteriorSD J vec of posterior (marginal) standard deviations
	// @return NegativeProb J vec of posterior (marginal) probability of being negative
	// @return ZeroProb J vec of posterior (marginal) probability of being zero
	arma::vec PosteriorMean() { return post_mean; }
	arma::vec PosteriorSD() { return arma::sqrt(post_var); }
	arma::vec PosteriorCov() { return post_var; }
	arma::vec NegativeProb() { return neg_prob; }
	arma::vec ZeroProb() { return zero_prob; }

private:
	// input of J vecs
	arma::vec b_vec;
	arma::vec s_vec;
	double v;
	arma::vec U_vec;
	// output of J vecs
	arma::vec post_mean;
	arma::vec post_var;
	arma::vec neg_prob;
	arma::vec zero_prob;
};

#endif
