// Gao Wang (c) 2017-2020 wang.gao@columbia.edu
#ifndef _MASH_H
#define _MASH_H
#include <cmath>
#include <armadillo>
#include <iostream>
#ifdef _OPENMP
# include <omp.h>
#endif

using std::log;
using std::exp;
using std::sqrt;

using arma::uword;
using arma::vec;
using arma::uvec;
using arma::rowvec;
using arma::colvec;
using arma::mat;
using arma::cube;
using arma::datum;
using arma::zeros;
using arma::eye;
using arma::size;
using arma::accu;
using arma::sum;
using arma::max;
using arma::abs;
using arma::sqrt;
using arma::pow;
using arma::exp;
using arma::log;
using arma::trace;
using arma::trans;
using arma::find;
using arma::inv;
using arma::trimatu;
using arma::chol;
using arma::dot;
using arma::intersect;
using arma::find;

// CONSTANTS
// ---------
const double LOG_2PI          = log(2.0 * M_PI);
const double INV_SQRT_2PI     = 1.0 / sqrt(2.0 * M_PI);
const double LOG_INV_SQRT_2PI = log(INV_SQRT_2PI);

// INLINE FUNCTION DEFINITONS
// --------------------------
inline vec
dnorm(const vec & x,
      const vec & mu,
      const vec & sigma2,
      bool              logd = false)
{
	vec res = LOG_INV_SQRT_2PI
	                - log(sqrt(sigma2))
	                - pow(x - mu, 2.0) / (2.0 * sigma2);

	if (logd) return res;
	else return exp(res);
}

inline vec
dmvnorm_mat(const mat & x,
            const vec & mean,
            const mat & sigma,
            bool              logd     = false,
            bool              inversed = false)
{
	double xdim = static_cast<double>(x.n_rows);

	vec out(x.n_cols);
	mat rooti;

	// we have previously computed rooti
	// in R eg rooti <- backsolve(chol(sigma), diag(ncol(x)))
	if (inversed) { rooti = sigma; } else {
		try {
			rooti = trans(inv(trimatu(chol(sigma))));
		} catch (const std::runtime_error & error) {
			if (logd) out.fill(-datum::inf);
			else out.fill(0.0);
			for (uword i = 0; i < x.n_cols; ++i)
				if (accu(abs(x.col(i) - mean)) < 1e-6) out.at(i) = datum::inf;
			return out;
		}
	}
	double rootisum  = sum(log(rooti.diag()));
	double constants = -(xdim / 2.0) * LOG_2PI;

	for (unsigned i = 0; i < x.n_cols; i++) {
		vec z = rooti * (x.col(i) - mean);
		out.at(i) = constants - 0.5 * sum(z % z) + rootisum;
	}

	if (logd == false) {
		out = exp(out);
	}
	return out;
}

inline double
dmvnorm(const vec & x,
        const vec & mean,
        const mat & sigma,
        bool              logd     = false,
        bool              inversed = false)
{
	mat rooti;

	if (inversed) { rooti = sigma; } else {
		try {
			rooti = trans(inv(trimatu(chol(sigma))));
		} catch (const std::runtime_error & error) {
			double diff = accu(abs(x - mean));
			if (logd) return (diff < 1e-6) ? datum::inf : -datum::inf;
			else return (diff < 1e-6) ? datum::inf : 0.0;
		}
	}
	double rootisum  = sum(log(rooti.diag()));
	double constants = -(static_cast<double>(x.n_elem) / 2.0) * LOG_2PI;

	vec z = rooti * (x - mean);
	double out  = constants - 0.5 * sum(z % z) + rootisum;

	if (logd == false) {
		out = exp(out);
	}
	return out;
}

template <class T, class U>
inline T
pnorm(const U & x, const T & m, const T & s,
      bool logd = false, bool lower_tail = true)
{
	// see `normalCDF` function at:
	// http://en.cppreference.com/w/cpp/numeric/math/erfc

	T res = 0.5 * arma::erfc((x - m) / s * M_SQRT1_2);

	// FIXME: unlike R::pnorm(0,0,0) = 1 and R::pnorm(-1,0,0) = 0, here it generates NaN
	// I manually fix it below.
	// "s == 0" check is not good enough to ensure that res doesn't have NaN due to division by zero
	uvec nan = arma::find_nonfinite(0 / s);
	if (nan.n_elem > 0) {
		res.elem(intersect(find(x >= m), nan)).ones();
		res.elem(intersect(find(x < m), nan)).zeros();
	}

	if (!lower_tail & !logd) {
		return 1.0 - res;
	} else if (lower_tail & !logd) {
		return res;
	} else if (!lower_tail & logd) {
		return log(1.0 - res);
	} else { // (lower_tail & logd)
		return log(res);
	}
}

// a quicker way to compute diag(s) %*% V %*% diag(s)
inline mat
get_cov(const vec & s, const mat & V, const mat & L)
{
	if (L.is_empty()) {
		/* return arma::diagmat(s) * V * arma::diagmat(s); */
		return (V.each_col() % s).each_row() % s.t();
	} else {
		mat svs = (V.each_col() % s).each_row() % s.t();
		return L * svs * L.t();
	}
}

inline mat
get_cov(const vec & s, const mat & V)
{
	/* return arma::diagmat(s) * V * arma::diagmat(s); */
	return (V.each_col() % s).each_row() % s.t();
}

// @title posterior_cov
// @param Vinv R x R inverse covariance matrix for the likelihood
// @param U R x R prior covariance matrix
// @return R x R posterior covariance matrix
// @description If bhat is N(b,V) and b is N(0,U) then b|bhat N(mu1,U1). This function returns U1.
inline mat
get_posterior_cov(const mat & Vinv, const mat & U)
{
	// U %*% solve(Vinv %*% U + diag(nrow(U)))
	mat S = Vinv * U;

	S.diag() += 1.0;
	return U * S.i();
}

// @title posterior_mean
// @param bhat R vector of observations
// @param Vinv R x R inverse covariance matrix for the likelihood
// @param U1 R x R posterior covariance matrix, computed using posterior_cov
// @return R vector of posterior mean
// @description If bhat is N(b,V) and b is N(0,U) then b|bhat N(mu1,U1). This function returns mu1.
inline vec
get_posterior_mean(const vec & bhat, const mat & Vinv,
                   const mat & U1)
{
	return U1 * Vinv * bhat;
}

inline mat
get_posterior_mean_mat(const mat & bhat, const mat & Vinv,
                       const mat & U1)
{
	return U1 * Vinv * bhat;
}

// SE CLASS
// --------
class SE
{
public:
SE(){
}

~SE(){
}

void
set(const mat & sbhat, const mat & sbhat_alpha)
{
	s = sbhat;
	if (sbhat_alpha.is_empty()) s_alpha.ones(sbhat.n_rows, sbhat.n_cols);
	else s_alpha = sbhat_alpha;
}

void
set(int J, int R)
{
	s.ones(J, R);
	s_alpha.ones(J, R);
}

void
set_original(const mat & value)
{
	s_orig        = value;
	is_orig_empty = s_orig.is_empty();
}

mat
get_original() const
{
	if (is_orig_empty) return (s);
	else return (s_orig);
}

mat
get() const
{
	return (s_alpha);
}

private:
mat s;
mat s_orig;
mat s_alpha;
bool is_orig_empty;
};

// FUNCTION DECLARATIONS
// ---------------------
int
mash_compute_posterior(const mat& b_mat, const SE& s_obj,
                       const mat& v_mat, const mat& l_mat,
                       const mat& a_mat, const cube& U_cube,
                       const cube& Vinv_cube,
                       const cube& U0_cube, mat& post_mean,
                       mat& post_var, mat& neg_prob,
                       mat& zero_prob, cube& post_cov,
                       const mat& posterior_weights,
                       const int& report_type);

int
mash_compute_posterior_comcov(const mat&   b_mat,
                              const SE &         s_obj,
                              const mat &  v_mat,
                              const mat &  l_mat,
                              const mat &  a_mat,
                              const cube & U_cube,
                              const cube & Vinv_cube,
                              const cube & U0_cube,
                              mat &        post_mean,
                              mat &        post_var,
                              mat &        neg_prob,
                              mat &        zero_prob,
                              cube &       post_cov,
                              const mat &  posterior_weights,
                              const int &        report_type);

int
mvsermix_compute_posterior(const mat&  b_mat,
                           const mat & s_mat,
                           mat &       v_mat,
                           cube &      U_cube,
                           cube &      Vinv_cube,
                           cube &      U0_cube,
                           cube &      Uinv_cube_drank,
                           mat &       post_mean,
                           mat &       post_var,
                           mat &       neg_prob,
                           mat &       zero_prob,
                           cube &      post_cov,
                           vec &       prior_scalar,
                           const mat & posterior_weights,
                           const mat & posterior_variable_weights);

int
mvsermix_compute_posterior_comcov(const mat&   b_mat,
                                  const mat &  s_mat,
                                  const mat &  v_mat,
                                  const cube & U_cube,
                                  const cube & Vinv_cube,
                                  const cube & U0_cube,
                                  const cube & Uinv_cube_drank,
                                  mat &        post_mean,
                                  mat &        post_var,
                                  mat &        neg_prob,
                                  mat &        zero_prob,
                                  cube &       post_cov,
                                  vec &        prior_scalar,
                                  const mat &  posterior_weights,
                                  const mat &  posterior_variable_weights);

// POSTERIORMASH CLASS
// -------------------
// @param b_mat R by J
// @param s_mat R by J
// @param s_orig_mat R by J
// @param s_alpha_mat R by J
// @param v_mat R by R
// @param l_mat R by R for the common baseline application (@Yuxin Zou)
// @param a_mat Q by R for the common baseline application (@Yuxin Zou)
// @param U_cube list of prior covariance matrices, for each mixture component P by R by R
class PosteriorMASH
{
public:
PosteriorMASH(const mat &  b_mat,
              const mat &  s_mat,
              const mat &  s_alpha_mat,
              const mat &  s_orig_mat,
              const mat &  v_mat,
              const mat &  l_mat,
              const mat &  a_mat,
              const cube & U_cube) :
	b_mat(b_mat), v_mat(v_mat), l_mat(l_mat), a_mat(a_mat), U_cube(U_cube)
{
	int J = b_mat.n_cols, R = b_mat.n_rows;

	if (s_mat.is_empty()) s_obj.set(R, J);
	else s_obj.set(s_mat, s_alpha_mat);
	s_obj.set_original(s_orig_mat);
	if (!a_mat.is_empty()) {
		R = a_mat.n_rows;
	}
	post_mean.set_size(R, J);
	post_var.set_size(R, J);
	post_cov.set_size(R, R, J);
	neg_prob.set_size(R, J);
	zero_prob.set_size(R, J);
	post_mean.zeros();
	post_var.zeros();
	post_cov.zeros();
	neg_prob.zeros();
	zero_prob.zeros();
	#ifdef _OPENMP
	omp_set_num_threads(1);
	#endif
}

~PosteriorMASH(){
}

// @title Compute posterior matrices
// @description More detailed description of function goes here.
// @param posterior_weights P X J matrix, the posterior probabilities of each mixture component for each effect
// @param report_type an integer: 1 for posterior mean only, 2 for posterior second moment, 3 for default mash output, 4 for additionally posterior covariance
int
compute_posterior(const mat & posterior_weights, const int & report_type)
{
	return mash_compute_posterior(b_mat, s_obj, v_mat, l_mat, a_mat, U_cube,
	                              Vinv_cube, U0_cube, post_mean, post_var,
	                              neg_prob, zero_prob, post_cov,
	                              posterior_weights, report_type);
}

// @title Compute posterior matrices when covariance SVS is the same for all J conditions
// @description More detailed description of function goes here.
// @param posterior_weights P X J matrix, the posterior probabilities of each mixture component for each effect
// @param report_type an integer: 1 for posterior mean only, 2 for posterior second moment, 3 for default mash output, 4 for additionally posterior covariance
int
compute_posterior_comcov(const mat & posterior_weights, const int & report_type)
{
	return mash_compute_posterior_comcov(b_mat, s_obj, v_mat, l_mat, a_mat,
	                                     U_cube, Vinv_cube, U0_cube, post_mean,
	                                     post_var, neg_prob, zero_prob,
	                                     post_cov, posterior_weights,
	                                     report_type);
}     // compute_posterior_comcov

// initializing some optinally precomputed quantities
int
set_vinv(const cube & value)
{
	Vinv_cube = value;
	return 0;
}

int
set_U0(const cube & value)
{
	U0_cube = value;
	return 0;
}

int
set_thread(const int & value)
{
	#ifdef _OPENMP
	omp_set_num_threads(value);
	#endif
	return 0;
}

// @return PosteriorMean JxR matrix of posterior means
// @return PosteriorSD JxR matrix of posterior (marginal) standard deviations
// @return NegativeProb JxR matrix of posterior (marginal) probability of being negative
// @return ZeroProb JxR matrix of posterior (marginal) probability of being zero
mat
PosteriorMean(){
	return post_mean.t();
}

mat
PosteriorSD(){
	return sqrt(post_var).t();
}

cube
PosteriorCov(){
	return post_cov;
}

mat
NegativeProb(){
	return neg_prob.t();
}

mat
ZeroProb(){
	return zero_prob.t();
}

private:
// input
mat b_mat;
SE s_obj;
mat v_mat;
mat l_mat;
mat a_mat;
cube U_cube;
cube Vinv_cube;
cube U0_cube;
// output
// all R X J mat
mat post_mean;
mat post_var;
mat neg_prob;
mat zero_prob;
// J X R X R cube
cube post_cov;
};

// POSTERIORASH CLASS
// ------------------
// @param b_vec of J
// @param s_vec of J
// @param s_alpha_vec of J
// @param v double
// @param U_vec of P
class PosteriorASH
{
public:
PosteriorASH(const vec & b_vec,
             const vec & s_vec,
             const vec & s_alpha,
             double            v,
             const vec & U_vec) :
	b_vec(b_vec), s_vec(s_vec), v(v), U_vec(U_vec)
{
	int J = b_vec.n_elem;

	if (s_alpha.is_empty()) s_alpha_vec.ones(J);
	else s_alpha_vec = s_alpha;

	post_mean.set_size(J);
	post_var.set_size(J);
	neg_prob.set_size(J);
	zero_prob.set_size(J);
}

~PosteriorASH(){
}

// @title Compute posterior matrices
// @description univariate version of PosteriorMASH::compute_posterior(), same logic
// @param posterior_weights P X J matrix, the posterior probabilities of each mixture component for each effect
int
compute_posterior(const mat & posterior_weights)
{
	vec vinv = 1 / (s_vec % s_vec * v);
	unsigned J     = b_vec.n_elem;
	unsigned P     = U_vec.n_elem;
	vec mean(J, arma::fill::zeros);
	// J X P matrices
	mat mu1_mat(J, P, arma::fill::zeros);
	mat mu2_mat(J, P, arma::fill::zeros);
	mat zero_mat(J, P, arma::fill::zeros);
	mat neg_mat(J, P, arma::fill::zeros);

	for (uword p = 0; p < P; ++p) {
		vec U1 = U_vec.at(p) / (vinv * U_vec.at(p) + 1.0);
		mu1_mat.col(p) = U1 % vinv % b_vec % s_alpha_vec;
		U1 = U1 % (s_alpha_vec % s_alpha_vec);
		mu2_mat.col(p) = pow(mu1_mat.col(p), 2.0) + U1;
		vec sigma = sqrt(U1);
		neg_mat.col(p) = pnorm(mu1_mat.col(p), mean, sigma);
		for (uword j = 0; j < J; ++j) {
			if (U1.at(j) == 0) {
				zero_mat.at(j, p) = 1.0;
				neg_mat.at(j, p)  = 0.0;
			}
		}
	}
	// compute weighted means of posterior arrays
	for (uword j = 0; j < J; ++j) {
		post_mean.at(j) = dot(mu1_mat.row(j), posterior_weights.col(j));
		post_var.at(j)  = dot(mu2_mat.row(j), posterior_weights.col(j));
		neg_prob.at(j)  = dot(neg_mat.row(j), posterior_weights.col(j));
		zero_prob.at(j) = dot(zero_mat.row(j), posterior_weights.col(j));
	}
	post_var -= pow(post_mean, 2.0);
	return 0;
}     // compute_posterior

// @return PosteriorMean J vec of posterior means
// @return PosteriorSD J vec of posterior (marginal) standard deviations
// @return NegativeProb J vec of posterior (marginal) probability of being negative
// @return ZeroProb J vec of posterior (marginal) probability of being zero
vec
PosteriorMean(){
	return post_mean;
}

vec
PosteriorSD(){
	return sqrt(post_var);
}

vec
PosteriorCov(){
	return post_var;
}

vec
NegativeProb(){
	return neg_prob;
}

vec
ZeroProb(){
	return zero_prob;
}

private:
// input of J vecs
vec b_vec;
vec s_vec;
vec s_alpha_vec;
double v;
vec U_vec;
// output of J vecs
vec post_mean;
vec post_var;
vec neg_prob;
vec zero_prob;
};

// MVSERMIX CLASS
// --------------
// @title Inferences for Multivariate Single Effect Regression with Mixture prior
// @param b_mat R by J
// @param s_mat R by J
// @param v_mat R by R
// @param U_cube list of prior covariance matrices, for each mixture component P by R by R
class MVSERMix
{
public:
MVSERMix(const mat &  b_mat,
         const mat &  s_mat,
         const mat &  v_mat,
         const cube & U_cube) :
	b_mat(b_mat), s_mat(s_mat), v_mat(v_mat), U_cube(U_cube)
{
	int J = b_mat.n_cols, R = b_mat.n_rows;

	post_mean.set_size(R, J);
	post_var.set_size(R, J);
	post_cov.set_size(R, R, J);
	neg_prob.set_size(R, J);
	zero_prob.set_size(R, J);
	post_mean.zeros();
	post_var.zeros();
	post_cov.zeros();
	neg_prob.zeros();
	zero_prob.zeros();
	prior_scalar.set_size(U_cube.n_slices);
	#ifdef _OPENMP
	omp_set_num_threads(1);
	#endif
}

~MVSERMix(){
}

// @title Compute posterior matrices and EM updates for prior scalar estimate
// @description Make posterior inferences, and also perform the EM update for prior scalar, for mvSuSiE model.
// @param posterior_weights P X J matrix, the posterior probabilities of each mixture component for each effect.
// @param posterior_variable_weights P X J matrix, the posterior inclusion probabilities of each effect in a single-effect model.
// posterior_variable_weights is only relevant when EM updates for prior scalar is needed.
int
compute_posterior(const mat & posterior_weights,
                  const mat & posterior_variable_weights)
{
	return mvsermix_compute_posterior(b_mat, s_mat, v_mat, U_cube, Vinv_cube,
	                                  U0_cube, Uinv_cube_drank, post_mean, post_var,
	                                  neg_prob, zero_prob, post_cov,
	                                  prior_scalar,
	                                  posterior_weights,
	                                  posterior_variable_weights);
}     // compute_posterior

// @title Compute posterior matrices when covariance SVS is the same for all J conditions
// @description More detailed description of function goes here.
// @param posterior_weights P X J matrix, the posterior probabilities of each mixture component for each effect
int
compute_posterior_comcov(const mat & posterior_weights,
                         const mat & posterior_variable_weights)
{
	return mvsermix_compute_posterior_comcov(b_mat, s_mat, v_mat, U_cube,
	                                         Vinv_cube, U0_cube, Uinv_cube_drank,
	                                         post_mean, post_var, neg_prob,
	                                         zero_prob, post_cov, prior_scalar,
	                                         posterior_weights,
	                                         posterior_variable_weights);
}     // compute_posterior_comcov

// initializing some optinally precomputed quantities
int
set_Vinv(const cube & value)
{
	Vinv_cube = value;
	return 0;
}

int
set_U0(const cube & value)
{
	U0_cube = value;
	return 0;
}

int
set_Uinv(const cube & value)
{
	Uinv_cube_drank        = value;
	return 0;
}

int
set_thread(const int & value)
{
	#ifdef _OPENMP
	omp_set_num_threads(value);
	#endif
	return 0;
}

// @return PosteriorMean JxR matrix of posterior means
// @return PosteriorSD JxR matrix of posterior (marginal) standard deviations
// @return NegativeProb JxR matrix of posterior (marginal) probability of being negative
// @return ZeroProb JxR matrix of posterior (marginal) probability of being zero
mat
PosteriorMean(){
	return post_mean.t();
}

mat
PosteriorSD(){
	return sqrt(post_var).t();
}

cube
PosteriorCov(){
	return post_cov;
}

mat
NegativeProb(){
	return neg_prob.t();
}

mat
ZeroProb(){
	return zero_prob.t();
}

vec
PriorScalar(){
	return prior_scalar;
}

private:
// input
mat b_mat;
mat s_mat;
mat v_mat;
cube U_cube;
cube Vinv_cube;
cube U0_cube;
cube Uinv_cube_drank;
// output
// all R X J mat
mat post_mean;
mat post_var;
mat neg_prob;
mat zero_prob;
// J X R X R cube
cube post_cov;
// P vector of scalars
vec prior_scalar;
};

// Softmax functions: yi = exp(xi) / sum(exp(xj))
inline vec
softmax(const vec & x)
{
	// Calculate exp()
	// Subtract the max - this prevents overflow, which happens for x ~ 1000
	vec y = exp(x - max(x));
	// Renormalise
	y /= sum(y);
	return y;
}

// function for "shrinking" the covariance matrix, to get $\hat U_k$.
inline mat
shrink_cov(const mat & V, const double & eps)
{
	vec eigval;
	mat eigvec;
	eig_sym(eigval, eigvec, V);
	for (uword i = 0; i < eigval.n_elem; ++i) {
		eigval(i) = (eigval(i) > 1.0) ? eigval(i) : (1.0 + eps);
	}
	return eigvec * diagmat(eigval) * trans(eigvec);
}

// TEEM CLASS
// ----------
// @title Truncated Eigenvalue Extreme deconvolution
// @description ...
// @param X
// @param w
// @param U
// @param maxiter
// @param tol
// @param verbose
class TEEM
{
public:
TEEM(const mat &  X_mat,
     const vec &  w_vec,
     const cube & U_cube) :
	X_mat(X_mat), w_vec(w_vec)
{
	T_cube = U_cube;
	for (unsigned j = 0; j < T_cube.n_slices; ++j) {
		T_cube.slice(j) += eye(size(T_cube.slice(j)));
	}
}

~TEEM(){
}

vec
get_objective() const {
	return objective;
}

vec
get_maxd() const {
	return maxd;
}

vec
get_w() const {
	return w_vec;
}

cube
get_U() const
{
	cube U_cube = T_cube;
	for (unsigned j = 0; j < U_cube.n_slices; ++j) {
		U_cube.slice(j) -= eye(size(U_cube.slice(j)));
	}
	return U_cube;
}

int
fit(const int & maxiter, const double & converge_tol, const double & eigen_tol, const bool & verbose)
{
	// initialize to store progress
	objective.zeros(maxiter);
	maxd.zeros(maxiter);
	int iter_out = 0;

	// Get the number of samples (n) and the number of mixture components (k)
	unsigned int n = X_mat.n_rows;
	unsigned int k = w_vec.size();

	for (unsigned int iter = 0; iter < (unsigned int) maxiter; ++iter) {
		// store parameters and likelihood in the previous step
		vec w0_vec = w_vec;

		// E-step: calculate posterior probabilities using the current mu and sigmas
		mat logP = zeros<mat>(n, k); // n by k matrix
		for (unsigned j = 0; j < k; ++j) {
			logP.col(j) = log(w_vec(j)) + dmvnorm_mat(trans(X_mat), zeros<vec>(
									  X_mat.n_cols), T_cube.slice(j), true); // ??
		}
		// softmax for renormalization
		mat P_mat = zeros<mat>(k, n); // k by n matrix. because of row/col vec converting

		for (uword i = 0; i < n; ++i) {
			colvec y = arma::conv_to<colvec>::from(logP.row(i));
			P_mat.col(i) = softmax(y);
		}
		P_mat = trans(P_mat); // n by k matrix

		// M-step:
		for (unsigned int j = 0; j < k; ++j) {
			T_cube.slice(j) = trans(X_mat) * (P_mat.col(j) % X_mat.each_col()) / accu(P_mat.col(j));
			T_cube.slice(j) = shrink_cov(T_cube.slice(j), eigen_tol);
		}
		// update mixture weights
		w_vec = arma::conv_to<colvec>::from(sum(P_mat, 0)) / n; // 0:sum by column;

		// Compute log-likelihood at the current estimates
		double f = compute_loglik();

		// Check stopping criterion
		double d = max(abs(w_vec - w0_vec));
		maxd(iter)      = d;
		objective(iter) = f;
		iter_out        = iter;

		if (d < converge_tol) {
			break;
		}
	}
	objective.resize(iter_out + 1);
	maxd.resize(iter_out + 1);
	return 0;
}     // fit

private:
mat X_mat;
vec w_vec;
cube T_cube;
vec objective;
vec maxd;
double
compute_loglik()
{
	unsigned int n = X_mat.n_rows;
	unsigned int k = w_vec.size();

	vec y = zeros<vec>(n);
	for (unsigned int j = 0; j < k; ++j) {
		y = y + w_vec(j) * dmvnorm_mat(trans(X_mat), zeros<vec>(X_mat.n_cols), T_cube.slice(j));
	}
	return (sum(log(y)));
}
};

// FUNCTION DEFINITIONS
// --------------------
// @title calc_lik
// @description computes matrix of likelihoods for each of J cols of Bhat for each of P prior covariances
// @param b_mat R by J
// @param s_mat R by J
// @param v_mat R by R
// @param l_mat R by R for the common baseline application (@Yuxin Zou)
// @param U_cube list of prior covariance matrices
// @param sigma_cube list of sigma which is result of get_cov(s_mat, v_mat, l_mat)
// @param logd if true computes log-likelihood
// @param common_cov if true use version for common covariance
// @return J x P matrix of multivariate normal likelihoods, p(bhat | U[p], V)
mat
calc_lik(const mat &  b_mat,
         const mat &  s_mat,
         const mat &  v_mat,
         const mat &  l_mat,
         const cube & U_cube,
         const cube & sigma_cube,
         bool               logd,
         bool               common_cov,
         int                n_thread = 1)
{
	// In armadillo data are stored with column-major ordering
	// slicing columns are therefore faster than rows
	// lik is a J by P matrix
	mat lik(b_mat.n_cols, U_cube.n_slices, arma::fill::zeros);
	vec mean(b_mat.n_rows, arma::fill::zeros);
	mat sigma;
    #ifdef _OPENMP
	omp_set_num_threads(n_thread);
    #endif
	if (common_cov) {
		if (!sigma_cube.is_empty()) sigma = sigma_cube.slice(0);
		else sigma = get_cov(s_mat.col(0), v_mat, l_mat);
	#pragma omp parallel for default(none) schedule(static) shared(lik, U_cube, mean, sigma, logd, b_mat)
		for (uword p = 0; p < lik.n_cols; ++p) {
			lik.col(p) = dmvnorm_mat(b_mat, mean, sigma + U_cube.slice(p), logd);
		}
	} else {
	#pragma \
		omp parallel for default(none) schedule(static) shared(lik, mean, logd, U_cube, b_mat, sigma_cube, l_mat, v_mat, s_mat) private(sigma)
		for (uword j = 0; j < lik.n_rows; ++j) {
			if (!sigma_cube.is_empty()) sigma = sigma_cube.slice(j);
			else sigma = get_cov(s_mat.col(j), v_mat, l_mat);
			for (uword p = 0; p < lik.n_cols; ++p) {
				lik.at(j, p) = dmvnorm(b_mat.col(j), mean, sigma + U_cube.slice(p), logd);
			}
		}
	}
	return lik;
}

// @title calc_lik multivariate common cov version with sigma inverse precomputed
// @description computes matrix of likelihoods for each of J cols of Bhat for each of P prior covariances
// @param b_mat R by J
// @param rooti_cube R by R by P, or R by R by J by P, if common_cov is False
// @param logd if true computes log-likelihood
// @param common_cov if true use version for common covariance
// @return J x P matrix of multivariate normal likelihoods, p(bhat | U[p], V)
mat
calc_lik(const mat & b_mat,
         const cube & rooti_cube,
         bool logd, bool common_cov, int n_thread = 1)
{
    #ifdef _OPENMP
	omp_set_num_threads(n_thread);
    #endif
	// In armadillo data are stored with column-major ordering
	// slicing columns are therefore faster than rows
	// lik is a J by P matrix
	int P;

	if (common_cov) P = rooti_cube.n_slices;
	else P = rooti_cube.n_slices / b_mat.n_cols;
	mat lik(b_mat.n_cols, P, arma::fill::zeros);
	vec mean(b_mat.n_rows, arma::fill::zeros);
	if (common_cov) {
	#pragma omp parallel for default(none) schedule(static) shared(lik, mean, logd, rooti_cube, b_mat)
		for (uword p = 0; p < lik.n_cols; ++p) {
			lik.col(p) = dmvnorm_mat(b_mat, mean, rooti_cube.slice(p), logd, true);
		}
	} else {
	#pragma omp parallel for default(none) schedule(static) shared(lik, mean, logd, rooti_cube, b_mat)
		for (uword j = 0; j < lik.n_rows; ++j) {
			for (uword p = 0; p < lik.n_cols; ++p) {
				lik.at(j, p) = dmvnorm(b_mat.col(j), mean, rooti_cube.slice(j * lik.n_cols + p), logd, true);
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
mat
calc_lik(const vec & b_vec,
         const vec & s_vec,
         double            v,
         const vec & U_vec,
         bool              logd)
{
	mat lik(b_vec.n_elem, U_vec.n_elem, arma::fill::zeros);
	vec sigma = s_vec % s_vec * v;
	vec mean(b_vec.n_elem, arma::fill::zeros);

	for (uword p = 0; p < lik.n_cols; ++p) {
		lik.col(p) = dnorm(b_vec, mean, sigma + U_vec.at(p), logd);
	}
	return lik;
}

// This implements the core part of the compute_posterior method in
// the PosteriorMASH class.
int
mash_compute_posterior(const mat& b_mat, const SE& s_obj,
                       const mat& v_mat, const mat& l_mat,
                       const mat& a_mat, const cube& U_cube,
                       const cube& Vinv_cube,
                       const cube& U0_cube, mat& post_mean,
                       mat& post_var, mat& neg_prob,
                       mat& zero_prob, cube& post_cov,
                       const mat& posterior_weights,
                       const int& report_type)
{
	vec mean(post_mean.n_rows);
	mean.fill(0);
	
    #pragma \
	omp parallel for schedule(static) default(none) shared(posterior_weights, report_type, mean, post_mean, post_var, neg_prob, zero_prob, post_cov, b_mat, s_obj, l_mat, v_mat, a_mat, U_cube, Vinv_cube, U0_cube)
	for (uword j = 0; j < post_mean.n_cols; ++j) {
		// FIXME: improved math may help here
		mat Vinv_j;
		if (Vinv_cube.is_empty())
			Vinv_j = inv_sympd(get_cov(s_obj.get_original().col(j),
			                                 v_mat, l_mat));
		else
			Vinv_j = Vinv_cube.slice(j);

		// R X P matrices
		mat mu1_mat(post_mean.n_rows, U_cube.n_slices);
		mat diag_mu2_mat(post_mean.n_rows, U_cube.n_slices);
		mat zero_mat(post_mean.n_rows, U_cube.n_slices);
		mat neg_mat(post_mean.n_rows, U_cube.n_slices);

		mu1_mat.fill(0);
		diag_mu2_mat.fill(0);
		zero_mat.fill(0);
		
		for (uword p = 0; p < U_cube.n_slices; ++p) {
			//
			mat U1(post_mean.n_rows, post_mean.n_rows);
			mat U0;
			U1.fill(0);
			
			if (U0_cube.is_empty())
				U0 = get_posterior_cov(Vinv_j, U_cube.slice(p));
			else
				U0 = U0_cube.slice(j * U_cube.n_slices + p);
			if (a_mat.is_empty()) {
				mu1_mat.col(p) = get_posterior_mean(b_mat.col(j), Vinv_j, U0) % s_obj.get().col(j);
				U1 = (U0.each_col() % s_obj.get().col(j)).each_row() % s_obj.get().col(j).t();
			} else {
				mu1_mat.col(p) = a_mat * (get_posterior_mean(b_mat.col(j), Vinv_j, U0) % s_obj.get().col(j));
				U1 = a_mat * (((U0.each_col() % s_obj.get().col(j)).each_row() % s_obj.get().col(j).t()) * a_mat.t());
			}

			if (report_type == 2 || report_type == 4) {
				post_cov.slice(j) += posterior_weights.at(p, j) * (U1 + mu1_mat.col(p) * mu1_mat.col(p).t());
			}

			vec sigma = sqrt(U1.diag()); // U1.diag() is the posterior covariance
			diag_mu2_mat.col(p) = pow(mu1_mat.col(p), 2.0) + U1.diag();
			neg_mat.col(p)      = pnorm(mu1_mat.col(p), mean, sigma);
			for (uword r = 0; r < sigma.n_elem; ++r) {
				if (sigma.at(r) == 0) {
					zero_mat.at(r, p) = 1.0;
					neg_mat.at(r, p)  = 0.0;
				}
			}
		}

		// compute weighted means of posterior arrays
		post_mean.col(j) = mu1_mat * posterior_weights.col(j);
		post_var.col(j)  = diag_mu2_mat * posterior_weights.col(j);
		neg_prob.col(j)  = neg_mat * posterior_weights.col(j);
		zero_prob.col(j) = zero_mat * posterior_weights.col(j);
		//
		if (report_type == 4)
			post_cov.slice(j) -= post_mean.col(j) * post_mean.col(j).t();
	}
	post_var -= pow(post_mean, 2.0);

	return 0;
} // mash_compute_posterior

// This implements the core part of the compute_posterior_comcov method in
// the PosteriorMASH class.
int
mash_compute_posterior_comcov(const mat&   b_mat,
                              const SE &         s_obj,
                              const mat &  v_mat,
                              const mat &  l_mat,
                              const mat &  a_mat,
                              const cube & U_cube,
                              const cube & Vinv_cube,
                              const cube & U0_cube,
                              mat &        post_mean,
                              mat &        post_var,
                              mat &        neg_prob,
                              mat &        zero_prob,
                              cube &       post_cov,
                              const mat &  posterior_weights,
                              const int &        report_type)
{
	mat mean(post_mean.n_rows, post_mean.n_cols);
	mean.fill(0);
	
	// R X R
	mat Vinv;

	if (Vinv_cube.is_empty())
		Vinv =
			inv_sympd(get_cov(s_obj.get_original().col(0), v_mat, l_mat));
	else Vinv = Vinv_cube.slice(0);

	rowvec ones(post_mean.n_cols);
	rowvec zeros(post_mean.n_cols);
	ones.fill(1);
	zeros.fill(0);
	
    #pragma \
	omp parallel for schedule(static) default(none) shared(posterior_weights, report_type, mean, Vinv, ones, zeros, post_mean, post_var, neg_prob, zero_prob, post_cov, b_mat, s_obj, a_mat, U_cube, U0_cube)
	for (uword p = 0; p < U_cube.n_slices; ++p) {
		mat zero_mat(post_mean.n_rows, post_mean.n_cols);
		// R X R
		mat U1(post_mean.n_rows, post_mean.n_rows);
		// R X J
		mat mu1_mat(post_mean.n_rows, post_mean.n_cols);
		mat U0;
		zero_mat.fill(0);
		U1.fill(0);
		mu1_mat.fill(0);
		
		if (U0_cube.is_empty()) U0 = get_posterior_cov(Vinv, U_cube.slice(p));
		else U0 = U0_cube.slice(p);
		if (a_mat.is_empty()) {
			mu1_mat = get_posterior_mean_mat(b_mat, Vinv, U0) % s_obj.get();
			U1      = (U0.each_col() % s_obj.get().col(0)).each_row() % s_obj.get().col(0).t();
		} else {
			mu1_mat = a_mat * (get_posterior_mean_mat(b_mat, Vinv, U0) % s_obj.get());
			U1      = a_mat
			          * (((U0.each_col() % s_obj.get().col(0)).each_row() % s_obj.get().col(0).t())
			             * a_mat.t());
		}
		// R X J
		mat diag_mu2_mat = pow(mu1_mat, 2.0);
		diag_mu2_mat.each_col() += U1.diag();
		// R X J
		// FIXME: any better way to init sigma?
		mat sigma(post_mean.n_rows, post_mean.n_cols);
		sigma.fill(0);
		
		sigma.each_col() += sqrt(U1.diag()); // U1.diag() is the posterior covariance
		mat neg_mat = pnorm(mu1_mat, mean, sigma);
		for (uword r = 0; r < sigma.n_rows; ++r) {
			if (sigma.at(r, 0) == 0) {
				zero_mat.row(r) = ones;
				neg_mat.row(r)  = zeros;
			}
		}
		// compute weighted means of posterior arrays
	#pragma omp critical
		{
			post_mean += mu1_mat.each_row() % posterior_weights.row(p);
			post_var  += diag_mu2_mat.each_row() % posterior_weights.row(p);
			neg_prob  += neg_mat.each_row() % posterior_weights.row(p);
			zero_prob += zero_mat.each_row() % posterior_weights.row(p);
			if (report_type == 2 || report_type == 4) {
				for (uword j = 0; j < post_mean.n_cols; ++j) {
					post_cov.slice(j) +=
						posterior_weights.at(p, j) * (U1 + mu1_mat.col(j) * mu1_mat.col(j).t());
				}
			}
		}
	}
	post_var -= pow(post_mean, 2.0);
	//
	if (report_type == 4) {
	#pragma omp parallel for schedule(static) default(none) shared(post_cov, post_mean)
		for (uword j = 0; j < post_mean.n_cols; ++j) {
			post_cov.slice(j) -= post_mean.col(j) * post_mean.col(j).t();
		}
	}
	return 0;
} // mash_compute_posterior_comcov

// This implements the core part of the compute_posterior method in
// the MVSERMix class.
int
mvsermix_compute_posterior(const mat&  b_mat,
                           const mat & s_mat,
                           mat &       v_mat,
                           cube &      U_cube,
                           cube &      Vinv_cube,
                           cube &      U0_cube,
                           cube &      Uinv_cube_drank,
                           mat &       post_mean,
                           mat &       post_var,
                           mat &       neg_prob,
                           mat &       zero_prob,
                           cube &      post_cov,
                           vec &       prior_scalar,
                           const mat & posterior_weights,
                           const mat & posterior_variable_weights)
{
	vec mean(post_mean.n_rows);
	mean.fill(0);
	
	// This is meant to store a length P of 2nd moment matrices,
	// each element is \sum_j posterior_{p,j} * mu2_{p,j}
	cube Eb2_cube;
	bool to_estimate_prior = !posterior_variable_weights.is_empty();
	if (to_estimate_prior) {
		// we will compute the EM update for prior scalar here
		// for use with mmbr package
		Eb2_cube.set_size(post_mean.n_rows, post_mean.n_rows, U_cube.n_slices);
		Eb2_cube.zeros();
	}
    #pragma \
	omp parallel for schedule(static) default(none) shared(posterior_weights, posterior_variable_weights, to_estimate_prior, mean, Eb2_cube, post_mean, post_var, neg_prob, zero_prob, post_cov, prior_scalar, b_mat, s_mat, v_mat, U_cube, Vinv_cube, U0_cube, Uinv_cube_drank)
	for (uword j = 0; j < post_mean.n_cols; ++j) {
		// FIXME: improved math may help here
		mat Vinv_j;
		if (Vinv_cube.is_empty())
			Vinv_j = inv_sympd(get_cov(s_mat.col(j), v_mat));
		else Vinv_j = Vinv_cube.slice(j);
		// R X P matrices
		mat mu1_mat(post_mean.n_rows, U_cube.n_slices);
		mat diag_mu2_mat(post_mean.n_rows, U_cube.n_slices);
		mat zero_mat(post_mean.n_rows, U_cube.n_slices);
		mat neg_mat(post_mean.n_rows, U_cube.n_slices);
		mu1_mat.fill(0);
		diag_mu2_mat.fill(0);
		zero_mat.fill(0);
		neg_mat.fill(0);
		
		// R X R X P
		cube mu2_cube;
		mu2_cube.set_size(post_mean.n_rows, post_mean.n_rows, U_cube.n_slices);
		for (uword p = 0; p < U_cube.n_slices; ++p) {
			mat U1;
			if (U0_cube.is_empty()) U1 = get_posterior_cov(Vinv_j, U_cube.slice(p));
			else U1 = U0_cube.slice(j * U_cube.n_slices + p);
			mu1_mat.col(p) = get_posterior_mean(b_mat.col(j), Vinv_j, U1);
			// this is posterior 2nd moment for the j-th variable and the p-th prior
			mu2_cube.slice(p) = U1 + mu1_mat.col(p) * mu1_mat.col(p).t();
			// add to posterior 2nd moment contribution of the p-th component
			post_cov.slice(j) += posterior_weights.at(p, j) * mu2_cube.slice(p);
			vec sigma = sqrt(U1.diag()); // U1.diag() is the posterior covariance
			diag_mu2_mat.col(p) = pow(mu1_mat.col(p), 2.0) + U1.diag();
			neg_mat.col(p)      = pnorm(mu1_mat.col(p), mean, sigma);
			for (uword r = 0; r < sigma.n_elem; ++r) {
				if (sigma.at(r) == 0) {
					zero_mat.at(r, p) = 1.0;
					neg_mat.at(r, p)  = 0.0;
				}
			}
		}
		// compute weighted means of posterior arrays
		post_mean.col(j)   = mu1_mat * posterior_weights.col(j);
		post_var.col(j)    = diag_mu2_mat * posterior_weights.col(j);
		neg_prob.col(j)    = neg_mat * posterior_weights.col(j);
		zero_prob.col(j)   = zero_mat * posterior_weights.col(j);
		post_cov.slice(j) -= post_mean.col(j) * post_mean.col(j).t();
		if (to_estimate_prior) {
	    #pragma omp critical
			{
				for (uword p = 0; p < U_cube.n_slices; ++p) {
					// we will compute some quantity to provide for
					// EM update for prior scalar in mmbr package
					// the M-step update is:
					// \sigma_0^2 = \sum_{p=1}^P p(\gamma_p) \mathrm{tr}(U_p^{-1} E[bb^T \,|\, \gamma_p])/r
					// where E[bb^T \,|\, \gamma_p] = \sum_j \alpha_{p,j} * mu2_mat_{p,j}
					Eb2_cube.slice(p) += posterior_variable_weights.at(p, j) * mu2_cube.slice(p);
				}
			}
		}
	}
	post_var -= pow(post_mean, 2.0);
	if (to_estimate_prior) {
		// now compute \mathrm{tr}(U_p^{-1} E[bb^T \,|\, \gamma_p])/r for each p
		for (uword p = 0; p < U_cube.n_slices; ++p) {
			prior_scalar.at(p) = trace(Uinv_cube_drank.slice(p) * Eb2_cube.slice(p));
		}
	}
	return 0;
} // mvsermix_compute_posterior

// This implements the core part of the compute_posterior_comcov method in
// the MVSERMix class.
int
mvsermix_compute_posterior_comcov(const mat&   b_mat,
                                  const mat &  s_mat,
                                  const mat &  v_mat,
                                  const cube & U_cube,
                                  const cube & Vinv_cube,
                                  const cube & U0_cube,
                                  const cube & Uinv_cube_drank,
                                  mat &        post_mean,
                                  mat &        post_var,
                                  mat &        neg_prob,
                                  mat &        zero_prob,
                                  cube &       post_cov,
                                  vec &        prior_scalar,
                                  const mat &  posterior_weights,
                                  const mat &  posterior_variable_weights)
{
	mat mean(post_mean.n_rows, post_mean.n_cols);
	mean.fill(0);
	
	// for Eb2_cube see compute_posterior() for detailed documentations.
	cube Eb2_cube;
	bool to_estimate_prior = !posterior_variable_weights.is_empty();
	if (to_estimate_prior) {
		Eb2_cube.set_size(post_mean.n_rows, post_mean.n_rows, U_cube.n_slices);
		Eb2_cube.zeros();
	}
	// R X R
	mat Vinv;
	if (Vinv_cube.is_empty()) Vinv =
			inv_sympd(get_cov(s_mat.col(0), v_mat));
	else Vinv = Vinv_cube.slice(0);

	rowvec ones(post_mean.n_cols);
	rowvec zeros(post_mean.n_cols);
	ones.fill(1);
	zeros.fill(0);
	
    #pragma \
	omp parallel for schedule(static) default(none) shared(posterior_weights, posterior_variable_weights, to_estimate_prior, mean, Vinv, zeros, ones, Eb2_cube, post_mean, post_var, neg_prob, zero_prob, post_cov, prior_scalar, b_mat, U_cube, U0_cube, Uinv_cube_drank)
	for (uword p = 0; p < U_cube.n_slices; ++p) {
		mat zero_mat(post_mean.n_rows, post_mean.n_cols);
		// R X R
		mat U1;
		// R X J
		mat mu1_mat;
		zero_mat.fill(0);
		
		if (U0_cube.is_empty()) U1 = get_posterior_cov(Vinv, U_cube.slice(p));
		else U1 = U0_cube.slice(p);
		mu1_mat = get_posterior_mean_mat(b_mat, Vinv, U1);
		cube mu2_cube;
		mu2_cube.set_size(post_mean.n_rows, post_mean.n_rows, post_mean.n_cols);
		for (uword j = 0; j < post_mean.n_cols; ++j) {
			mu2_cube.slice(j) = U1 + mu1_mat.col(j) * mu1_mat.col(j).t();
			if (to_estimate_prior) Eb2_cube.slice(p) += posterior_variable_weights.at(p, j) * mu2_cube.slice(j);
		}
		if (to_estimate_prior) prior_scalar.at(p) = trace(Uinv_cube_drank.slice(p) * Eb2_cube.slice(p));
		// R X J
		mat diag_mu2_mat = pow(mu1_mat, 2.0);
		diag_mu2_mat.each_col() += U1.diag();
		// R X J
		// FIXME: any better way to init sigma?
		mat sigma(post_mean.n_rows, post_mean.n_cols);
		sigma.fill(0);
		sigma.each_col() += sqrt(U1.diag()); // U1.diag() is the posterior covariance
		mat neg_mat = pnorm(mu1_mat, mean, sigma);
		for (uword r = 0; r < sigma.n_rows; ++r) {
			if (sigma.at(r, 0) == 0) {
				zero_mat.row(r) = ones;
				neg_mat.row(r)  = zeros;
			}
		}
	#pragma omp critical
		{
			// compute weighted means of posterior arrays
			post_mean += mu1_mat.each_row() % posterior_weights.row(p);
			post_var  += diag_mu2_mat.each_row() % posterior_weights.row(p);
			neg_prob  += neg_mat.each_row() % posterior_weights.row(p);
			zero_prob += zero_mat.each_row() % posterior_weights.row(p);
			for (uword j = 0; j < post_mean.n_cols; ++j) {
				post_cov.slice(j) += posterior_weights.at(p, j) * mu2_cube.slice(j);
			}
		}
	}
	post_var -= pow(post_mean, 2.0);
    #pragma omp parallel for schedule(static) default(none) shared(post_cov, post_mean)
	for (uword j = 0; j < post_mean.n_cols; ++j) {
		post_cov.slice(j) -= post_mean.col(j) * post_mean.col(j).t();
	}
	return 0;
} // mvsermix_compute_posterior_comcov

#endif // ifndef _MASH_H
