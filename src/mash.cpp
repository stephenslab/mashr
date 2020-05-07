// Wrapper to various C++ functions/objects for inference in MASH
// Gao Wang (c) 2017-2019 wang.gao@columbia.edu
#include <iostream>
#include <stdexcept>
#include "RcppArmadillo.h"
#include "mash.h"
// [[Rcpp::depends(RcppArmadillo)]]

// Enable C++11.
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
Rcpp::List
inv_chol_tri_rcpp(Rcpp::NumericMatrix x_mat)
{
    arma::mat res = arma::trans(arma::inv(arma::trimatu(arma::chol(Rcpp::as<arma::mat>(x_mat)))));

    return Rcpp::List::create(Rcpp::Named("data") = res,
             Rcpp::Named("status") = 0);
}

// [[Rcpp::export]]
Rcpp::List
calc_lik_rcpp(Rcpp::NumericMatrix b_mat,
  Rcpp::NumericMatrix             s_mat,
  Rcpp::NumericMatrix             v_mat,
  Rcpp::NumericMatrix             l_mat,
  Rcpp::NumericVector             U_3d,
  Rcpp::NumericVector             sigma_3d,
  bool                            logd,
  bool                            common_cov)
{
    // hide armadillo warning / error messages
    // std::ostream nullstream(0);
    // arma::set_stream_err2(nullstream);
    arma::mat res;
    if (!Rf_isNull(U_3d.attr("dim"))) {
        // matrix version
        // set cube data from R 3D array
        Rcpp::IntegerVector dimU = U_3d.attr("dim");
        arma::cube U_cube(U_3d.begin(), dimU[0], dimU[1], dimU[2]);
        arma::cube sigma_cube;
        if (!Rf_isNull(sigma_3d.attr("dim"))) {
          Rcpp::IntegerVector dimSigma = sigma_3d.attr("dim");
          arma::cube tmp_cube(sigma_3d.begin(), dimSigma[0], dimSigma[1], dimSigma[2]);
          sigma_cube = tmp_cube;
        }
        res = calc_lik(Rcpp::as<arma::mat>(b_mat),
            Rcpp::as<arma::mat>(s_mat),
            Rcpp::as<arma::mat>(v_mat),
            Rcpp::as<arma::mat>(l_mat),
            U_cube,
            sigma_cube,
            logd,
            common_cov);
    } else {
        // vector version
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
Rcpp::List
calc_lik_precomputed_rcpp(Rcpp::NumericMatrix b_mat,
  Rcpp::NumericVector                         rooti_3d,
  bool                                        logd,
  bool                                        common_cov)
{
    // hide armadillo warning / error messages
    // std::ostream nullstream(0);
    // arma::set_stream_err2(nullstream);
    arma::mat res;
    // set cube data from R 3D array
    Rcpp::IntegerVector dimR = rooti_3d.attr("dim");
    arma::cube rooti_cube(rooti_3d.begin(), dimR[0], dimR[1], dimR[2]);

    res = calc_lik(Rcpp::as<arma::mat>(b_mat),
        rooti_cube,
        logd, common_cov);
    return Rcpp::List::create(Rcpp::Named("data") = res,
             Rcpp::Named("status") = 0);
}

// [[Rcpp::export]]
Rcpp::List
calc_post_rcpp(Rcpp::NumericMatrix b_mat,
  Rcpp::NumericMatrix              s_mat,
  Rcpp::NumericMatrix              s_alpha_mat,
  Rcpp::NumericMatrix              s_orig_mat,
  Rcpp::NumericMatrix              v_mat,
  Rcpp::NumericMatrix              l_mat,
  Rcpp::NumericMatrix              a_mat,
  Rcpp::NumericVector              U_3d,
  Rcpp::NumericMatrix              posterior_weights,
  bool                             common_cov,
  int                              report_type)
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
            Rcpp::Named("post_sd")   = pc.PosteriorSD(),
            Rcpp::Named("post_cov")  = pc.PosteriorCov(),
            Rcpp::Named("post_zero") = pc.ZeroProb(),
            Rcpp::Named("post_neg")  = pc.NegativeProb());
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
            Rcpp::Named("post_cov")  = pc.PosteriorCov(),
            Rcpp::Named("post_sd")   = pc.PosteriorSD(),
            Rcpp::Named("post_zero") = pc.ZeroProb(),
            Rcpp::Named("post_neg")  = pc.NegativeProb());
    }
} // calc_post_rcpp

// [[Rcpp::export]]
Rcpp::List
calc_sermix_rcpp(Rcpp::NumericMatrix b_mat,
  Rcpp::NumericMatrix                s_mat,
  Rcpp::NumericMatrix                v_mat,
  Rcpp::NumericVector                vinv_3d,
  Rcpp::NumericVector                U_3d,
  Rcpp::NumericVector                Uinv_3d,
  Rcpp::NumericVector                U0_3d,
  Rcpp::NumericMatrix                posterior_mixture_weights,
  Rcpp::NumericMatrix                posterior_variable_weights,
  double                             sigma0,
  bool                               common_cov)
{
    // hide armadillo warning / error messages
    // std::ostream nullstream(0);
    // arma::set_stream_err2(nullstream);
    if (Rf_isNull(U_3d.attr("dim")) && Rf_isNull(U0_3d.attr("dim"))) {
        throw std::invalid_argument(
                  "Either U_3d (prior matrices) or U0_3d (precomputed prior quantaties) has to be specified");
    }
    // set cube data from R 3D array
    arma::cube U_cube;
    Rcpp::IntegerVector dimU = (Rf_isNull(U_3d.attr("dim"))) ? U0_3d.attr("dim") : U_3d.attr("dim");
    if (Rf_isNull(U_3d.attr("dim"))) {
        if (!common_cov) U_cube.resize(dimU[0], dimU[1], dimU[2] / b_mat.ncol());
        else U_cube.resize(dimU[0], dimU[1], dimU[2]);
    } else {
        U_cube = arma::cube(U_3d.begin(), dimU[0], dimU[1], dimU[2]);
    }
    MVSERMix pc(Rcpp::as<arma::mat>(b_mat),
      Rcpp::as<arma::mat>(s_mat),
      Rcpp::as<arma::mat>(v_mat),
      U_cube);
    if (!Rf_isNull(U0_3d.attr("dim"))) {
        Rcpp::IntegerVector dimU0 = U0_3d.attr("dim");
        arma::cube U0_cube(U0_3d.begin(), dimU0[0], dimU0[1], dimU0[2]);
        pc.set_U0(U0_cube);
    }
    if (!Rf_isNull(vinv_3d.attr("dim"))) {
        // should be dimV = dimU
        Rcpp::IntegerVector dimV = vinv_3d.attr("dim");
        arma::cube vinv_cube(vinv_3d.begin(), dimV[0], dimV[1], dimV[2]);
        pc.set_Vinv(vinv_cube);
    }
    if (!Rf_isNull(Uinv_3d.attr("dim"))) {
        // inverse of prior matrices
        // relevent to computing updated scalar for prior variances
        // a feature used in mmbr package
        arma::cube Uinv_cube(Uinv_3d.begin(), dimU[0], dimU[1], dimU[2]);
        pc.set_Uinv(Uinv_cube);
    }
    if (!common_cov) pc.compute_posterior(Rcpp::as<arma::mat>(posterior_mixture_weights),
          Rcpp::as<arma::mat>(posterior_variable_weights), sigma0);
    else pc.compute_posterior_comcov(Rcpp::as<arma::mat>(posterior_mixture_weights),
          Rcpp::as<arma::mat>(posterior_variable_weights), sigma0);
    Rcpp::List res = Rcpp::List::create(
        Rcpp::Named("post_mean") = pc.PosteriorMean(),
        Rcpp::Named("post_sd")   = pc.PosteriorSD(),
        Rcpp::Named("post_cov")  = pc.PosteriorCov(),
        Rcpp::Named("post_zero") = pc.ZeroProb(),
        Rcpp::Named("post_neg")  = pc.NegativeProb());
    if (posterior_variable_weights.nrow() > 0) res.push_back(pc.PriorScalar(), "prior_scale_em_update");
    return res;
} // calc_sermix_rcpp

// [[Rcpp::export]]
Rcpp::List
fit_teem_rcpp(Rcpp::NumericMatrix X_mat,
  Rcpp::NumericVector             w_vec,
  Rcpp::NumericVector             U_3d,
  int                             maxiter,
  double                          converge_tol,
  double                          eigen_tol,
  bool                            verbose)
{
    // Convert R 3d array to Rcpp cube
    if (Rf_isNull(U_3d.attr("dim"))) {
        throw std::invalid_argument(
                  "U_3d has to be a 3D array");
    }
    Rcpp::IntegerVector dimU = U_3d.attr("dim");
    arma::cube U_cube(U_3d.begin(), dimU[0], dimU[1], dimU[2]);
    //
    TEEM teem(Rcpp::as<arma::mat>(X_mat),
      Rcpp::as<arma::vec>(w_vec),
      U_cube);
    teem.fit(maxiter, converge_tol, eigen_tol, verbose);
    Rcpp::List res = Rcpp::List::create(
        Rcpp::Named("w")         = teem.get_w(),
        Rcpp::Named("U")         = teem.get_U(),
        Rcpp::Named("objective") = teem.get_objective(),
        Rcpp::Named("maxd")      = teem.get_maxd());
    return res;
}
