#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_randist.h>
#include <cmath>
using namespace Rcpp;

// --------------------------------------------------------------------------- //
// General Helper Functions
// --------------------------------------------------------------------------- //

// [[Rcpp::export]]
void rcpp_rprintf(arma::vec v){
    // printing values of all the elements of Rcpp vector  
    int n = v.n_elem;
    for(int i=0; i < n; ++i){
        Rprintf("the value of v[%i] : %f \n", i, v(i));
    }
}

//[[Rcpp::export]]
IntegerVector csample_int(IntegerVector x, int size, bool replace, NumericVector prob = NumericVector::create()) {
    return sample(x, size, true, prob);
}

//[[Rcpp::export]]
arma::vec v_fmod(arma::vec v, double d) {
    int n = v.n_elem;
    arma::vec v1(n);
    
    for (int i = 0; i < n; i++) {
        v1(i) = fmod(v(i), d);
        if (v1(i) < 0) {
            v1(i) = v1(i) + d;
        }
    }
    
    return v1;
}

// [[Rcpp::export]]
int my_sign(double num) {
    if (num > 0) {
        return 1;
    } else if (num < 0){
        return -1;
    } else {
        return 0;
    }
}

//[[Rcpp::export]]
double find_max_norm(double a, double b, double mu) {
    if (mu > a & mu < b) {
        return mu;
    } else if (mu < a) {
        return a;
    } else {
        return b;
    }
}

// --------------------------------------------------------------------------- //
// Distribution Functions
// --------------------------------------------------------------------------- //

double my_pnorm(double x, double mu, double sd, bool lower_tail = true, bool take_log = false) {
    NumericVector X(1), f(1);
    X(0) = x;
    f = pnorm(X, mu, sd, lower_tail, take_log);
    
    return f(0);
}

// --------------------------------------------------------------------------- //
// Density Functions
// --------------------------------------------------------------------------- //

//[[Rcpp::export]]
double my_dnorm(double x, double mu = 0.0, double sd = 1.0, bool take_log = false) {
    NumericVector X(1), f(1);
    X(0) = x;
    f = dnorm(X, mu, sd, take_log);
    
    return f(0);
}

//[[Rcpp::export]]
double my_dgamma(double x, double shape, double scale, bool take_log = false) {
    NumericVector X(1), f(1);
    X(0) = x;
    f = dgamma(X, shape, scale, take_log);
    
    return f(0);
}

// [[Rcpp::export]]
double dvonmises_cpp(double theta, double mu, double kappa, bool take_log = false) {
    double log_dens;
    
    log_dens = -log(2*M_PI) - log(gsl_sf_bessel_I0(kappa)) + kappa * cos(theta - mu);
    if (take_log) {
        return log_dens;
    }
    return exp(log_dens);
}

// [[Rcpp::export]]
double wn_approx_small_sigma(double theta, double mu, double sd, int order, bool take_log = false) {
    double p, dens = 0.0, sigma_sq = pow(sd,2);
    
    if (order == 0) {
        dens = 1/(2*M_PI);
    } else {
        for (int k = -order; k <= order; k++) {
            p = exp(-pow(theta + 2*M_PI*k - mu, 2) / (2*sigma_sq));
            dens += p;
        }
        dens = dens / sqrt(2*M_PI*sigma_sq);
    }
    
    if (take_log) {
        return log(dens);
    }
    return dens;
}

//[[Rcpp::export]]
double wn_approx_large_sigma(double theta, double mu, double sd, int order, bool take_log = false) {
    double p, dens = 0.0, sigma_sq = pow(sd, 2);
    
    if (order == 0) {
        dens = 1/(2*M_PI);
    } else {
        for (int k = 0; k <= order; k++) {
            p = pow(exp(-sigma_sq*.5), pow(order,2)) * cos(k * (theta - mu));
            dens += p;
        }
        dens = (1 + dens) / (2*M_PI);
    }
    
    if (take_log) {
        return log(dens);
    }
    return dens;
}

// [[Rcpp::export]]
double dwrappednorm_cpp(double theta, double mu, double sd, bool take_log = false) {
    double log_dens;
    
    if (sd < 0.76) {
        log_dens = wn_approx_small_sigma(theta, mu, sd, 1, true);
    } else if (sd < 1.53) {
        log_dens = wn_approx_small_sigma(theta, mu, sd, 1, true);
    } else if (sd < 2.31) {
        log_dens = wn_approx_small_sigma(theta, mu, sd, 2, true);
    } else if (sd < 2.73) {
        log_dens = wn_approx_large_sigma(theta, mu, sd, 3, true);
    } else if (sd < 4.09) {
        log_dens = wn_approx_large_sigma(theta, mu, sd, 2, true);
    } else {
        log_dens = wn_approx_large_sigma(theta, mu, sd, 1, true);
    }
    
    if (take_log) {
        return log_dens;
    }
    return exp(log_dens);
}

//[[Rcpp::export]]
double dtruncnorm_cpp(double x, double a, double b, double mu, double sd, bool take_log = false) {
    NumericVector rats = NumericVector::create(x,a,b);
    rats = (rats - mu)/sd;
    double log_dens = -log(sd) + my_dnorm(rats(0),0.0,1.0,true) - log(my_pnorm(rats(2),0.0,1.0) - my_pnorm(rats(1),0.0,1.0));
    if (take_log) {
        return log_dens;
    }
    return exp(log_dens);
}

//[[Rcpp::export]]
double dabeley_cpp(double theta, double x, double alpha, double beta, double kappa, double mu, double lambda) {
    
    double c = log(alpha) + alpha * log(beta)  - log(2*M_PI)  - log(cosh(kappa));
    double wc = log(1 + lambda * sin(theta - mu));
    double weib = (alpha - 1) * log(x) - pow(beta * x, alpha) * (1 - tanh(kappa) * cos(theta - mu));
    return exp(c + wc + weib);
}

//[[Rcpp::export]]
double my_dbeta(double x, double shape1, double shape2, bool take_log = false) {
    NumericVector X(1), f(1);
    X(0) = x;
    f = dbeta(X, shape1, shape2, take_log);
    
    return f(0);
}

//[[Rcpp::export]]
double dbeta_rescale_cpp(double x, double shape1, double shape2, bool take_log = false) {
    // Rescales the Beta(shape1, shape2) distribution to lie on (-1,1) instead of (0,1)
    return 0.5 * my_dbeta((x/2 + .5), shape1, shape2, take_log);
}

// --------------------------------------------------------------------------- //
// RNG Functions
// --------------------------------------------------------------------------- //

//[[Rcpp::export]]
arma::vec onesided_rtruncnorm(int n, double a) {
    // Assumes left-sided truncation at a
    arma::vec smpl(n);
    double z, q_z, u, alpha_star = .5 * (a + sqrt(pow(a,2) + 4));
    
    for (int i = 0; i < n; i++) {
        z = rexp(1, alpha_star)(0);
        z += a;
        q_z = exp(-0.5 * pow(z - alpha_star,2));
        u = runif(1)(0);
        while (u > q_z) {
            z = rexp(1, alpha_star)(0);
            z += a;
            q_z = exp(-0.5 * pow(z - alpha_star,2));
            u = runif(1)(0);
        }
        smpl(i) = z;
    }
    
    return smpl;
}

//[[Rcpp::export]]
arma::vec twosided_rtruncnorm(int n, double a, double b) {
    arma::vec smpl(n);
    double z, u, q_z;
    
    for (int i = 0; i < n; i++) {
        z = runif(1, a, b)(0);
        if (0 > a & 0 < b) {
            q_z = exp(-pow(z,2) * .5);
        } else if (b < 0) {
            q_z = exp((pow(b,2) - pow(z,2)) * .5);
        } else if (a > 0) {
            q_z = exp((pow(a,2) - pow(z,2)) * .5);
        }
        u = runif(1)(0);
        while (u > q_z) {
            z = runif(1, a, b)(0);
            if (0 > a & 0 < b) {
                q_z = exp(-pow(z,2) * .5);
            } else if (b < 0) {
                q_z = exp((pow(b,2) - pow(z,2)) * .5);
            } else if (a > 0) {
                q_z = exp((pow(a,2) - pow(z,2)) * .5);
            }
        }
        smpl(i) = z;
    }
    return smpl;
}

//[[Rcpp::export]]
arma::vec rtruncnorm_cpp(int n, double a, double b, double mu, double sd) {
    arma::vec smpl(n);
    
    if (not arma::is_finite(a) & not arma::is_finite(b)) {
        smpl = rnorm(n, mu, sd);
    } else if (not arma::is_finite(a) & arma::is_finite(b)) {
        double tmp = a;
        a = -b;
        b = -tmp;
    }
    
    if (not arma::is_finite(b)) {
        a = (a - mu) / sd;
        smpl = onesided_rtruncnorm(n, a);
    } else {
        a = (a - mu) / sd;
        b = (b - mu) / sd;
        smpl = twosided_rtruncnorm(n, a, b);
    }
    for (int i = 0; i < n; i++){
        smpl(i) = sd * smpl(i) + mu;
    }
    
    if (not arma::is_finite(a) & arma::is_finite(b)) {
        smpl = -1.0 * smpl;
    }
    
    return smpl;
}

//[[Rcpp::export]]
arma::vec rbeta_rescale_cpp(int n, double shape1, double shape2) {
    arma::vec U = rbeta(n, shape1, shape2);
    arma::vec V = 2 * U - 1.0;
    return V;
}

//[[Rcpp::export]]
arma::vec rwrappednorm_cpp(int n, double mu, double sd) {
    return v_fmod(rnorm(n, mu, sd), 2.0*M_PI);
}

//[[Rcpp::export]]
arma::vec rwrappedcauchy_cpp(int n, double mu, double kappa) {
    return v_fmod(rcauchy(n, mu, kappa), 2.0*M_PI);
}

// [[Rcpp::export]]
arma::vec rvonmises_cpp(int n, double mu, double kappa) {
    arma::vec theta(n), U(3);
    double t_star, f, tau, rho, r, z, c;
    int i;
    
    for (i = 0; i < n; i++) {
        t_star = -1;
        while (t_star < 0) {
            U =  runif(3,0,1);
            tau = 1 + sqrt(1 + 4*pow(kappa, 2));
            rho = (tau - sqrt(2*tau)) / (2*kappa);
            r = (1 + pow(rho,2))/(2*rho);
            z = cos(M_PI* U(0));
            f = (1 + r*z) / (r + z);
            c = kappa * (r - f);
            t_star = log(c / U(1)) + 1 - c;
        }
        theta(i) = my_sign(U(2) - .5) * acos(f) + mu;
    }
    
    return theta;
}

// [[Rcpp::export]]
arma::vec rssvm_cpp(int n, double mu, double kappa, double lambda) {
    arma::vec t(n), t1 = rvonmises_cpp(n, 0, kappa);
    double U;
    for (int i = 0; i < n; i++) {
        U = runif(1,0,1)(0);
        if (U < .5*(1 + lambda * sin(t1(i)))) {
            t(i) = t1(i);
        } else {
            t(i) = -t1(i);
        }
    }
    t = v_fmod(t + mu, 2.0*M_PI);
    return(t);
}

// [[Rcpp::export]]
arma::mat rdirichlet_cpp(int size, arma::vec alpha) {
    int K = alpha.n_elem;
    arma::mat smp(size, K);
    arma::vec X(K);
    double tot;
    
    for (int i = 0; i < size; i++) {
        for (int k = 0; k < K; k++){
            X(k) = rgamma(1, alpha(k), 1)(0);
        }
        tot = sum(X);
        smp.row(i) = X.t() / tot;
    }
    
    return smp;
}

//[[Rcpp::export]]
arma::mat rabeley_cpp(int n, double alpha, double beta, double kappa, double mu, double lambda) {
    arma::mat dat(n,2);
    arma::vec drawn_obs(2);
    NumericVector U = runif(n, 0,1);
    double t1, t, scale_par, x;
    
    for (int i = 0; i < n; i++) {
        t1 = rwrappedcauchy_cpp(1, 0, tanh(kappa / 2))(0);
        t = t1;
        if (U(i) > (1 + lambda * sin(t1)) / 2) {
            t = t1;
        }
        t += M_PI;
        t = fmod(t, 2*M_PI);
        scale_par = beta * pow(1 - tanh(kappa) * cos(mu - t), 1/alpha);
        x = rweibull(1, alpha, 1/scale_par)(0);
        drawn_obs = {t,x};
        dat.row(i) = drawn_obs.t();
    }
    
    return dat;
}

// [[Rcpp::export]]
arma::vec rssvm_mixture_cpp(int n, arma::vec mix_prop, arma::vec mu, arma::vec kappa, arma::vec lambda) {
    int k, K = mix_prop.n_elem;
    arma::vec t(n);
    IntegerVector labs(K);
    
    for (int k = 0; k < K; k++) {
        labs(k) = k;
    }
    
    for (int i = 0; i < n; i++) {
        k = csample_int(labs, 1, false, as<NumericVector>(wrap(mix_prop)))(0);
        t(i) = rssvm_cpp(1, mu(k), kappa(k), lambda(k))(0);
    }
    return v_fmod(t, 2*M_PI);
}

// --------------------------------------------------------------------------- //
// Sampler Helper Functions
// --------------------------------------------------------------------------- //

// [[Rcpp::export]]
arma::ivec count_cluster_size(arma::ivec class_labels, int K) {
    arma::ivec counts(K);
    int sum = 0, n = class_labels.n_elem;
    
    for (int k = 0; k < K; k++) {
        sum = 0;
        for (int i = 0; i < n; i++) {
            if (class_labels(i) == k+1) {
                sum++;
            }
        }
        counts(k) = sum;
    }
    
    
    return counts;
}

// [[Rcpp::export]]
arma::vec sample_proposal_dist_cpp(arma::vec mu_prev, arma::vec sd_prev) {
    // Calls:
    //    : rtruncnorm_cpp(int n, double a, double b, double mu, double sd)
    //    : rwrappednorm_cpp(int n, double mu, double sd)
    arma::vec nu_sample(5);
    nu_sample(0) = rtruncnorm_cpp(1, 0, R_PosInf, mu_prev(0), sd_prev(0))(0);
    nu_sample(1) = rtruncnorm_cpp(1, 0, R_PosInf, mu_prev(1), sd_prev(1))(0);
    nu_sample(2) = rtruncnorm_cpp(1, 0, R_PosInf, mu_prev(2), sd_prev(2))(0);
    nu_sample(3) = rwrappednorm_cpp(1, mu_prev(3), sd_prev(3))(0);
    nu_sample(4) = rtruncnorm_cpp(1, -1, 1, mu_prev(4), sd_prev(4))(0);
    return nu_sample;
}

// [[Rcpp::export]]
arma::vec accept_ratio_cpp(arma::mat dat, arma::vec y_star, arma::vec y_prev, arma::vec sd_prev, 
                           Rcpp::List hyperpar = Rcpp::List::create(_["alpha_shape"] = 1000, _["alpha_scale"] = 0.001,
                                                              _["beta_shape"] = 1000, _["beta_scale"] = 0.001,
                                                              _["kappa_shape"] = 1000, _["kappa_scale"] = 0.001,
                                                              _["mu_mean"] = 0, _["mu_conc"] = 0.001, 
                                                              _["lambda_a"] = 1, _["lambda_b"] = 1, _["alpha0"] = 1)) {
    // Calls:
    //   x: dabeley_cpp(double t, double x, double alpha, double beta, double kappa, double mu, double lambda)
    //    : dtruncnorm_cpp(double x, double a, double b, double mu, double sd)
    //    : dbeta_rescale_cpp(x, double a, double b)
    // parameter order: alpha, beta, kappa, mu, lambda
    
    arma::vec t = dat.col(0), x = dat.col(1);
    int n = t.n_elem;
    double p_ystar, p_yprev, J_prev, J_star, r, logL_prev = 0.0, logL_prop = 0.0;
    arma::vec y_tmp(5), shape = {hyperpar["alpha_shape"], hyperpar["beta_shape"], hyperpar["kappa_shape"]}, scale = {hyperpar["alpha_scale"], hyperpar["beta_scale"], hyperpar["kappa_scale"]};
    for (int i = 0; i < n; i++) {
        logL_prev += log(dabeley_cpp(t(i), x(i), y_prev(0), y_prev(1), y_prev(2), y_prev(3), y_prev(4)));
    }
    arma::vec A(5);
    
    for (int j = 0; j < 5; j++) {
        y_tmp = y_prev;
        y_tmp(j) = y_star(j);
        logL_prop = 0.0;
        for (int i = 0; i < n; i++) {
            logL_prop += log(dabeley_cpp(t(i), x(i), y_tmp(0), y_tmp(1), y_tmp(2), y_tmp(3), y_tmp(4)));
        }
        if (j <= 2) {
            p_ystar = my_dgamma(y_star(j), shape(j), scale(j), true);
            p_yprev = my_dgamma(y_prev(j), shape(j), scale(j), true);
            J_prev = dtruncnorm_cpp(y_prev(j), 0, R_PosInf, y_star(j), sd_prev(j), true);
            J_star = dtruncnorm_cpp(y_star(j), 0, R_PosInf, y_prev(j), sd_prev(j), true);
        } else if (j == 3) {
            p_ystar = dvonmises_cpp(y_star(j), hyperpar["mu_mean"], hyperpar["mu_conc"], true);
            p_yprev = dvonmises_cpp(y_prev(j), hyperpar["mu_mean"], hyperpar["mu_conc"], true);
            J_prev = dwrappednorm_cpp(y_prev(j), y_star(j), sd_prev(j), true);
            J_star = dwrappednorm_cpp(y_star(j), y_prev(j), sd_prev(j), true);
        } else {
            p_ystar = dbeta_rescale_cpp(y_star(j), hyperpar["lambda_a"], hyperpar["lambda_b"], true);
            p_yprev = dbeta_rescale_cpp(y_prev(j), hyperpar["lambda_a"], hyperpar["lambda_b"], true);
            J_prev = dtruncnorm_cpp(y_prev(j), -1, 1, y_star(j), sd_prev(j), true);
            J_prev = dtruncnorm_cpp(y_star(j), -1, 1, y_prev(j), sd_prev(j), true);
        }
        r = logL_prop + p_ystar + J_prev - logL_prev - p_yprev - J_star;
        A(j) = min(NumericVector::create(0,r));
    }
    return A;
    }

// [[Rcpp::export]]
arma::ivec sample_class_labels_cpp(arma::mat dat, arma::mat nu_t, int K) {
    // Calls:
    //   x: csample_int(IntegerVector K, int size, bool replace, NumericVector probs)
    //   x: dabeley_cpp(double t, double x, double alpha, double beta, double kappa, double mu, double lambda)
    
    int n = dat.n_rows;
    double tot_prob;
    IntegerVector labs(K);
    NumericVector probs(K);
    arma::ivec label_sample(n);
    
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < K; k++) {
            labs(k) = k+1;
            probs(k) = nu_t(k,5) * dabeley_cpp(dat(i,0), dat(i,1), nu_t(k,0), nu_t(k,1), nu_t(k,2), nu_t(k,3), nu_t(k,4));
        }
        tot_prob = sum(probs);
        label_sample(i) = csample_int(labs, 1, true, probs/tot_prob)(0);
        // Rprintf("obs_%i assigned %i \n", i, label_sample(i));
    }
    return label_sample;
}

// [[Rcpp::export]]
arma::vec adjust_sd_cpp(int iter, arma::vec sd_prev, int R, arma::vec accept_probs){
    // Calls:
    //   None
    int K = sd_prev.n_elem;
    double s;
    arma::vec sd_new = sd_prev;
    
    if ((iter % 50) == 0) {
        for (int k = 0; k < 5; k++) {
            s = min(NumericVector::create(0.01, sqrt((double)R / iter)));
            if (accept_probs(k) > 0.44) {
                sd_new(k) = sqrt(pow(sd_prev(k), 2) + s);
            } else{
                if (sqrt(pow(sd_prev(k), 2) - s) > 0) {
                    sd_new(k) = sqrt(pow(sd_prev(k), 2) - s);
                } else {
                    sd_new(k) = sd_prev(k);
                }
            }
        }
    }
    return sd_new;
}

// --------------------------------------------------------------------------- //
// Sampler
// --------------------------------------------------------------------------- //

// [[Rcpp::export]]
Rcpp::List abeley_MM_sampler_MH_cpp(arma::mat dat, int K,
                                    Rcpp::List control = List::create(_["Q"] = 2000, _["burnin"] = 1000, _["sd_init"] = 1, _["thin"] = 1),
                                    Rcpp::List hyperpar = List::create(_["alpha_shape"] = 1000, _["alpha_scale"] = 0.001,
                                                                       _["beta_shape"] = 1000, _["beta_scale"] = 0.001,
                                                                       _["kappa_shape"] = 1000, _["kappa_scale"] = 0.001,
                                                                       _["mu_mean"] = 0, _["mu_conc"] = 0.001,
                                                                       _["lambda_a"] = 1, _["lambda_b"] = 1, _["alpha0"] = 1),
                                    IntegerVector obs_lbls = IntegerVector::create(),
                                    IntegerVector R_ind = IntegerVector::create()) {
    // data should be:
    // theta     x
    // t1        x1
    // t2        x2
    // ...       ...

    // Calls:
    //   x: rvonmises_cpp(int n, double mu, double kappa)
    //   x: rbeta_rescale_cpp(int n, double a, double b)
    //   x: rdirichlet_cpp(int n, NumericVector alpha)
    //   x: csample_int(IntegerVector x, int size, bool replace, NumericVector prob)
    //   x: rssvm_mixture_cpp(int n, NumericVector mix_prop, NumericVector mu, NumericVector kappa, NumericVector lambda)
    //   x: adjust_sd_cpp(int iter, NumericVector sd_prev, int R, NumericVector accept_probs)
    //   x: sample_proposal_dist_cpp(NumericVector mu_prev, NumericVector sd_prev)
    //   x: accept_ratio_cpp(NumericMatrix dat, NumericVector y_star, NumericVector y_prev, NumericVector sd_prev, Rcpp::List hyperpar)
    //   x: sample_class_labels_cpp(NumericVector dat, arma::mat nu_t, int K)

    // MCMC Draws will be stored in a QxKx6 array
    // The Q dim is the rows; for each iteration of MCMC
    // The K dim is the columns; for each mixture component
    // The J dim = 6 is the depth; for each parameter
    //             J = 1: alpha
    //             J = 2: beta
    //             J = 3: kappa
    //             J = 4: mu
    //             J = 5: lambda
    //             J = 6: tau (mixing props)

    int Q = control["Q"], n = dat.n_rows;
    double sd_init = control["sd_init"], alpha0 = hyperpar["alpha0"];
    arma::vec tau(K), shape = {hyperpar["alpha_shape"], hyperpar["beta_shape"], hyperpar["kappa_shape"]}, scale = {hyperpar["alpha_scale"], hyperpar["beta_scale"], hyperpar["kappa_scale"]};
    IntegerVector labs(K);
    arma::mat nu_t(K, 6);
    arma::cube chains(Q, K, 6);
    arma::imat class_labels(Q,n);

    // Setting initial random class assignments
    arma::vec alpha_0(K);
    for (int k = 0; k < K; k++) {
        labs(k)= k+1;
        alpha_0(k) = hyperpar["alpha0"];
    }
    tau = rdirichlet_cpp(1, alpha_0).row(0).t();
    arma::ivec s = csample_int(labs, n, true, as<NumericVector>(wrap(tau)));
    class_labels.row(0) = s.t();
    
    // // Checking for observed labels
    // if (obs_lbls.length() == n) {
    //     for (int i = 0; i < n; i++) {
    //         if (R_ind(i) == 1) {
    //             class_labels(0,i) = obs_lbls(i);
    //         }
    //     }
    // }
    
    // Initialization
    arma::mat sd_prev(K,5);
    sd_prev.fill(sd_init);
    arma::vec s_p(5), mu_prev = {1,1,1,0,0};
    arma::ivec cluster_sizes = count_cluster_size(class_labels.row(0).t(), K);
    
    // for (int i = 0; i < K; i++) {
    //     Rprintf("Inital Clusters: \n n_%i = %i \n", i+1, cluster_sizes(i));
    // }
    
    for (int k = 0; k < K; k++) {
        s_p = sd_prev.row(k).t();
        nu_t(k, arma::span(0,4)) = sample_proposal_dist_cpp(mu_prev, s_p).t();
    }
    // double conc_pred, theta_new, theta_pred;
    // if (obs_lbls.length() == n) {
    //     for (int i = 0; i < n; i++) {
    // 
    //     }
    // }
    
    nu_t.col(5) = tau;
    chains.row(0) = nu_t;
    arma::mat a_probs(K, 5);
    arma::mat nu_prev = chains.row(0);
    arma::vec sd_prev_k(5), a_probs_k(5), nu_prev_k(5), nu_star_k(5), A(5), U(5);
    
    for (int q = 1; q < Q; q++) {
    arma::vec alpha_n = as<arma::vec>(wrap(cluster_sizes)) + alpha0;
    tau = rdirichlet_cpp(1, alpha_n).row(0).t();
    chains.subcube(arma::span(q), arma::span(0,K-1), arma::span(5)) = tau;

    for (int k = 0; k < K; k++) {
        // subset the data
        arma::mat data_k = dat.rows(arma::find(class_labels.row(q-1) == k+1));
        int n_k = data_k.n_rows;
        // Rprintf("n_%i = %i , ", k+1, n_k);
        nu_prev = chains.row(q-1);
        nu_prev_k = nu_prev(k, arma::span(0,4)).t();
        sd_prev_k = sd_prev.row(k).t();
        a_probs_k = a_probs.row(k).t();
        sd_prev.row(k) = adjust_sd_cpp(q, sd_prev_k, 50, a_probs_k/50.0).t();

        nu_star_k = sample_proposal_dist_cpp(nu_prev_k, sd_prev_k);

        A = accept_ratio_cpp(data_k, nu_star_k, nu_prev_k, sd_prev_k, hyperpar);
        
        U = runif(5);
        for (int j = 0; j < 5; j++) {
            if (log(U(j)) >= A(j)) {
                nu_star_k(j) = nu_prev_k(j);
            } else {
                a_probs_k(j) += 1;
                nu_star_k(j) = nu_star_k(j);
            }
        }

        if ((q % 50) == 0) {
            a_probs.zeros(K,5);
        }

        nu_t(k, arma::span(0,4)) = nu_star_k.t();

        cluster_sizes(k) = n_k;
    }
    // Rprintf("\n");
    nu_t.col(5) = tau;
    chains.subcube(arma::span(q),arma::span(0,K-1),arma::span(0,5)) = nu_t;
    
    arma::ivec lbls_t(n);
    if (K > 1) {
        lbls_t = sample_class_labels_cpp(dat, nu_t, K);
    }

    class_labels.row(q) = lbls_t.t();
    if (obs_lbls.length() == n) {
        for (int i = 0; i < n; i++) {
            if (R_ind(i) == 1) {
                class_labels(q,i) = obs_lbls(i);
            }
        }
    }
    
    }

return Rcpp::List::create(Rcpp::Named("chains") = chains,
                          Rcpp::Named("class_labels") = class_labels);
}


