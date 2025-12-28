#include "basics.h"
#include "mathfunc_probit.h"

// [[Rcpp::export]]

NumericMatrix Compute_logl_probit_with_theta(const List & X,  const NumericMatrix & Y, 
                                             const NumericVector & beta, const NumericMatrix & Theta ){
    
    const Map<MatrixXd> Y_ (as <Map<MatrixXd>> (Y));
    
    double epsilon = 1e-10;
    
    NumericMatrix logl_est = Compute_index_with_theta(X, beta, Theta);
    Map<MatrixXd> logl_est_ (as <Map<MatrixXd>> (logl_est));
    
    for (int i = 0; i < logl_est_.rows(); i++) {
        for (int j = 0; j < logl_est_.cols(); j++) {
            logl_est_(i, j) = R::pnorm5(logl_est_(i, j), 0.0, 1.0, 1, 0);
        }
    }
    
    logl_est_ = Y_.array() * (logl_est_.array() + epsilon).log()
        + (1 - Y_.array() ) * (1 - logl_est_.array() + epsilon).log();
    return wrap(logl_est_);
}

// [[Rcpp::export]]

NumericMatrix Compute_logl_probit_with_theta_FE(const List & X,  const NumericMatrix & Y, 
                                               const NumericVector & beta, 
                                               const NumericVector & fe_N, const NumericVector & fe_T,
                                               const NumericMatrix & Theta ){
    
    const Map<MatrixXd> Y_ (as <Map<MatrixXd>> (Y));
    
    double epsilon = 1e-10;
    
    NumericMatrix logl_est = Compute_index_with_theta_FE(X, beta, fe_N, fe_T, Theta);
    Map<MatrixXd> logl_est_ (as <Map<MatrixXd>> (logl_est));
    
    for (int i = 0; i < logl_est_.rows(); i++) {
        for (int j = 0; j < logl_est_.cols(); j++) {
            logl_est_(i, j) = R::pnorm5(logl_est_(i, j), 0.0, 1.0, 1, 0);
        }
    }
    
    logl_est_ = Y_.array() * (logl_est_.array() + epsilon).log()
        + (1 - Y_.array() ) * (1 - logl_est_.array() + epsilon).log();
    
    return wrap(logl_est_);
}

// [[Rcpp::export]]

NumericMatrix Compute_logl_probit_with_LR(const List & X,  const NumericMatrix & Y, 
                                          const NumericVector & beta, 
                                          const NumericMatrix & L, const NumericMatrix & R ){
    
    const Map<MatrixXd> Y_ (as <Map<MatrixXd>> (Y));
    double epsilon = 1e-10;
    
    NumericMatrix logl_est = Compute_index_with_LR(X, beta, L, R);
    Map<MatrixXd> logl_est_ (as <Map<MatrixXd>> (logl_est));
    
    for (int i = 0; i < logl_est_.rows(); i++) {
        for (int j = 0; j < logl_est_.cols(); j++) {
            logl_est_(i, j) = R::pnorm5(logl_est_(i, j), 0.0, 1.0, 1, 0);
        }
    }
    
    logl_est_ = Y_.array() * (logl_est_.array() + epsilon).log()
        + (1 - Y_.array() ) * (1 - logl_est_.array() + epsilon).log();
    
    return wrap(logl_est_);
}


// [[Rcpp::export]]

NumericMatrix Compute_logl_probit_with_LR_FE(const List & X,  const NumericMatrix & Y, 
                                            const NumericVector & beta, 
                                            const NumericVector & fe_N, const NumericVector & fe_T,
                                            const NumericMatrix & L, const NumericMatrix & R ){
    
    const Map<MatrixXd> Y_ (as <Map<MatrixXd>> (Y));
    double epsilon = 1e-10;
    
    NumericMatrix logl_est = Compute_index_with_LR_FE(X, beta,fe_N, fe_T, L, R);
    Map<MatrixXd> logl_est_ (as <Map<MatrixXd>> (logl_est));
    
    for (int i = 0; i < logl_est_.rows(); i++) {
        for (int j = 0; j < logl_est_.cols(); j++) {
            logl_est_(i, j) = R::pnorm5(logl_est_(i, j), 0.0, 1.0, 1, 0);
        }
    }
    
    logl_est_ = Y_.array() * (logl_est_.array() + epsilon).log()
        + (1 - Y_.array() ) * (1 - logl_est_.array() + epsilon).log();
    
    return wrap(logl_est_);
}

// [[Rcpp::export]]

List Compute_score_probit_with_theta(const List & X, const NumericMatrix & Y, 
                                    const NumericVector & beta, const NumericMatrix & Theta){
    const Map<MatrixXd> Y_ (as <Map<MatrixXd>> (Y));
    const int N = Y_.rows();
    const int T = Y_.cols();
    double epsilon = 1e-10;
    
    
    NumericMatrix score = Compute_index_with_theta(X, beta, Theta);
    Map<MatrixXd> score_ (as <Map<MatrixXd>> (score));
    NumericVector score_beta = wrap(VectorXd::Constant(X.length(), 1.0));
    
    MatrixXd pdf_ = MatrixXd::Constant(N, T, 0);
    
    
    for (int i = 0; i!=N; i++) {
        for (int j = 0; j!=T; j++) {
            pdf_(i, j) = R::dnorm4(score_(i, j), 0.0, 1.0, 0);
        }
    }
    for (int i = 0; i!=N; i++) {
        for (int j = 0; j!=T; j++) {
            score_(i, j) = R::pnorm5(score_(i, j), 0.0, 1.0, 1, 0);
        }
    }
    
    score_ = (score_.array() * (1.0 - score_.array()) + epsilon ).inverse() * (Y_.array() - score_.array() );
    score_ = score_.array() * pdf_.array();
    
    for(int i=0; i!=X.length(); i++){
        const Map<MatrixXd> X_ (as <Map<MatrixXd>> (X[i]));
        score_beta(i) = (score_.array() * X_.array()).sum();
    }
    return List::create(Named("score_beta") = score_beta, 
                        Named("score_theta") =  score_);
}

// [[Rcpp::export]]

List Compute_score_probit_with_theta_FE(const List & X, const NumericMatrix & Y, 
                                       const NumericVector & beta, 
                                       const NumericVector & fe_N, const NumericVector & fe_T,
                                       const NumericMatrix & Theta){
    const Map<MatrixXd> Y_ (as <Map<MatrixXd>> (Y));
    const int N = Y_.rows();
    const int T = Y_.cols();
    double epsilon = 1e-10;
    
    NumericMatrix score = Compute_index_with_theta_FE(X, beta, fe_N, fe_T, Theta);
    Map<MatrixXd> score_ (as <Map<MatrixXd>> (score));
    NumericVector score_beta = wrap(VectorXd::Constant(X.length(), 1.0));
    
    MatrixXd pdf_ = MatrixXd::Constant(N, T, 0);
    
    
    for (int i = 0; i!=N; i++) {
        for (int j = 0; j!=T; j++) {
            pdf_(i, j) = R::dnorm4(score_(i, j), 0.0, 1.0, 0);
        }
    }
    for (int i = 0; i!=N; i++) {
        for (int j = 0; j!=T; j++) {
            score_(i, j) = R::pnorm5(score_(i, j), 0.0, 1.0, 1, 0);
        }
    }
    
    score_ = (score_.array() * (1.0 - score_.array()) + epsilon ).inverse() * (Y_.array() - score_.array() );
    score_ = score_.array() * pdf_.array();
    
    for(int i=0; i!=X.length(); i++){
        const Map<MatrixXd> X_ (as <Map<MatrixXd>> (X[i]));
        score_beta(i) = (score_.array() * X_.array()).sum();
    }
    
    return List::create(Named("score_beta") = score_beta, 
                        Named("score_theta") =  score_);
}

// [[Rcpp::export]]

List Compute_score_probit_with_LR(const List & X, const NumericMatrix & Y, 
                                 const NumericVector & beta, 
                                 const NumericMatrix & L, const NumericMatrix & R){
    
    const Map<MatrixXd> Y_ (as <Map<MatrixXd>> (Y));
    const Map<MatrixXd> L_ (as <Map<MatrixXd>> (L));
    const Map<MatrixXd> R_ (as <Map<MatrixXd>> (R));
    const int N = Y_.rows();
    const int T = Y_.cols();
    double epsilon = 1e-10;
    
    
    NumericMatrix score = Compute_index_with_LR(X, beta, L, R);
    Map<MatrixXd> score_ (as <Map<MatrixXd>> (score));
    NumericVector score_beta = wrap(VectorXd::Constant(X.length(), 1.0));
    
    MatrixXd pdf_ = MatrixXd::Constant(N, T, 0);
    
    for (int i = 0; i!=N; i++) {
        for (int j = 0; j!=T; j++) {
            pdf_(i, j) = R::dnorm4(score_(i, j), 0.0, 1.0, 0);
        }
    }
    for (int i = 0; i!=N; i++) {
        for (int j = 0; j!=T; j++) {
            score_(i, j) = R::pnorm5(score_(i, j), 0.0, 1.0, 1, 0);
        }
    }
    
    score_ = (score_.array() * (1.0 - score_.array()) + epsilon ).inverse() * (Y_.array() - score_.array() );
    score_ = score_.array() * pdf_.array();
    
    for(int i=0; i!=X.length(); i++){
        const Map<MatrixXd> X_ (as <Map<MatrixXd>> (X[i]));
        score_beta(i) = (score_.array() * X_.array()).sum();
    }
    MatrixXd score_L_ = score_ * R_;
    MatrixXd score_R_ = score_.transpose() * L_;
    
    return List::create(Named("score_beta") = score_beta, 
                        Named("score_L") =  wrap(score_L_), 
                        Named("score_R") =  wrap(score_R_));
}    


// [[Rcpp::export]]

List Compute_score_probit_with_LR_FE(const List & X, const NumericMatrix & Y, 
                                    const NumericVector & beta, 
                                    const NumericVector & fe_N, const NumericVector & fe_T,
                                    const NumericMatrix & L, const NumericMatrix & R){
    
    const Map<MatrixXd> Y_ (as <Map<MatrixXd>> (Y));
    const Map<MatrixXd> L_ (as <Map<MatrixXd>> (L));
    const Map<MatrixXd> R_ (as <Map<MatrixXd>> (R));
    const int N = Y_.rows();
    const int T = Y_.cols();
    double epsilon = 1e-10;
    
    
    NumericMatrix score = Compute_index_with_LR_FE(X, beta, fe_N, fe_T, L, R);
    Map<MatrixXd> score_ (as <Map<MatrixXd>> (score));
    NumericVector score_beta = wrap(VectorXd::Constant(X.length(), 1.0));
    
    MatrixXd pdf_ = MatrixXd::Constant(N, T, 0);
    
    for (int i = 0; i!=N; i++) {
        for (int j = 0; j!=T; j++) {
            pdf_(i, j) = R::dnorm4(score_(i, j), 0.0, 1.0, 0);
        }
    }
    for (int i = 0; i!=N; i++) {
        for (int j = 0; j!=T; j++) {
            score_(i, j) = R::pnorm5(score_(i, j), 0.0, 1.0, 1, 0);
        }
    }
    
    score_ = (score_.array() * (1.0 - score_.array()) + epsilon ).inverse() * (Y_.array() - score_.array() );
    score_ = score_.array() * pdf_.array();
    
    for(int i=0; i!=X.length(); i++){
        const Map<MatrixXd> X_ (as <Map<MatrixXd>> (X[i]));
        score_beta(i) = (score_.array() * X_.array()).sum();
    }
    MatrixXd score_L_ = score_ * R_;
    MatrixXd score_R_ = score_.transpose() * L_;
    
    return List::create(Named("score_beta") = score_beta, 
                        Named("score_L") =  wrap(score_L_), 
                        Named("score_R") =  wrap(score_R_));
}    

// [[Rcpp::export]]

double Compute_obj_probit_with_theta(const List & X,  const NumericMatrix & Y, 
                                    const NumericVector & beta, const NumericMatrix & Theta, 
                                    const double sigma_sum, const double phi){
    NumericMatrix logl_probit = Compute_logl_probit_with_theta( X,  Y, beta, Theta );
    Map<MatrixXd> logl_probit_ (as <Map<MatrixXd>> (logl_probit));
    double loss = logl_probit_.sum();
    double obj  = -loss + phi * sigma_sum;
    return obj;
}

// [[Rcpp::export]]

double Compute_obj_probit_with_theta_FE(const List & X,  const NumericMatrix & Y, 
                                       const NumericVector & beta, 
                                       const NumericVector & fe_N, const NumericVector & fe_T,
                                       const NumericMatrix & Theta, 
                                       const double sigma_sum, const double phi){
    NumericMatrix logl_probit = Compute_logl_probit_with_theta_FE( X,  Y, beta, fe_N, fe_T, Theta );
    Map<MatrixXd> logl_probit_ (as <Map<MatrixXd>> (logl_probit));
    double loss = logl_probit_.sum();
    double obj  = -loss + phi * sigma_sum;
    return obj;
}

// [[Rcpp::export]]

double Compute_obj_probit_with_LR(const List & X,  const NumericMatrix & Y, 
                                 const NumericVector & beta, 
                                 const NumericMatrix & L, const NumericMatrix & R){
    NumericMatrix logl_probit = Compute_logl_probit_with_LR( X,  Y, beta, L, R );
    Map<MatrixXd> logl_probit_ (as <Map<MatrixXd>> (logl_probit));
    double loss = logl_probit_.sum();
    double obj  = -loss;
    return obj;
}

// [[Rcpp::export]]

double Compute_obj_probit_with_LR_FE(const List & X,  const NumericMatrix & Y, 
                                    const NumericVector & beta, 
                                    const NumericVector & fe_N, const NumericVector & fe_T,
                                    const NumericMatrix & L, const NumericMatrix & R){
    NumericMatrix logl_probit = Compute_logl_probit_with_LR_FE( X,  Y, beta,fe_N, fe_T, L, R );
    Map<MatrixXd> logl_probit_ (as <Map<MatrixXd>> (logl_probit));
    double loss = logl_probit_.sum();
    double obj  = -loss;
    return obj;
}

// [[Rcpp::export]]

NumericMatrix Compute_first_order_probit_with_LR(const List & X, const NumericMatrix & Y,
                                                const NumericVector & beta, 
                                                const NumericMatrix & L, const NumericMatrix & R){
    const Map<MatrixXd> Y_ (as <Map<MatrixXd>> (Y));
    const int N = Y_.rows();
    const int T = Y_.cols();
    double epsilon = 1e-10;

    NumericMatrix first_score = Compute_index_with_LR(X, beta, L, R);
    Map<MatrixXd> first_score_ (as <Map<MatrixXd>> (first_score));
    
    MatrixXd pdf_ = MatrixXd::Constant(N, T, 0);
    
    for (int i = 0; i!=N; i++) {
        for (int j = 0; j!=T; j++) {
            pdf_(i, j) = R::dnorm4(first_score_(i, j), 0.0, 1.0, 0);
        }
    }
    for (int i = 0; i!=N; i++) {
        for (int j = 0; j!=T; j++) {
            first_score_(i, j) = R::pnorm5(first_score_(i, j), 0.0, 1.0, 1, 0);
        }
    }
    
    first_score_ = (first_score_.array() * (1.0 - first_score_.array()) + epsilon ).inverse() * (Y_.array() - first_score_.array() );
    first_score_ = first_score_.array() * pdf_.array();
    
    return wrap(first_score_);
}

// [[Rcpp::export]]

NumericMatrix Compute_second_order_probit_with_LR(const List & X, const NumericMatrix & Y,
                                                  const NumericVector & beta, 
                                                  const NumericMatrix & L, const NumericMatrix & R){
    const Map<MatrixXd> Y_ (as <Map<MatrixXd>> (Y));
    NumericMatrix index = Compute_index_with_LR(X, beta, L, R);
    Map<MatrixXd> index_ (as <Map<MatrixXd>> (index));
    int N = Y_.rows();
    int T = Y_.cols();
    
    MatrixXd second_score_ = MatrixXd::Constant(N, T, 0);
    MatrixXd cdf_ = MatrixXd::Constant(N, T, 0);
    MatrixXd pdf_ = MatrixXd::Constant(N, T, 0);
    
    for (int i = 0; i!=N; i++) {
        for (int j = 0; j!=T; j++) {
            pdf_(i, j) = R::dnorm4(index_(i, j), 0.0, 1.0, 0);
        }
    }
    for (int i = 0; i!=N; i++) {
        for (int j = 0; j!=T; j++) {
            cdf_(i, j) = R::pnorm5(index_(i, j), 0.0, 1.0, 1, 0);
        }
    }
    
    
    double epsilon = 1e-10;
    
    second_score_ = index_.array() * pdf_.array() * (Y_.array() - cdf_.array()) * (cdf_.array() * (1-cdf_.array()) + epsilon).inverse();
    second_score_ = second_score_.array() + pdf_.array() *pdf_.array() * (cdf_.array() * cdf_.array() + Y_.array()* (1 - 2* cdf_.array())) * ( cdf_.array() * cdf_.array() * (1 - cdf_.array()) * (1- cdf_.array()) + epsilon).inverse();
    for(int i=0;i!=N; i++){
        for(int t=0; t!=T; t++){
            if(second_score_(i, t) <= 1e-10){
                second_score_(i, t) = 1e-10;
            }
        }
    }
    return wrap(second_score_);
}





