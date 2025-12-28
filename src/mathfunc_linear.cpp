#include "basics.h"
#include "mathfunc_linear.h"

// [[Rcpp::export]]

NumericMatrix Compute_logl_linear_with_theta(const List & X,  const NumericMatrix & Y, 
                                             const NumericVector & beta, const NumericMatrix & Theta ){
    
    const Map<MatrixXd> Y_ (as <Map<MatrixXd>> (Y));
    

    NumericMatrix logl_est = Compute_index_with_theta(X, beta, Theta);
    Map<MatrixXd> logl_est_ (as <Map<MatrixXd>> (logl_est));
    
    logl_est_ = -1.0 /2.0 * (Y_.array() - logl_est_.array()).array() * (Y_.array() - logl_est_.array()).array();
    
    return wrap(logl_est_);
}

// [[Rcpp::export]]

NumericMatrix Compute_logl_linear_with_theta_FE(const List & X,  const NumericMatrix & Y, 
                                                const NumericVector & beta, 
                                                const NumericVector & fe_N, const NumericVector & fe_T,
                                                const NumericMatrix & Theta ){
    
    const Map<MatrixXd> Y_ (as <Map<MatrixXd>> (Y));
    

    NumericMatrix logl_est = Compute_index_with_theta_FE(X, beta, fe_N, fe_T, Theta);
    Map<MatrixXd> logl_est_ (as <Map<MatrixXd>> (logl_est));
    
    logl_est_ = -1.0 /2.0 * (Y_.array() - logl_est_.array()).array() * (Y_.array() - logl_est_.array()).array();
    
    return wrap(logl_est_);
}

// [[Rcpp::export]]

NumericMatrix Compute_logl_linear_with_LR(const List & X,  const NumericMatrix & Y, 
                                          const NumericVector & beta, 
                                          const NumericMatrix & L, const NumericMatrix & R ){
    
    const Map<MatrixXd> Y_ (as <Map<MatrixXd>> (Y));

    NumericMatrix logl_est = Compute_index_with_LR(X, beta, L, R);
    Map<MatrixXd> logl_est_ (as <Map<MatrixXd>> (logl_est));
    
    logl_est_ = -1.0 /2.0 * (Y_.array() - logl_est_.array()).array() * (Y_.array() - logl_est_.array()).array();
    
    return wrap(logl_est_);
}


// [[Rcpp::export]]

NumericMatrix Compute_logl_linear_with_LR_FE(const List & X,  const NumericMatrix & Y, 
                                             const NumericVector & beta, 
                                             const NumericVector & fe_N, const NumericVector & fe_T,
                                             const NumericMatrix & L, const NumericMatrix & R ){
    
    const Map<MatrixXd> Y_ (as <Map<MatrixXd>> (Y));

    NumericMatrix logl_est = Compute_index_with_LR_FE(X, beta,fe_N, fe_T, L, R);
    Map<MatrixXd> logl_est_ (as <Map<MatrixXd>> (logl_est));
    
    logl_est_ = -1.0 /2.0 * (Y_.array() - logl_est_.array()).array() * (Y_.array() - logl_est_.array()).array();
    
    return wrap(logl_est_);
}

// [[Rcpp::export]]

List Compute_score_linear_with_theta(const List & X, const NumericMatrix & Y, 
                                     const NumericVector & beta, const NumericMatrix & Theta){
    const Map<MatrixXd> Y_ (as <Map<MatrixXd>> (Y));
    const int N = Y_.rows();
    const int T = Y_.cols();

    
    NumericMatrix score = Compute_index_with_theta(X, beta, Theta);
    Map<MatrixXd> score_ (as <Map<MatrixXd>> (score));
    NumericVector score_beta = wrap(VectorXd::Constant(X.length(), 1.0));
    

    score_= Y_ -  score_; 
    
    for(int i=0; i!=X.length(); i++){
        const Map<MatrixXd> X_ (as <Map<MatrixXd>> (X[i]));
        score_beta(i) = (score_.array() * X_.array()).sum();
    }
    return List::create(Named("score_beta") = score_beta, 
                        Named("score_theta") =  score_);
}

// [[Rcpp::export]]

List Compute_score_linear_with_theta_FE(const List & X, const NumericMatrix & Y, 
                                        const NumericVector & beta, 
                                        const NumericVector & fe_N, const NumericVector & fe_T,
                                        const NumericMatrix & Theta){
    const Map<MatrixXd> Y_ (as <Map<MatrixXd>> (Y));
    const int N = Y_.rows();
    const int T = Y_.cols();

    NumericMatrix score = Compute_index_with_theta_FE(X, beta, fe_N, fe_T, Theta);
    Map<MatrixXd> score_ (as <Map<MatrixXd>> (score));
    NumericVector score_beta = wrap(VectorXd::Constant(X.length(), 1.0));
    
    score_= Y_ -  score_; 
    
    for(int i=0; i!=X.length(); i++){
        const Map<MatrixXd> X_ (as <Map<MatrixXd>> (X[i]));
        score_beta(i) = (score_.array() * X_.array()).sum();
    }
    
    return List::create(Named("score_beta") = score_beta, 
                        Named("score_theta") =  score_);
}

// [[Rcpp::export]]

List Compute_score_linear_with_LR(const List & X, const NumericMatrix & Y, 
                                  const NumericVector & beta, 
                                  const NumericMatrix & L, const NumericMatrix & R){
    
    const Map<MatrixXd> Y_ (as <Map<MatrixXd>> (Y));
    const Map<MatrixXd> L_ (as <Map<MatrixXd>> (L));
    const Map<MatrixXd> R_ (as <Map<MatrixXd>> (R));
    const int N = Y_.rows();
    const int T = Y_.cols();

    
    NumericMatrix score = Compute_index_with_LR(X, beta, L, R);
    Map<MatrixXd> score_ (as <Map<MatrixXd>> (score));
    NumericVector score_beta = wrap(VectorXd::Constant(X.length(), 1.0));
    
    score_= Y_ -  score_; 
    
    
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

List Compute_score_linear_with_LR_FE(const List & X, const NumericMatrix & Y, 
                                     const NumericVector & beta, 
                                     const NumericVector & fe_N, const NumericVector & fe_T,
                                     const NumericMatrix & L, const NumericMatrix & R){
    
    const Map<MatrixXd> Y_ (as <Map<MatrixXd>> (Y));
    const Map<MatrixXd> L_ (as <Map<MatrixXd>> (L));
    const Map<MatrixXd> R_ (as <Map<MatrixXd>> (R));
    const int N = Y_.rows();
    const int T = Y_.cols();

    
    NumericMatrix score = Compute_index_with_LR_FE(X, beta, fe_N, fe_T, L, R);
    Map<MatrixXd> score_ (as <Map<MatrixXd>> (score));
    NumericVector score_beta = wrap(VectorXd::Constant(X.length(), 1.0));
    
    score_= Y_ -  score_; 
    
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

double Compute_obj_linear_with_theta(const List & X,  const NumericMatrix & Y, 
                                     const NumericVector & beta, const NumericMatrix & Theta, 
                                     const double sigma_sum, const double phi){
    NumericMatrix logl_linear = Compute_logl_linear_with_theta( X,  Y, beta, Theta );
    Map<MatrixXd> logl_linear_ (as <Map<MatrixXd>> (logl_linear));
    double loss = logl_linear_.sum();
    double obj  = -loss + phi * sigma_sum;
    return obj;
}

// [[Rcpp::export]]

double Compute_obj_linear_with_theta_FE(const List & X,  const NumericMatrix & Y, 
                                        const NumericVector & beta, 
                                        const NumericVector & fe_N, const NumericVector & fe_T,
                                        const NumericMatrix & Theta, 
                                        const double sigma_sum, const double phi){
    NumericMatrix logl_linear = Compute_logl_linear_with_theta_FE( X,  Y, beta, fe_N, fe_T, Theta );
    Map<MatrixXd> logl_linear_ (as <Map<MatrixXd>> (logl_linear));
    double loss = logl_linear_.sum();
    double obj  = -loss + phi * sigma_sum;
    return obj;
}

// [[Rcpp::export]]

double Compute_obj_linear_with_LR(const List & X,  const NumericMatrix & Y, 
                                  const NumericVector & beta, 
                                  const NumericMatrix & L, const NumericMatrix & R){
    NumericMatrix logl_linear = Compute_logl_linear_with_LR( X,  Y, beta, L, R );
    Map<MatrixXd> logl_linear_ (as <Map<MatrixXd>> (logl_linear));
    double loss = logl_linear_.sum();
    double obj  = -loss;
    return obj;
}

// [[Rcpp::export]]

double Compute_obj_linear_with_LR_FE(const List & X,  const NumericMatrix & Y, 
                                     const NumericVector & beta, 
                                     const NumericVector & fe_N, const NumericVector & fe_T,
                                     const NumericMatrix & L, const NumericMatrix & R){
    NumericMatrix logl_linear = Compute_logl_linear_with_LR_FE( X,  Y, beta,fe_N, fe_T, L, R );
    Map<MatrixXd> logl_linear_ (as <Map<MatrixXd>> (logl_linear));
    double loss = logl_linear_.sum();
    double obj  = -loss;
    return obj;
}

// [[Rcpp::export]]

NumericMatrix Compute_first_order_linear_with_LR(const List & X, const NumericMatrix & Y,
                                                 const NumericVector & beta, 
                                                 const NumericMatrix & L, const NumericMatrix & R){
    const Map<MatrixXd> Y_ (as <Map<MatrixXd>> (Y));
    const int N = Y_.rows();
    const int T = Y_.cols();

    NumericMatrix first_score = Compute_index_with_LR(X, beta, L, R);
    Map<MatrixXd> first_score_ (as <Map<MatrixXd>> (first_score));
    
    
    first_score_= Y_ -  first_score_; 
    
    return wrap(first_score_);
}



