#include "basics.h"
#include "mathfunc_poisson.h"


// [[Rcpp::export]]

NumericMatrix Compute_logl_poisson_with_theta(const List & X,  const NumericMatrix & Y, 
                                              const NumericVector & beta, 
                                              const NumericMatrix & Theta ){
    
    const Map<MatrixXd> Y_ (as <Map<MatrixXd>> (Y));
    
    NumericMatrix logl_est = Compute_index_with_theta(X, beta, Theta);
    Map<MatrixXd> logl_est_ (as <Map<MatrixXd>> (logl_est));
    
    logl_est_ = logl_est_.array() * Y_.array() - logl_est_.array().exp();
    return wrap(logl_est_);
}

// [[Rcpp::export]]

NumericMatrix Compute_logl_poisson_with_theta_FE(const List & X,  const NumericMatrix & Y, 
                                                 const NumericVector & beta, 
                                                 const NumericVector & fe_N, const NumericVector & fe_T,
                                                 const NumericMatrix & Theta ){
    
    const Map<MatrixXd> Y_ (as <Map<MatrixXd>> (Y));
    
    NumericMatrix logl_est = Compute_index_with_theta_FE(X, beta, fe_N, fe_T, Theta);
    Map<MatrixXd> logl_est_ (as <Map<MatrixXd>> (logl_est));
    
    logl_est_ = logl_est_.array() * Y_.array() - logl_est_.array().exp();
    return wrap(logl_est_);
}


// [[Rcpp::export]]

NumericMatrix Compute_logl_poisson_with_LR(const List & X,  const NumericMatrix & Y, 
                                           const NumericVector & beta, 
                                           const NumericMatrix & L, const NumericMatrix & R ){
    
    const Map<MatrixXd> Y_ (as <Map<MatrixXd>> (Y));

    NumericMatrix logl_est = Compute_index_with_LR(X, beta, L, R);
    Map<MatrixXd> logl_est_ (as <Map<MatrixXd>> (logl_est));
    
    logl_est_ = logl_est_.array() * Y_.array() - logl_est_.array().exp();
    
    return wrap(logl_est_);
}


// [[Rcpp::export]]

NumericMatrix Compute_logl_poisson_with_LR_FE(const List & X,  const NumericMatrix & Y, 
                                              const NumericVector & beta, 
                                              const NumericVector & fe_N, const NumericVector & fe_T,
                                              const NumericMatrix & L, const NumericMatrix & R ){
    
    const Map<MatrixXd> Y_ (as <Map<MatrixXd>> (Y));
    
    NumericMatrix logl_est = Compute_index_with_LR_FE(X, beta, fe_N, fe_T, L, R);
    Map<MatrixXd> logl_est_ (as <Map<MatrixXd>> (logl_est));
    
    logl_est_ = logl_est_.array() * Y_.array() - logl_est_.array().exp();
    
    return wrap(logl_est_);
}

// [[Rcpp::export]]

List Compute_score_poisson_with_theta(const List & X, const NumericMatrix Y, 
                                    const NumericVector & beta, const NumericMatrix & Theta){
    const Map<MatrixXd> Y_ (as <Map<MatrixXd>> (Y));
    
    NumericMatrix score = Compute_index_with_theta(X, beta, Theta);
    Map<MatrixXd> score_ (as <Map<MatrixXd>> (score));
    
    NumericVector score_beta = wrap(VectorXd::Constant(X.length(), 1.0));
    
    score_ = Y_.array() - score_.array().exp();
    
    for(int i=0; i!=X.length(); i++){
        const Map<MatrixXd> X_ (as <Map<MatrixXd>> (X[i]));
        score_beta(i) = (score_.array() * X_.array()).sum();
    }
    return List::create(Named("score_beta") = score_beta, 
                        Named("score_theta") =  wrap(score_));
}

// [[Rcpp::export]]

List Compute_score_poisson_with_theta_FE(const List & X, const NumericMatrix Y, 
                                         const NumericVector & beta, 
                                         const NumericVector & fe_N, const NumericVector & fe_T,
                                         const NumericMatrix & Theta){
    const Map<MatrixXd> Y_ (as <Map<MatrixXd>> (Y));
    
    NumericMatrix score = Compute_index_with_theta_FE(X, beta, fe_N, fe_T, Theta);
    Map<MatrixXd> score_ (as <Map<MatrixXd>> (score));
    
    NumericVector score_beta = wrap(VectorXd::Constant(X.length(), 1.0));
    
    score_ = Y_.array() - score_.array().exp();
    
    for(int i=0; i!=X.length(); i++){
        const Map<MatrixXd> X_ (as <Map<MatrixXd>> (X[i]));
        score_beta(i) = (score_.array() * X_.array()).sum();
    }
    VectorXd score_fe_N_ = score_.rowwise().sum();
    VectorXd score_fe_T_ = score_.colwise().sum();
    
    return List::create(Named("score_beta") = score_beta, 
                        Named("score_fe_N") = wrap(score_fe_N_),
                        Named("score_fe_T") = wrap(score_fe_T_),
                        Named("score_theta") =  wrap(score_));
}

// [[Rcpp::export]]

List Compute_score_poisson_with_LR(const List & X, const NumericMatrix Y, 
                                 const NumericVector & beta, 
                                 const NumericMatrix & L, const NumericMatrix & R){
    
    const Map<MatrixXd> Y_ (as <Map<MatrixXd>> (Y));
    const Map<MatrixXd> L_ (as <Map<MatrixXd>> (L));
    const Map<MatrixXd> R_ (as <Map<MatrixXd>> (R));
    
    NumericMatrix score = Compute_index_with_LR(X, beta, L, R);
    Map<MatrixXd> score_ (as <Map<MatrixXd>> (score));
    
    NumericVector score_beta = wrap(VectorXd::Constant(X.length(), 1.0));
    
    score_ = Y_.array() - score_.array().exp();
    
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

List Compute_score_poisson_with_LR_FE(const List & X, const NumericMatrix Y, 
                                      const NumericVector & beta, 
                                      const NumericVector & fe_N, const NumericVector & fe_T,
                                      const NumericMatrix & L, const NumericMatrix & R){
    
    const Map<MatrixXd> Y_ (as <Map<MatrixXd>> (Y));
    const Map<MatrixXd> L_ (as <Map<MatrixXd>> (L));
    const Map<MatrixXd> R_ (as <Map<MatrixXd>> (R));
    
    //int N = Y_.rows();
    //int T = Y_.cols();
    
    NumericMatrix score = Compute_index_with_LR_FE(X, beta, fe_N, fe_T, L, R);
    Map<MatrixXd> score_ (as <Map<MatrixXd>> (score));
    
    NumericVector score_beta = wrap(VectorXd::Constant(X.length(), 1.0));
    
    score_ = Y_.array() - score_.array().exp();
    
    for(int i=0; i!=X.length(); i++){
        const Map<MatrixXd> X_ (as <Map<MatrixXd>> (X[i]));
        score_beta(i) = (score_.array() * X_.array()).sum();
    }
    MatrixXd score_L_ = score_ * R_;
    MatrixXd score_R_ = score_.transpose() * L_;
    
    VectorXd score_fe_N_ = score_.rowwise().sum();
    VectorXd score_fe_T_ = score_.colwise().sum();
    
    return List::create(Named("score_beta") = score_beta, 
                        Named("score_fe_N") =  wrap(score_fe_N_),
                        Named("score_fe_T") =  wrap(score_fe_T_),
                        Named("score_L") =  wrap(score_L_), 
                        Named("score_R") =  wrap(score_R_));
}    

// [[Rcpp::export]]

double Compute_obj_poisson_with_theta(const List & X,  const NumericMatrix & Y, 
                                      const NumericVector & beta, const NumericMatrix & Theta, 
                                      const double sigma_sum, const double phi){
    NumericMatrix logl_poisson = Compute_logl_poisson_with_theta( X,  Y, beta, Theta );
    Map<MatrixXd> logl_poisson_ (as <Map<MatrixXd>> (logl_poisson));
    double loss = logl_poisson_.sum();
    double obj  = -loss + phi * sigma_sum;
    return obj;
}

// [[Rcpp::export]]

double Compute_obj_poisson_with_theta_FE(const List & X,  const NumericMatrix & Y, 
                                         const NumericVector & beta, 
                                         const NumericVector & fe_N, const NumericVector & fe_T,
                                         const NumericMatrix & Theta, 
                                         const double sigma_sum, const double phi){
    NumericMatrix logl_poisson = Compute_logl_poisson_with_theta_FE( X,  Y, beta, fe_N, fe_T, Theta );
    Map<MatrixXd> logl_poisson_ (as <Map<MatrixXd>> (logl_poisson));
    double loss = logl_poisson_.sum();
    double obj  = -loss + phi * sigma_sum;
    return obj;
}


// [[Rcpp::export]]

double Compute_obj_poisson_with_LR(const List & X,  const NumericMatrix & Y, 
                                 const NumericVector & beta, 
                                 const NumericMatrix & L, const NumericMatrix & R){
    NumericMatrix logl_poisson = Compute_logl_poisson_with_LR( X,  Y, beta, L, R );
    Map<MatrixXd> logl_poisson_ (as <Map<MatrixXd>> (logl_poisson));
    double loss = logl_poisson_.sum();
    double obj  = -loss;
    return obj;
}

// [[Rcpp::export]]

double Compute_obj_poisson_with_LR_FE(const List & X,  const NumericMatrix & Y, 
                                      const NumericVector & beta,
                                      const NumericVector & fe_N, const NumericVector & fe_T,
                                      const NumericMatrix & L, const NumericMatrix & R){
    NumericMatrix logl_poisson = Compute_logl_poisson_with_LR_FE( X,  Y, beta, fe_N, fe_T, L, R );
    Map<MatrixXd> logl_poisson_ (as <Map<MatrixXd>> (logl_poisson));
    double loss = logl_poisson_.sum();
    double obj  = -loss;
    return obj;
}


// [[Rcpp::export]]

NumericMatrix Compute_first_order_poisson_with_LR(const List & X, const NumericMatrix & Y,
                                                  const NumericVector & beta, 
                                                  const NumericMatrix & L, const NumericMatrix & R){
    const Map<MatrixXd> Y_ (as <Map<MatrixXd>> (Y));
    
    NumericMatrix first_score = Compute_index_with_LR(X, beta, L, R);
    Map<MatrixXd> first_score_ (as <Map<MatrixXd>> (first_score));
    
    first_score_ = Y_.array() - first_score_.array().exp();

    return wrap(first_score_);
}


// [[Rcpp::export]]

NumericMatrix Compute_first_order_poisson_with_LR_FE(const List & X, const NumericMatrix & Y,
                                                     const NumericVector & beta, 
                                                     const NumericVector & fe_N, const NumericVector & fe_T,
                                                     const NumericMatrix & L, const NumericMatrix & R){
    const Map<MatrixXd> Y_ (as <Map<MatrixXd>> (Y));
    
    NumericMatrix first_score = Compute_index_with_LR_FE(X, beta, fe_N, fe_T, L, R);
    Map<MatrixXd> first_score_ (as <Map<MatrixXd>> (first_score));
    
    first_score_ = Y_.array() - first_score_.array().exp();
    
    return wrap(first_score_);
}


// [[Rcpp::export]]

NumericMatrix Compute_second_order_poisson_with_LR(const List & X,
                                                  const NumericVector & beta, 
                                                  const NumericMatrix & L, const NumericMatrix & R){
    
    NumericMatrix second_score = Compute_index_with_LR(X, beta, L, R);
    Map<MatrixXd> second_score_ (as <Map<MatrixXd>> (second_score));
    Map<MatrixXd> L_ (as <Map<MatrixXd>> (L));
    Map<MatrixXd> R_ (as <Map<MatrixXd>> (R));
    
    int N = L_.rows();
    int T = R_.rows();
    
    second_score_ = second_score_.array().exp();
    for(int i=0; i!= N; i++){
        for(int t=0; t!=T; t++ ){
            if(second_score_(i, t)>=1e10){
                second_score_(i, t) = 1e10;
            }
        }
    }
    
    return wrap(second_score_);
}


// [[Rcpp::export]]

NumericMatrix Compute_second_order_poisson_with_LR_FE(const List & X,
                                                      const NumericVector & beta, 
                                                      const NumericVector & fe_N, const NumericVector & fe_T,
                                                      const NumericMatrix & L, const NumericMatrix & R){
    
    NumericMatrix second_score = Compute_index_with_LR_FE(X, beta, fe_N, fe_T, L, R);
    Map<MatrixXd> second_score_ (as <Map<MatrixXd>> (second_score));
    //Map<MatrixXd> L_ (as <Map<MatrixXd>> (L));
    //Map<MatrixXd> R_ (as <Map<MatrixXd>> (R));
    
    //int N = L_.rows();
    //int T = R_.rows();
    
    second_score_ = second_score_.array().exp();
    
    return wrap(second_score_);
}


// [[Rcpp::export]]

NumericMatrix Sim_dynamic_poisson(const List Z, 
                                  const double beta_W, const NumericVector & beta_Z, 
                                  const NumericMatrix & Theta){
    
    Map<MatrixXd> Theta_ (as <Map<MatrixXd>> (Theta));
    
    int N = Theta_.rows();
    int T = Theta_.cols();
    
    NumericMatrix index = Compute_index_with_theta(Z, beta_Z, Theta);
    Map<MatrixXd> index_ (as <Map<MatrixXd>> (index));
    
    MatrixXd W_ = MatrixXd::Constant(N, T, 0.0);
    
    for(int i = 0; i!=N; i++){
        W_(i, 0) = index_(i, 0);
        W_(i,0) = R::rpois(std::exp(W_(i,0)));
        for(int t = 1; t!=T; t++){
            W_(i, t) = beta_W * std::log(1+  W_(i, t-1)) +  index_(i, t);
            W_(i, t) = R::rpois(std::exp(W_(i,t)));
        }
    }
    return wrap(W_);
}





