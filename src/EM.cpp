#include "basics.h"
#include "mathfunc_logit.h"
#include "mathfunc_poisson.h"

// [[Rcpp::export]]

NumericMatrix Compute_E_step(const NumericMatrix & index,
                             const NumericMatrix & first_order, 
                             const NumericMatrix & second_order){
    
    const Map<MatrixXd> index_ (as<Map<MatrixXd>> (index));
    const Map<MatrixXd> first_order_ (as<Map<MatrixXd>> (first_order));
    const Map<MatrixXd> second_order_ (as<Map<MatrixXd>> (second_order));
    
    int N = index_.rows();
    int T = index_.cols();
    
    //MatrixXd second_order_stable_(second_order_);
    //second_order_stable_ = second_order_stable_.unaryExpr([](double x) { return x < 1e-3 ? 1e-3 : x; });
    
    //MatrixXd Mu_ = index_.array() + first_order_.array()/second_order_stable_.array();
    MatrixXd Mu_ = index_.array() + first_order_.array()/second_order_.array();
    for(int i=0; i!=N; i++){
        for(int t=0; t!=T; t++){
            if(Mu_(i, t) >=100.0){
                Mu_(i, t) = 100.0;
            }
            else if(Mu_(i, t) <= -100.0){
                Mu_(i, t) = -100.0;
            }
        }
    }
    return(wrap(Mu_));
}


// [[Rcpp::export]]

List Compute_M_step(NumericMatrix & Mu, 
                    const List & X, 
                    const NumericMatrix & L, 
                    const NumericMatrix & R, 
                    int num_factor){
    
    const Map<MatrixXd> Mu_ (as<Map<MatrixXd>> (Mu));
    const Map<MatrixXd> L_ (as<Map<MatrixXd>> (L));
    const Map<MatrixXd> R_ (as<Map<MatrixXd>> (R));
    
    //int N = Mu_.rows();
    //int T = Mu.cols();
    
    NumericVector beta_est(X.length());
    NumericMatrix X_vec = Convert_List_to_vector(X);
    Map<MatrixXd> X_vec_ (as<Map<MatrixXd>> (X_vec));
    
    MatrixXd theta_ = L_ * R_.transpose();
    MatrixXd theta_vec_ = Map<MatrixXd>(theta_.data(), theta_.size(), 1);
    
    MatrixXd beta_est_ = (X_vec_.transpose() * X_vec_).inverse() * X_vec_.transpose() * (Mu_ - theta_vec_);
    for(int i=0; i!=X.length(); i++){
        beta_est[i] = beta_est_(i,0);
    }
    
    NumericMatrix index = clone(Mu);
    Map<MatrixXd> index_ (as<Map<MatrixXd>> (index));
    for(int i=0; i!=X.length(); i++){
        Map<MatrixXd> X_ (as <Map<MatrixXd>> (X[i]));
        index_ = index_ - beta_est[i] * X_;
    }
    
    
    List list_low_rank = Low_rank_appro(index, num_factor);
    NumericMatrix L_est = as<NumericMatrix>(list_low_rank["L"]);
    NumericMatrix R_est = as<NumericMatrix>(list_low_rank["R"]);
    
    return List::create(Named("beta") = beta_est,
                        Named("L") = L_est,
                        Named("R") = R_est);
}


// [[Rcpp::export]]

NumericVector EM_logit(const NumericMatrix & Y, 
                       const List & X, 
                       const NumericVector  beta_0, 
                       const NumericMatrix & L_0, const NumericMatrix & R_0, 
                       const int num_factor, 
                       const int iter_max, const double tol){
    
    //const Map<MatrixXd> Y_(as<Map<MatrixXd>> (Y));
    //const Map<MatrixXd> L_0_(as<Map<MatrixXd>> (L_0));
    //const Map<MatrixXd> R_0_(as<Map<MatrixXd>> (R_0));
    
    NumericVector beta_new = clone(beta_0);
    NumericVector beta_old =clone(beta_0);
    NumericMatrix L = clone(L_0);
    NumericMatrix R = clone(R_0);
    
    NumericMatrix index = Compute_index_with_LR(X, beta_0, L_0, R_0);
    NumericMatrix first_order = Compute_first_order_logit_with_LR(X, Y, beta_0, L_0, R_0);
    NumericMatrix second_order = Compute_second_order_logit_with_LR(X, beta_0, L_0, R_0);
    
    int num_iter = 0;
    
    while(num_iter !=iter_max){
        //Rcout << num_iter <<std::endl;
        first_order = Compute_first_order_logit_with_LR(X, Y, beta_old, L, R);
        second_order = Compute_second_order_logit_with_LR(X, beta_old, L, R);
        index = Compute_E_step(index, first_order, second_order);
        Rcout << index(0,0) << std::endl;
        Rcout << index(0,5) << std::endl;
        
        //Rcout << index(0,0) << std::endl;
        for(int i=0; i!=1; i++){
            List list_M_step = Compute_M_step(index, X, L, R, num_factor);
            beta_new = as<NumericVector>(list_M_step["beta"]);
            L = as<NumericMatrix>(list_M_step["L"]);
            R = as<NumericMatrix>(list_M_step["R"]);
        }
        
        if(sqrt(sum((beta_new - beta_old) * (beta_new - beta_old))) > tol){
            beta_old = beta_new;
            Rcout << num_iter << std::endl;
            Rcout << beta_new << std::endl;
            num_iter = num_iter + 1;
        }
        else{
            break;
        }
    }
    
    return (beta_new);
}

    
    
    
    
    
    
    
    
    
    
