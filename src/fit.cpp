#include "basics.h"
#include "mathfunc_logit.h"
#include "mathfunc_poisson.h"
#include "mathfunc_probit.h"
#include "mathfunc_linear.h"



// [[Rcpp::export]]

NumericVector Update_beta(const NumericVector & beta, const NumericVector & score_beta, 
                          const double s){
    NumericVector beta_est = clone(beta);
    for(int i=0; i!=beta.size(); i++){
        beta_est[i] = beta[i] + s * score_beta[i];
    }
    return beta_est;
}


// [[Rcpp::export]]

List Update_theta(const NumericMatrix & Theta, const NumericMatrix & score_theta,
                          const double s, const double phi){
    const Map<MatrixXd> Theta_ (as<Map<MatrixXd>> (Theta));
    const Map<MatrixXd> score_theta_ (as<Map<MatrixXd>> (score_theta));
    
    MatrixXd A_est_ = Theta_ + s * score_theta_;
    NumericMatrix A_est = wrap(A_est_);
    
    List list_svd = SVD(A_est);
    MatrixXd U_ =  list_svd["U_"];
    MatrixXd V_ =  list_svd["V_"];
    VectorXd sigma_ = list_svd["sigma_"];
    NumericMatrix U = wrap(U_);
    NumericMatrix V = wrap(V_);
    NumericVector sigma = wrap(sigma_);
    
    List list_SVT = SVT(U, V, sigma, s * phi);
    
    return List::create(Named("Theta_est") = list_SVT["A"],
                        Named("sigma") = list_SVT["sigma"]);
}

// [[Rcpp::export]]

List Update_LR(const NumericMatrix & L, const NumericMatrix & R, 
               const NumericMatrix & score_L, const NumericMatrix & score_R,
               const double s_L, const double s_R){
    const Map<MatrixXd> L_ (as<Map<MatrixXd>> (L));
    const Map<MatrixXd> R_ (as<Map<MatrixXd>> (R));
    const Map<MatrixXd> score_L_ (as<Map<MatrixXd>> (score_L));
    const Map<MatrixXd> score_R_ (as<Map<MatrixXd>> (score_R));
    

    MatrixXd L_est_ = L_ + s_L * score_L_;
    MatrixXd R_est_ = R_ + s_R * score_R_;

    return List::create(Named("L_est") = wrap(L_est_),
                        Named("R_est") = wrap(R_est_));
}

// [[Rcpp::export]]

List Update_FE(const NumericVector & fe_N, const NumericVector & fe_T, 
               const NumericVector & score_fe_N, const NumericVector & score_fe_T,
               const double s_N, const double s_T){
    const Map<MatrixXd> fe_N_ (as<Map<MatrixXd>> (fe_N));
    const Map<MatrixXd> fe_T_ (as<Map<MatrixXd>> (fe_T));
    const Map<MatrixXd> score_fe_N_ (as<Map<MatrixXd>> (score_fe_N));
    const Map<MatrixXd> score_fe_T_ (as<Map<MatrixXd>> (score_fe_T));
    
    
    VectorXd fe_N_est_ = fe_N_ + s_N * score_fe_N_;
    MatrixXd fe_T_est_ = fe_T_ + s_T * score_fe_T_;
    
    return List::create(Named("fe_N_est") = wrap(fe_N_est_),
                        Named("fe_T_est") = wrap(fe_T_est_));
}


// [[Rcpp::export]]

List fit_logit(const List & X,  const NumericMatrix & Y, const double phi,  
               const double s, const int iter_max, const double tol){
    
    const Map<MatrixXd> Y_(as<Map<MatrixXd>> (Y));
    
    NumericMatrix Theta_new =  wrap(MatrixXd::Constant(Y_.rows(), Y_.cols(), 0));
    NumericMatrix Theta_old =  wrap(MatrixXd::Constant(Y_.rows(), Y_.cols(), 0));
    NumericVector beta_new = wrap(VectorXd::Constant(X.length(),0));
    NumericVector beta_old = wrap(VectorXd::Constant(X.length(),0));
    NumericMatrix score_theta = wrap(MatrixXd::Constant(Y_.rows(), Y_.cols(), 0));
    NumericVector score_beta =  wrap(VectorXd::Constant(X.length(),0));
    
    double sigma_sum = 0;
    
    double obj_old = Compute_obj_logit_with_theta(X, Y, beta_old, Theta_old, 
                                                  sigma_sum, phi);
    double obj_new = obj_old;
    
    int num_iter = 0;
    
    while(num_iter !=iter_max){
        double s_t = s;
        List list_score = Compute_score_logit_with_theta(X, Y, beta_old, Theta_old);
        score_beta = as<NumericVector>(list_score["score_beta"]);
        score_theta = as<NumericMatrix>(list_score["score_theta"]);
        while(1){
            beta_new = Update_beta(beta_old, score_beta, s_t/(Y_.size()));
            List list_theta = Update_theta(Theta_old, score_theta, s_t, phi);
            Theta_new = as<NumericMatrix>(list_theta["Theta_est"]);
            const Map<VectorXd> sigma_(as<Map<VectorXd>> 
                                       (as<NumericVector>(list_theta["sigma"])));
            sigma_sum = sigma_.sum();
            obj_new = Compute_obj_logit_with_theta(X, Y, beta_new, Theta_new, 
                                                   sigma_sum, phi);
            //Rcout << obj_old<< std::endl;
            //Rcout << obj_new<< std::endl;
            if (obj_new > obj_old){
                s_t = s_t/2;
                //Rcout << s_t << std::endl;
            } 
            else{
                beta_old =  beta_new;
                Theta_old =  Theta_new;
                num_iter = num_iter + 1;
                Rcout << "NNR - Iter " << num_iter <<": Coefficiet = " << beta_new << std::endl;
                //Rcout << "Hi\n" << std::endl;
                break;
            }
            
        }
        if(obj_new < (1-tol) * obj_old){
            obj_old = obj_new;
            //Rcout << obj_new << std::endl;
        }
        else{
            break;
        }
        
    }
    return List::create(Named("beta_est") = beta_new,
                        Named("Theta_est") = Theta_new );
}


// [[Rcpp::export]]

List MLE_logit(const List & X,  const NumericMatrix & Y, 
               const NumericVector  beta_0, 
               const NumericMatrix & L_0, const NumericMatrix & R_0, 
               const double s, const int iter_max, const double tol){
    
    const Map<MatrixXd> Y_(as<Map<MatrixXd>> (Y));
    const Map<MatrixXd> L_0_(as<Map<MatrixXd>> (L_0));
    const Map<MatrixXd> R_0_(as<Map<MatrixXd>> (R_0));
    
    NumericVector beta_new = clone(beta_0);
    NumericVector beta_old =clone(beta_0);
    NumericMatrix L_new = clone(L_0);
    NumericMatrix L_old = clone(L_0);
    NumericMatrix R_new = clone(R_0);
    NumericMatrix R_old = clone(R_0);
    
    NumericMatrix score_L = wrap(MatrixXd::Constant(L_0_.rows(), L_0_.cols(), 0));
    NumericMatrix score_R = wrap(MatrixXd::Constant(R_0_.rows(), R_0_.cols(), 0));
    NumericVector score_beta =  wrap(VectorXd::Constant(X.length(),0));
    
    double obj_old = Compute_obj_logit_with_LR(X, Y, beta_old, L_old, R_old);
    double obj_new = obj_old;
    
    int num_iter = 0;
    
    while(num_iter !=iter_max){
        double s_t = s;
        List list_score = Compute_score_logit_with_LR(X, Y, beta_old, L_old, R_old);
        score_beta = as<NumericVector>(list_score["score_beta"]);
        score_L = as<NumericMatrix>(list_score["score_L"]);
        score_R = as<NumericMatrix>(list_score["score_R"]);
        
        while(1){
            beta_new = Update_beta(beta_old, score_beta, s_t/(Y_.size()));
            List list_LR = Update_LR(L_old, R_old, score_L, score_R, s_t/(Y_.cols()), s_t/(Y_.rows()));
            L_new = as<NumericMatrix>(list_LR["L_est"]);
            R_new = as<NumericMatrix>(list_LR["R_est"]);
            //Rcout << 1 << std::endl;
            //Rcout << 2 << std::endl;
            
            obj_new = Compute_obj_logit_with_LR(X, Y, beta_new, L_new, R_new);
            //Rcout << obj_old<< std::endl;
            //Rcout << obj_new<< std::endl;
            if (obj_new > obj_old && s_t >=1e-5){
                s_t = s_t/2;
                //Rcout << s_t << std::endl;
            } 
            else{
                beta_old =  beta_new;
                L_old = L_new;
                R_old = R_new;
                num_iter = num_iter + 1;
                Rcout << "MLE - Iter " << num_iter <<": Coefficiet = " << beta_new << std::endl;
                //Rcout << "Hi\n" << std::endl;
                break;
            }
        }

        if(obj_new < (1-tol) * obj_old){
            //Rcout << obj_old << std::endl;
            obj_old = obj_new;
        }
        else{
            break;
        }
    }
    return List::create(Named("beta_est") = beta_new,
                        Named("L_est") = L_new, 
                        Named("R_est") = R_new);
}

// [[Rcpp::export]]

List fit_poisson(const List & X,  const NumericMatrix & Y, const double phi,  
               const double s, const int iter_max, const double tol, 
               NumericVector beta_init = NumericVector(0)){
    
    if (beta_init.size()==0) {
        beta_init = NumericVector(X.length());
    }
    
    const Map<MatrixXd> Y_(as<Map<MatrixXd>> (Y));
    
    NumericMatrix Theta_new =  wrap(MatrixXd::Constant(Y_.rows(), Y_.cols(), 0));
    NumericMatrix Theta_old =  wrap(MatrixXd::Constant(Y_.rows(), Y_.cols(), 0));
    NumericVector beta_new = clone(beta_init);
    NumericVector beta_old = clone(beta_init);
    NumericMatrix score_theta = wrap(MatrixXd::Constant(Y_.rows(), Y_.cols(), 0));
    NumericVector score_beta =  wrap(VectorXd::Constant(X.length(),0));
    
    double sigma_sum = 0;
    
    double obj_old = Compute_obj_poisson_with_theta(X, Y, beta_old, Theta_old, 
                                                  sigma_sum, phi);
    double obj_new = obj_old;
    
    int num_iter = 0;
    
    while(num_iter !=iter_max){
        double s_t = s;
        List list_score = Compute_score_poisson_with_theta(X, Y, beta_old, Theta_old);
        score_beta = as<NumericVector>(list_score["score_beta"]);
        score_theta = as<NumericMatrix>(list_score["score_theta"]);
        while(1){
            beta_new = Update_beta(beta_old, score_beta, s_t/(Y_.size()));
            List list_theta = Update_theta(Theta_old, score_theta, s_t, phi);
            Theta_new = as<NumericMatrix>(list_theta["Theta_est"]);
            const Map<VectorXd> sigma_(as<Map<VectorXd>> 
                                       (as<NumericVector>(list_theta["sigma"])));
            sigma_sum = sigma_.sum();
            obj_new = Compute_obj_poisson_with_theta(X, Y, beta_new, Theta_new, 
                                                   sigma_sum, phi);
            //Rcout << obj_old<< std::endl;
            //Rcout << obj_new<< std::endl;
            if (obj_new > obj_old && s_t >=1e-10){
                s_t = s_t/2;
                //Rcout << s_t << std::endl;
            } 
            else{
                beta_old =  beta_new;
                Theta_old =  Theta_new;
                num_iter = num_iter + 1;
                Rcout << "NNR - Iter " << num_iter <<": Coefficiet = " << beta_new << std::endl;
                //Rcout << "Hi\n" << std::endl;
                break;
            }
            
        }
        //Rcout << obj_old << std::endl;
        //Rcout << obj_new << std::endl;
        //Rcout << abs(obj_new - obj_old)  /  abs(obj_old) << std::endl;
        //Rcout << abs(obj_new - obj_old)  / abs(obj_new)  <<std::endl;
        if( abs(obj_new - obj_old) > tol *  abs(obj_old) && abs(obj_new - obj_old) > tol *  abs(obj_new) ){
            obj_old = obj_new;
        }
        else{
            break;
        }
        
    }
    return List::create(Named("beta_est") = beta_new,
                        Named("Theta_est") = Theta_new );
}



// [[Rcpp::export]]

List fit_poisson_FE(const List & X,  const NumericMatrix & Y, const double phi,  
                    const double s, const int iter_max, const double tol, 
                    NumericVector beta_init = NumericVector(0)){
    
    if (beta_init.size()==0) {
        beta_init = NumericVector(X.length());
    }
    
    const Map<MatrixXd> Y_(as<Map<MatrixXd>> (Y));
    
    int N = Y.rows();
    int T = Y.cols();
    
    NumericMatrix Theta_new =  wrap(MatrixXd::Constant(Y_.rows(), Y_.cols(), 0));
    NumericMatrix Theta_old =  wrap(MatrixXd::Constant(Y_.rows(), Y_.cols(), 0));
    NumericVector beta_new = clone(beta_init);
    NumericVector beta_old = clone(beta_init);
    NumericVector fe_N_old =  wrap(VectorXd::Constant(Y_.rows(),0));
    NumericVector fe_N_new =  wrap(VectorXd::Constant(Y_.rows(),0));
    NumericVector fe_T_old =  wrap(VectorXd::Constant(Y_.cols(),0));
    NumericVector fe_T_new =  wrap(VectorXd::Constant(Y_.cols(),0));
    NumericMatrix score_theta = wrap(MatrixXd::Constant(Y_.rows(), Y_.cols(), 0));
    NumericVector score_beta =  wrap(VectorXd::Constant(X.length(),0));
    NumericVector score_fe_N =  wrap(VectorXd::Constant(Y_.rows(),0));
    NumericVector score_fe_T =  wrap(VectorXd::Constant(Y_.cols(),0));
    
    
    double sigma_sum = 0;
    
    double obj_old = Compute_obj_poisson_with_theta_FE(X, Y, beta_old, 
                                                       fe_N_old, fe_T_old,
                                                       Theta_old, 
                                                       sigma_sum, phi);
    double obj_new = obj_old;
    
    int num_iter = 0;
    
    while(num_iter !=iter_max){
        double s_t = s;
        List list_score = Compute_score_poisson_with_theta_FE(X, Y, beta_old, fe_N_old, fe_T_old, Theta_old);
        score_beta = as<NumericVector>(list_score["score_beta"]);
        score_fe_N = as<NumericVector>(list_score["score_fe_N"]);
        score_fe_T = as<NumericVector>(list_score["score_fe_T"]);
        score_theta = as<NumericMatrix>(list_score["score_theta"]);
        while(1){
            beta_new = Update_beta(beta_old, score_beta, s_t/(Y_.size()));
            List list_fe = Update_FE(fe_N_old, fe_T_old, score_fe_N, score_fe_T, s_t/T, s_t/N);
            fe_N_new = as<NumericVector>(list_fe["fe_N_est"]);
            fe_T_new = as<NumericVector>(list_fe["fe_T_est"]);
            List list_theta = Update_theta(Theta_old, score_theta, s_t, phi);
            Theta_new = as<NumericMatrix>(list_theta["Theta_est"]);
            const Map<VectorXd> sigma_(as<Map<VectorXd>> 
                                       (as<NumericVector>(list_theta["sigma"])));
            sigma_sum = sigma_.sum();
            obj_new = Compute_obj_poisson_with_theta_FE(X, Y, beta_new, 
                                                        fe_N_new, fe_T_new, Theta_new, 
                                                        sigma_sum, phi);
            //Rcout << obj_old<< std::endl;
            //Rcout << obj_new<< std::endl;
            if (obj_new > obj_old && s_t >=1e-10){
                s_t = s_t/2;
                //Rcout << s_t << std::endl;
            } 
            else{
                beta_old =  beta_new;
                Theta_old =  Theta_new;
                fe_N_old = fe_N_new;
                fe_T_old = fe_T_new;
                num_iter = num_iter + 1;
                Rcout << beta_new << std::endl;
                //Rcout << "Hi\n" << std::endl;
                break;
            }
            
        }
        if( abs(obj_new - obj_old) > tol *  abs(obj_old) && abs(obj_new - obj_old) > tol *  abs(obj_new) ){
            obj_old = obj_new;
        }
        else{
            break;
        }
        
    }
    return List::create(Named("beta_est") = beta_new,
                        Named("fe_N_est") = fe_N_new,
                        Named("fe_T_est") = fe_T_new,
                        Named("Theta_est") = Theta_new );
}


// [[Rcpp::export]]

List fit_poisson_FE_w(const List & X,  const NumericMatrix & Y, const double phi,  
                      const double s, const int iter_max, const double tol,
                      const NumericVector fe_N_0, const NumericVector fe_T_0, 
                      NumericVector beta_init = NumericVector(0)){
    
    if (beta_init.size()==0) {
        beta_init = NumericVector(X.length());
    }
    
    const Map<MatrixXd> Y_(as<Map<MatrixXd>> (Y));
    
    int N = Y.rows();
    int T = Y.cols();
    
    NumericMatrix Theta_new =  wrap(MatrixXd::Constant(Y_.rows(), Y_.cols(), 0));
    NumericMatrix Theta_old =  wrap(MatrixXd::Constant(Y_.rows(), Y_.cols(), 0));
    NumericVector beta_new = clone(beta_init);
    NumericVector beta_old = clone(beta_init);
    NumericVector fe_N_old =  clone(fe_N_0);
    NumericVector fe_N_new =  clone(fe_N_0);
    NumericVector fe_T_old =  clone(fe_T_0);
    NumericVector fe_T_new =  clone(fe_T_0);

    //NumericVector fe_N_old =  wrap(VectorXd::Constant(Y_.rows(),0));
    //NumericVector fe_N_new =  wrap(VectorXd::Constant(Y_.rows(),0));
    //NumericVector fe_T_old =  wrap(VectorXd::Constant(Y_.cols(),0));
    //NumericVector fe_T_new =  wrap(VectorXd::Constant(Y_.cols(),0));
    NumericMatrix score_theta = wrap(MatrixXd::Constant(Y_.rows(), Y_.cols(), 0));
    NumericVector score_beta =  wrap(VectorXd::Constant(X.length(),0));
    NumericVector score_fe_N =  wrap(VectorXd::Constant(Y_.rows(),0));
    NumericVector score_fe_T =  wrap(VectorXd::Constant(Y_.cols(),0));
    
    
    double sigma_sum = 0;
    
    double obj_old = Compute_obj_poisson_with_theta_FE(X, Y, beta_old, 
                                                       fe_N_old, fe_T_old,
                                                       Theta_old, 
                                                       sigma_sum, phi);
    double obj_new = obj_old;
    
    int num_iter = 0;
    
    while(num_iter !=iter_max){
        double s_t = s;
        List list_score = Compute_score_poisson_with_theta_FE(X, Y, beta_old, fe_N_old, fe_T_old, Theta_old);
        score_beta = as<NumericVector>(list_score["score_beta"]);
        score_fe_N = as<NumericVector>(list_score["score_fe_N"]);
        score_fe_T = as<NumericVector>(list_score["score_fe_T"]);
        score_theta = as<NumericMatrix>(list_score["score_theta"]);
        while(1){
            beta_new = Update_beta(beta_old, score_beta, s_t/(Y_.size()));
            List list_fe = Update_FE(fe_N_old, fe_T_old, score_fe_N, score_fe_T, s_t/T, s_t/N);
            fe_N_new = as<NumericVector>(list_fe["fe_N_est"]);
            fe_T_new = as<NumericVector>(list_fe["fe_T_est"]);
            List list_theta = Update_theta(Theta_old, score_theta, s_t, phi);
            Theta_new = as<NumericMatrix>(list_theta["Theta_est"]);
            const Map<VectorXd> sigma_(as<Map<VectorXd>> 
                                       (as<NumericVector>(list_theta["sigma"])));
            sigma_sum = sigma_.sum();
            obj_new = Compute_obj_poisson_with_theta_FE(X, Y, beta_new, 
                                                        fe_N_new, fe_T_new, Theta_new, 
                                                        sigma_sum, phi);
            //Rcout << obj_old<< std::endl;
            //Rcout << obj_new<< std::endl;
            if (obj_new > obj_old && s_t >=1e-8){
                s_t = s_t/2;
                //Rcout << s_t << std::endl;
            } 
            else{
                beta_old =  beta_new;
                //Theta_old =  Theta_new;
                fe_N_old = fe_N_new;
                fe_T_old = fe_T_new;
                num_iter = num_iter + 1;
                Rcout << beta_new << std::endl;
                //Rcout << "Hi\n" << std::endl;
                break;
            }
            
        }
        if( abs(obj_new - obj_old) > tol *  abs(obj_old) && abs(obj_new - obj_old) > tol *  abs(obj_new) ){
            obj_old = obj_new;
        }
        else{
            break;
        }
        
    }
    return List::create(Named("beta_est") = beta_new,
                        Named("fe_N_est") = fe_N_new,
                        Named("fe_T_est") = fe_T_new,
                        Named("Theta_est") = Theta_new );
}


// [[Rcpp::export]]

List MLE_poisson(const List & X,  const NumericMatrix & Y, 
               const NumericVector  beta_0, 
               const NumericMatrix & L_0, const NumericMatrix & R_0, 
               const double s, const int iter_max, const double tol){
    
    const Map<MatrixXd> Y_(as<Map<MatrixXd>> (Y));
    const Map<MatrixXd> L_0_(as<Map<MatrixXd>> (L_0));
    const Map<MatrixXd> R_0_(as<Map<MatrixXd>> (R_0));
    
    NumericVector beta_new = clone(beta_0);
    NumericVector beta_old =clone(beta_0);
    NumericMatrix L_new = clone(L_0);
    NumericMatrix L_old = clone(L_0);
    NumericMatrix R_new = clone(R_0);
    NumericMatrix R_old = clone(R_0);
    
    NumericMatrix score_L = wrap(MatrixXd::Constant(L_0_.rows(), L_0_.cols(), 0));
    NumericMatrix score_R = wrap(MatrixXd::Constant(R_0_.rows(), R_0_.cols(), 0));
    NumericVector score_beta =  wrap(VectorXd::Constant(X.length(),0));
    
    double obj_old = Compute_obj_poisson_with_LR(X, Y, beta_old, L_old, R_old);
    double obj_new = obj_old;
    
    int num_iter = 0;
    
    while(num_iter !=iter_max){
        double s_t = s;
        List list_score = Compute_score_poisson_with_LR(X, Y, beta_old, L_old, R_old);
        score_beta = as<NumericVector>(list_score["score_beta"]);
        score_L = as<NumericMatrix>(list_score["score_L"]);
        score_R = as<NumericMatrix>(list_score["score_R"]);
        
        while(1){
            beta_new = Update_beta(beta_old, score_beta, s_t/(Y_.size()));
            List list_LR = Update_LR(L_old, R_old, score_L, score_R, s_t/(Y_.cols()), s_t/(Y_.rows()));
            L_new = as<NumericMatrix>(list_LR["L_est"]);
            R_new = as<NumericMatrix>(list_LR["R_est"]);
            //Rcout << L_old(0,0) << std::endl;
            //Rcout << L_new(0,0) << std::endl;
            
            obj_new = Compute_obj_poisson_with_LR(X, Y, beta_new, L_new, R_new);
            //Rcout << obj_old<< std::endl;
            //Rcout << obj_new<< std::endl;
            if (obj_new > obj_old && s_t >=1e-10){
                s_t = s_t/2;
                //Rcout << obj_old << std::endl;
                //Rcout << obj_new << std::endl;
            } 
            else{
                beta_old =  beta_new;
                L_old = L_new;
                R_old = R_new;
                num_iter = num_iter + 1;
                Rcout << "MLE - Iter " << num_iter <<": Coefficiet = " << beta_new << std::endl;
                //Rcout << "Hi\n" << std::endl;
                break;
            }
        }
        
        if( abs(obj_new - obj_old) > tol *  abs(obj_old) && abs(obj_new - obj_old) > tol *  abs(obj_new) ){
            obj_old = obj_new;
        }
        else{
            break;
        }
    }
    return List::create(Named("beta_est") = beta_new,
                        Named("L_est") = L_new, 
                        Named("R_est") = R_new);
}





// [[Rcpp::export]]

List MLE_poisson_FE(const List & X,  const NumericMatrix & Y, 
                    const NumericVector  beta_0, 
                    const NumericVector & fe_N_0, const NumericVector & fe_T_0,
                    const NumericMatrix & L_0, const NumericMatrix & R_0, 
                    const double s, const int iter_max, const double tol){
    
    const Map<MatrixXd> Y_(as<Map<MatrixXd>> (Y));
    const Map<MatrixXd> L_0_(as<Map<MatrixXd>> (L_0));
    const Map<MatrixXd> R_0_(as<Map<MatrixXd>> (R_0));
    //const Map<VectorXd> fe_N_0_(as<Map<VectorXd>> (fe_N_0));
    //const Map<VectorXd> fe_T_0_(as<Map<VectorXd>> (fe_T_0));
    
    int N = Y_.rows();
    int T = Y_.cols();
    
    NumericVector beta_new = clone(beta_0);
    NumericVector beta_old =clone(beta_0);
    NumericMatrix L_new = clone(L_0);
    NumericMatrix L_old = clone(L_0);
    NumericMatrix R_new = clone(R_0);
    NumericMatrix R_old = clone(R_0);
    NumericVector fe_N_old = clone(fe_N_0);
    NumericVector fe_N_new = clone(fe_N_0);
    NumericVector fe_T_old = clone(fe_T_0);
    NumericVector fe_T_new = clone(fe_T_0);
    
    
    NumericMatrix score_L = wrap(MatrixXd::Constant(L_0_.rows(), L_0_.cols(), 0));
    NumericMatrix score_R = wrap(MatrixXd::Constant(R_0_.rows(), R_0_.cols(), 0));
    NumericVector score_beta =  wrap(VectorXd::Constant(X.length(),0));
    NumericVector score_fe_N =  wrap(VectorXd::Constant(Y_.rows(),0));
    NumericVector score_fe_T =  wrap(VectorXd::Constant(Y_.cols(),0));
    
    
    double obj_old = Compute_obj_poisson_with_LR_FE(X, Y, beta_old, fe_N_old, fe_T_old, L_old, R_old);
    double obj_new = obj_old;
    
    int num_iter = 0;
    
    while(num_iter !=iter_max){
        double s_t = s;
        List list_score = Compute_score_poisson_with_LR_FE(X, Y, beta_old, fe_N_old, fe_T_old, L_old, R_old);
        score_beta = as<NumericVector>(list_score["score_beta"]);
        score_L = as<NumericMatrix>(list_score["score_L"]);
        score_R = as<NumericMatrix>(list_score["score_R"]);
        score_fe_N = as<NumericVector>(list_score["score_fe_N"]);
        score_fe_T = as<NumericVector>(list_score["score_fe_T"]);
        
        while(1){
            beta_new = Update_beta(beta_old, score_beta, s_t/(Y_.size()));
            List list_fe = Update_FE(fe_N_old, fe_T_old, score_fe_N, score_fe_T, s_t/T, s_t/N);
            fe_N_new = as<NumericVector>(list_fe["fe_N_est"]);
            fe_T_new = as<NumericVector>(list_fe["fe_T_est"]);
            List list_LR = Update_LR(L_old, R_old, score_L, score_R, s_t/(Y_.cols()), s_t/(Y_.rows()));
            L_new = as<NumericMatrix>(list_LR["L_est"]);
            R_new = as<NumericMatrix>(list_LR["R_est"]);
            //Rcout << L_old(0,0) << std::endl;
            //Rcout << L_new(0,0) << std::endl;
            
            obj_new = Compute_obj_poisson_with_LR_FE(X, Y, beta_new, fe_N_new, fe_T_new, L_new, R_new);
            //Rcout << obj_old<< std::endl;
            //Rcout << obj_new<< std::endl;
            if (obj_new > obj_old && s_t >=1e-12){
                s_t = s_t/2;
                //Rcout << obj_old << std::endl;
                //Rcout << obj_new << std::endl;
            } 
            else{
                beta_old =  beta_new;
                L_old = L_new;
                R_old = R_new;
                fe_N_old = fe_N_new;
                fe_T_old = fe_T_new;
                num_iter = num_iter + 1;
                Rcout << beta_new << std::endl;
                //Rcout << "Hi\n" << std::endl;
                break;
            }
        }
        
        if( abs(obj_new - obj_old) > tol *  abs(obj_old) && abs(obj_new - obj_old) > tol *  abs(obj_new) ){
            obj_old = obj_new;
        }
        else{
            break;
        }
    }
    return List::create(Named("beta_est") = beta_new,
                        Named("fe_N_est") = fe_N_new, 
                        Named("fe_T_est") = fe_T_new, 
                        Named("L_est") = L_new, 
                        Named("R_est") = R_new);
}


// [[Rcpp::export]]

List MLE_poisson_LR(const List & X,  const NumericMatrix & Y, 
                 const NumericVector  beta_0, 
                 const NumericMatrix & L_0, const NumericMatrix & R_0, 
                 const double s, const int iter_max, const double tol){
    
    const Map<MatrixXd> Y_(as<Map<MatrixXd>> (Y));
    const Map<MatrixXd> L_0_(as<Map<MatrixXd>> (L_0));
    const Map<MatrixXd> R_0_(as<Map<MatrixXd>> (R_0));
    
    NumericMatrix L_new = clone(L_0);
    NumericMatrix L_old = clone(L_0);
    NumericMatrix R_new = clone(R_0);
    NumericMatrix R_old = clone(R_0);
    
    NumericMatrix score_L = wrap(MatrixXd::Constant(L_0_.rows(), L_0_.cols(), 0));
    NumericMatrix score_R = wrap(MatrixXd::Constant(R_0_.rows(), R_0_.cols(), 0));

    double obj_old = Compute_obj_poisson_with_LR(X, Y, beta_0, L_old, R_old);
    double obj_new = obj_old;
    
    int num_iter = 0;
    
    while(num_iter !=iter_max){
        double s_t = s;
        List list_score = Compute_score_poisson_with_LR(X, Y, beta_0, L_old, R_old);
        score_L = as<NumericMatrix>(list_score["score_L"]);
        score_R = as<NumericMatrix>(list_score["score_R"]);
        
        while(1){
            List list_LR = Update_LR(L_old, R_old, score_L, score_R, s_t/(Y_.cols()), s_t/(Y_.rows()));
            L_new = as<NumericMatrix>(list_LR["L_est"]);
            R_new = as<NumericMatrix>(list_LR["R_est"]);
            //Rcout << L_old(0,0) << std::endl;
            //Rcout << L_new(0,0) << std::endl;
            
            obj_new = Compute_obj_poisson_with_LR(X, Y, beta_0, L_new, R_new);
            //Rcout << obj_old<< std::endl;
            //Rcout << obj_new<< std::endl;
            if (obj_new > obj_old && s_t >=1e-5){
                s_t = s_t/2;
                //Rcout << obj_old << std::endl;
                //Rcout << obj_new << std::endl;
            } 
            else{
                L_old = L_new;
                R_old = R_new;
                num_iter = num_iter + 1;
                Rcout << beta_0 << std::endl;
                //Rcout << "Hi\n" << std::endl;
                break;
            }
        }
        
        if( abs(obj_new - obj_old) > tol *  abs(obj_old) && abs(obj_new - obj_old) > tol *  abs(obj_new) ){
            obj_old = obj_new;
        }
        else{
            break;
        }
    }
    return List::create(Named("beta_est") = beta_0,
                        Named("L_est") = L_new, 
                        Named("R_est") = R_new);
}



// [[Rcpp::export]]

List fit_probit(const List & X,  const NumericMatrix & Y, const double phi,  
                const double s, const int iter_max, const double tol){
    
    const Map<MatrixXd> Y_(as<Map<MatrixXd>> (Y));
    
    NumericMatrix Theta_new =  wrap(MatrixXd::Constant(Y_.rows(), Y_.cols(), 0));
    NumericMatrix Theta_old =  wrap(MatrixXd::Constant(Y_.rows(), Y_.cols(), 0));
    NumericVector beta_new = wrap(VectorXd::Constant(X.length(),0));
    NumericVector beta_old = wrap(VectorXd::Constant(X.length(),0));
    NumericMatrix score_theta = wrap(MatrixXd::Constant(Y_.rows(), Y_.cols(), 0));
    NumericVector score_beta =  wrap(VectorXd::Constant(X.length(),0));
    
    double sigma_sum = 0;
    
    double obj_old = Compute_obj_probit_with_theta(X, Y, beta_old, Theta_old, 
                                                  sigma_sum, phi);
    double obj_new = obj_old;
    
    int num_iter = 0;
    
    while(num_iter !=iter_max){
        double s_t = s;
        List list_score = Compute_score_probit_with_theta(X, Y, beta_old, Theta_old);
        score_beta = as<NumericVector>(list_score["score_beta"]);
        score_theta = as<NumericMatrix>(list_score["score_theta"]);
        while(1){
            beta_new = Update_beta(beta_old, score_beta, s_t/(Y_.size()));
            List list_theta = Update_theta(Theta_old, score_theta, s_t, phi);
            Theta_new = as<NumericMatrix>(list_theta["Theta_est"]);
            const Map<VectorXd> sigma_(as<Map<VectorXd>> 
                                       (as<NumericVector>(list_theta["sigma"])));
            sigma_sum = sigma_.sum();
            obj_new = Compute_obj_probit_with_theta(X, Y, beta_new, Theta_new, 
                                                    sigma_sum, phi);
            //Rcout << obj_old<< std::endl;
            //Rcout << obj_new<< std::endl;
            if (obj_new > obj_old){
                s_t = s_t/2;
                //Rcout << s_t << std::endl;
            } 
            else{
                beta_old =  beta_new;
                Theta_old =  Theta_new;
                num_iter = num_iter + 1;
                Rcout << "NNR - Iter " << num_iter <<": Coefficiet = " << beta_new << std::endl;
                //Rcout << "Hi\n" << std::endl;
                break;
            }
            
        }
        if(obj_new < (1-tol) * obj_old){
            obj_old = obj_new;
            //Rcout << obj_new << std::endl;
        }
        else{
            break;
        }
        
    }
    return List::create(Named("beta_est") = beta_new,
                        Named("Theta_est") = Theta_new );
}



// [[Rcpp::export]]

List MLE_probit(const List & X,  const NumericMatrix & Y, 
                const NumericVector  beta_0, 
                const NumericMatrix & L_0, const NumericMatrix & R_0, 
                const double s, const int iter_max, const double tol){
    
    const Map<MatrixXd> Y_(as<Map<MatrixXd>> (Y));
    const Map<MatrixXd> L_0_(as<Map<MatrixXd>> (L_0));
    const Map<MatrixXd> R_0_(as<Map<MatrixXd>> (R_0));
    
    NumericVector beta_new = clone(beta_0);
    NumericVector beta_old =clone(beta_0);
    NumericMatrix L_new = clone(L_0);
    NumericMatrix L_old = clone(L_0);
    NumericMatrix R_new = clone(R_0);
    NumericMatrix R_old = clone(R_0);
    
    NumericMatrix score_L = wrap(MatrixXd::Constant(L_0_.rows(), L_0_.cols(), 0));
    NumericMatrix score_R = wrap(MatrixXd::Constant(R_0_.rows(), R_0_.cols(), 0));
    NumericVector score_beta =  wrap(VectorXd::Constant(X.length(),0));
    
    double obj_old = Compute_obj_probit_with_LR(X, Y, beta_old, L_old, R_old);
    double obj_new = obj_old;
    
    int num_iter = 0;
    
    while(num_iter !=iter_max){
        double s_t = s;
        List list_score = Compute_score_probit_with_LR(X, Y, beta_old, L_old, R_old);
        score_beta = as<NumericVector>(list_score["score_beta"]);
        score_L = as<NumericMatrix>(list_score["score_L"]);
        score_R = as<NumericMatrix>(list_score["score_R"]);
        
        while(1){
            beta_new = Update_beta(beta_old, score_beta, s_t/(Y_.size()));
            List list_LR = Update_LR(L_old, R_old, score_L, score_R, s_t/(Y_.cols()), s_t/(Y_.rows()));
            L_new = as<NumericMatrix>(list_LR["L_est"]);
            R_new = as<NumericMatrix>(list_LR["R_est"]);
            //Rcout << 1 << std::endl;
            //Rcout << 2 << std::endl;
            
            obj_new = Compute_obj_probit_with_LR(X, Y, beta_new, L_new, R_new);
            //Rcout << obj_old<< std::endl;
            //Rcout << obj_new<< std::endl;
            if (obj_new > obj_old && s_t >=1e-5){
                s_t = s_t/2;
                //Rcout << s_t << std::endl;
            } 
            else{
                beta_old =  beta_new;
                L_old = L_new;
                R_old = R_new;
                num_iter = num_iter + 1;
                Rcout << "MLE - Iter " << num_iter <<": Coefficiet = " << beta_new << std::endl;
                //Rcout << "Hi\n" << std::endl;
                break;
            }
        }
        
        if(obj_new < (1-tol) * obj_old){
            //Rcout << obj_old << std::endl;
            obj_old = obj_new;
        }
        else{
            break;
        }
    }
    return List::create(Named("beta_est") = beta_new,
                        Named("L_est") = L_new, 
                        Named("R_est") = R_new);
}


// [[Rcpp::export]]

List fit_linear(const List & X,  const NumericMatrix & Y, const double phi,  
                const double s, const int iter_max, const double tol){
    
    const Map<MatrixXd> Y_(as<Map<MatrixXd>> (Y));
    
    NumericMatrix Theta_new =  wrap(MatrixXd::Constant(Y_.rows(), Y_.cols(), 0));
    NumericMatrix Theta_old =  wrap(MatrixXd::Constant(Y_.rows(), Y_.cols(), 0));
    NumericVector beta_new = wrap(VectorXd::Constant(X.length(),0));
    NumericVector beta_old = wrap(VectorXd::Constant(X.length(),0));
    NumericMatrix score_theta = wrap(MatrixXd::Constant(Y_.rows(), Y_.cols(), 0));
    NumericVector score_beta =  wrap(VectorXd::Constant(X.length(),0));
    
    double sigma_sum = 0;
    
    double obj_old = Compute_obj_linear_with_theta(X, Y, beta_old, Theta_old, 
                                                   sigma_sum, phi);
    double obj_new = obj_old;
    
    int num_iter = 0;
    
    while(num_iter !=iter_max){
        double s_t = s;
        List list_score = Compute_score_linear_with_theta(X, Y, beta_old, Theta_old);
        score_beta = as<NumericVector>(list_score["score_beta"]);
        score_theta = as<NumericMatrix>(list_score["score_theta"]);
        while(1){
            beta_new = Update_beta(beta_old, score_beta, s_t/(Y_.size()));
            List list_theta = Update_theta(Theta_old, score_theta, s_t, phi);
            Theta_new = as<NumericMatrix>(list_theta["Theta_est"]);
            const Map<VectorXd> sigma_(as<Map<VectorXd>> 
                                       (as<NumericVector>(list_theta["sigma"])));
            sigma_sum = sigma_.sum();
            obj_new = Compute_obj_linear_with_theta(X, Y, beta_new, Theta_new, 
                                                    sigma_sum, phi);
            //Rcout << obj_old<< std::endl;
            //Rcout << obj_new<< std::endl;
            if (obj_new > obj_old){
                s_t = s_t/2;
                //Rcout << s_t << std::endl;
            } 
            else{
                beta_old =  beta_new;
                Theta_old =  Theta_new;
                num_iter = num_iter + 1;
                Rcout << "NNR - Iter " << num_iter <<": Coefficiet = " << beta_new << std::endl;
                //Rcout << "Hi\n" << std::endl;
                break;
            }
            
        }
        if(obj_new < (1-tol) * obj_old){
            obj_old = obj_new;
            //Rcout << obj_new << std::endl;
        }
        else{
            break;
        }
        
    }
    return List::create(Named("beta_est") = beta_new,
                        Named("Theta_est") = Theta_new );
}


// [[Rcpp::export]]

List MLE_linear(const List & X,  const NumericMatrix & Y, 
                const NumericVector  beta_0, 
                const NumericMatrix & L_0, const NumericMatrix & R_0, 
                const double s, const int iter_max, const double tol){
    
    const Map<MatrixXd> Y_(as<Map<MatrixXd>> (Y));
    const Map<MatrixXd> L_0_(as<Map<MatrixXd>> (L_0));
    const Map<MatrixXd> R_0_(as<Map<MatrixXd>> (R_0));
    
    NumericVector beta_new = clone(beta_0);
    NumericVector beta_old =clone(beta_0);
    NumericMatrix L_new = clone(L_0);
    NumericMatrix L_old = clone(L_0);
    NumericMatrix R_new = clone(R_0);
    NumericMatrix R_old = clone(R_0);
    
    NumericMatrix score_L = wrap(MatrixXd::Constant(L_0_.rows(), L_0_.cols(), 0));
    NumericMatrix score_R = wrap(MatrixXd::Constant(R_0_.rows(), R_0_.cols(), 0));
    NumericVector score_beta =  wrap(VectorXd::Constant(X.length(),0));
    
    double obj_old = Compute_obj_linear_with_LR(X, Y, beta_old, L_old, R_old);
    double obj_new = obj_old;
    
    int num_iter = 0;
    
    while(num_iter !=iter_max){
        double s_t = s;
        List list_score = Compute_score_linear_with_LR(X, Y, beta_old, L_old, R_old);
        score_beta = as<NumericVector>(list_score["score_beta"]);
        score_L = as<NumericMatrix>(list_score["score_L"]);
        score_R = as<NumericMatrix>(list_score["score_R"]);
        
        while(1){
            beta_new = Update_beta(beta_old, score_beta, s_t/(Y_.size()));
            List list_LR = Update_LR(L_old, R_old, score_L, score_R, s_t/(Y_.cols()), s_t/(Y_.rows()));
            L_new = as<NumericMatrix>(list_LR["L_est"]);
            R_new = as<NumericMatrix>(list_LR["R_est"]);
            //Rcout << 1 << std::endl;
            //Rcout << 2 << std::endl;
            
            obj_new = Compute_obj_linear_with_LR(X, Y, beta_new, L_new, R_new);
            //Rcout << obj_old<< std::endl;
            //Rcout << obj_new<< std::endl;
            if (obj_new > obj_old && s_t >=1e-5){
                s_t = s_t/2;
                //Rcout << s_t << std::endl;
            } 
            else{
                beta_old =  beta_new;
                L_old = L_new;
                R_old = R_new;
                num_iter = num_iter + 1;
                Rcout << "MLE - Iter " << num_iter <<": Coefficiet = " << beta_new << std::endl;
                //Rcout << "Hi\n" << std::endl;
                break;
            }
        }
        
        if(obj_new < (1-tol) * obj_old){
            //Rcout << obj_old << std::endl;
            obj_old = obj_new;
        }
        else{
            break;
        }
    }
    return List::create(Named("beta_est") = beta_new,
                        Named("L_est") = L_new, 
                        Named("R_est") = R_new);
}



/*
// [[Rcpp::export]]

List fit_logit_pool(const List & X,  const NumericMatrix & Y, const double phi,  
               const double s, const int iter_max, const double tol){
    
    const Map<MatrixXd> Y_(as<Map<MatrixXd>> (Y));
    
    NumericMatrix Theta_new =  wrap(MatrixXd::Constant(Y_.rows(), Y_.cols(), 0));
    NumericMatrix Theta_old =  wrap(MatrixXd::Constant(Y_.rows(), Y_.cols(), 0));
    NumericVector beta_new = wrap(VectorXd::Constant(X.length(),0));
    NumericVector beta_old = wrap(VectorXd::Constant(X.length(),0));
    NumericMatrix score_theta = wrap(MatrixXd::Constant(Y_.rows(), Y_.cols(), 0));
    NumericVector score_beta =  wrap(VectorXd::Constant(X.length(),0));
    
    double sigma_sum = 0;
    
    double obj_old = Compute_obj_logit_with_theta(X, Y, beta_old, Theta_old, 
                                                  sigma_sum, phi);
    double obj_new = obj_old;
    
    int num_iter = 0;
    
    while(num_iter !=iter_max){
        double s_t = s;
        List list_score = Compute_score_logit_with_theta(X, Y, beta_old, Theta_old);
        score_beta = as<NumericVector>(list_score["score_beta"]);
        score_theta = as<NumericMatrix>(list_score["score_theta"]);
        while(1){
            beta_new = Update_beta(beta_old, score_beta, s_t/(Y_.size()));
            List list_theta = Update_theta(Theta_old, score_theta, s_t, phi);
            //Theta_new = as<NumericMatrix>(list_theta["Theta_est"]);
            const Map<VectorXd> sigma_(as<Map<VectorXd>> 
                                       (as<NumericVector>(list_theta["sigma"])));
            sigma_sum = sigma_.sum();
            obj_new = Compute_obj_logit_with_theta(X, Y, beta_new, Theta_new, 
                                                   sigma_sum, phi);
            Rcout << obj_old<< std::endl;
            Rcout << obj_new<< std::endl;
            if (obj_new > obj_old){
                s_t = s_t/2;
                Rcout << s_t << std::endl;
            } 
            else{
                beta_old =  beta_new;
                Theta_old =  Theta_new;
                num_iter = num_iter + 1;
                //Rcout << "Hi\n" << std::endl;
                break;
            }
            
        }
        if(obj_new < (1-tol) * obj_old){
            obj_old = obj_new;
        }
        else{
            break;
        }
        
    }
    
    
    
    return List::create(Named("beta_est") = beta_new,
                        Named("Theta_est") = Theta_new );
    
    
    
}
*/






