#include "basics.h"
#include "mathfunc_logit.h"
#include "mathfunc_poisson.h"
#include "mathfunc_probit.h"

using std::vector;


// [[Rcpp::export]]

NumericMatrix Compute_weighted(const NumericMatrix & X, 
                               const NumericMatrix & Lambda, const NumericMatrix & Gamma, 
                               const NumericMatrix W){
    
    const Eigen::Map<MatrixXd> X_(as<Eigen::Map<MatrixXd>> (X));
    const Eigen::Map<MatrixXd> Lambda_(as<Eigen::Map<MatrixXd>> (Lambda));
    const Eigen::Map<MatrixXd> Gamma_(as<Eigen::Map<MatrixXd>> (Gamma));
    NumericMatrix W_truc = clone(W);
    Eigen::Map<MatrixXd> W_truc_(as<Eigen::Map<MatrixXd>> (W_truc));
    
    int N = X_.rows();
    int T = X_.cols();
    int num_factor = Lambda_.cols();
    
    for(int i=0; i!= N; i++){
        for(int t=0; t!=T; t++ ){
            if(W_truc(i, t)<=1e-4){
                W_truc(i, t) = 1e-4;
            }
        }
    }
    for(int i=0; i!= N; i++){
        for(int t=0; t!=T; t++ ){
            if(W_truc(i, t)>=1e10){
                W_truc(i, t) = 1e10;
            }
        }
    }
    
    IntegerMatrix id(N, T);
    IntegerMatrix time(N, T);
    for(int i=0; i!= N; i++){
        for(int t=0; t!=T; t++ ){
            id(i, t) = i;
            time(i, t) = t;
        }
    }
    
    IntegerMatrix id_vec = Convert_integer_matrix_to_vector(id);
    IntegerMatrix time_vec = Convert_integer_matrix_to_vector(time);
    NumericMatrix X_vec = Convert_numeric_matrix_to_vector(X);
    NumericVector W_vec = wrap(Map<VectorXd>(W_truc_.data(), W_truc_.size()));
    
    Eigen::MatrixXd Lambda_rep_ = Lambda_.replicate(T, 1);
    Eigen::MatrixXd Gamma_rep_ = Eigen::MatrixXd::Constant(N*T, num_factor, 1);
    for (int t=0; t!=T; t++){
        Gamma_rep_.block(t*N, 0, N, num_factor) = Gamma_.row(t).replicate(N, 1);
    }
    
    
    DataFrame df = DataFrame::create(
        Named("X") = X_vec,
        Named("id") = id_vec,
        Named("time") = time_vec
    );
    
    
    for (int r = 0; r!=num_factor; r++) {
        df.push_back(Lambda_rep_.col(r), "Lambda" + std::to_string(r + 1)); 
    }
    for (int r = 0; r!=num_factor; r++) {
        df.push_back(Gamma_rep_.col(r), "Gamma" + std::to_string(r + 1)); 
    }
    
    std::string formula_str = "X ~ 0 |";
    for (int r = 0; r!=num_factor; r++) {
        if(r==0){
            formula_str += "time[Lambda" + std::to_string(r + 1) + "]";
            formula_str += " + id[Gamma" + std::to_string(r + 1) + "]"; 
        }
        else{
            formula_str += " + time[Lambda" + std::to_string(r + 1) + "]";
            formula_str += " + id[Gamma" + std::to_string(r + 1) + "]"; 
        }
    }
    Function as_formula("as.formula");
    SEXP formula = as_formula(formula_str);
    Environment fixest_env = Environment::namespace_env("fixest");
    Function feols = fixest_env["feols"];
    
    List result = feols(_["fml"] = formula, _["data"] = df, _["weights"] = W_vec);
    
    NumericVector residuals_vec = as<NumericVector>(result["residuals"]);
    NumericMatrix residuals(N, T);
    for (int i = 0; i != N; i++) {
        for (int t = 0; t!= T; t++) {
            residuals(i, t) = residuals_vec[t * N + i];
        }
    }
    return(residuals);
}


// [[Rcpp::export]]

NumericMatrix Compute_weighted_FE(const NumericMatrix & X, 
                                  const NumericMatrix & Lambda, const NumericMatrix & Gamma, 
                                  const NumericMatrix W){
    
    const Eigen::Map<MatrixXd> X_(as<Eigen::Map<MatrixXd>> (X));
    const Eigen::Map<MatrixXd> Lambda_(as<Eigen::Map<MatrixXd>> (Lambda));
    const Eigen::Map<MatrixXd> Gamma_(as<Eigen::Map<MatrixXd>> (Gamma));
    NumericMatrix W_truc = clone(W);
    Eigen::Map<MatrixXd> W_truc_(as<Eigen::Map<MatrixXd>> (W_truc));
    
    int N = X_.rows();
    int T = X_.cols();
    int num_factor = Lambda_.cols();
    
    for(int i=0; i!= N; i++){
        for(int t=0; t!=T; t++ ){
            if(W_truc(i, t)<=1e-4){
                W_truc(i, t) = 1e-4;
            }
        }
    }
    for(int i=0; i!= N; i++){
        for(int t=0; t!=T; t++ ){
            if(W_truc(i, t)>=1e10){
                W_truc(i, t) = 1e10;
            }
        }
    }
    
    IntegerMatrix id(N, T);
    IntegerMatrix time(N, T);
    for(int i=0; i!= N; i++){
        for(int t=0; t!=T; t++ ){
            id(i, t) = i;
            time(i, t) = t;
        }
    }
    
    IntegerMatrix id_vec = Convert_integer_matrix_to_vector(id);
    IntegerMatrix time_vec = Convert_integer_matrix_to_vector(time);
    NumericMatrix X_vec = Convert_numeric_matrix_to_vector(X);
    NumericVector W_vec = wrap(Map<VectorXd>(W_truc_.data(), W_truc_.size()));
    
    Eigen::MatrixXd Lambda_rep_ = Lambda_.replicate(T, 1);
    Eigen::MatrixXd Gamma_rep_ = Eigen::MatrixXd::Constant(N*T, num_factor, 1);
    for (int t=0; t!=T; t++){
        Gamma_rep_.block(t*N, 0, N, num_factor) = Gamma_.row(t).replicate(N, 1);
    }
    
    
    DataFrame df = DataFrame::create(
        Named("X") = X_vec,
        Named("id") = id_vec,
        Named("time") = time_vec
    );
    
    
    for (int r = 0; r!=num_factor; r++) {
        df.push_back(Lambda_rep_.col(r), "Lambda" + std::to_string(r + 1)); 
    }
    for (int r = 0; r!=num_factor; r++) {
        df.push_back(Gamma_rep_.col(r), "Gamma" + std::to_string(r + 1)); 
    }
    
    std::string formula_str = "X ~ 0 |";
    for (int r = 0; r!=num_factor; r++) {
        if(r==0){
            formula_str += "time[Lambda" + std::to_string(r + 1) + "]";
            formula_str += " + id[Gamma" + std::to_string(r + 1) + "]"; 
        }
        else{
            formula_str += " + time[Lambda" + std::to_string(r + 1) + "]";
            formula_str += " + id[Gamma" + std::to_string(r + 1) + "]"; 
        }
    }
    formula_str += " + time"; 
    formula_str += " + id"; 
    Function as_formula("as.formula");
    SEXP formula = as_formula(formula_str);
    Environment fixest_env = Environment::namespace_env("fixest");
    Function feols = fixest_env["feols"];
    
    List result = feols(_["fml"] = formula, _["data"] = df, _["weights"] = W_vec);
    
    NumericVector residuals_vec = as<NumericVector>(result["residuals"]);
    NumericMatrix residuals(N, T);
    for (int i = 0; i != N; i++) {
        for (int t = 0; t!= T; t++) {
            residuals(i, t) = residuals_vec[t * N + i];
        }
    }
    return(residuals);
}



/*
// [[Rcpp::export]]

NumericMatrix Compute_weighted_FE(const NumericMatrix & X, 
                         const NumericMatrix & Lambda, const NumericMatrix & Gamma, 
                         const NumericMatrix W){
    
    const Eigen::Map<MatrixXd> X_(as<Eigen::Map<MatrixXd>> (X));
    const Eigen::Map<MatrixXd> Lambda_(as<Eigen::Map<MatrixXd>> (Lambda));
    const Eigen::Map<MatrixXd> Gamma_(as<Eigen::Map<MatrixXd>> (Gamma));
    NumericMatrix W_truc = clone(W);
    Eigen::Map<MatrixXd> W_truc_(as<Eigen::Map<MatrixXd>> (W_truc));
    
    int N = X_.rows();
    int T = X_.cols();
    int num_factor = Lambda_.cols();
    
    for(int i=0; i!= N; i++){
        for(int t=0; t!=T; t++ ){
            if(W_truc(i, t)<=1e-4){
                W_truc(i, t) = 1e-4;
            }
        }
    }
    
    IntegerMatrix id(N, T);
    IntegerMatrix time(N, T);
    for(int i=0; i!= N; i++){
        for(int t=0; t!=T; t++ ){
            id(i, t) = i;
            time(i, t) = t;
        }
    }
    
    IntegerMatrix id_vec = Convert_integer_matrix_to_vector(id);
    IntegerMatrix time_vec = Convert_integer_matrix_to_vector(time);
    NumericMatrix X_vec = Convert_numeric_matrix_to_vector(X);
    NumericVector W_vec = wrap(Map<VectorXd>(W_truc_.data(), W_truc_.size()));

    Eigen::MatrixXd Lambda_rep_ = Lambda_.replicate(T, 1);
    Eigen::MatrixXd Gamma_rep_ = Eigen::MatrixXd::Constant(N*T, num_factor, 1);
    for (int t=0; t!=T; t++){
        Gamma_rep_.block(t*N, 0, N, num_factor) = Gamma_.row(t).replicate(N, 1);
    }
    
    
    DataFrame df = DataFrame::create(
        Named("X") = X_vec,
        Named("id") = id_vec,
        Named("time") = time_vec
    );
    
    for (int r = 0; r!=num_factor; r++) {
        df.push_back(Lambda_rep_.col(r), "Lambda" + std::to_string(r + 1)); 
    }
    for (int r = 0; r!=num_factor; r++) {
        df.push_back(Gamma_rep_.col(r), "Gamma" + std::to_string(r + 1)); 
    }
    
    std::string formula_str = "X ~ 0";
    for (int r = 0; r!=num_factor; r++) {
        formula_str += " + Lambda" + std::to_string(r + 1) + "*factor(time)";
        formula_str += " + Gamma" + std::to_string(r + 1) + "*factor(id)";
    }
    Function as_formula("as.formula");
    SEXP formula = as_formula(formula_str);
    Environment fixest_env = Environment::namespace_env("fixest");
    Function feols = fixest_env["feols"];

    List result = feols(_["fml"] = formula, _["data"] = df, _["weights"] = W_vec);

    NumericVector residuals_vec = as<NumericVector>(result["residuals"]);
    NumericMatrix residuals(N, T);
    for (int i = 0; i != N; i++) {
        for (int t = 0; t!= T; t++) {
            residuals(i, t) = residuals_vec[t * N + i];
        }
    }
    return(residuals);
}
*/

// [[Rcpp::export]]

std::vector<std::vector<Eigen::MatrixXd>> Create_X_res  (const List & X, 
                                                         const NumericMatrix & Lambda, 
                                                         const NumericMatrix & Gamma, 
                                                         const NumericMatrix & W) {
    
    const Map<Eigen::MatrixXd> X_ (as <Map<Eigen::MatrixXd>> (X[0]));
    int dx = X.length();
    int N = X_.rows();
    int T = X_.cols();
    
    std::vector<std::vector<Eigen::MatrixXd>> X_res_(N, std::vector<Eigen::MatrixXd>(T));
    for (int i = 0; i < N; i++) {
        for (int t = 0; t < T; t++) {
            Eigen::MatrixXd x_res_ = Eigen::MatrixXd::Constant(dx, 1, 0);
            X_res_[i][t] = x_res_;
        }
    }
    for(int d=0; d!=X.length() ; d++){
        NumericMatrix X_res_dx = Compute_weighted(as<NumericMatrix>(X[d]), Lambda, Gamma, W);
        for (int i = 0; i < N; i++){
            for (int t = 0; t < T; t++){
                X_res_[i][t](d, 0) = X_res_dx(i, t);
            }
        }
    }
    return X_res_;
} 


// [[Rcpp::export]]

NumericMatrix Compute_bias_corr_W(const NumericMatrix & second_order, 
                                  const std::vector<std::vector<Eigen::MatrixXd>> & X_res_){
    const Map<MatrixXd> second_order_ (as <Map<MatrixXd>> (second_order));
    int N = second_order_.rows();
    int T = second_order_.cols();
    int dx = X_res_[0][0].rows();
    
    MatrixXd W_ = MatrixXd::Constant(dx, dx, 0);
    for(int i=0; i!=N; i++){
        for(int t=0; t!=T; t++){
            W_ = W_ + second_order(i, t) * (X_res_[i][t] * X_res_[i][t].transpose());
        }
    }
    W_ = W_ / N;
    W_ = W_ /T;
    
    return wrap(W_);
}



// [[Rcpp::export]]

NumericMatrix Compute_bias_corr_V(const NumericMatrix & first_order, 
                                  const std::vector<std::vector<Eigen::MatrixXd>> & X_res_){
    // const Map<MatrixXd> first_order_ (as <Map<MatrixXd>> (first_order));
    int N = first_order.rows();
    int T = first_order.cols();
    int dx = X_res_[0][0].rows();
    
    MatrixXd W_ = MatrixXd::Constant(dx, dx, 0);
    for(int i=0; i!=N; i++){
        for(int t=0; t!=T; t++){
            W_ = W_ + first_order(i, t) * first_order(i, t) * (X_res_[i][t] * X_res_[i][t].transpose());
        }
    }
    for(int i=0; i!=N; i++){
        for(int t=0; t!=T; t++){
            W_ = W_ + first_order(i, t) * first_order(t, i) * (X_res_[i][t] * X_res_[t][i].transpose());
        }
    }
    
    W_ = W_ / N;
    W_ = W_ /T;
    
    return wrap(W_);
}




// [[Rcpp::export]]

NumericMatrix Compute_bias_corr_D(const NumericMatrix & first_order, 
                                  const NumericMatrix & second_order, 
                                  const NumericMatrix & third_order,
                                  const NumericMatrix & Lambda,
                                  const std::vector<std::vector<Eigen::MatrixXd>> & X_res_){
    
    const Map<MatrixXd> first_order_ (as <Map<MatrixXd>> (first_order));
    const Map<MatrixXd> second_order_ (as <Map<MatrixXd>> (second_order));
    const Map<MatrixXd> third_order_ (as <Map<MatrixXd>> (third_order));
    const Map<MatrixXd> Lambda_ (as <Map<MatrixXd>> (Lambda));
    
    int N = second_order_.rows();
    int T = second_order_.cols();
    int dx = X_res_[0][0].rows();
    int num_factor = Lambda_.cols();
    
    MatrixXd D_ = MatrixXd::Constant(dx, 1, 0);
    VectorXd lambda_i_ = VectorXd::Constant(num_factor, 0);
    for(int d = 0; d!=dx; d++){
        for(int t = 0; t!=T; t++){
            MatrixXd D_3_ = MatrixXd::Constant(num_factor, num_factor, 0);
            for (int i=0; i!=N; i++){
                lambda_i_ = Lambda_.row(i);
                D_3_ = D_3_ + second_order_(i, t) * lambda_i_* lambda_i_.transpose();
            }
            D_3_ = D_3_.inverse();
            double D_1 = 0.0;
            double D_2 = 0.0;
            for (int i=0; i!=N; i++){
                lambda_i_ = Lambda_.row(i);
                D_1 = D_1 +  first_order_(i, t) * second_order_(i, t) * X_res_[i][t](d, 0) 
                      * lambda_i_.transpose() * D_3_ * lambda_i_;
                D_2 = D_2 +  1.0/2.0 * third_order_(i, t)  *  X_res_[i][t](d, 0) 
                      * lambda_i_.transpose() * D_3_ * lambda_i_;
            }
            D_(d, 0) = D_(d, 0) + (D_1 + D_2 );
        }
        D_(d, 0) = -1.0 * D_(d, 0)/T;
    }
    
    return wrap(D_);
}



// [[Rcpp::export]]

NumericMatrix Compute_bias_corr_D_probit(const NumericVector & beta, 
                                         const NumericMatrix & first_order, 
                                         const NumericMatrix & second_order, 
                                         const NumericMatrix & Lambda,
                                         const std::vector<std::vector<Eigen::MatrixXd>> & X_res_){
    
    const Map<MatrixXd> first_order_ (as <Map<MatrixXd>> (first_order));
    const Map<MatrixXd> second_order_ (as <Map<MatrixXd>> (second_order));
    const Map<MatrixXd> Lambda_ (as <Map<MatrixXd>> (Lambda));
    const Map<VectorXd> beta_ (as <Map<VectorXd>> (beta));
    
    int N = second_order_.rows();
    int T = second_order_.cols();
    int dx = X_res_[0][0].rows();
    int num_factor = Lambda_.cols();
    
    MatrixXd D_ = MatrixXd::Constant(dx, 1, 0);
    VectorXd lambda_i_ = VectorXd::Constant(num_factor, 0);
    for(int d = 0; d!=dx; d++){
        for(int t = 0; t!=T; t++){
            MatrixXd D_3_ = MatrixXd::Constant(num_factor, num_factor, 0);
            for (int i=0; i!=N; i++){
                lambda_i_ = Lambda_.row(i);
                D_3_ = D_3_ + second_order_(i, t) * lambda_i_* lambda_i_.transpose();
            }
            D_3_ = D_3_.inverse();
            double D_1 = 0.0;
            for (int i=0; i!=N; i++){
                lambda_i_ = Lambda_.row(i);
                D_1 = D_1 +  first_order_(i, t) * second_order_(i, t) * X_res_[i][t](d, 0) * X_res_[i][t](d, 0) * beta_(d) * lambda_i_.transpose() * D_3_ * lambda_i_;
            }
            D_(d, 0) = D_(d, 0) + (D_1  );
        }
        D_(d, 0) = -1.0 * D_(d, 0)/T;
    }
    
    return wrap(D_);
}

// [[Rcpp::export]]

NumericMatrix Compute_bias_corr_B(const NumericMatrix & first_order, 
                                  const NumericMatrix & second_order, 
                                  const NumericMatrix & third_order,
                                  const NumericMatrix & Gamma,
                                  const std::vector<std::vector<Eigen::MatrixXd>> & X_res_, 
                                  const int truc){
    
    const Map<MatrixXd> first_order_ (as <Map<MatrixXd>> (first_order));
    const Map<MatrixXd> second_order_ (as <Map<MatrixXd>> (second_order));
    const Map<MatrixXd> third_order_ (as <Map<MatrixXd>> (third_order));
    const Map<MatrixXd> Gamma_ (as <Map<MatrixXd>> (Gamma));
    
    int N = second_order_.rows();
    int T = second_order_.cols();
    int dx = X_res_[0][0].rows();
    int num_factor = Gamma_.cols();
    
    MatrixXd B_ = MatrixXd::Constant(dx, 1, 0);
    VectorXd gamma_t_ = VectorXd::Constant(num_factor, 0);
    VectorXd gamma_tl_ = VectorXd::Constant(num_factor, 0);
    
    for(int d = 0; d!=dx; d++){
        for(int i = 0; i!=N; i++){
            MatrixXd B_3_ = MatrixXd::Constant(num_factor, num_factor, 0);
            for (int t=0; t!=T; t++){
                gamma_t_ = Gamma_.row(t);
                B_3_ = B_3_ + second_order_(i, t) * gamma_t_* gamma_t_.transpose();
            }
            B_3_ = B_3_.inverse();
            double B_2 = 0;
            VectorXd B_cov_ = VectorXd::Constant(truc + 1, 0);
            for (int t=0; t!=T; t++){
                gamma_t_ = Gamma_.row(t);
                for (int l = 0; l!=truc +1; l++){
                    if (t + l <T){
                        gamma_tl_ = Gamma_.row(t + l);
                        B_cov_(l) = B_cov_(l) +  first_order_(i, t + l) * second_order_(i, t) 
                            *  X_res_[i][t](d, 0) * gamma_tl_.transpose() * B_3_ * gamma_tl_;
                    }
                }
                B_2 = B_2 + 1.0/2.0 * third_order_(i, t)  * X_res_[i][t](d, 0) 
                            * gamma_tl_.transpose() * B_3_ * gamma_tl_;;
            }
            for (int l = 0; l!=truc +1; l++){
                B_cov_(l) = (B_cov_(l) * T) /(T-l);
            }
            B_(d, 0) = B_(d, 0) + (B_cov_.sum() + B_2);
            
        }
        B_(d, 0) = -1.0 * B_(d, 0)/N;
    }
    
    return wrap(B_);
}


// [[Rcpp::export]]

NumericMatrix Compute_bias_corr_B_probit(const NumericVector & beta,
                                       const NumericMatrix & first_order, 
                                       const NumericMatrix & second_order, 
                                       const NumericMatrix & Gamma,
                                       const std::vector<std::vector<Eigen::MatrixXd>> & X_res_, 
                                       const int truc){
    
    const Map<MatrixXd> first_order_ (as <Map<MatrixXd>> (first_order));
    const Map<MatrixXd> second_order_ (as <Map<MatrixXd>> (second_order));
    const Map<MatrixXd> Gamma_ (as <Map<MatrixXd>> (Gamma));
    const Map<VectorXd> beta_ (as <Map<VectorXd>> (beta));
    
    
    int N = second_order_.rows();
    int T = second_order_.cols();
    int dx = X_res_[0][0].rows();
    int num_factor = Gamma_.cols();
    
    MatrixXd B_ = MatrixXd::Constant(dx, 1, 0);
    VectorXd gamma_t_ = VectorXd::Constant(num_factor, 0);
    VectorXd gamma_tl_ = VectorXd::Constant(num_factor, 0);
    
    for(int d = 0; d!=dx; d++){
        for(int i = 0; i!=N; i++){
            MatrixXd B_3_ = MatrixXd::Constant(num_factor, num_factor, 0);
            for (int t=0; t!=T; t++){
                gamma_t_ = Gamma_.row(t);
                B_3_ = B_3_ + second_order_(i, t) * gamma_t_* gamma_t_.transpose();
            }
            B_3_ = B_3_.inverse();
            VectorXd B_cov_ = VectorXd::Constant(truc + 1, 0);
            for (int t=0; t!=T; t++){
                gamma_t_ = Gamma_.row(t);
                for (int l = 0; l!=truc +1; l++){
                    if (t + l <T){
                        gamma_tl_ = Gamma_.row(t + l);
                        B_cov_(l) = B_cov_(l) +  second_order_(i, t) *  X_res_[i][t](d, 0) *  X_res_[i][t](d, 0) * beta_(d) * gamma_tl_.transpose() * B_3_ * gamma_tl_;
                    }
                }
            }
            for (int l = 0; l!=truc +1; l++){
                B_cov_(l) = (B_cov_(l) * T) /(T-l);
            }
            B_(d, 0) = B_(d, 0) + (B_cov_.sum());
            
        }
        B_(d, 0) = -1.0 * B_(d, 0)/N;
    }
    
    return wrap(B_);
}


// [[Rcpp::export]]

List Compute_bias_corr_logit(const NumericMatrix & Y, const List & X, 
                             const NumericVector & beta, 
                             const NumericMatrix & L, const NumericMatrix & R, 
                             const int truc){
    
    const Map<MatrixXd> Y_ (as <Map<MatrixXd>> (Y));
    const Map<MatrixXd> L_ (as <Map<MatrixXd>> (L));
    //const Map<MatrixXd> R_ (as <Map<MatrixXd>> (R));
    const Map<VectorXd> beta_ (as <Map<VectorXd>> (beta));
    
    int N = Y_.rows();
    int T = Y_.cols();
    //int num_factor = L_.cols();
    //int dx = X.length();
    
    NumericMatrix first_order = Compute_first_order_logit_with_LR(X, Y, beta, L, R);
    NumericMatrix second_order = Compute_second_order_logit_with_LR(X, beta, L, R);
    NumericMatrix third_order = Compute_third_order_logit_with_LR(X, beta, L, R);
    
    
    std::vector<std::vector<Eigen::MatrixXd>> X_res_ = Create_X_res(X, L, R, second_order);

    NumericMatrix B = Compute_bias_corr_B(first_order, second_order, third_order,
                                          R, X_res_, truc);
    NumericMatrix D = Compute_bias_corr_D(first_order, second_order, third_order,
                                          L, X_res_);
    NumericMatrix W = Compute_bias_corr_W(second_order, X_res_);
    
    Map<MatrixXd> B_ (as <Map<MatrixXd>> (B));
    Map<MatrixXd> D_ (as <Map<MatrixXd>> (D));
    Map<MatrixXd> W_ (as <Map<MatrixXd>> (W));
    
    
    MatrixXd Cov_matrixX_ = W_.inverse();
    
    VectorXd beta_corr_ = beta_ + 1.0/T * Cov_matrixX_ * B_  + 1.0/N * Cov_matrixX_ * D_;

    
    VectorXd std_est_ = Cov_matrixX_.diagonal();
    std_est_ = std_est_/(N*T);
    std_est_ = std_est_.array().sqrt();
    
    
    return List::create(Named("beta_corr") = wrap(beta_corr_), 
                        Named("B_est") = B, 
                        Named("D_est") = D,
                        Named("W_est") = W, 
                        Named("Cov_est") = wrap(Cov_matrixX_), 
                        Named("std_est") = wrap(std_est_));
}


// [[Rcpp::export]]

List Compute_bias_corr_logit_robust(const NumericMatrix & Y, const List & X, 
                                    const NumericVector & beta, 
                                    const NumericMatrix & L, const NumericMatrix & R, 
                                    const int truc){
    
    const Map<MatrixXd> Y_ (as <Map<MatrixXd>> (Y));
    const Map<MatrixXd> L_ (as <Map<MatrixXd>> (L));
    //const Map<MatrixXd> R_ (as <Map<MatrixXd>> (R));
    const Map<VectorXd> beta_ (as <Map<VectorXd>> (beta));
    
    int N = Y_.rows();
    int T = Y_.cols();
    //int num_factor = L_.cols();
    //int dx = X.length();
    
    NumericMatrix first_order = Compute_first_order_logit_with_LR(X, Y, beta, L, R);
    NumericMatrix second_order = Compute_second_order_logit_with_LR(X, beta, L, R);
    NumericMatrix third_order = Compute_third_order_logit_with_LR(X, beta, L, R);
    
    
    std::vector<std::vector<Eigen::MatrixXd>> X_res_ = Create_X_res(X, L, R, second_order);
    
    NumericMatrix B = Compute_bias_corr_B(first_order, second_order, third_order,
                                          R, X_res_, truc);
    NumericMatrix D = Compute_bias_corr_D(first_order, second_order, third_order,
                                          L, X_res_);
    NumericMatrix W = Compute_bias_corr_W(second_order, X_res_);
    
    NumericMatrix V = Compute_bias_corr_V(first_order, X_res_);
    
    
    Map<MatrixXd> B_ (as <Map<MatrixXd>> (B));
    Map<MatrixXd> D_ (as <Map<MatrixXd>> (D));
    Map<MatrixXd> W_ (as <Map<MatrixXd>> (W));
    Map<MatrixXd> V_ (as <Map<MatrixXd>> (V));
    
    
    MatrixXd Cov_matrixX_ = W_.inverse();
    
    MatrixXd Cov_matrixX_robust_ = W_.inverse() * V_ * W_.inverse();
    
    VectorXd beta_corr_ = beta_ + 1.0/T * Cov_matrixX_ * B_  + 1.0/N * Cov_matrixX_ * D_;
    
    
    VectorXd std_est_ = Cov_matrixX_.diagonal();
    std_est_ = std_est_/(N*T);
    std_est_ = std_est_.array().sqrt();
    
    VectorXd std_est_robust_ = Cov_matrixX_robust_.diagonal();
    std_est_robust_ = std_est_robust_/(N*T);
    std_est_robust_ = std_est_robust_.array().sqrt();
    
    
    
    return List::create(Named("beta_corr") = wrap(beta_corr_), 
                        Named("B_est") = B, 
                        Named("D_est") = D,
                        Named("W_est") = W,
                        Named("V_est") = V,
                        Named("Cov_est") = wrap(Cov_matrixX_), 
                        Named("std_est") = wrap(std_est_), 
                        Named("Cov_est_robust") = wrap(Cov_matrixX_robust_), 
                        Named("std_est_robust") = wrap(std_est_robust_));
    
}


// [[Rcpp::export]]

List Compute_bias_corr_poisson(const NumericMatrix & Y, const List & X, 
                               const NumericVector & beta, 
                               const NumericMatrix & L, const NumericMatrix & R, 
                               const int truc){
    
    const Map<MatrixXd> Y_ (as <Map<MatrixXd>> (Y));
    const Map<MatrixXd> L_ (as <Map<MatrixXd>> (L));
    //const Map<MatrixXd> R_ (as <Map<MatrixXd>> (R));
    const Map<VectorXd> beta_ (as <Map<VectorXd>> (beta));
    
    int N = Y_.rows();
    int T = Y_.cols();
    //int num_factor = L_.cols();
    //int dx = X.length();
    
    NumericMatrix first_order = Compute_first_order_poisson_with_LR(X, Y, beta, L, R);
    NumericMatrix second_order = Compute_second_order_poisson_with_LR(X, beta, L, R);
    NumericMatrix third_order = Compute_second_order_poisson_with_LR(X, beta, L, R);
    
    
    std::vector<std::vector<Eigen::MatrixXd>> X_res_ = Create_X_res(X, L, R, second_order);
    
    NumericMatrix B = Compute_bias_corr_B(first_order, second_order, third_order,
                                          R, X_res_, truc);
    NumericMatrix D = Compute_bias_corr_D(first_order, second_order, third_order,
                                          L, X_res_);
    NumericMatrix W = Compute_bias_corr_W(second_order, X_res_);
    

    
    Map<MatrixXd> B_ (as <Map<MatrixXd>> (B));
    //Map<MatrixXd> D_ (as <Map<MatrixXd>> (D));
    Map<MatrixXd> W_ (as <Map<MatrixXd>> (W));

    
    MatrixXd Cov_matrixX_ = W_.inverse();
    

    VectorXd beta_corr_ = beta_ + 1.0/T * Cov_matrixX_ * B_;
    
    
    VectorXd std_est_ = Cov_matrixX_.diagonal();
    std_est_ = std_est_/(N*T);
    std_est_ = std_est_.array().sqrt();
    
    
    return List::create(Named("beta_corr") = wrap(beta_corr_), 
                        Named("B_est") = B, 
                        Named("D_est") = D,
                        Named("W_est") = W, 
                        Named("Cov_est") = wrap(Cov_matrixX_), 
                        Named("std_est") = wrap(std_est_));
}

// [[Rcpp::export]]

List Compute_bias_corr_poisson_robust(const NumericMatrix & Y, const List & X, 
                                      const NumericVector & beta, 
                                      const NumericMatrix & L, const NumericMatrix & R, 
                                      const int truc){
    
    const Map<MatrixXd> Y_ (as <Map<MatrixXd>> (Y));
    const Map<MatrixXd> L_ (as <Map<MatrixXd>> (L));
    //const Map<MatrixXd> R_ (as <Map<MatrixXd>> (R));
    const Map<VectorXd> beta_ (as <Map<VectorXd>> (beta));
    
    int N = Y_.rows();
    int T = Y_.cols();
    //int num_factor = L_.cols();
    int dx = X.length();
    
    NumericMatrix first_order = Compute_first_order_poisson_with_LR(X, Y, beta, L, R);
    NumericMatrix second_order = Compute_second_order_poisson_with_LR(X, beta, L, R);
    NumericMatrix third_order = Compute_second_order_poisson_with_LR(X, beta, L, R);
    
    
    std::vector<std::vector<Eigen::MatrixXd>> X_res_ = Create_X_res(X, L, R, second_order);
    
    NumericMatrix B = Compute_bias_corr_B(first_order, second_order, third_order,
                                          R, X_res_, truc);
    NumericMatrix D = Compute_bias_corr_D(first_order, second_order, third_order,
                                          L, X_res_);
    NumericMatrix W = Compute_bias_corr_W(second_order, X_res_);
    
    NumericMatrix V = Compute_bias_corr_V(first_order, X_res_);
    
    
    Map<MatrixXd> B_ (as <Map<MatrixXd>> (B));
    Map<MatrixXd> D_ (as <Map<MatrixXd>> (D));
    Map<MatrixXd> W_ (as <Map<MatrixXd>> (W));
    Map<MatrixXd> V_ (as <Map<MatrixXd>> (V));
    
    
    MatrixXd Cov_matrixX_ = W_.inverse();
    
    MatrixXd Cov_matrixX_robust_ = W_.inverse() * V_ * W_.inverse();
    
    VectorXd beta_corr_ = beta_ + 1.0/T * Cov_matrixX_ * B_  + 1.0/N * Cov_matrixX_ * D_;

    
    VectorXd std_est_ = Cov_matrixX_.diagonal();
    std_est_ = std_est_/(N*T);
    std_est_ = std_est_.array().sqrt();
    
    VectorXd std_est_robust_ = Cov_matrixX_robust_.diagonal();
    std_est_robust_ = std_est_robust_/(N*T);
    std_est_robust_ = std_est_robust_.array().sqrt();
    

    
    return List::create(Named("beta_corr") = wrap(beta_corr_), 
                        Named("B_est") = B, 
                        Named("D_est") = D,
                        Named("W_est") = W,
                        Named("V_est") = V,
                        Named("Cov_est") = wrap(Cov_matrixX_), 
                        Named("std_est") = wrap(std_est_), 
                        Named("Cov_est_robust") = wrap(Cov_matrixX_robust_), 
                        Named("std_est_robust") = wrap(std_est_robust_));
}


// [[Rcpp::export]]

List Compute_bias_corr_probit(const NumericMatrix & Y, const List & X, 
                              const NumericVector & beta, 
                              const NumericMatrix & L, const NumericMatrix & R, 
                              const int truc){
    
    const Map<MatrixXd> Y_ (as <Map<MatrixXd>> (Y));
    const Map<MatrixXd> L_ (as <Map<MatrixXd>> (L));
    //const Map<MatrixXd> R_ (as <Map<MatrixXd>> (R));
    const Map<VectorXd> beta_ (as <Map<VectorXd>> (beta));
    
    int N = Y_.rows();
    int T = Y_.cols();
    //int num_factor = L_.cols();
    //int dx = X.length();
    
    NumericMatrix first_order = Compute_first_order_probit_with_LR(X, Y, beta, L, R);
    NumericMatrix second_order = Compute_second_order_probit_with_LR(X, Y, beta, L, R);
    
    
    std::vector<std::vector<Eigen::MatrixXd>> X_res_ = Create_X_res(X, L, R, second_order);
    
    NumericMatrix B = Compute_bias_corr_B_probit(beta, first_order, second_order,
                                                 R, X_res_, truc);
    NumericMatrix D = Compute_bias_corr_D_probit(beta, first_order, second_order,
                                                 L, X_res_);
    NumericMatrix W = Compute_bias_corr_W(second_order, X_res_);
    

    
    Map<MatrixXd> B_ (as <Map<MatrixXd>> (B));
    Map<MatrixXd> D_ (as <Map<MatrixXd>> (D));
    Map<MatrixXd> W_ (as <Map<MatrixXd>> (W));

    
    MatrixXd Cov_matrixX_ = W_.inverse();
    

    VectorXd beta_corr_ = beta_ + 1.0/T * Cov_matrixX_ * B_  + 1.0/N * Cov_matrixX_ * D_;
    
    
    VectorXd std_est_ = Cov_matrixX_.diagonal();
    std_est_ = std_est_/(N*T);
    std_est_ = std_est_.array().sqrt();
    
    
    
    return List::create(Named("beta_corr") = wrap(beta_corr_), 
                        Named("B_est") = B, 
                        Named("D_est") = D,
                        Named("W_est") = W,
                        Named("Cov_est") = wrap(Cov_matrixX_), 
                        Named("std_est") = wrap(std_est_));
    
}

// [[Rcpp::export]]

List Compute_bias_corr_probit_robust(const NumericMatrix & Y, const List & X, 
                                     const NumericVector & beta, 
                                     const NumericMatrix & L, const NumericMatrix & R, 
                                     const int truc){
    
    const Map<MatrixXd> Y_ (as <Map<MatrixXd>> (Y));
    const Map<MatrixXd> L_ (as <Map<MatrixXd>> (L));
    //const Map<MatrixXd> R_ (as <Map<MatrixXd>> (R));
    const Map<VectorXd> beta_ (as <Map<VectorXd>> (beta));
    
    int N = Y_.rows();
    int T = Y_.cols();
    //int num_factor = L_.cols();
    //int dx = X.length();
    
    NumericMatrix first_order = Compute_first_order_probit_with_LR(X, Y, beta, L, R);
    NumericMatrix second_order = Compute_second_order_probit_with_LR(X, Y, beta, L, R);

    
    std::vector<std::vector<Eigen::MatrixXd>> X_res_ = Create_X_res(X, L, R, second_order);
    
    NumericMatrix B = Compute_bias_corr_B_probit(beta, first_order, second_order,
                                          R, X_res_, truc);
    NumericMatrix D = Compute_bias_corr_D_probit(beta, first_order, second_order,
                                          L, X_res_);
    NumericMatrix W = Compute_bias_corr_W(second_order, X_res_);
    
    NumericMatrix V = Compute_bias_corr_V(first_order, X_res_);
    
    
    Map<MatrixXd> B_ (as <Map<MatrixXd>> (B));
    Map<MatrixXd> D_ (as <Map<MatrixXd>> (D));
    Map<MatrixXd> W_ (as <Map<MatrixXd>> (W));
    Map<MatrixXd> V_ (as <Map<MatrixXd>> (V));
    
    
    MatrixXd Cov_matrixX_ = W_.inverse();
    
    MatrixXd Cov_matrixX_robust_ = W_.inverse() * V_ * W_.inverse();
    
    VectorXd beta_corr_ = beta_ + 1.0/T * Cov_matrixX_ * B_  + 1.0/N * Cov_matrixX_ * D_;
    
    
    VectorXd std_est_ = Cov_matrixX_.diagonal();
    std_est_ = std_est_/(N*T);
    std_est_ = std_est_.array().sqrt();
    
    VectorXd std_est_robust_ = Cov_matrixX_robust_.diagonal();
    std_est_robust_ = std_est_robust_/(N*T);
    std_est_robust_ = std_est_robust_.array().sqrt();
    
    
    
    return List::create(Named("beta_corr") = wrap(beta_corr_), 
                        Named("B_est") = B, 
                        Named("D_est") = D,
                        Named("W_est") = W,
                        Named("V_est") = V,
                        Named("Cov_est") = wrap(Cov_matrixX_), 
                        Named("std_est") = wrap(std_est_), 
                        Named("Cov_est_robust") = wrap(Cov_matrixX_robust_), 
                        Named("std_est_robust") = wrap(std_est_robust_));
    
}
