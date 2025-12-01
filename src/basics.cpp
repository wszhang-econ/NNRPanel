#include "basics.h"

// [[Rcpp::export]]

List SVD(NumericMatrix & A){
    
    Map<MatrixXd> A_(as<Map<MatrixXd>> (A));
    for(int j=0; j!=A_.cols(); j++){
        for(int i=0; i!=A_.rows(); i++){
            if(abs(A(i,j))<=1e-15){
                A(i,j) = 0;
            }
        }
    }
    
    A_ = A_ +  MatrixXd::Identity(A_.rows(), A_.cols()) * 1e-10;
    
    BDCSVD<MatrixXd> svd(A_.rows(), A_.cols(), ComputeThinV | ComputeThinU);
    svd.compute(A_);
    if(isnan(svd.matrixU()(0,0))==TRUE){
        Rcout << "BDCSVD does not success\n" << std::endl;
        const Map<MatrixXd> A_(as<Map<MatrixXd>> (A));
        JacobiSVD<MatrixXd> svd(A_.rows(), A_.cols(), ComputeThinV | ComputeThinU);
        svd.compute(A_);
    }
    return List::create(Named("U_") = svd.matrixU(),
                        Named("V_") = svd.matrixV(),
                        Named("sigma_") = svd.singularValues());
}


// [[Rcpp::export]]

List SVT(const NumericMatrix & U, const NumericMatrix & V, 
         const NumericVector & sigma, double lambda){
    
    const Map<MatrixXd> U_ (as<Map<MatrixXd>> (U));
    const Map<MatrixXd> V_ (as<Map<MatrixXd>> (V));
    const Map<VectorXd> sigma_ (as<Map<VectorXd>> (sigma));
    
    VectorXd sigma_trunc_ = sigma_ - VectorXd::Constant(sigma_.size(), lambda);
    sigma_trunc_ = sigma_trunc_.cwiseMax(0);
    int num_positive = 0;
    for(int i=0; i!=sigma_.size(); i++){
        if(sigma_trunc_(i) > 1e-10){
            num_positive = num_positive + 1;
        }
        else{
            break;
        }
    }
    if(num_positive ==0){
        num_positive = num_positive + 1;
    }
    MatrixXd A_truc_ = U_.leftCols(num_positive) 
        * sigma_trunc_.head(num_positive).asDiagonal()
        * V_.leftCols(num_positive).transpose();
        NumericMatrix A_truc = wrap(A_truc_);
        NumericVector sigma_trunc = wrap(sigma_trunc_);
        return List::create(Named("A") = A_truc, 
                            Named("sigma") = sigma_trunc);
}

// [[Rcpp::export]]

List Low_rank_appro(NumericMatrix & Theta, const int num_factor){
    
    const Map<MatrixXd> Theta_ (as<Map<MatrixXd>> (Theta));
    int N = Theta_.rows();
    int T = Theta_.cols();
    
    MatrixXd L_truc_;
    MatrixXd R_truc_;
    
    if (num_factor > 0){
        List svd_list = SVD(Theta);
        
        MatrixXd U_ =  svd_list["U_"];
        MatrixXd V_ =  svd_list["V_"];
        VectorXd sigma_ = svd_list["sigma_"];
        
        L_truc_ = U_.leftCols(num_factor);
        VectorXd sigma_truc_ = sigma_.head(num_factor);
        R_truc_ = V_.leftCols(num_factor);
        
        
        L_truc_ = L_truc_ * sigma_truc_.asDiagonal() / std::sqrt(Theta_.cols());
        R_truc_ = std::sqrt(Theta_.cols()) * R_truc_;
        
        if(sigma_truc_(0) <=2e-10){
            L_truc_.setZero();
            R_truc_.setZero();
        }
    }
    else{
        L_truc_  = MatrixXd::Constant(N, 1, 0);
        R_truc_  = MatrixXd::Constant(T, 1, 0);
    }
    
    return List::create(Named("L") = wrap(L_truc_),
                        Named("R") = wrap(R_truc_));
} 


// [[Rcpp::export]]
NumericMatrix Compute_generalized_inverse(NumericMatrix & A){
    
    const Map<MatrixXd> A_ (as <Map<MatrixXd>> (A));
    
    List list_svd = SVD(A);
    MatrixXd U_ =  list_svd["U_"];
    MatrixXd V_ =  list_svd["V_"];
    VectorXd sigma_ = list_svd["sigma_"];
    
    MatrixXd sigma_inverse_ = MatrixXd::Constant(A_.cols(), A_.rows(), 0);
    for (int i = 0; i < sigma_.size(); i++){
        if (sigma_(i) > 1e-9){ 
            sigma_inverse_(i, i) = 1.0 / sigma_(i);
        }
    }
    Rcout << sigma_inverse_ << std::endl;
    MatrixXd A_ginv_ = V_ * sigma_inverse_ * U_.transpose();
    return wrap(A_ginv_);
}
    
    
// [[Rcpp::export]]

NumericMatrix Compute_index_with_theta(const List & X, const NumericVector & beta, 
                                       const NumericMatrix & Theta ){
    
    const Map<MatrixXd> Theta_ (as <Map<MatrixXd>> (Theta));
    
    MatrixXd index_est_ = Theta_;
    for(int i=0; i!=X.length() ; i++){
        const Map<MatrixXd> X_ (as <Map<MatrixXd>> (X[i]));
        index_est_ = index_est_ + beta[i] * X_;
    }
    return wrap(index_est_);
}

// [[Rcpp::export]]

NumericMatrix Compute_index_with_theta_FE(const List & X, const NumericVector & beta, 
                                          const NumericVector & fe_N, const NumericVector & fe_T,
                                          const NumericMatrix & Theta ){
    
    const Map<MatrixXd> Theta_ (as <Map<MatrixXd>> (Theta));
    const Map<VectorXd> fe_N_ (as <Map<VectorXd>> (fe_N));
    const Map<VectorXd> fe_T_ (as <Map<VectorXd>> (fe_T));
    
    int N = Theta_.rows();
    int T = Theta_.cols();
    
    MatrixXd index_est_ = Theta_;
    for(int i=0; i!=X.length() ; i++){
        const Map<MatrixXd> X_ (as <Map<MatrixXd>> (X[i]));
        index_est_ = index_est_ + beta[i] * X_;
    }
    index_est_ = index_est_ + fe_N_ * MatrixXd::Constant(1, T, 1);
    index_est_ = index_est_ + MatrixXd::Constant(N, 1, 1) * fe_T_.transpose();
    
    return wrap(index_est_);
}


// [[Rcpp::export]]

NumericMatrix Compute_index_with_LR(const List & X, const NumericVector & beta, 
                                    const NumericMatrix & L,  const NumericMatrix & R){
    
    const Map<MatrixXd> L_ (as <Map<MatrixXd>> (L));
    const Map<MatrixXd> R_ (as <Map<MatrixXd>> (R));
    
    MatrixXd index_est_ = L_ * R_.transpose();
    for(int i=0; i!=X.length() ; i++){
        const Map<MatrixXd> X_ (as <Map<MatrixXd>> (X[i]));
        index_est_ = index_est_ + beta[i] * X_;
    }
    return wrap(index_est_);
}


// [[Rcpp::export]]

NumericMatrix Compute_index_with_LR_FE(const List & X, const NumericVector & beta, 
                                       const NumericVector & fe_N, const NumericVector & fe_T,
                                       const NumericMatrix & L,  const NumericMatrix & R){
    
    const Map<MatrixXd> L_ (as <Map<MatrixXd>> (L));
    const Map<MatrixXd> R_ (as <Map<MatrixXd>> (R));
    const Map<VectorXd> fe_N_ (as <Map<VectorXd>> (fe_N));
    const Map<VectorXd> fe_T_ (as <Map<VectorXd>> (fe_T));
    
    int N = L_.rows();
    int T = R_.rows();
    
    MatrixXd index_est_ = L_ * R_.transpose();
    for(int i=0; i!=X.length() ; i++){
        const Map<MatrixXd> X_ (as <Map<MatrixXd>> (X[i]));
        index_est_ = index_est_ + beta[i] * X_;
    }
    index_est_ = index_est_ + fe_N_ * MatrixXd::Constant(1, T, 1);
    index_est_ = index_est_ + MatrixXd::Constant(N, 1, 1) * fe_T_.transpose();
    
    return wrap(index_est_);
}


// [[Rcpp::export]]

NumericMatrix Convert_numeric_matrix_to_vector(const NumericMatrix & A){
    const Map<MatrixXd> A_ (as <Map<MatrixXd>> (A));
    MatrixXd A_vec_ = Map<const MatrixXd>(A_.data(), A_.size(), 1);
    return wrap(A_vec_);
}

// [[Rcpp::export]]

IntegerMatrix Convert_integer_matrix_to_vector(const IntegerMatrix & A){
    const Map<MatrixXi> A_ (as <Map<MatrixXi>> (A));
    MatrixXi A_vec_ = Map<const MatrixXi>(A_.data(), A_.size(), 1);
    return wrap(A_vec_);
}

// [[Rcpp::export]]

NumericMatrix Convert_List_to_vector(const List & X){
    const Map<MatrixXd> X0_ (as <Map<MatrixXd>> (X[0]));
    MatrixXd A_ = MatrixXd::Constant(X0_.size(), X.length() ,0);
    for(int i=0; i!=X.length(); i++){
        Map<MatrixXd> X_vec_ (as <Map<MatrixXd>> (Convert_numeric_matrix_to_vector(X[i])));
        A_.col(i) = X_vec_;
    }
    return wrap(A_);
}


// [[Rcpp::export]]

Eigen::SparseMatrix<double> Convert_FE_to_regressor(const NumericMatrix & L, const NumericMatrix & R){
    const Map<MatrixXd> L_ (as <Map<MatrixXd>> (L));
    const Map<MatrixXd> R_ (as <Map<MatrixXd>> (R));
    int N = L_.rows();
    int T = R_.rows();
    int num_factor = L_.cols();
    
    Eigen::SparseMatrix<double> A_(N*T, (N + T) * num_factor);
    
    for(int i = 0; i != N; i++){
        for(int t = 0; t!=T; t++){
            for(int r = 0; r!= num_factor; r++){
                A_.insert(t*N + i, i*num_factor + r) = L_(i, r);
            } 
        }
    }
    for(int t = 0; t!=T; t++){
        for(int i = 0; i!= N; i++){
            for(int r = 0; r!= num_factor; r++){
                A_.insert(i + t * N, N*num_factor + t*num_factor + r) = R_(t, r);
            }  
        }
    }
  
    return A_;
}

/*
// [[Rcpp::export]]

Eigen::SparseMatrix<double> compute_FEtWFE(const Eigen::SparseMatrix<double> & FE_, 
                                          const Eigen::VectorXd & W_vec_) {

    Eigen::SparseMatrix<double> WFE_ = FE_;
    
    for (int k = 0; k < WFE_.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(WFE_, k); it; ++it) {
            it.valueRef() *= W_vec_(it.row());
        }
    }
    
    Eigen::SparseMatrix<double> FEtWFE_ = FE_.transpose() * WFE_;
    
    return FEtWFE_;
} 
 
*/

/*
// [[Rcpp::export]]

NumericVector Compute_porjection_residual(const NumericMatrix & X, 
                                          const NumericMatrix & L, const NumericMatrix & R,
                                          const NumericMatrix & W){
    NumericVector W_vec = Convert_matrix_to_vector(W);
    NumericVector X_vec = Convert_matrix_to_vector(X);
    const Map<VectorXd> W_vec_ (as <Map<VectorXd>> (W_vec));
    const Map<VectorXd> X_vec_ (as <Map<VectorXd>> (X_vec));
    Eigen::SparseMatrix<double> FE_ = Convert_FE_to_regressor(L, R);
    
    Eigen::SparseMatrix<double> WFE_ = FE_;
    
    for (int k = 0; k < WFE_.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(WFE_, k); it; ++it) {
            it.valueRef() *= W_vec_(it.row());
        }
    }
    
    MatrixXd FEtWFE_ = FE_.transpose() * WFE_;
    NumericMatrix FEtWFE = wrap(FEtWFE_);
    NumericMatrix FEtWFE_ginv = Compute_generalized_inverse(FEtWFE);
    Map<MatrixXd> FEtWFE_ginv_ (as <Map<MatrixXd>> (FEtWFE_ginv));
    VectorXd R_X_ = W_vec_.array() * X_vec_.array();
    R_X_ = FE_ * FEtWFE_ginv_ * FE_.transpose() * R_X_;
    
    return wrap(R_X_);
}
 */

