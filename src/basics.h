#ifndef BASICS_H
#define BASICS_H

#include <Rcpp.h>
#include <iostream>
#include <cmath>
#include "Eigen/Dense"
#include "Eigen/SVD"
#include "Eigen/Core"
#include <Eigen/Sparse>
#include <RcppEigen.h>
#include <random>
#include <stdlib.h>
#include <Rmath.h>

using namespace Eigen;
using namespace Rcpp;
using Eigen::Map;


List SVD(NumericMatrix & A);

List SVT(const NumericMatrix & U, const NumericMatrix & V, 
         const NumericVector & sigma, double lambda);

List Low_rank_appro(NumericMatrix & Theta, const int num_factor);

NumericMatrix Compute_generalized_inverse(NumericMatrix & A);

NumericMatrix Compute_index_with_theta(const List & X, const NumericVector & beta, 
                                       const NumericMatrix & Theta );

NumericMatrix Compute_index_with_theta_FE(const List & X, const NumericVector & beta, 
                                          const NumericVector & fe_N, const NumericVector & fe_T,
                                          const NumericMatrix & Theta );


NumericMatrix Compute_index_with_LR(const List & X, const NumericVector & beta, 
                                    const NumericMatrix & L,  const NumericMatrix & R);

NumericMatrix Compute_index_with_LR_FE(const List & X, const NumericVector & beta, 
                                       const NumericVector & fe_N, const NumericVector & fe_T,
                                       const NumericMatrix & L,  const NumericMatrix & R);

NumericMatrix Convert_numeric_matrix_to_vector(const NumericMatrix & A);

IntegerMatrix Convert_integer_matrix_to_vector(const IntegerMatrix & A);

NumericMatrix Convert_List_to_vector(const List & X);

Eigen::SparseMatrix<double> Convert_FE_to_regressor(const NumericMatrix & L, const NumericMatrix & R);

#endif 
