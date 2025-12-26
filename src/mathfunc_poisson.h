#ifndef MATHFUNC_POISSON_H
#define MATHFUNC_POISSON_H

#include <Rcpp.h>
#include <iostream>
#include <cmath>
#include <Rcpp.h>
#include "Eigen/Dense"
#include "Eigen/SVD"
#include "Eigen/Core"
#include "Eigen/Sparse"
#include <RcppEigen.h>
#include <random>
#include <stdlib.h>

#include "mathfunc_logit.h"


using namespace Eigen;
using namespace Rcpp;
using Eigen::Map;


NumericMatrix Compute_logl_poisson_with_theta(const List & X,  const NumericMatrix & Y, 
                                              const NumericVector & beta, const NumericMatrix & Theta );

NumericMatrix Compute_logl_poisson_with_theta_FE(const List & X,  const NumericMatrix & Y, 
                                                 const NumericVector & beta, 
                                                 const NumericVector & fe_N, const NumericVector & fe_T,
                                                 const NumericMatrix & Theta );

NumericMatrix Compute_logl_poisson_with_LR(const List & X,  const NumericMatrix & Y, 
                                           const NumericVector & beta, 
                                           const NumericMatrix & L, const NumericMatrix & R );

NumericMatrix Compute_logl_poisson_with_LR_FE(const List & X,  const NumericMatrix & Y, 
                                              const NumericVector & beta, 
                                              const NumericVector & fe_N, const NumericVector & fe_T,
                                              const NumericMatrix & L, const NumericMatrix & R );

List Compute_score_poisson_with_theta(const List & X, const NumericMatrix Y, 
                                    const NumericVector & beta, const NumericMatrix & Theta);

List Compute_score_poisson_with_theta_FE(const List & X, const NumericMatrix Y, 
                                         const NumericVector & beta, 
                                         const NumericVector & fe_N, const NumericVector & fe_T,
                                         const NumericMatrix & Theta);

List Compute_score_poisson_with_LR(const List & X, const NumericMatrix Y, 
                                 const NumericVector & beta, 
                                 const NumericMatrix & L, const NumericMatrix & R);

List Compute_score_poisson_with_LR_FE(const List & X, const NumericMatrix Y, 
                                      const NumericVector & beta, 
                                      const NumericVector & fe_N, const NumericVector & fe_T,
                                      const NumericMatrix & L, const NumericMatrix & R);

double Compute_obj_poisson_with_theta(const List & X,  const NumericMatrix & Y, 
                                    const NumericVector & beta, const NumericMatrix & Theta, 
                                    const double sigma_sum, const double phi);

double Compute_obj_poisson_with_theta_FE(const List & X,  const NumericMatrix & Y, 
                                         const NumericVector & beta, 
                                         const NumericVector & fe_N, const NumericVector & fe_T,
                                         const NumericMatrix & Theta, 
                                         const double sigma_sum, const double phi);

double Compute_obj_poisson_with_LR(const List & X,  const NumericMatrix & Y, 
                                 const NumericVector & beta, 
                                 const NumericMatrix & L, const NumericMatrix & R);

double Compute_obj_poisson_with_LR_FE(const List & X,  const NumericMatrix & Y, 
                                      const NumericVector & beta,
                                      const NumericVector & fe_N, const NumericVector & fe_T,
                                      const NumericMatrix & L, const NumericMatrix & R);

NumericMatrix Compute_first_order_poisson_with_LR(const List & X, const NumericMatrix & Y,
                                                  const NumericVector & beta, 
                                                  const NumericMatrix & L, const NumericMatrix & R);

NumericMatrix Compute_first_order_poisson_with_LR_FE(const List & X, const NumericMatrix & Y,
                                                     const NumericVector & beta, 
                                                     const NumericVector & fe_N, const NumericVector & fe_T,
                                                     const NumericMatrix & L, const NumericMatrix & R);

NumericMatrix Compute_second_order_poisson_with_LR(const List & X,
                                                   const NumericVector & beta, 
                                                   const NumericMatrix & L, const NumericMatrix & R);

NumericMatrix Compute_second_order_poisson_with_LR_FE(const List & X,
                                                      const NumericVector & beta, 
                                                      const NumericVector & fe_N, const NumericVector & fe_T,
                                                      const NumericMatrix & L, const NumericMatrix & R);

NumericMatrix Sim_dynamic_poisson(const List Z, 
                                  const double beta_W, const NumericVector & beta_Z, 
                                  const NumericMatrix & Theta);

#endif 
