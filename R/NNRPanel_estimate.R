#' The proposed estimator with bias corrections. 
#'
#' Implements the nonlinear nuclear norm–regularized (NNR) estimator with
#' interactive fixed effects for panel data, followed by a second-step MLE and
#' analytical / sample-splitting bias corrections. The function also reports the standard error.
#' The procedure supports both logit and Poisson outcome models. It accommodates 
#' either a prespecified tuning parameter \eqn{\varphi} or a data-driven one 
#' obtained via \code{calculate_tuning_parameter()}. 
#' In addition, the function computes results using both the data-driven estimate 
#' of the number of factors and, when provided by the user, a prespecified \code{R_true}. 
#'
#' @param data_frame Panel data in long format (data frame or matrix). The first
#'   column must contain individual IDs, the second column the time index, and
#'   the third column the outcome variable \code{Y}; all remaining columns are
#'   treated as covariates in the regression model. The panel is required to be
#'   balanced.
#' @param phi Optional tuning parameter for the NNR estimator. If omitted, a
#'   data-driven value is computed internally via
#'   \code{calculate_tuning_parameter()} using the same \code{data_frame},
#'   \code{func}, \code{delta}, \code{R_max}, \code{s}, \code{iter_max}, and
#'   \code{tol}.
#' @param func A character string indicating the outcome model: \code{"linear"} for linear regression,
#'   \code{"logit"} for logistic regression, \code{"probit"} for Probit regression, and \code{"poisson"} for Poisson
#'   regression.
#' @param delta Inflation factor used when computing a data-driven \eqn{\varphi}.
#'   The tuning parameter is set proportional to \eqn{(1 + \delta)} times the
#'   leading singular value of an appropriate residual matrix.
#' @param R_max Maximum number of factors allowed when selecting the factor
#'   dimension in the data-driven rank selection step.
#' @param R_true Optional integer giving the “true” number of factors. If
#'   supplied, the function reports both the FE and bias-corrected estimates
#'   using \code{R_true}, in addition to the data-driven results based on
#'   \code{num_factor_est}. If missing, only the data-driven versions are
#'   returned.
#' @param s Initial step size used in the inner optimization routines for both
#'   the NNR and MLE steps. This value serves as the starting step size for the line-search procedure.  
#' @param iter_max Maximum number of iterations used in the MLE step.
#' @param tol Convergence tolerance for the optimization routines. Iterations
#'   stop when the change in parameters or objective falls below this threshold.
#'
#' @return A list containing at least
#'   \describe{
#'     \item{\code{beta_nnr}}{NNR estimator of the regression coefficients.}
#'     \item{\code{num_factor_est}}{Data-driven estimate of the number of factors.}
#'     \item{\code{beta_fe_data}}{Second-step MLE coefficient estimates based
#'       on \code{num_factor_est}.}
#'     \item{\code{beta_corr_data}}{Analytically bias-corrected estimates based
#'       on \code{num_factor_est}.}
#'     \item{\code{beta_corr_sp_data}}{Sample-splitting bias-corrected estimates
#'       based on \code{num_factor_est}.}
#'     \item{\code{std_beta_corr_data}}{Estimated standard errors  based on
#'       based on \code{num_factor_est}.}
#'   }
#'   If \code{R_true} is provided, the list additionally includes
#'   \describe{
#'     \item{\code{beta_fe}}{Second step MLE coefficient estimates using \code{R_true}.}
#'     \item{\code{beta_corr}}{Analytically bias-corrected estimates using
#'       \code{R_true}.}
#'     \item{\code{beta_corr_sp}}{Sample-splitting bias-corrected estimates
#'       using \code{R_true}.}
#'     \item{\code{std_beta_corr}}{Estimated standard errors based on 
#'       \code{R_true}.}
#'   }
#'
#'
#' @examples
#' \dontrun{
#' data_list <- gen_data_logit(N = 200, T = 200, num_factor = 2,
#'                               num_Z = 1, beta = c(0.2))
#' data_frame <- list_to_data_frame(data_list)
#'
#' # Logit model with data-driven phi and rank
#' est_logit <- NNRPanel_estimate(
#'   data_frame = data_frame,
#'   func       = "logit",
#'   delta      = 0.05,
#'   R_max      = 5
#' )
#'
#' # Poisson model with user-specified phi and known R_true
#' est_logit <- NNRPanel_estimate(
#'   data_frame = data_frame,
#'   phi        = 5,
#'   func       = "logit",
#'   delta      = 0.05,
#'   R_max      = 5,
#'   R_true     = 2
#' )
#' }
#'
#' @export
NNRPanel_estimate <- function(data_frame, phi, func, delta = 0.05, R_max = 5, R_true, s = 10, iter_max = 10000, tol = 1e-8){
    
    
    data_list <- data_frame_to_list(data_frame)
    
    Y <- data_list$Y
    X <- data_list$X
    
    N <- nrow(Y)
    T <- ncol(Y)
    
    if (missing(phi)){
        phi <- calculate_tuning_parameter(data_frame,  func, delta, R_max = R_max, s = s, iter_max = iter_max, tol = tol)
    }
    
    if(func =="poisson"){
        nnr_fit <- fit_poisson(X, Y,  phi = phi, s = s, iter_max = 300, tol = tol)
    }
    if(func =="logit"){
        nnr_fit <- fit_logit(X, Y,  phi = phi, s = s, iter_max = 300, tol = tol)
    }
    if(func =="probit"){
        nnr_fit <- fit_probit(X, Y,  phi = phi, s = s, iter_max = 300, tol = tol)
    }
    if(func =="linear"){
        nnr_fit <- fit_linear(X, Y,  phi = phi, s = s, iter_max = 300, tol = tol)
    }
    
    
    beta_nnr <-  nnr_fit$beta_est
    
    # determine the number of factors 
    
    sigma <-  SVD(nnr_fit$Theta_est)$sigma
    num_factor_est <- determine_num_factor(sigma, R_max)

    
    # beta_fe with correct number of estimator R_ture
    
    if(!missing(R_true)){
        if(func == "poisson"){
            
            # Second step estimation 
            list_LR = Low_rank_appro(nnr_fit$Theta_est, R_true);
            L_0 = list_LR$L
            R_0 = list_LR$R
            beta_0 = beta_nnr  
            fe_fit = MLE_poisson(X, Y, beta_0, L_0, R_0, s = s, iter_max = iter_max, tol = tol)
            beta_fe <- fe_fit$beta_est
            L_fe = fe_fit$L_est
            R_fe = fe_fit$R_est
            
            # Analytical bias correction 
            
            bias_corr <- Compute_bias_corr_poisson(Y, X, beta_fe, L_fe, R_fe, 0)
            print(bias_corr)
            beta_corr <- bias_corr$beta_corr
            if(is.na(beta_corr[1])){
                beta_corr <- beta_fe
            }
            std_beta_corr <- bias_corr$std_est
            
            # Sample splitting bias correction with correct number of factor
            
            data_sample_split <-  sample_split(data_list, L_0, R_0)
            data_left <- data_sample_split$data_left
            data_right <- data_sample_split$data_right
            data_upper <- data_sample_split$data_upper
            data_bottom <- data_sample_split$data_bottom
            
            fe_fit_left = MLE_poisson(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, s = s, iter_max = iter_max, tol = tol)
            beta_fe_left <- fe_fit_left$beta_est
            fe_fit_right = MLE_poisson(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, s = s, iter_max = iter_max, tol = tol)
            beta_fe_right <- fe_fit_right$beta_est
            fe_fit_upper = MLE_poisson(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, s = s, iter_max = iter_max, tol = tol)
            beta_fe_upper <- fe_fit_upper$beta_est
            fe_fit_bottom = MLE_poisson(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, s = s, iter_max = iter_max, tol = tol)
            beta_fe_bottom <- fe_fit_bottom$beta_est
            
            beta_corr_sp <- 3* beta_fe - (beta_fe_left + beta_fe_right + beta_fe_upper + beta_fe_bottom)/2
            
        }
        if(func == "logit"){
            
            # Second step estimation 
            list_LR = Low_rank_appro(nnr_fit$Theta_est, R_true);
            L_0 = list_LR$L
            R_0 = list_LR$R
            beta_0 = beta_nnr  
            fe_fit = MLE_logit(X, Y, beta_0, L_0, R_0, s = s, iter_max = iter_max, tol = tol)
            beta_fe <- fe_fit$beta_est
            L_fe = fe_fit$L_est
            R_fe = fe_fit$R_est
            
            # Analytical bias correction 
            
            bias_corr <- Compute_bias_corr_logit(Y, X, beta_fe, L_fe, R_fe, 0)
            print(bias_corr)
            beta_corr <- bias_corr$beta_corr
            if(is.na(beta_corr[1])){
                beta_corr <- beta_fe
            }
            std_beta_corr <- bias_corr$std_est
            
            # Sample splitting bias correction with correct number of factor
            
            data_sample_split <-  sample_split(data_list, L_0, R_0)
            data_left <- data_sample_split$data_left
            data_right <- data_sample_split$data_right
            data_upper <- data_sample_split$data_upper
            data_bottom <- data_sample_split$data_bottom
            
            fe_fit_left = MLE_logit(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, s = s, iter_max = iter_max, tol = tol)
            beta_fe_left <- fe_fit_left$beta_est
            fe_fit_right = MLE_logit(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, s = s, iter_max = iter_max, tol = tol)
            beta_fe_right <- fe_fit_right$beta_est
            fe_fit_upper = MLE_logit(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, s = s, iter_max = iter_max, tol = tol)
            beta_fe_upper <- fe_fit_upper$beta_est
            fe_fit_bottom = MLE_logit(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, s = s, iter_max = iter_max, tol = tol)
            beta_fe_bottom <- fe_fit_bottom$beta_est
            
            beta_corr_sp <- 3* beta_fe - (beta_fe_left + beta_fe_right + beta_fe_upper + beta_fe_bottom)/2
            
        }
        if(func == "probit"){
            
            # Second step estimation 
            list_LR = Low_rank_appro(nnr_fit$Theta_est, R_true);
            L_0 = list_LR$L
            R_0 = list_LR$R
            beta_0 = beta_nnr  
            fe_fit = MLE_probit(X, Y, beta_0, L_0, R_0, s = s, iter_max = iter_max, tol = tol)
            beta_fe <- fe_fit$beta_est
            L_fe = fe_fit$L_est
            R_fe = fe_fit$R_est
            
            # Analytical bias correction 
            
            bias_corr <- Compute_bias_corr_probit(Y, X, beta_fe, L_fe, R_fe, 0)
            print(bias_corr)
            beta_corr <- bias_corr$beta_corr
            if(is.na(beta_corr[1])){
                beta_corr <- beta_fe
            }
            std_beta_corr <- bias_corr$std_est
            
            # Sample splitting bias correction with correct number of factor
            
            data_sample_split <-  sample_split(data_list, L_0, R_0)
            data_left <- data_sample_split$data_left
            data_right <- data_sample_split$data_right
            data_upper <- data_sample_split$data_upper
            data_bottom <- data_sample_split$data_bottom
            
            fe_fit_left = MLE_probit(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, s = s, iter_max = iter_max, tol = tol)
            beta_fe_left <- fe_fit_left$beta_est
            fe_fit_right = MLE_probit(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, s = s, iter_max = iter_max, tol = tol)
            beta_fe_right <- fe_fit_right$beta_est
            fe_fit_upper = MLE_probit(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, s = s, iter_max = iter_max, tol = tol)
            beta_fe_upper <- fe_fit_upper$beta_est
            fe_fit_bottom = MLE_probit(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, s = s, iter_max = iter_max, tol = tol)
            beta_fe_bottom <- fe_fit_bottom$beta_est
            
            beta_corr_sp <- 3* beta_fe - (beta_fe_left + beta_fe_right + beta_fe_upper + beta_fe_bottom)/2
            
        }
        if(func == "linear"){
            
            # Second step estimation 
            list_LR = Low_rank_appro(nnr_fit$Theta_est, R_true);
            L_0 = list_LR$L
            R_0 = list_LR$R
            beta_0 = beta_nnr  
            fe_fit = MLE_linear(X, Y, beta_0, L_0, R_0, s = s, iter_max = iter_max, tol = tol)
            beta_fe <- fe_fit$beta_est
            L_fe = fe_fit$L_est
            R_fe = fe_fit$R_est
            
            # Analytical bias correction 
            
            bias_corr <- Compute_bias_corr_linear(Y, X, beta_fe, L_fe, R_fe, 0)
            print(bias_corr)
            beta_corr <- bias_corr$beta_corr
            if(is.na(beta_corr[1])){
                beta_corr <- beta_fe
            }
            std_beta_corr <- bias_corr$std_est
            
            # Sample splitting bias correction with correct number of factor
            
            data_sample_split <-  sample_split(data_list, L_0, R_0)
            data_left <- data_sample_split$data_left
            data_right <- data_sample_split$data_right
            data_upper <- data_sample_split$data_upper
            data_bottom <- data_sample_split$data_bottom
            
            fe_fit_left = MLE_linear(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, s = s, iter_max = iter_max, tol = tol)
            beta_fe_left <- fe_fit_left$beta_est
            fe_fit_right = MLE_linear(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, s = s, iter_max = iter_max, tol = tol)
            beta_fe_right <- fe_fit_right$beta_est
            fe_fit_upper = MLE_linear(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, s = s, iter_max = iter_max, tol = tol)
            beta_fe_upper <- fe_fit_upper$beta_est
            fe_fit_bottom = MLE_linear(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, s = s, iter_max = iter_max, tol = tol)
            beta_fe_bottom <- fe_fit_bottom$beta_est
            
            beta_corr_sp <- 3* beta_fe - (beta_fe_left + beta_fe_right + beta_fe_upper + beta_fe_bottom)/2
            
        }
    }
    
    # beta_fe with number of factor (data_driven)
    
    if(func=="poisson"){
        
        # Second step estimation 
        
        list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor_est);
        L_0 = list_LR$L
        R_0 = list_LR$R
        beta_0 = beta_nnr
        if (num_factor_est==0){
            L_0 <- matrix(0, N, 1)
            R_0 <- matrix(0, T, 1)
        }
        
        fe_fit_data = MLE_poisson(X, Y, beta_0, L_0, R_0, s = s, iter_max = iter_max, tol = tol)
        beta_fe_data <- fe_fit_data$beta_est
        L_fe_data = fe_fit_data$L_est
        R_fe_data = fe_fit_data$R_est
        
        
        # Analytical bias correction with number of factor (data_driven)
        bias_corr_data <- Compute_bias_corr_poisson(Y, X, beta_fe_data, L_fe_data, R_fe_data, 0)
        beta_corr_data <- bias_corr_data$beta_corr
        if(is.na(beta_corr_data[1])){
            beta_corr_data <- beta_fe_data
        }
        std_beta_corr_data <- bias_corr_data$std_est
        
        
        # Sample splitting bias correction with number of factor (data_driven)
        
        data_sample_split <-  sample_split(data_list, L_0, R_0)
        data_left <- data_sample_split$data_left
        data_right <- data_sample_split$data_right
        data_upper <- data_sample_split$data_upper
        data_bottom <- data_sample_split$data_bottom
        
        fe_fit_left = MLE_poisson(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, s = s, iter_max = iter_max, tol = tol)
        beta_fe_left_data <- fe_fit_left$beta_est
        fe_fit_right = MLE_poisson(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, s = s,  iter_max = iter_max, tol = tol)
        beta_fe_right_data <- fe_fit_right$beta_est
        fe_fit_upper = MLE_poisson(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, s = s,  iter_max = iter_max, tol = tol)
        beta_fe_upper_data <- fe_fit_upper$beta_est
        fe_fit_bottom = MLE_poisson(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, s = s,  iter_max = iter_max, tol = tol)
        beta_fe_bottom_data <- fe_fit_bottom$beta_est
        
        beta_corr_sp_data <- 3* beta_fe_data - (beta_fe_left_data + beta_fe_right_data + beta_fe_upper_data + beta_fe_bottom_data)/2
        
    }
    
    
    if(func=="logit"){
        
        # Second step estimation 
        
        list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor_est);
        L_0 = list_LR$L
        R_0 = list_LR$R
        beta_0 = beta_nnr
        if (num_factor_est==0){
            L_0 <- matrix(0, N, 1)
            R_0 <- matrix(0, T, 1)
        }
        
        fe_fit_data = MLE_logit(X, Y, beta_0, L_0, R_0, s = s, iter_max = iter_max, tol = tol)
        beta_fe_data <- fe_fit_data$beta_est
        L_fe_data = fe_fit_data$L_est
        R_fe_data = fe_fit_data$R_est
        
        
        # Analytical bias correction with number of factor (data_driven)
        bias_corr_data <- Compute_bias_corr_logit(Y, X, beta_fe_data, L_fe_data, R_fe_data, 0)
        beta_corr_data <- bias_corr_data$beta_corr
        if(is.na(beta_corr_data[1])){
            beta_corr_data <- beta_fe_data
        }
        std_beta_corr_data <- bias_corr_data$std_est
        
        
        # Sample splitting bias correction with number of factor (data_driven)
        
        data_sample_split <-  sample_split(data_list, L_0, R_0)
        data_left <- data_sample_split$data_left
        data_right <- data_sample_split$data_right
        data_upper <- data_sample_split$data_upper
        data_bottom <- data_sample_split$data_bottom
        
        fe_fit_left = MLE_logit(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, s = s, iter_max = iter_max, tol = tol)
        beta_fe_left_data <- fe_fit_left$beta_est
        fe_fit_right = MLE_logit(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, s = s,  iter_max = iter_max, tol = tol)
        beta_fe_right_data <- fe_fit_right$beta_est
        fe_fit_upper = MLE_logit(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, s = s,  iter_max = iter_max, tol = tol)
        beta_fe_upper_data <- fe_fit_upper$beta_est
        fe_fit_bottom = MLE_logit(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, s = s,  iter_max = iter_max, tol = tol)
        beta_fe_bottom_data <- fe_fit_bottom$beta_est
        
        beta_corr_sp_data <- 3* beta_fe_data - (beta_fe_left_data + beta_fe_right_data + beta_fe_upper_data + beta_fe_bottom_data)/2
        
    }
    if(func=="probit"){
        
        # Second step estimation 
        
        list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor_est);
        L_0 = list_LR$L
        R_0 = list_LR$R
        beta_0 = beta_nnr
        if (num_factor_est==0){
            L_0 <- matrix(0, N, 1)
            R_0 <- matrix(0, T, 1)
        }
        
        fe_fit_data = MLE_probit(X, Y, beta_0, L_0, R_0, s = s, iter_max = iter_max, tol = tol)
        beta_fe_data <- fe_fit_data$beta_est
        L_fe_data = fe_fit_data$L_est
        R_fe_data = fe_fit_data$R_est
        
        
        # Analytical bias correction with number of factor (data_driven)
        bias_corr_data <- Compute_bias_corr_probit(Y, X, beta_fe_data, L_fe_data, R_fe_data, 0)
        beta_corr_data <- bias_corr_data$beta_corr
        if(is.na(beta_corr_data[1])){
            beta_corr_data <- beta_fe_data
        }
        std_beta_corr_data <- bias_corr_data$std_est
        
        
        # Sample splitting bias correction with number of factor (data_driven)
        
        data_sample_split <-  sample_split(data_list, L_0, R_0)
        data_left <- data_sample_split$data_left
        data_right <- data_sample_split$data_right
        data_upper <- data_sample_split$data_upper
        data_bottom <- data_sample_split$data_bottom
        
        fe_fit_left = MLE_probit(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, s = s, iter_max = iter_max, tol = tol)
        beta_fe_left_data <- fe_fit_left$beta_est
        fe_fit_right = MLE_probit(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, s = s,  iter_max = iter_max, tol = tol)
        beta_fe_right_data <- fe_fit_right$beta_est
        fe_fit_upper = MLE_probit(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, s = s,  iter_max = iter_max, tol = tol)
        beta_fe_upper_data <- fe_fit_upper$beta_est
        fe_fit_bottom = MLE_probit(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, s = s,  iter_max = iter_max, tol = tol)
        beta_fe_bottom_data <- fe_fit_bottom$beta_est
        
        beta_corr_sp_data <- 3* beta_fe_data - (beta_fe_left_data + beta_fe_right_data + beta_fe_upper_data + beta_fe_bottom_data)/2
        
    }
    if(func=="linear"){
        
        # Second step estimation 
        
        list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor_est);
        L_0 = list_LR$L
        R_0 = list_LR$R
        beta_0 = beta_nnr
        if (num_factor_est==0){
            L_0 <- matrix(0, N, 1)
            R_0 <- matrix(0, T, 1)
        }
        
        fe_fit_data = MLE_linear(X, Y, beta_0, L_0, R_0, s = s, iter_max = iter_max, tol = tol)
        beta_fe_data <- fe_fit_data$beta_est
        L_fe_data = fe_fit_data$L_est
        R_fe_data = fe_fit_data$R_est
        
        
        # Analytical bias correction with number of factor (data_driven)
        bias_corr_data <- Compute_bias_corr_linear(Y, X, beta_fe_data, L_fe_data, R_fe_data, 0)
        beta_corr_data <- bias_corr_data$beta_corr
        if(is.na(beta_corr_data[1])){
            beta_corr_data <- beta_fe_data
        }
        std_beta_corr_data <- bias_corr_data$std_est
        
        
        # Sample splitting bias correction with number of factor (data_driven)
        
        data_sample_split <-  sample_split(data_list, L_0, R_0)
        data_left <- data_sample_split$data_left
        data_right <- data_sample_split$data_right
        data_upper <- data_sample_split$data_upper
        data_bottom <- data_sample_split$data_bottom
        
        fe_fit_left = MLE_linear(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, s = s, iter_max = iter_max, tol = tol)
        beta_fe_left_data <- fe_fit_left$beta_est
        fe_fit_right = MLE_linear(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, s = s,  iter_max = iter_max, tol = tol)
        beta_fe_right_data <- fe_fit_right$beta_est
        fe_fit_upper = MLE_linear(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, s = s,  iter_max = iter_max, tol = tol)
        beta_fe_upper_data <- fe_fit_upper$beta_est
        fe_fit_bottom = MLE_linear(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, s = s,  iter_max = iter_max, tol = tol)
        beta_fe_bottom_data <- fe_fit_bottom$beta_est
        
        beta_corr_sp_data <- 3* beta_fe_data - (beta_fe_left_data + beta_fe_right_data + beta_fe_upper_data + beta_fe_bottom_data)/2
        
    }
    
    
    if(!missing(R_true)){
        result <- list(
            beta_nnr = beta_nnr,
            num_factor_est = num_factor_est,
            beta_fe = beta_fe,
            beta_corr = beta_corr,
            beta_corr_sp = beta_corr_sp,
            std_beta_corr = std_beta_corr,
            beta_fe_data = beta_fe_data, 
            beta_corr_data = beta_corr_data, 
            beta_corr_sp_data = beta_corr_sp_data,
            std_beta_corr_data = std_beta_corr_data
        ) 
    }
    if(missing(R_true)){
        result <- list(
            beta_nnr = beta_nnr,
            num_factor_est = num_factor_est,
            beta_fe_data = beta_fe_data, 
            beta_corr_data = beta_corr_data, 
            beta_corr_sp_data = beta_corr_sp_data,
            std_beta_corr_data = std_beta_corr_data
        ) 
    }
    
    
    
    return(result)
}