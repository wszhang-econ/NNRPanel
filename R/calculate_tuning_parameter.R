#' Calculate the tuning parameter for nuclear norm regularied estimators
#'
#' This function provides a data-driven way to calculate tuning parameter \eqn{\varphi} 
#' for non-linear nuclear-norm regularized (NNR) estimator with interactive fixed
#' effects. 
#' It supports both logistic and Poisson models.
#'
#' The function first estimates a two–way fixed-effects model, obtains the 
#' fitted values, and sets the initial \eqn{\varphi} to be the  (1 + \eqn{\alpha}) * leading 
#' singular value of the residual matrix. It then runs the NNR step, 
#' selects the number of factors, refines the estimates via MLE, and 
#' updates \eqn{\varphi} based on the new residuals.
#'
#' @param data_frame Panel data in long format (could be a data frame or 
#'   matrix). The first column must contain individual IDs, the second column is
#'   time, and the third column is the  outcome variable \(Y\); all remaining 
#'   columns are treated as covariates in the regression model. The panel is 
#'   required to be balanced. 
#' @param func A character string indicating the model:
#'   \code{"linear"}, \code{"logit"},  \code{"probit"}, and \code{"poisson"}.
#' @param delta Positive scalar inflation factor used to set
#'   \eqn{\varphi = (1+\delta)\sigma_1}, where \eqn{\sigma_1} is the leading
#'   singular value of the residual matrix.
#' @param R_max Maximum number of factors allowed when determining the
#'   number of factors.
#' @param s Initial step size  for both the NNR and MLE steps. 
#' This value serves as the starting step size for the line-search procedure. 
#' @param iter_max Maximum number of iterations for optimizations. 
#' @param tol Convergence tolerance passed to the optimization routines.
#'
#' @return A single numeric value giving the \eqn{\varphi}.
#'
#'
#' @examples
#' \dontrun{
#' data_list <- gen_data_poisson(N = 200, T = 200, num_factor = 2, num_Z = 1, beta = c(0.2))
#' data_frame <- list_to_data_frame(data_list) 
#' phi <- calculate_tuning_parameter(
#'   data_frame = data_frame,
#'   func       = "poisson",
#'   delta      = 0.05,
#'   R_max      = 5
#' )
#' }
#'
#' @export
#' 
calculate_tuning_parameter <- function(data_frame, func, delta, R_max = 5, s = 10, iter_max = 10000, tol = 1e-8){
    
    data_list <- data_frame_to_list(data_frame)
    
    Y <- data_list$Y
    X <- data_list$X
    
    N <- nrow(Y)
    T <- ncol(Y)
    
    
    colnames(data_frame) <- c("id", "time", "Y", paste0("X", seq_len(ncol(data_frame) - 3)))
    data_frame <- as.data.frame(data_frame)
    
    X_cols <- grep("^X", names(data_frame), value = TRUE)
    formula <- as.formula(paste("Y ~", paste(X_cols, collapse = " + "), "| id + time"))
    if(func=="poisson"){
        TW_model <- fixest::fepois(formula, data = data_frame)
        Y_fit <- fitted(TW_model)
        Y_fit <- t(matrix(Y_fit, T, N))
    }
    if(func=="logit"){
        TW_model <- fixest::feglm(formula , family = binomial(link = "logit"), data = data_frame)
        Y_fit <- fitted(TW_model)
        Y_fit <- t(matrix(Y_fit, T, N))
    }
    if(func=="probit"){
        TW_model <- fixest::feglm(formula , family = binomial(link = "probit"), data = data_frame)
        Y_fit <- fitted(TW_model)
        Y_fit <- t(matrix(Y_fit, T, N))
    }
    if(func=="linear"){
        TW_model <- fixest::felos(formula, data = data_frame)
        Y_fit <- fitted(TW_model)
        Y_fit <- t(matrix(Y_fit, T, N))
    }
    
    
    
    
    phi <- SVD(Y -Y_fit)$sigma[1] * (1 + delta)
    
    
    if(func=="poisson"){
        
        nnr_fit = fit_poisson(X, Y,  phi, s, 1000, tol)
        beta_nnr <-  nnr_fit$beta_est
        
        sigma <-  SVD(nnr_fit$Theta_est)$sigma
        num_factor_est <- determine_num_factor(sigma, R_max)
        
        list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor_est);
        L_0 = list_LR$L
        R_0 = list_LR$R
        beta_0 = beta_nnr
        
        fe_fit = MLE_poisson(X, Y, beta_0, L_0, R_0, s, iter_max, tol)
        beta_fe <- fe_fit$beta_est
        L_fe = fe_fit$L_est
        R_fe = fe_fit$R_est
        
        
        index <- Compute_index_with_LR(X, beta_fe, L_fe, R_fe)
        index <- exp(index)
        phi <- SVD(Y -index)$sigma[1] * (1 + delta)
    }
    if(func=="logit"){
        nnr_fit = fit_logit(X, Y,  phi, s, 1000, tol)
        beta_nnr <-  nnr_fit$beta_est
        
        sigma <-  SVD(nnr_fit$Theta_est)$sigma
        num_factor_est <- determine_num_factor(sigma, R_max)
        
        list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor_est);
        L_0 = list_LR$L
        R_0 = list_LR$R
        beta_0 = beta_nnr
        
        fe_fit = MLE_logit(X, Y, beta_0, L_0, R_0, s, iter_max, tol)
        beta_fe <- fe_fit$beta_est
        L_fe = fe_fit$L_est
        R_fe = fe_fit$R_est
        
        
        index <- Compute_index_with_LR(X, beta_fe, L_fe, R_fe)
        P <-  1 - 1/(exp(index) + 1) 
        phi <- SVD(Y -P)$sigma[1] * (1 + delta)
    }
    if(func=="probit"){
        nnr_fit = fit_probit(X, Y,  phi, s, 1000, tol)
        beta_nnr <-  nnr_fit$beta_est
        
        sigma <-  SVD(nnr_fit$Theta_est)$sigma
        num_factor_est <- determine_num_factor(sigma, R_max)
        
        list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor_est);
        L_0 = list_LR$L
        R_0 = list_LR$R
        beta_0 = beta_nnr
        
        fe_fit = MLE_probit(X, Y, beta_0, L_0, R_0, s, iter_max, tol)
        beta_fe <- fe_fit$beta_est
        L_fe = fe_fit$L_est
        R_fe = fe_fit$R_est
        
        
        index <- Compute_index_with_LR(X, beta_fe, L_fe, R_fe)
        P <-  pnorm(index)
        phi <- SVD(Y -P)$sigma[1] * (1 + delta)
    }
    if(func=="linear"){
        nnr_fit = fit_linear(X, Y,  phi, s, 1000, tol)
        beta_nnr <-  nnr_fit$beta_est
        
        sigma <-  SVD(nnr_fit$Theta_est)$sigma
        num_factor_est <- determine_num_factor(sigma, R_max)
        
        list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor_est);
        L_0 = list_LR$L
        R_0 = list_LR$R
        beta_0 = beta_nnr
        
        fe_fit = MLE_linear(X, Y, beta_0, L_0, R_0, s, iter_max, tol)
        beta_fe <- fe_fit$beta_est
        L_fe = fe_fit$L_est
        R_fe = fe_fit$R_est
        
        
        index <- Compute_index_with_LR(X, beta_fe, L_fe, R_fe)
        phi <- SVD(Y -index)$sigma[1] * (1 + delta)
    }
    
    
    
    return(phi)
}