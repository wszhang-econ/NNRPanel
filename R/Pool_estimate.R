#' Pooled Estimator for Panel Data
#'
#' Computes the pooled estimator for a panel dataset by fitting a generalized 
#' linear model without fixed effects. 
#'
#' @param data_frame Panel data in long format (could be a data frame or 
#'   matrix). The first column must contain individual IDs, the second column is
#'   time, and the third column is the  outcome variable \(Y\); all remaining 
#'   columns are treated as covariates in the regression model. The panel is 
#'   required to be balanced.
#' @param func A character string indicating the model:
#'   \code{"logit"},  \code{"probit"}, or \code{"poisson"}.
#' @return A numeric vector of pooled coefficient estimates.
#'
#' @details
#' The function fits a pooled GLM without individual or time effects. For 
#' \code{func="poisson"}, it uses Poisson regression; for 
#' \code{func="logit"}, it uses logistic regression. 
#'
#' @examples
#' \dontrun{
#' data_list <- gen_data_poisson(N = 100, T = 50, num_factor = 2, num_Z = 1, beta = c(0.2))
#' data_frame <- list_to_data_frame(data_list)
#' 
#' beta <- Pool_estimate(data_frame, "poisson")
#' }
#' 
#' @export
Pool_estimate <- function(data_frame, func){
    colnames(data_frame) <- c("id", "time", "Y", paste0("X", seq_len(ncol(data_frame) - 3)))
    data_frame <- as.data.frame(data_frame)
    
    X_cols <- grep("^X", names(data_frame), value = TRUE)
    formula <- as.formula(paste("Y ~", paste(X_cols, collapse = " + ")))
    
    if(func == "poisson"){
        beta_pool <- glm(formula, data = data_frame, family = poisson(link = "log"))
    }
    if(func =="logit"){
        beta_pool <- glm(formula, data = data_frame, family = binomial(link = "logit"))
    }
    if(func =="probit"){
        beta_pool <- glm(formula, data = data_frame, family = binomial(link = "probit"))
    }
   
    
    beta_pool <- beta_pool$coefficients
    
    return(beta_pool)
}