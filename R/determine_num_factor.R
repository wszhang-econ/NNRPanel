#' Determine the Number of Factors using eigenvalue ratio test. 
#'
#' @description
#' Computes the number of factors using eigenvalue ratio test 
#' (see Ahn, S. C., & Horenstein, A. R. (2013)) 
#'
#' @references
#' Ahn, S. C., & Horenstein, A. R. (2013). Eigenvalue ratio test for the number of factors. Econometrica, 81(3), 1203-1227.
#' 
#' @param sigma A numeric vector of non-negative values, typically singular 
#'   values (or eigenvalues) sorted in non-increasing order. Must have length 
#'   at least \code{R_max + 1}.
#' @param R_max An integer indicating the maximum number of factors to consider.
#'
#' @details
#' The algorithm follows two steps:
#' \itemize{
#'   \item If any element \code{sigma[i] <= 1e-6}, the estimated number of 
#'   factors is set to \code{i - 1}.
#'   \item Otherwise, it identifies the index \code{i} that maximizes 
#'   \code{sigma[i] / sigma[i + 1]}, capturing a sharp drop in the first \code{R_max + 1} 
#'   elements of \code{sigma}. 
#' }
#'
#'
#' @return An integer giving the estimated number of factor.
#'
#' @examples
#' sigma <- c(5.0, 3.1, 1.9, 0.2, 1e-8, 1e-9)
#' determine_num_factor(sigma, R_max = 4) 
#'
#' @export


determine_num_factor <- function(sigma, R_max){
    num_factor = 0
    zero_sigma = FALSE
    for (i in 1:(R_max + 1)) {
        if (sigma[i] <= 1e-6) {
            num_factor = i - 1
            zero_sigma = TRUE
            break
        }
    }
    ratio = 1
    if (zero_sigma == FALSE) {
        for (i in 1:R_max) {
            new_ratio = sigma[i] / sigma[i + 1]
            if (ratio < new_ratio) {
                ratio = new_ratio
                num_factor = i
            }
        }
    }
    return(num_factor)
}


# determine_num_factor <- function(sigma, R_max){
#     num_factor = 0;
#     zero_sigma = FALSE;
#     for (i in 1:(R_max+1)){
#         if (sigma[i] <=1e-6) {
#             num_factor = i-1;
#             zero_sigma = TRUE;
#             break;
#         }
#     }
#     ratio = 1
#     if (zero_sigma==FALSE){
#         for (i in 1:R_max){
#             if (ratio < sigma[i] / sigma[i+1]){
#                 ratio = sigma[i] / sigma[i+1];
#                 num_factor = i;
#             }
#         }
#     }
#     return(num_factor);
# }