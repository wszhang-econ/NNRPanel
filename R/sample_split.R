sample_split <- function(data, L, R){
    N <- nrow(data$Y)
    T <- ncol(data$Y)
    num_X <- length(data$X)
    num_factor <- ncol(L)
    
    N_upper <- 1:round(N/2)
    N_bottom <- (1+round(N/2)):N
    T_left <- 1:round(T/2)
    T_right <- (1+round(T/2)):T
    
    data_left <- list(Y = 1, X = vector("list", num_X), L = 0, R = 0)
    data_left$Y <- data$Y[, T_left]
    data_left$L <- L
    data_left$R <- R[T_left, , drop = FALSE]
    for(i in 1:num_X){
        data_left$X[[i]] <- data$X[[i]][, T_left]
    }
    
    data_right <- list(Y = 1, X = vector("list", num_X), L = 0, R = 0)
    data_right$Y <- data$Y[, T_right]
    data_right$L <- L
    data_right$R <- R[T_right, , drop = FALSE]
    for(i in 1:num_X){
        data_right$X[[i]] <- data$X[[i]][, T_right]
    }
    
    data_upper <- list(Y = 1, X = vector("list", num_X), L = 0, R = 0)
    data_upper$Y <- data$Y[N_upper, ]
    data_upper$L <- L[N_upper, , drop = FALSE]
    data_upper$R <- R
    for(i in 1:num_X){
        data_upper$X[[i]] <- data$X[[i]][N_upper, ]
    }
    
    data_bottom <- list(Y = 1, X = vector("list", num_X), L = 0, R = 0)
    data_bottom$Y <- data$Y[N_bottom, ]
    data_bottom$L <- L[N_bottom, , drop = FALSE]
    data_bottom$R <- R
    for(i in 1:num_X){
        data_bottom$X[[i]] <- data$X[[i]][N_bottom, ]
    }
    return(list(
        data_left = data_left,
        data_right = data_right,
        data_upper = data_upper,
        data_bottom = data_bottom
    ))
    
}