gen_factor <- function(N, T, num_factor, sigma = 1){
    L <- rnorm(N * num_factor,  0, sigma) %>% matrix(N, num_factor)
    R <- rnorm(T * num_factor,  0, sigma) %>% matrix(T, num_factor)
    return(
        list(
            L = L,
            R = R
        )
    )
}

gen_Z <- function(N, T, num_Z, list_factor, sigma = 1){
    L <- list_factor$L
    R <- list_factor$R
    num_factor <- ncol(L)
    Z <- vector("list", length = num_Z)
    for (i in 1:num_Z){
        Z[[i]] <- rnorm(N * T, 0, 2)  %>% matrix(N, T)
        Z[[i]] <- Z[[i]]  +  L %*% matrix(1, num_factor, T)
        Z[[i]] <- Z[[i]]  +  matrix(1, N, num_factor) %*% t(R)
        Z[[i]] <- Z[[i]] + (L %*% t(R))
        list_factor_Z <- gen_factor(N, T, 1);
        L_Z <- list_factor_Z$L
        R_Z <- list_factor_Z$R
        Z[[i]] <- Z[[i]] + (L_Z ) %*% t(R_Z )
        #Z[[i]] <- -1 * Z[[i]]
    }
    return(Z)
}

gen_data_logit <- function(N, T, num_factor, num_Z, beta, sigma = 1){
    index <- matrix(0, N, T)
    list_factor <- gen_factor(N, T, num_factor, sigma)
    L <- list_factor$L
    R <- list_factor$R
    Z <- gen_Z(N, T, num_Z, list_factor, sigma)
    for (i in 1:num_Z){
        index <- index + beta[i] * Z[[i]]
    }
    index <- index + L %*% t(R)
    P <-  1 - 1/(exp(index) + 1) 
    epsilon <- runif(N*T) %>% matrix(N, T)
    Y <- matrix(0, N, T)
    Y[epsilon < P] <- 1
    return(
        list(
            X = Z,
            list_factor = list_factor, 
            index = index, 
            interact_f = L %*% t(R), 
            Y = Y, 
            P = P
        )
    )
}

gen_data_dynamic_logit <- function(N, T, num_factor, num_Z, num_W, beta_W, beta_Z, sigma = 1){
    index <- matrix(0, N, T)
    list_factor <- gen_factor(N, T, num_factor, sigma)
    L <- list_factor$L
    R <- list_factor$R
    Z <- gen_Z(N, T, num_Z, list_factor, sigma)
    for (i in 1:num_Z){
        index <- index + beta_Z[i] * Z[[i]]
    }
    Y = Sim_dynamic_logit(Z, beta_W, beta_Z, L %*% t(R))
    X = matrix(0, N, T)
    X[, -1] = Y[, -T]
    X <- c(list(X), Z)
    return(
        list(
            list_factor = list_factor, 
            interact_f = L %*% t(R), 
            Y = Y, 
            X = X
        )
    )
}

gen_data_poisson <- function(N, T, num_factor, num_Z, beta, sigma = 1){
    index <- matrix(0, N, T)
    list_factor <- gen_factor(N, T, num_factor, sigma)
    L <- list_factor$L
    R <- list_factor$R
    Z <- gen_Z(N, T, num_Z, list_factor, sigma)
    for (i in 1:num_Z){
        index <- index + beta[i] * Z[[i]]
    }
    index <- index + L %*% t(R)/2
    index <- exp(index)
    Y <- matrix(rpois(N*T, index), N, T)
    return(
        list(
            X = Z,
            list_factor = list_factor, 
            index = index, 
            interact_f = L %*% t(R), 
            Y = Y
        )
    )
}

gen_data_dynamic_poisson <- function(N, T, num_factor, num_Z, num_W, beta_W, beta_Z, sigma = 1){
    index <- matrix(0, N, T)
    list_factor <- gen_factor(N, T, num_factor, sigma)
    L <- list_factor$L
    R <- list_factor$R
    Z <- gen_Z(N, T, num_Z, list_factor,sigma)
    for (i in 1:num_Z){
        index <- index + beta_Z[i] * Z[[i]]
    }
    Y = Sim_dynamic_poisson(Z, beta_W, beta_Z, L %*% t(R)/2)
    X = matrix(0, N, T)
    X[, -1] = log(1 + Y[, -T])
    X <- c(list(X), Z)
    return(
        list(
            list_factor = list_factor, 
            interact_f = L %*% t(R), 
            Y = Y, 
            X = X
        )
    )
}


gen_data_probit <- function(N, T, num_factor, num_Z, beta, sigma = 1){
    index <- matrix(0, N, T)
    list_factor <- gen_factor(N, T, num_factor, sigma)
    L <- list_factor$L
    R <- list_factor$R
    Z <- gen_Z(N, T, num_Z, list_factor, sigma)
    for (i in 1:num_Z){
        index <- index + beta[i] * Z[[i]]
    }
    index <- index + L %*% t(R)
    P <-  pnorm(index)  
    epsilon <- runif(N*T) %>% matrix(N, T)
    Y <- matrix(0, N, T)
    Y[epsilon < P] <- 1
    return(
        list(
            X = Z,
            list_factor = list_factor, 
            index = index, 
            interact_f = L %*% t(R), 
            Y = Y, 
            P = P
        )
    )
}

gen_data_linear <- function(N, T, num_factor, num_Z, beta, sigma = 1){
    index <- matrix(0, N, T)
    list_factor <- gen_factor(N, T, num_factor, sigma)
    L <- list_factor$L
    R <- list_factor$R
    Z <- gen_Z(N, T, num_Z, list_factor, sigma)
    for (i in 1:num_Z){
        index <- index + beta[i] * Z[[i]]
    }
    index <- index + L %*% t(R)
    epsilon <- rnorm(N*T, 0, 1) %>% matrix(N, T)
    Y <- index + epsilon
    return(
        list(
            X = Z,
            list_factor = list_factor, 
            index = index, 
            interact_f = L %*% t(R), 
            Y = Y
        )
    )
}