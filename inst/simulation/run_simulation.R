library(NNRPanel)


gen_factor <- function(N, T, num_factor, sigma = 1){
    L <- rnorm(N * num_factor,  0, sigma) 
    L <- matrix(L, N, num_factor)
    R <- rnorm(T * num_factor,  0, sigma) 
    R <- matrix(R, T, num_factor)
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
        Z[[i]] <- rnorm(N * T, 0, 2)  
        Z[[i]] <- matrix(Z[[i]] , N, T)
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
    epsilon <- runif(N*T) 
    epsilon <- matrix(epsilon, N, T)
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


rearrange_result <- function(result, beta, num_X){
    sim = length(result)
    
    # number of factors 
    num_factor_mean <- 0
    for(i in 1 : sim){
        num_factor_mean <- num_factor_mean + result[[i]]$num_factor_est;
    }
    num_factor_mean = 1.0 / sim * num_factor_mean
    
    # MSE of beta_pool
    
    mean_beta_pool <- rep(0, num_X)
    for(i in 1 : sim){
        mean_beta_pool <- mean_beta_pool + result[[i]]$beta_pool;
    }
    mean_beta_pool <-  1.0 / sim * mean_beta_pool
    bias_beta_pool <- mean_beta_pool - beta
    std_beta_pool <-  rep(0, num_X)
    for(i in 1 : sim){
        std_beta_pool <- std_beta_pool + (result[[i]]$beta_pool - mean_beta_pool)^2 ;
    }
    std_beta_pool <-  1.0 / sim * std_beta_pool
    std_beta_pool <- sqrt(std_beta_pool)
    
    # MSE of beta_nnr
    
    mean_beta_nnr <- rep(0, num_X)
    for(i in 1 : sim){
        mean_beta_nnr <- mean_beta_nnr + result[[i]]$beta_nnr;
    }
    mean_beta_nnr <-  1.0 / sim * mean_beta_nnr
    bias_beta_nnr <- mean_beta_nnr - beta
    std_beta_nnr <-  rep(0, num_X)
    for(i in 1 : sim){
        std_beta_nnr <- std_beta_nnr + (result[[i]]$beta_nnr - mean_beta_nnr)^2 ;
    }
    std_beta_nnr <-  1.0 / sim * std_beta_nnr
    std_beta_nnr <- sqrt(std_beta_nnr)
    
    # MSE of beta_fe with correct number of factors
    
    mean_beta_fe <- rep(0, num_X)
    for(i in 1 : sim){
        mean_beta_fe <- mean_beta_fe + result[[i]]$beta_fe;
    }
    mean_beta_fe <-  1.0 / sim * mean_beta_fe
    bias_beta_fe <- mean_beta_fe - beta
    std_beta_fe <-  rep(0, num_X)
    for(i in 1 : sim){
        std_beta_fe <- std_beta_fe + (result[[i]]$beta_fe - mean_beta_fe)^2 ;
    }
    std_beta_fe <-  1.0 / sim * std_beta_fe
    std_beta_fe <- sqrt(std_beta_fe)
    
    # MSE of beta_corr with correct number of factors
    
    mean_beta_corr <- rep(0, num_X)
    for(i in 1 : sim){
        mean_beta_corr <- mean_beta_corr + result[[i]]$beta_corr;
    }
    mean_beta_corr <-  1.0 / sim * mean_beta_corr
    bias_beta_corr <- mean_beta_corr - beta
    std_beta_corr <-  rep(0, num_X)
    for(i in 1 : sim){
        std_beta_corr <- std_beta_corr + (result[[i]]$beta_corr - mean_beta_corr)^2 ;
    }
    std_beta_corr <-  1.0 / sim * std_beta_corr
    std_beta_corr <- sqrt(std_beta_corr)
    
    # MSE of beta_corr sample splitting with correct number of factors
    
    mean_beta_corr_sp <- rep(0, num_X)
    for(i in 1 : sim){
        mean_beta_corr_sp <- mean_beta_corr_sp + result[[i]]$beta_corr_sp;
    }
    mean_beta_corr_sp <-  1.0 / sim * mean_beta_corr_sp
    bias_beta_corr_sp <- mean_beta_corr_sp - beta
    std_beta_corr_sp <-  rep(0, num_X)
    for(i in 1 : sim){
        std_beta_corr_sp <- std_beta_corr_sp + (result[[i]]$beta_corr_sp - mean_beta_corr_sp)^2 ;
    }
    std_beta_corr_sp <-  1.0 / sim * std_beta_corr_sp
    std_beta_corr_sp <- sqrt(std_beta_corr_sp)
    
    # MSE of beta_fe_data with data driven number of factors 
    
    mean_beta_fe_data <- rep(0, num_X)
    for(i in 1 : sim){
        mean_beta_fe_data <- mean_beta_fe_data + result[[i]]$beta_fe_data;
    }
    mean_beta_fe_data <-  1.0 / sim * mean_beta_fe_data
    bias_beta_fe_data <- mean_beta_fe_data - beta
    std_beta_fe_data <-  rep(0, num_X)
    for(i in 1 : sim){
        std_beta_fe_data <- std_beta_fe_data + (result[[i]]$beta_fe_data - mean_beta_fe_data)^2 ;
    }
    std_beta_fe_data <-  1.0 / sim * std_beta_fe_data
    std_beta_fe_data <- sqrt(std_beta_fe_data)
    
    
    # MSE of beta_corr_data with data driven number of factors 
    
    mean_beta_corr_data <- rep(0, num_X)
    for(i in 1 : sim){
        mean_beta_corr_data <- mean_beta_corr_data + result[[i]]$beta_corr_data;
    }
    mean_beta_corr_data <-  1.0 / sim * mean_beta_corr_data
    bias_beta_corr_data <- mean_beta_corr_data - beta
    std_beta_corr_data <-  rep(0, num_X)
    for(i in 1 : sim){
        std_beta_corr_data <- std_beta_corr_data + (result[[i]]$beta_corr_data - mean_beta_corr_data)^2 ;
    }
    std_beta_corr_data <-  1.0 / sim * std_beta_corr_data
    std_beta_corr_data <- sqrt(std_beta_corr_data)
    
    # MSE of beta_corr_sp_data sample splitting with data driven number of factors 
    
    mean_beta_corr_sp_data <- rep(0, num_X)
    for(i in 1 : sim){
        mean_beta_corr_sp_data <- mean_beta_corr_sp_data + result[[i]]$beta_corr_sp_data;
    }
    mean_beta_corr_sp_data <-  1.0 / sim * mean_beta_corr_sp_data
    bias_beta_corr_sp_data <- mean_beta_corr_sp_data - beta
    std_beta_corr_sp_data <-  rep(0, num_X)
    for(i in 1 : sim){
        std_beta_corr_sp_data <- std_beta_corr_sp_data + (result[[i]]$beta_corr_sp_data - mean_beta_corr_sp_data)^2 ;
    }
    std_beta_corr_sp_data <-  1.0 / sim * std_beta_corr_sp_data
    std_beta_corr_sp_data <- sqrt(std_beta_corr_sp_data)
    
    rearrange_list <-  list(
                            num_factor_mean = num_factor_mean,
                            bias_beta_pool = bias_beta_pool,
                            std_beta_pool = std_beta_pool,
                            bias_beta_nnr = bias_beta_nnr,
                            std_beta_nnr = std_beta_nnr,
                            bias_beta_fe = bias_beta_fe,
                            std_beta_fe = std_beta_fe,
                            bias_beta_corr = bias_beta_corr,
                            std_beta_corr = std_beta_corr,
                            bias_beta_corr_sp = bias_beta_corr_sp,
                            std_beta_corr_sp = std_beta_corr_sp,
                            bias_beta_fe_data = bias_beta_fe_data,
                            std_beta_fe_data = std_beta_fe_data,
                            bias_beta_corr_data = bias_beta_corr_data,
                            std_beta_corr_data = std_beta_corr_data, 
                            bias_beta_corr_sp_data = bias_beta_corr_sp_data, 
                            std_beta_corr_sp_data = std_beta_corr_sp_data)
    
    return(rearrange_list)

}



















library(dplyr)
library(fixest)
library(doParallel)

set.seed(199505)

delta = 0.05

num_Z = 1
num_W = 0
num_X = num_Z + num_W
num_factor = 2
beta = c(0.2)
R_max = 5
sim = 1000


seed <- sample.int(100000, sim, replace = FALSE)

N = 50
T = 40

myCluster <- makeCluster(detectCores()-1, type = "FORK")
registerDoParallel(myCluster)
raw_data <- foreach (k=1:sim) %dopar% {
    set.seed(seed[k])
    data <- gen_data_logit(N, T, num_factor, num_Z, beta)
    
    # Pool
    Yvec <- data$Y %>% c() 
    X1vec <- data$X[[1]]%>% c()
    id <- rep(1:N, T) %>% matrix(N, T) %>% c()
    time <- rep(1:T, N) %>% matrix(T, N) %>% t() %>% c()
    df <- data.frame(Yvec, X1vec, id, time)
    
    X <- data$X
    Y <- data$Y
    
    pool <- glm(Yvec ~ X1vec, data = df, family = binomial(link = "logit"))
    
    beta_pool <- pool$coefficients[-1]
    beta_pool <- c(unlist(beta_pool))
    print(beta_pool)
    
    # Compute tuning parameter
    
    TW_model <- feglm(Yvec ~ X1vec | id + time, family = binomial(link = "logit"), data = df)
    Y_fit <- fitted(TW_model) %>% matrix(N, T)
    phi <- SVD(Y -Y_fit)$sigma[1] * (1 + delta)
    print(phi)
    
    nnr_fit = fit_logit(X, Y,  phi, 10, 300, 1e-8)
    beta_nnr <-  nnr_fit$beta_est
    
    sigma <-  SVD(nnr_fit$Theta_est)$sigma
    num_factor_est <- determine_num_factor(sigma, R_max)
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    
    fe_fit = MLE_logit(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe <- fe_fit$beta_est
    L_fe = fe_fit$L_est
    R_fe = fe_fit$R_est
    
    
    index <- Compute_index_with_LR(X, beta_fe, L_fe, R_fe)
    P <-  1 - 1/(exp(index) + 1) 
    phi <- SVD(Y -P)$sigma[1] * (1 + delta)
    
    
    
    # NNR 
    
    nnr_fit = fit_logit(X, Y,  phi, 10, 300, 1e-8)
    beta_nnr <-  nnr_fit$beta_est
    print(beta_pool)
    
    # determine the number of factor
    
    sigma <-  SVD(nnr_fit$Theta_est)$sigma
    num_factor_est <- determine_num_factor(sigma, R_max)
    print(num_factor_est)
    
    # beta_fe with correct number of factor
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    
    fe_fit = MLE_logit(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe <- fe_fit$beta_est
    L_fe = fe_fit$L_est
    R_fe = fe_fit$R_est
    
    # beta_fe analytical bias correction with correct number of factor
    
    bias_corr <- Compute_bias_corr_logit_robust(Y, X, beta_fe, L_fe, R_fe, 0)
    print(bias_corr)
    beta_corr <- bias_corr$beta_corr
    if(is.na(beta_corr[1])){
        beta_corr <- beta_fe
    }
    std_beta_corr <- bias_corr$std_est
    
    # beta_fe sample splitting bias correction with correct number of factor
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_logit(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left <- fe_fit_left$beta_est
    fe_fit_right = MLE_logit(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right <- fe_fit_right$beta_est
    fe_fit_upper = MLE_logit(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_logit(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom <- fe_fit_bottom$beta_est
    
    beta_corr_sp <- 3* beta_fe - (beta_fe_left + beta_fe_right + beta_fe_upper + beta_fe_bottom)/2
    
    
    # beta_fe with number of factor (data_driven)
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor_est);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    if (num_factor_est==0){
        L_0 <- matrix(0, N, 1)
        R_0 <- matrix(0, T, 1)
    }
    
    fe_fit_data = MLE_logit(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe_data <- fe_fit_data$beta_est
    L_fe_data = fe_fit_data$L_est
    R_fe_data = fe_fit_data$R_est
    
    # beta_fe analytical bias correction with number of factor (data_driven)
    bias_corr_data <- Compute_bias_corr_logit(Y, X, beta_fe_data, L_fe_data, R_fe_data, 0)
    beta_corr_data <- bias_corr_data$beta_corr
    if(is.na(beta_corr_data[1])){
        beta_corr_data <- beta_fe_data
    }
    std_beta_corr_data <- bias_corr_data$std_est
    
    # beta_fe sample splitting bias correction with number of factor (data_driven)
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_logit(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left_data <- fe_fit_left$beta_est
    fe_fit_right = MLE_logit(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right_data <- fe_fit_right$beta_est
    fe_fit_upper = MLE_logit(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper_data <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_logit(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom_data <- fe_fit_bottom$beta_est
    
    beta_corr_sp_data <- 3* beta_fe_data - (beta_fe_left_data + beta_fe_right_data + beta_fe_upper_data + beta_fe_bottom_data)/2
    
    
    
    result <- list(
        beta_pool = beta_pool,
        beta_nnr = beta_nnr,
        num_factor_est = num_factor_est,
        beta_fe = beta_fe,
        beta_corr = beta_corr,
        beta_corr_sp = beta_corr_sp,
        std_beta_corr = std_beta_corr,
        beta_fe_data = beta_fe_data, 
        beta_corr_data = beta_corr_data, 
        beta_corr_sp_data = beta_corr_sp_data,
        std_beta_corr_data = std_beta_corr_data)
}
stopCluster(myCluster)    
summary_le_50_40 <- rearrange_result(raw_data, beta, num_X)




seed <- sample.int(100000, sim, replace = FALSE)

N = 100
T = 40

myCluster <- makeCluster(detectCores()-1, type = "FORK")
registerDoParallel(myCluster)
raw_data <- foreach (k=1:sim) %dopar% {
    set.seed(seed[k])
    data <- gen_data_logit(N, T, num_factor, num_Z, beta)
    
    # Pool
    Yvec <- data$Y %>% c() 
    X1vec <- data$X[[1]]%>% c()
    id <- rep(1:N, T) %>% matrix(N, T) %>% c()
    time <- rep(1:T, N) %>% matrix(T, N) %>% t() %>% c()
    df <- data.frame(Yvec, X1vec, id, time)
    
    X <- data$X
    Y <- data$Y
    
    pool <- glm(Yvec ~ X1vec, data = df, family = binomial(link = "logit"))
    
    beta_pool <- pool$coefficients[-1]
    beta_pool <- c(unlist(beta_pool))
    print(beta_pool)
    
    # Compute tuning parameter
    
    TW_model <- feglm(Yvec ~ X1vec | id + time, family = binomial(link = "logit"), data = df)
    Y_fit <- fitted(TW_model) %>% matrix(N, T)
    phi <- SVD(Y -Y_fit)$sigma[1] * (1 + delta)
    print(phi)
    
    nnr_fit = fit_logit(X, Y,  phi, 10, 300, 1e-8)
    beta_nnr <-  nnr_fit$beta_est
    
    sigma <-  SVD(nnr_fit$Theta_est)$sigma
    num_factor_est <- determine_num_factor(sigma, R_max)
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    
    fe_fit = MLE_logit(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe <- fe_fit$beta_est
    L_fe = fe_fit$L_est
    R_fe = fe_fit$R_est
    
    
    index <- Compute_index_with_LR(X, beta_fe, L_fe, R_fe)
    P <-  1 - 1/(exp(index) + 1) 
    phi <- SVD(Y -P)$sigma[1] * (1 + delta)
    
    
    
    # NNR 
    
    nnr_fit = fit_logit(X, Y,  phi, 10, 300, 1e-8)
    beta_nnr <-  nnr_fit$beta_est
    print(beta_pool)
    
    # determine the number of factor
    
    sigma <-  SVD(nnr_fit$Theta_est)$sigma
    num_factor_est <- determine_num_factor(sigma, R_max)
    print(num_factor_est)
    
    # beta_fe with correct number of factor
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    
    fe_fit = MLE_logit(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe <- fe_fit$beta_est
    L_fe = fe_fit$L_est
    R_fe = fe_fit$R_est
    
    # beta_fe analytical bias correction with correct number of factor
    
    bias_corr <- Compute_bias_corr_logit(Y, X, beta_fe, L_fe, R_fe, 0)
    print(bias_corr)
    beta_corr <- bias_corr$beta_corr
    if(is.na(beta_corr[1])){
        beta_corr <- beta_fe
    }
    std_beta_corr <- bias_corr$std_est
    
    # beta_fe sample splitting bias correction with correct number of factor
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_logit(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left <- fe_fit_left$beta_est
    fe_fit_right = MLE_logit(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right <- fe_fit_right$beta_est
    fe_fit_upper = MLE_logit(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_logit(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom <- fe_fit_bottom$beta_est
    
    beta_corr_sp <- 3* beta_fe - (beta_fe_left + beta_fe_right + beta_fe_upper + beta_fe_bottom)/2
    
    
    # beta_fe with number of factor (data_driven)
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor_est);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    if (num_factor_est==0){
        L_0 <- matrix(0, N, 1)
        R_0 <- matrix(0, T, 1)
    }
    
    fe_fit_data = MLE_logit(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe_data <- fe_fit_data$beta_est
    L_fe_data = fe_fit_data$L_est
    R_fe_data = fe_fit_data$R_est
    
    # beta_fe analytical bias correction with number of factor (data_driven)
    bias_corr_data <- Compute_bias_corr_logit(Y, X, beta_fe_data, L_fe_data, R_fe_data, 0)
    beta_corr_data <- bias_corr_data$beta_corr
    if(is.na(beta_corr_data[1])){
        beta_corr_data <- beta_fe_data
    }
    std_beta_corr_data <- bias_corr_data$std_est
    
    # beta_fe sample splitting bias correction with number of factor (data_driven)
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_logit(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left_data <- fe_fit_left$beta_est
    fe_fit_right = MLE_logit(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right_data <- fe_fit_right$beta_est
    fe_fit_upper = MLE_logit(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper_data <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_logit(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom_data <- fe_fit_bottom$beta_est
    
    beta_corr_sp_data <- 3* beta_fe_data - (beta_fe_left_data + beta_fe_right_data + beta_fe_upper_data + beta_fe_bottom_data)/2
    
    
    
    result <- list(
        beta_pool = beta_pool,
        beta_nnr = beta_nnr,
        num_factor_est = num_factor_est,
        beta_fe = beta_fe,
        beta_corr = beta_corr,
        beta_corr_sp = beta_corr_sp,
        std_beta_corr = std_beta_corr,
        beta_fe_data = beta_fe_data, 
        beta_corr_data = beta_corr_data, 
        beta_corr_sp_data = beta_corr_sp_data,
        std_beta_corr_data = std_beta_corr_data)
}
stopCluster(myCluster)    
summary_le_100_40 <- rearrange_result(raw_data, beta, num_X)



seed <- sample.int(100000, sim, replace = FALSE)

N = 200
T = 40

myCluster <- makeCluster(detectCores()-1, type = "FORK")
registerDoParallel(myCluster)
raw_data <- foreach (k=1:sim) %dopar% {
    set.seed(seed[k])
    data <- gen_data_logit(N, T, num_factor, num_Z, beta)
    
    # Pool
    Yvec <- data$Y %>% c() 
    X1vec <- data$X[[1]]%>% c()
    id <- rep(1:N, T) %>% matrix(N, T) %>% c()
    time <- rep(1:T, N) %>% matrix(T, N) %>% t() %>% c()
    df <- data.frame(Yvec, X1vec, id, time)
    
    X <- data$X
    Y <- data$Y
    
    pool <- glm(Yvec ~ X1vec, data = df, family = binomial(link = "logit"))
    
    beta_pool <- pool$coefficients[-1]
    beta_pool <- c(unlist(beta_pool))
    print(beta_pool)
    
    # Compute tuning parameter
    
    TW_model <- feglm(Yvec ~ X1vec | id + time, family = binomial(link = "logit"), data = df)
    Y_fit <- fitted(TW_model) %>% matrix(N, T)
    phi <- SVD(Y -Y_fit)$sigma[1] * (1 + delta)
    print(phi)
    
    nnr_fit = fit_logit(X, Y,  phi, 10, 300, 1e-8)
    beta_nnr <-  nnr_fit$beta_est
    
    sigma <-  SVD(nnr_fit$Theta_est)$sigma
    num_factor_est <- determine_num_factor(sigma, R_max)
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    
    fe_fit = MLE_logit(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe <- fe_fit$beta_est
    L_fe = fe_fit$L_est
    R_fe = fe_fit$R_est
    
    
    index <- Compute_index_with_LR(X, beta_fe, L_fe, R_fe)
    P <-  1 - 1/(exp(index) + 1) 
    phi <- SVD(Y -P)$sigma[1] * (1 + delta)
    
    
    
    # NNR 
    
    nnr_fit = fit_logit(X, Y,  phi, 10, 300, 1e-8)
    beta_nnr <-  nnr_fit$beta_est
    print(beta_pool)
    
    # determine the number of factor
    
    sigma <-  SVD(nnr_fit$Theta_est)$sigma
    num_factor_est <- determine_num_factor(sigma, R_max)
    print(num_factor_est)
    
    # beta_fe with correct number of factor
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    
    fe_fit = MLE_logit(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe <- fe_fit$beta_est
    L_fe = fe_fit$L_est
    R_fe = fe_fit$R_est
    
    # beta_fe analytical bias correction with correct number of factor
    
    bias_corr <- Compute_bias_corr_logit(Y, X, beta_fe, L_fe, R_fe, 0)
    print(bias_corr)
    beta_corr <- bias_corr$beta_corr
    if(is.na(beta_corr[1])){
        beta_corr <- beta_fe
    }
    std_beta_corr <- bias_corr$std_est
    
    # beta_fe sample splitting bias correction with correct number of factor
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_logit(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left <- fe_fit_left$beta_est
    fe_fit_right = MLE_logit(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right <- fe_fit_right$beta_est
    fe_fit_upper = MLE_logit(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_logit(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom <- fe_fit_bottom$beta_est
    
    beta_corr_sp <- 3* beta_fe - (beta_fe_left + beta_fe_right + beta_fe_upper + beta_fe_bottom)/2
    
    
    # beta_fe with number of factor (data_driven)
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor_est);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    if (num_factor_est==0){
        L_0 <- matrix(0, N, 1)
        R_0 <- matrix(0, T, 1)
    }
    
    fe_fit_data = MLE_logit(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe_data <- fe_fit_data$beta_est
    L_fe_data = fe_fit_data$L_est
    R_fe_data = fe_fit_data$R_est
    
    # beta_fe analytical bias correction with number of factor (data_driven)
    bias_corr_data <- Compute_bias_corr_logit(Y, X, beta_fe_data, L_fe_data, R_fe_data, 0)
    beta_corr_data <- bias_corr_data$beta_corr
    if(is.na(beta_corr_data[1])){
        beta_corr_data <- beta_fe_data
    }
    std_beta_corr_data <- bias_corr_data$std_est
    
    # beta_fe sample splitting bias correction with number of factor (data_driven)
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_logit(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left_data <- fe_fit_left$beta_est
    fe_fit_right = MLE_logit(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right_data <- fe_fit_right$beta_est
    fe_fit_upper = MLE_logit(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper_data <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_logit(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom_data <- fe_fit_bottom$beta_est
    
    beta_corr_sp_data <- 3* beta_fe_data - (beta_fe_left_data + beta_fe_right_data + beta_fe_upper_data + beta_fe_bottom_data)/2
    
    
    
    result <- list(
        beta_pool = beta_pool,
        beta_nnr = beta_nnr,
        num_factor_est = num_factor_est,
        beta_fe = beta_fe,
        beta_corr = beta_corr,
        beta_corr_sp = beta_corr_sp,
        std_beta_corr = std_beta_corr,
        beta_fe_data = beta_fe_data, 
        beta_corr_data = beta_corr_data, 
        beta_corr_sp_data = beta_corr_sp_data,
        std_beta_corr_data = std_beta_corr_data)
}
stopCluster(myCluster)         
summary_le_200_40 <- rearrange_result(raw_data, beta, num_X)



seed <- sample.int(100000, sim, replace = FALSE)

N = 100
T = 100

myCluster <- makeCluster(detectCores()-1, type = "FORK")
registerDoParallel(myCluster)
raw_data <- foreach (k=1:sim) %dopar% {
    set.seed(seed[k])
    data <- gen_data_logit(N, T, num_factor, num_Z, beta)
    
    # Pool
    Yvec <- data$Y %>% c() 
    X1vec <- data$X[[1]]%>% c()
    id <- rep(1:N, T) %>% matrix(N, T) %>% c()
    time <- rep(1:T, N) %>% matrix(T, N) %>% t() %>% c()
    df <- data.frame(Yvec, X1vec, id, time)
    
    X <- data$X
    Y <- data$Y
    
    pool <- glm(Yvec ~ X1vec, data = df, family = binomial(link = "logit"))
    
    beta_pool <- pool$coefficients[-1]
    beta_pool <- c(unlist(beta_pool))
    print(beta_pool)
    
    # Compute tuning parameter
    
    TW_model <- feglm(Yvec ~ X1vec | id + time, family = binomial(link = "logit"), data = df)
    Y_fit <- fitted(TW_model) %>% matrix(N, T)
    phi <- SVD(Y -Y_fit)$sigma[1] * (1 + delta)
    print(phi)
    
    nnr_fit = fit_logit(X, Y,  phi, 10, 300, 1e-8)
    beta_nnr <-  nnr_fit$beta_est
    
    sigma <-  SVD(nnr_fit$Theta_est)$sigma
    num_factor_est <- determine_num_factor(sigma, R_max)
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    
    fe_fit = MLE_logit(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe <- fe_fit$beta_est
    L_fe = fe_fit$L_est
    R_fe = fe_fit$R_est
    
    
    index <- Compute_index_with_LR(X, beta_fe, L_fe, R_fe)
    P <-  1 - 1/(exp(index) + 1) 
    phi <- SVD(Y -P)$sigma[1] * (1 + delta)
    
    
    
    # NNR 
    
    nnr_fit = fit_logit(X, Y,  phi, 10, 300, 1e-8)
    beta_nnr <-  nnr_fit$beta_est
    print(beta_pool)
    
    # determine the number of factor
    
    sigma <-  SVD(nnr_fit$Theta_est)$sigma
    num_factor_est <- determine_num_factor(sigma, R_max)
    print(num_factor_est)
    
    # beta_fe with correct number of factor
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    
    fe_fit = MLE_logit(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe <- fe_fit$beta_est
    L_fe = fe_fit$L_est
    R_fe = fe_fit$R_est
    
    # beta_fe analytical bias correction with correct number of factor
    
    bias_corr <- Compute_bias_corr_logit(Y, X, beta_fe, L_fe, R_fe, 0)
    print(bias_corr)
    beta_corr <- bias_corr$beta_corr
    if(is.na(beta_corr[1])){
        beta_corr <- beta_fe
    }
    std_beta_corr <- bias_corr$std_est
    
    # beta_fe sample splitting bias correction with correct number of factor
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_logit(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left <- fe_fit_left$beta_est
    fe_fit_right = MLE_logit(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right <- fe_fit_right$beta_est
    fe_fit_upper = MLE_logit(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_logit(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom <- fe_fit_bottom$beta_est
    
    beta_corr_sp <- 3* beta_fe - (beta_fe_left + beta_fe_right + beta_fe_upper + beta_fe_bottom)/2
    
    
    # beta_fe with number of factor (data_driven)
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor_est);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    if (num_factor_est==0){
        L_0 <- matrix(0, N, 1)
        R_0 <- matrix(0, T, 1)
    }
    
    fe_fit_data = MLE_logit(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe_data <- fe_fit_data$beta_est
    L_fe_data = fe_fit_data$L_est
    R_fe_data = fe_fit_data$R_est
    
    # beta_fe analytical bias correction with number of factor (data_driven)
    bias_corr_data <- Compute_bias_corr_logit(Y, X, beta_fe_data, L_fe_data, R_fe_data, 0)
    beta_corr_data <- bias_corr_data$beta_corr
    if(is.na(beta_corr_data[1])){
        beta_corr_data <- beta_fe_data
    }
    std_beta_corr_data <- bias_corr_data$std_est
    
    # beta_fe sample splitting bias correction with number of factor (data_driven)
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_logit(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left_data <- fe_fit_left$beta_est
    fe_fit_right = MLE_logit(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right_data <- fe_fit_right$beta_est
    fe_fit_upper = MLE_logit(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper_data <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_logit(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom_data <- fe_fit_bottom$beta_est
    
    beta_corr_sp_data <- 3* beta_fe_data - (beta_fe_left_data + beta_fe_right_data + beta_fe_upper_data + beta_fe_bottom_data)/2
    
    
    
    result <- list(
        beta_pool = beta_pool,
        beta_nnr = beta_nnr,
        num_factor_est = num_factor_est,
        beta_fe = beta_fe,
        beta_corr = beta_corr,
        beta_corr_sp = beta_corr_sp,
        std_beta_corr = std_beta_corr,
        beta_fe_data = beta_fe_data, 
        beta_corr_data = beta_corr_data, 
        beta_corr_sp_data = beta_corr_sp_data,
        std_beta_corr_data = std_beta_corr_data)
}
stopCluster(myCluster)   
summary_le_100_100 <- rearrange_result(raw_data, beta, num_X)


N = 200
T = 100

myCluster <- makeCluster(detectCores()-1, type = "FORK")
registerDoParallel(myCluster)
raw_data <- foreach (k=1:sim) %dopar% {
    set.seed(seed[k])
    data <- gen_data_logit(N, T, num_factor, num_Z, beta)
    
    # Pool
    Yvec <- data$Y %>% c() 
    X1vec <- data$X[[1]]%>% c()
    id <- rep(1:N, T) %>% matrix(N, T) %>% c()
    time <- rep(1:T, N) %>% matrix(T, N) %>% t() %>% c()
    df <- data.frame(Yvec, X1vec, id, time)
    
    X <- data$X
    Y <- data$Y
    
    pool <- glm(Yvec ~ X1vec, data = df, family = binomial(link = "logit"))
    
    beta_pool <- pool$coefficients[-1]
    beta_pool <- c(unlist(beta_pool))
    print(beta_pool)
    
    # Compute tuning parameter
    
    TW_model <- feglm(Yvec ~ X1vec | id + time, family = binomial(link = "logit"), data = df)
    Y_fit <- fitted(TW_model) %>% matrix(N, T)
    phi <- SVD(Y -Y_fit)$sigma[1] * (1 + delta)
    print(phi)
    
    nnr_fit = fit_logit(X, Y,  phi, 10, 300, 1e-8)
    beta_nnr <-  nnr_fit$beta_est
    
    sigma <-  SVD(nnr_fit$Theta_est)$sigma
    num_factor_est <- determine_num_factor(sigma, R_max)
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    
    fe_fit = MLE_logit(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe <- fe_fit$beta_est
    L_fe = fe_fit$L_est
    R_fe = fe_fit$R_est
    
    
    index <- Compute_index_with_LR(X, beta_fe, L_fe, R_fe)
    P <-  1 - 1/(exp(index) + 1) 
    phi <- SVD(Y -P)$sigma[1] * (1 + delta)
    
    
    
    # NNR 
    
    nnr_fit = fit_logit(X, Y,  phi, 10, 300, 1e-8)
    beta_nnr <-  nnr_fit$beta_est
    print(beta_pool)
    
    # determine the number of factor
    
    sigma <-  SVD(nnr_fit$Theta_est)$sigma
    num_factor_est <- determine_num_factor(sigma, R_max)
    print(num_factor_est)
    
    # beta_fe with correct number of factor
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    
    fe_fit = MLE_logit(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe <- fe_fit$beta_est
    L_fe = fe_fit$L_est
    R_fe = fe_fit$R_est
    
    # beta_fe analytical bias correction with correct number of factor
    
    bias_corr <- Compute_bias_corr_logit(Y, X, beta_fe, L_fe, R_fe, 0)
    print(bias_corr)
    beta_corr <- bias_corr$beta_corr
    if(is.na(beta_corr[1])){
        beta_corr <- beta_fe
    }
    std_beta_corr <- bias_corr$std_est
    
    # beta_fe sample splitting bias correction with correct number of factor
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_logit(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left <- fe_fit_left$beta_est
    fe_fit_right = MLE_logit(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right <- fe_fit_right$beta_est
    fe_fit_upper = MLE_logit(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_logit(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom <- fe_fit_bottom$beta_est
    
    beta_corr_sp <- 3* beta_fe - (beta_fe_left + beta_fe_right + beta_fe_upper + beta_fe_bottom)/2
    
    
    # beta_fe with number of factor (data_driven)
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor_est);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    if (num_factor_est==0){
        L_0 <- matrix(0, N, 1)
        R_0 <- matrix(0, T, 1)
    }
    
    fe_fit_data = MLE_logit(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe_data <- fe_fit_data$beta_est
    L_fe_data = fe_fit_data$L_est
    R_fe_data = fe_fit_data$R_est
    
    # beta_fe analytical bias correction with number of factor (data_driven)
    bias_corr_data <- Compute_bias_corr_logit(Y, X, beta_fe_data, L_fe_data, R_fe_data, 0)
    beta_corr_data <- bias_corr_data$beta_corr
    if(is.na(beta_corr_data[1])){
        beta_corr_data <- beta_fe_data
    }
    std_beta_corr_data <- bias_corr_data$std_est
    
    # beta_fe sample splitting bias correction with number of factor (data_driven)
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_logit(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left_data <- fe_fit_left$beta_est
    fe_fit_right = MLE_logit(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right_data <- fe_fit_right$beta_est
    fe_fit_upper = MLE_logit(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper_data <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_logit(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom_data <- fe_fit_bottom$beta_est
    
    beta_corr_sp_data <- 3* beta_fe_data - (beta_fe_left_data + beta_fe_right_data + beta_fe_upper_data + beta_fe_bottom_data)/2
    
    
    
    result <- list(
        beta_pool = beta_pool,
        beta_nnr = beta_nnr,
        num_factor_est = num_factor_est,
        beta_fe = beta_fe,
        beta_corr = beta_corr,
        beta_corr_sp = beta_corr_sp,
        std_beta_corr = std_beta_corr,
        beta_fe_data = beta_fe_data, 
        beta_corr_data = beta_corr_data, 
        beta_corr_sp_data = beta_corr_sp_data,
        std_beta_corr_data = std_beta_corr_data)
}
stopCluster(myCluster)      
summary_le_200_100 <- rearrange_result(raw_data, beta, num_X)

N = 200
T = 200

myCluster <- makeCluster(detectCores()-1, type = "FORK")
registerDoParallel(myCluster)
raw_data <- foreach (k=1:sim) %dopar% {
    set.seed(seed[k])
    data <- gen_data_logit(N, T, num_factor, num_Z, beta)
    
    # Pool
    Yvec <- data$Y %>% c() 
    X1vec <- data$X[[1]]%>% c()
    id <- rep(1:N, T) %>% matrix(N, T) %>% c()
    time <- rep(1:T, N) %>% matrix(T, N) %>% t() %>% c()
    df <- data.frame(Yvec, X1vec, id, time)
    
    X <- data$X
    Y <- data$Y
    
    pool <- glm(Yvec ~ X1vec, data = df, family = binomial(link = "logit"))
    
    beta_pool <- pool$coefficients[-1]
    beta_pool <- c(unlist(beta_pool))
    print(beta_pool)
    
    # Compute tuning parameter
    
    TW_model <- feglm(Yvec ~ X1vec | id + time, family = binomial(link = "logit"), data = df)
    Y_fit <- fitted(TW_model) %>% matrix(N, T)
    phi <- SVD(Y -Y_fit)$sigma[1] * (1 + delta)
    print(phi)
    
    nnr_fit = fit_logit(X, Y,  phi, 10, 300, 1e-8)
    beta_nnr <-  nnr_fit$beta_est
    
    sigma <-  SVD(nnr_fit$Theta_est)$sigma
    num_factor_est <- determine_num_factor(sigma, R_max)
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    
    fe_fit = MLE_logit(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe <- fe_fit$beta_est
    L_fe = fe_fit$L_est
    R_fe = fe_fit$R_est
    
    
    index <- Compute_index_with_LR(X, beta_fe, L_fe, R_fe)
    P <-  1 - 1/(exp(index) + 1) 
    phi <- SVD(Y -P)$sigma[1] * (1 + delta)
    
    
    
    # NNR 
    
    nnr_fit = fit_logit(X, Y,  phi, 10, 300, 1e-8)
    beta_nnr <-  nnr_fit$beta_est
    print(beta_pool)
    
    # determine the number of factor
    
    sigma <-  SVD(nnr_fit$Theta_est)$sigma
    num_factor_est <- determine_num_factor(sigma, R_max)
    print(num_factor_est)
    
    # beta_fe with correct number of factor
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    
    fe_fit = MLE_logit(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe <- fe_fit$beta_est
    L_fe = fe_fit$L_est
    R_fe = fe_fit$R_est
    
    # beta_fe analytical bias correction with correct number of factor
    
    bias_corr <- Compute_bias_corr_logit(Y, X, beta_fe, L_fe, R_fe, 0)
    print(bias_corr)
    beta_corr <- bias_corr$beta_corr
    if(is.na(beta_corr[1])){
        beta_corr <- beta_fe
    }
    std_beta_corr <- bias_corr$std_est
    
    # beta_fe sample splitting bias correction with correct number of factor
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_logit(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left <- fe_fit_left$beta_est
    fe_fit_right = MLE_logit(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right <- fe_fit_right$beta_est
    fe_fit_upper = MLE_logit(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_logit(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom <- fe_fit_bottom$beta_est
    
    beta_corr_sp <- 3* beta_fe - (beta_fe_left + beta_fe_right + beta_fe_upper + beta_fe_bottom)/2
    
    
    # beta_fe with number of factor (data_driven)
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor_est);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    if (num_factor_est==0){
        L_0 <- matrix(0, N, 1)
        R_0 <- matrix(0, T, 1)
    }
    
    fe_fit_data = MLE_logit(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe_data <- fe_fit_data$beta_est
    L_fe_data = fe_fit_data$L_est
    R_fe_data = fe_fit_data$R_est
    
    # beta_fe analytical bias correction with number of factor (data_driven)
    bias_corr_data <- Compute_bias_corr_logit(Y, X, beta_fe_data, L_fe_data, R_fe_data, 0)
    beta_corr_data <- bias_corr_data$beta_corr
    if(is.na(beta_corr_data[1])){
        beta_corr_data <- beta_fe_data
    }
    std_beta_corr_data <- bias_corr_data$std_est
    
    # beta_fe sample splitting bias correction with number of factor (data_driven)
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_logit(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left_data <- fe_fit_left$beta_est
    fe_fit_right = MLE_logit(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right_data <- fe_fit_right$beta_est
    fe_fit_upper = MLE_logit(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper_data <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_logit(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom_data <- fe_fit_bottom$beta_est
    
    beta_corr_sp_data <- 3* beta_fe_data - (beta_fe_left_data + beta_fe_right_data + beta_fe_upper_data + beta_fe_bottom_data)/2
    
    
    
    result <- list(
        beta_pool = beta_pool,
        beta_nnr = beta_nnr,
        num_factor_est = num_factor_est,
        beta_fe = beta_fe,
        beta_corr = beta_corr,
        beta_corr_sp = beta_corr_sp,
        std_beta_corr = std_beta_corr,
        beta_fe_data = beta_fe_data, 
        beta_corr_data = beta_corr_data, 
        beta_corr_sp_data = beta_corr_sp_data,
        std_beta_corr_data = std_beta_corr_data)
}
stopCluster(myCluster)      
summary_le_200_200 <- rearrange_result(raw_data, beta, num_X)

N = 1000
T = 200

myCluster <- makeCluster(detectCores()-1, type = "FORK")
registerDoParallel(myCluster)
raw_data <- foreach (k=1:sim) %dopar% {
    set.seed(seed[k])
    data <- gen_data_logit(N, T, num_factor, num_Z, beta)
    
    # Pool
    Yvec <- data$Y %>% c() 
    X1vec <- data$X[[1]]%>% c()
    id <- rep(1:N, T) %>% matrix(N, T) %>% c()
    time <- rep(1:T, N) %>% matrix(T, N) %>% t() %>% c()
    df <- data.frame(Yvec, X1vec, id, time)
    
    X <- data$X
    Y <- data$Y
    
    pool <- glm(Yvec ~ X1vec, data = df, family = binomial(link = "logit"))
    
    beta_pool <- pool$coefficients[-1]
    beta_pool <- c(unlist(beta_pool))
    print(beta_pool)
    
    # Compute tuning parameter
    
    TW_model <- feglm(Yvec ~ X1vec | id + time, family = binomial(link = "logit"), data = df)
    Y_fit <- fitted(TW_model) %>% matrix(N, T)
    phi <- SVD(Y -Y_fit)$sigma[1] * (1 + delta)
    print(phi)
    
    nnr_fit = fit_logit(X, Y,  phi, 10, 300, 1e-8)
    beta_nnr <-  nnr_fit$beta_est
    
    sigma <-  SVD(nnr_fit$Theta_est)$sigma
    num_factor_est <- determine_num_factor(sigma, R_max)
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    
    fe_fit = MLE_logit(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe <- fe_fit$beta_est
    L_fe = fe_fit$L_est
    R_fe = fe_fit$R_est
    
    
    index <- Compute_index_with_LR(X, beta_fe, L_fe, R_fe)
    P <-  1 - 1/(exp(index) + 1) 
    phi <- SVD(Y -P)$sigma[1] * (1 + delta)
    
    
    
    # NNR 
    
    nnr_fit = fit_logit(X, Y,  phi, 10, 300, 1e-8)
    beta_nnr <-  nnr_fit$beta_est
    print(beta_pool)
    
    # determine the number of factor
    
    sigma <-  SVD(nnr_fit$Theta_est)$sigma
    num_factor_est <- determine_num_factor(sigma, R_max)
    print(num_factor_est)
    
    # beta_fe with correct number of factor
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    
    fe_fit = MLE_logit(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe <- fe_fit$beta_est
    L_fe = fe_fit$L_est
    R_fe = fe_fit$R_est
    
    # beta_fe analytical bias correction with correct number of factor
    
    bias_corr <- Compute_bias_corr_logit(Y, X, beta_fe, L_fe, R_fe, 0)
    print(bias_corr)
    beta_corr <- bias_corr$beta_corr
    if(is.na(beta_corr[1])){
        beta_corr <- beta_fe
    }
    std_beta_corr <- bias_corr$std_est
    
    # beta_fe sample splitting bias correction with correct number of factor
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_logit(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left <- fe_fit_left$beta_est
    fe_fit_right = MLE_logit(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right <- fe_fit_right$beta_est
    fe_fit_upper = MLE_logit(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_logit(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom <- fe_fit_bottom$beta_est
    
    beta_corr_sp <- 3* beta_fe - (beta_fe_left + beta_fe_right + beta_fe_upper + beta_fe_bottom)/2
    
    
    # beta_fe with number of factor (data_driven)
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor_est);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    if (num_factor_est==0){
        L_0 <- matrix(0, N, 1)
        R_0 <- matrix(0, T, 1)
    }
    
    fe_fit_data = MLE_logit(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe_data <- fe_fit_data$beta_est
    L_fe_data = fe_fit_data$L_est
    R_fe_data = fe_fit_data$R_est
    
    # beta_fe analytical bias correction with number of factor (data_driven)
    bias_corr_data <- Compute_bias_corr_logit(Y, X, beta_fe_data, L_fe_data, R_fe_data, 0)
    beta_corr_data <- bias_corr_data$beta_corr
    if(is.na(beta_corr_data[1])){
        beta_corr_data <- beta_fe_data
    }
    std_beta_corr_data <- bias_corr_data$std_est
    
    # beta_fe sample splitting bias correction with number of factor (data_driven)
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_logit(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left_data <- fe_fit_left$beta_est
    fe_fit_right = MLE_logit(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right_data <- fe_fit_right$beta_est
    fe_fit_upper = MLE_logit(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper_data <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_logit(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom_data <- fe_fit_bottom$beta_est
    
    beta_corr_sp_data <- 3* beta_fe_data - (beta_fe_left_data + beta_fe_right_data + beta_fe_upper_data + beta_fe_bottom_data)/2
    
    
    
    result <- list(
        beta_pool = beta_pool,
        beta_nnr = beta_nnr,
        num_factor_est = num_factor_est,
        beta_fe = beta_fe,
        beta_corr = beta_corr,
        beta_corr_sp = beta_corr_sp,
        std_beta_corr = std_beta_corr,
        beta_fe_data = beta_fe_data, 
        beta_corr_data = beta_corr_data, 
        beta_corr_sp_data = beta_corr_sp_data,
        std_beta_corr_data = std_beta_corr_data)
}
stopCluster(myCluster)      
summary_le_1000_200 <- rearrange_result(raw_data, beta, num_X)





#=================================================

beta = c(0.0)


seed <- sample.int(100000, sim, replace = FALSE)

N = 50
T = 40

delta = 0.2

myCluster <- makeCluster(detectCores()-1, type = "FORK")
registerDoParallel(myCluster)
raw_data <- foreach (k=1:sim) %dopar% {
    set.seed(seed[k])
    data <- gen_data_poisson(N, T, num_factor, num_Z, beta)
    
    Yvec <- data$Y %>% c() 
    X1vec <- data$X[[1]]%>% c()
    id <- rep(1:N, T) %>% matrix(N, T) %>% c()
    time <- rep(1:T, N) %>% matrix(T, N) %>% t() %>% c()
    df <- data.frame(Yvec, X1vec, id, time)
    
    X <- data$X
    Y <- data$Y
    
    pool <- glm(Yvec ~ X1vec, data = df, family = poisson(link = "log"))
    
    beta_pool <- pool$coefficients[-1]
    beta_pool <- c(unlist(beta_pool))
    print(beta_pool)
    
    # Compute tuning parameter
    
    TW_model <- fepois(Yvec ~ X1vec | id + time,  data = df)
    Y_fit <- fitted(TW_model) %>% matrix(N, T)
    phi <- SVD(Y -Y_fit)$sigma[1] * (1 + delta)
    print(sqrt(N) * mean(Y))
    print(phi)
    
    nnr_fit = fit_poisson(X, Y,  phi, 1, 1000, 1e-8)
    beta_nnr <-  nnr_fit$beta_est
    
    sigma <-  SVD(nnr_fit$Theta_est)$sigma
    num_factor_est <- determine_num_factor(sigma, R_max)
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    
    fe_fit = MLE_poisson(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe <- fe_fit$beta_est
    L_fe = fe_fit$L_est
    R_fe = fe_fit$R_est
    
    
    index <- Compute_index_with_LR(X, beta_fe, L_fe, R_fe)
    index <- exp(index)
    phi <- SVD(Y -index)$sigma[1] * (1 + delta)
    
    
    
    
    
    # NNR 
    
    nnr_fit = fit_poisson(X, Y,  phi, 1, 1000, 1e-8)
    beta_nnr <-  nnr_fit$beta_est
    print(beta_pool)
    
    # determine the number of factor
    
    sigma <-  SVD(nnr_fit$Theta_est)$sigma
    num_factor_est <- determine_num_factor(sigma, R_max)
    print(num_factor_est)
    
    # beta_fe with correct number of factor
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    
    fe_fit = MLE_poisson(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe <- fe_fit$beta_est
    L_fe = fe_fit$L_est
    R_fe = fe_fit$R_est
    
    # beta_fe analytical bias correction with correct number of factor
    
    
    bias_corr <- Compute_bias_corr_poisson(Y, X, beta_fe, L_fe, R_fe, 0)
    beta_corr <- bias_corr$beta_corr
    if(is.na(beta_corr[1])){
        beta_corr <- beta_fe
    }
    std_beta_corr <- bias_corr$std_est
    
    # beta_fe sample splitting bias correction with correct number of factor
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_poisson(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left <- fe_fit_left$beta_est
    fe_fit_right = MLE_poisson(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right <- fe_fit_right$beta_est
    fe_fit_upper = MLE_poisson(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_poisson(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom <- fe_fit_bottom$beta_est
    
    beta_corr_sp <- 3* beta_fe - (beta_fe_left + beta_fe_right + beta_fe_upper + beta_fe_bottom)/2
    
    
    # beta_fe with number of factor (data_driven)
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor_est);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    if (num_factor_est==0){
        L_0 <- matrix(0, N, 1)
        R_0 <- matrix(0, T, 1)
    }
    
    fe_fit_data = MLE_poisson(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe_data <- fe_fit_data$beta_est
    L_fe_data = fe_fit_data$L_est
    R_fe_data = fe_fit_data$R_est
    
    # beta_fe analytical bias correction with number of factor (data_driven)
    bias_corr_data <- Compute_bias_corr_poisson(Y, X, beta_fe_data, L_fe_data, R_fe_data, 0)
    beta_corr_data <- bias_corr_data$beta_corr
    if(is.na(beta_corr_data[1])){
        beta_corr_data <- beta_fe_data
    }
    std_beta_corr_data <- bias_corr_data$std_est
    
    # beta_fe sample splitting bias correction with number of factor (data_driven)
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_poisson(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left_data <- fe_fit_left$beta_est
    fe_fit_right = MLE_poisson(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right_data <- fe_fit_right$beta_est
    fe_fit_upper = MLE_poisson(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper_data <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_poisson(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom_data <- fe_fit_bottom$beta_est
    
    beta_corr_sp_data <- 3* beta_fe_data - (beta_fe_left_data + beta_fe_right_data + beta_fe_upper_data + beta_fe_bottom_data)/2
    
    
    
    result <- list(
        beta_pool = beta_pool,
        beta_nnr = beta_nnr,
        num_factor_est = num_factor_est,
        beta_fe = beta_fe,
        beta_corr = beta_corr,
        beta_corr_sp = beta_corr_sp,
        std_beta_corr = std_beta_corr,
        beta_fe_data = beta_fe_data, 
        beta_corr_data = beta_corr_data, 
        beta_corr_sp_data = beta_corr_sp_data,
        std_beta_corr_data = std_beta_corr_data)
    
}
stopCluster(myCluster)    
summary_pe_50_40 <- rearrange_result(raw_data, beta, num_X)



seed <- sample.int(100000, sim, replace = FALSE)

N = 100
T = 40

myCluster <- makeCluster(detectCores()-1, type = "FORK")
registerDoParallel(myCluster)
raw_data <- foreach (k=1:sim) %dopar% {
    set.seed(seed[k])
    data <- gen_data_poisson(N, T, num_factor, num_Z, beta)
    
    Yvec <- data$Y %>% c() 
    X1vec <- data$X[[1]]%>% c()
    id <- rep(1:N, T) %>% matrix(N, T) %>% c()
    time <- rep(1:T, N) %>% matrix(T, N) %>% t() %>% c()
    df <- data.frame(Yvec, X1vec, id, time)
    
    X <- data$X
    Y <- data$Y
    
    pool <- glm(Yvec ~ X1vec, data = df, family = poisson(link = "log"))
    
    beta_pool <- pool$coefficients[-1]
    beta_pool <- c(unlist(beta_pool))
    print(beta_pool)
    
    # Compute tuning parameter
    
    TW_model <- fepois(Yvec ~ X1vec | id + time,  data = df)
    Y_fit <- fitted(TW_model) %>% matrix(N, T)
    phi <- SVD(Y -Y_fit)$sigma[1] * (1 + delta)
    print(sqrt(N) * mean(Y))
    print(phi)
    
    nnr_fit = fit_poisson(X, Y,  phi, 1, 1000, 1e-8)
    beta_nnr <-  nnr_fit$beta_est
    
    sigma <-  SVD(nnr_fit$Theta_est)$sigma
    num_factor_est <- determine_num_factor(sigma, R_max)
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    
    fe_fit = MLE_poisson(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe <- fe_fit$beta_est
    L_fe = fe_fit$L_est
    R_fe = fe_fit$R_est
    
    
    index <- Compute_index_with_LR(X, beta_fe, L_fe, R_fe)
    index <- exp(index)
    phi <- SVD(Y -index)$sigma[1] * (1 + delta)
    
    
    
    
    
    # NNR 
    
    nnr_fit = fit_poisson(X, Y,  phi, 1, 1000, 1e-8)
    beta_nnr <-  nnr_fit$beta_est
    print(beta_pool)
    
    # determine the number of factor
    
    sigma <-  SVD(nnr_fit$Theta_est)$sigma
    num_factor_est <- determine_num_factor(sigma, R_max)
    print(num_factor_est)
    
    # beta_fe with correct number of factor
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    
    fe_fit = MLE_poisson(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe <- fe_fit$beta_est
    L_fe = fe_fit$L_est
    R_fe = fe_fit$R_est
    
    # beta_fe analytical bias correction with correct number of factor
    
    
    bias_corr <- Compute_bias_corr_poisson(Y, X, beta_fe, L_fe, R_fe, 0)
    beta_corr <- bias_corr$beta_corr
    if(is.na(beta_corr[1])){
        beta_corr <- beta_fe
    }
    std_beta_corr <- bias_corr$std_est
    
    # beta_fe sample splitting bias correction with correct number of factor
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_poisson(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left <- fe_fit_left$beta_est
    fe_fit_right = MLE_poisson(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right <- fe_fit_right$beta_est
    fe_fit_upper = MLE_poisson(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_poisson(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom <- fe_fit_bottom$beta_est
    
    beta_corr_sp <- 3* beta_fe - (beta_fe_left + beta_fe_right + beta_fe_upper + beta_fe_bottom)/2
    
    
    # beta_fe with number of factor (data_driven)
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor_est);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    if (num_factor_est==0){
        L_0 <- matrix(0, N, 1)
        R_0 <- matrix(0, T, 1)
    }
    
    fe_fit_data = MLE_poisson(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe_data <- fe_fit_data$beta_est
    L_fe_data = fe_fit_data$L_est
    R_fe_data = fe_fit_data$R_est
    
    # beta_fe analytical bias correction with number of factor (data_driven)
    bias_corr_data <- Compute_bias_corr_poisson(Y, X, beta_fe_data, L_fe_data, R_fe_data, 0)
    beta_corr_data <- bias_corr_data$beta_corr
    if(is.na(beta_corr_data[1])){
        beta_corr_data <- beta_fe_data
    }
    std_beta_corr_data <- bias_corr_data$std_est
    
    # beta_fe sample splitting bias correction with number of factor (data_driven)
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_poisson(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left_data <- fe_fit_left$beta_est
    fe_fit_right = MLE_poisson(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right_data <- fe_fit_right$beta_est
    fe_fit_upper = MLE_poisson(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper_data <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_poisson(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom_data <- fe_fit_bottom$beta_est
    
    beta_corr_sp_data <- 3* beta_fe_data - (beta_fe_left_data + beta_fe_right_data + beta_fe_upper_data + beta_fe_bottom_data)/2
    
    
    
    result <- list(
        beta_pool = beta_pool,
        beta_nnr = beta_nnr,
        num_factor_est = num_factor_est,
        beta_fe = beta_fe,
        beta_corr = beta_corr,
        beta_corr_sp = beta_corr_sp,
        std_beta_corr = std_beta_corr,
        beta_fe_data = beta_fe_data, 
        beta_corr_data = beta_corr_data, 
        beta_corr_sp_data = beta_corr_sp_data,
        std_beta_corr_data = std_beta_corr_data)
    
}
stopCluster(myCluster)            
summary_pe_100_40 <- rearrange_result(raw_data, beta, num_X)




seed <- sample.int(100000, sim, replace = FALSE)

N = 200
T = 40

myCluster <- makeCluster(detectCores()-1, type = "FORK")
registerDoParallel(myCluster)
raw_data <- foreach (k=1:sim) %dopar% {
    set.seed(seed[k])
    data <- gen_data_poisson(N, T, num_factor, num_Z, beta)
    
    Yvec <- data$Y %>% c() 
    X1vec <- data$X[[1]]%>% c()
    id <- rep(1:N, T) %>% matrix(N, T) %>% c()
    time <- rep(1:T, N) %>% matrix(T, N) %>% t() %>% c()
    df <- data.frame(Yvec, X1vec, id, time)
    
    X <- data$X
    Y <- data$Y
    
    pool <- glm(Yvec ~ X1vec, data = df, family = poisson(link = "log"))
    
    beta_pool <- pool$coefficients[-1]
    beta_pool <- c(unlist(beta_pool))
    print(beta_pool)
    
    # Compute tuning parameter
    
    TW_model <- fepois(Yvec ~ X1vec | id + time,  data = df)
    Y_fit <- fitted(TW_model) %>% matrix(N, T)
    phi <- SVD(Y -Y_fit)$sigma[1] * (1 + delta)
    print(sqrt(N) * mean(Y))
    print(phi)
    
    nnr_fit = fit_poisson(X, Y,  phi, 10, 300, 1e-8)
    beta_nnr <-  nnr_fit$beta_est
    
    sigma <-  SVD(nnr_fit$Theta_est)$sigma
    num_factor_est <- determine_num_factor(sigma, R_max)
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    
    fe_fit = MLE_poisson(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe <- fe_fit$beta_est
    L_fe = fe_fit$L_est
    R_fe = fe_fit$R_est
    
    
    index <- Compute_index_with_LR(X, beta_fe, L_fe, R_fe)
    index <- exp(index)
    phi <- SVD(Y -index)$sigma[1] * (1 + delta)
    
    
    
    
    
    # NNR 
    
    nnr_fit = fit_poisson(X, Y,  phi, 10, 300, 1e-8)
    beta_nnr <-  nnr_fit$beta_est
    print(beta_pool)
    
    # determine the number of factor
    
    sigma <-  SVD(nnr_fit$Theta_est)$sigma
    num_factor_est <- determine_num_factor(sigma, R_max)
    print(num_factor_est)
    
    # beta_fe with correct number of factor
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    
    fe_fit = MLE_poisson(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe <- fe_fit$beta_est
    L_fe = fe_fit$L_est
    R_fe = fe_fit$R_est
    
    # beta_fe analytical bias correction with correct number of factor
    
    
    bias_corr <- Compute_bias_corr_poisson(Y, X, beta_fe, L_fe, R_fe, 0)
    beta_corr <- bias_corr$beta_corr
    std_beta_corr <- bias_corr$std_est
    
    # beta_fe sample splitting bias correction with correct number of factor
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_poisson(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left <- fe_fit_left$beta_est
    fe_fit_right = MLE_poisson(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right <- fe_fit_right$beta_est
    fe_fit_upper = MLE_poisson(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_poisson(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom <- fe_fit_bottom$beta_est
    
    beta_corr_sp <- 3* beta_fe - (beta_fe_left + beta_fe_right + beta_fe_upper + beta_fe_bottom)/2
    
    
    # beta_fe with number of factor (data_driven)
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor_est);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    if (num_factor_est==0){
        L_0 <- matrix(0, N, 1)
        R_0 <- matrix(0, T, 1)
    }
    
    fe_fit_data = MLE_poisson(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe_data <- fe_fit_data$beta_est
    L_fe_data = fe_fit_data$L_est
    R_fe_data = fe_fit_data$R_est
    
    # beta_fe analytical bias correction with number of factor (data_driven)
    bias_corr_data <- Compute_bias_corr_poisson(Y, X, beta_fe_data, L_fe_data, R_fe_data, 0)
    beta_corr_data <- bias_corr_data$beta_corr
    std_beta_corr_data <- bias_corr_data$std_est
    
    # beta_fe sample splitting bias correction with number of factor (data_driven)
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_poisson(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left_data <- fe_fit_left$beta_est
    fe_fit_right = MLE_poisson(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right_data <- fe_fit_right$beta_est
    fe_fit_upper = MLE_poisson(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper_data <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_poisson(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom_data <- fe_fit_bottom$beta_est
    
    beta_corr_sp_data <- 3* beta_fe_data - (beta_fe_left_data + beta_fe_right_data + beta_fe_upper_data + beta_fe_bottom_data)/2
    
    
    
    result <- list(
        beta_pool = beta_pool,
        beta_nnr = beta_nnr,
        num_factor_est = num_factor_est,
        beta_fe = beta_fe,
        beta_corr = beta_corr,
        beta_corr_sp = beta_corr_sp,
        std_beta_corr = std_beta_corr,
        beta_fe_data = beta_fe_data, 
        beta_corr_data = beta_corr_data, 
        beta_corr_sp_data = beta_corr_sp_data,
        std_beta_corr_data = std_beta_corr_data)
    
}
stopCluster(myCluster)     
summary_pe_200_40 <- rearrange_result(raw_data, beta, num_X)


seed <- sample.int(100000, sim, replace = FALSE)

N = 100
T = 100

myCluster <- makeCluster(detectCores()-1, type = "FORK")
registerDoParallel(myCluster)
raw_data <- foreach (k=1:sim) %dopar% {
    set.seed(seed[k])
    data <- gen_data_poisson(N, T, num_factor, num_Z, beta)
    
    Yvec <- data$Y %>% c() 
    X1vec <- data$X[[1]]%>% c()
    id <- rep(1:N, T) %>% matrix(N, T) %>% c()
    time <- rep(1:T, N) %>% matrix(T, N) %>% t() %>% c()
    df <- data.frame(Yvec, X1vec, id, time)
    
    X <- data$X
    Y <- data$Y
    
    pool <- glm(Yvec ~ X1vec, data = df, family = poisson(link = "log"))
    
    beta_pool <- pool$coefficients[-1]
    beta_pool <- c(unlist(beta_pool))
    print(beta_pool)
    
    # Compute tuning parameter
    
    TW_model <- fepois(Yvec ~ X1vec | id + time,  data = df)
    Y_fit <- fitted(TW_model) %>% matrix(N, T)
    phi <- SVD(Y -Y_fit)$sigma[1] * (1 + delta)
    print(sqrt(N) * mean(Y))
    print(phi)
    
    nnr_fit = fit_poisson(X, Y,  phi, 1, 1000, 1e-8)
    beta_nnr <-  nnr_fit$beta_est
    
    sigma <-  SVD(nnr_fit$Theta_est)$sigma
    num_factor_est <- determine_num_factor(sigma, R_max)
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    
    fe_fit = MLE_poisson(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe <- fe_fit$beta_est
    L_fe = fe_fit$L_est
    R_fe = fe_fit$R_est
    
    
    index <- Compute_index_with_LR(X, beta_fe, L_fe, R_fe)
    index <- exp(index)
    phi <- SVD(Y -index)$sigma[1] * (1 + delta)
    
    
    
    
    
    # NNR 
    
    nnr_fit = fit_poisson(X, Y,  phi, 1, 1000, 1e-8)
    beta_nnr <-  nnr_fit$beta_est
    print(beta_pool)
    
    # determine the number of factor
    
    sigma <-  SVD(nnr_fit$Theta_est)$sigma
    num_factor_est <- determine_num_factor(sigma, R_max)
    print(num_factor_est)
    
    # beta_fe with correct number of factor
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    
    fe_fit = MLE_poisson(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe <- fe_fit$beta_est
    L_fe = fe_fit$L_est
    R_fe = fe_fit$R_est
    
    # beta_fe analytical bias correction with correct number of factor
    
    
    bias_corr <- Compute_bias_corr_poisson(Y, X, beta_fe, L_fe, R_fe, 0)
    beta_corr <- bias_corr$beta_corr
    if(is.na(beta_corr[1])){
        beta_corr <- beta_fe
    }
    std_beta_corr <- bias_corr$std_est
    
    # beta_fe sample splitting bias correction with correct number of factor
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_poisson(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left <- fe_fit_left$beta_est
    fe_fit_right = MLE_poisson(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right <- fe_fit_right$beta_est
    fe_fit_upper = MLE_poisson(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_poisson(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom <- fe_fit_bottom$beta_est
    
    beta_corr_sp <- 3* beta_fe - (beta_fe_left + beta_fe_right + beta_fe_upper + beta_fe_bottom)/2
    
    
    # beta_fe with number of factor (data_driven)
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor_est);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    if (num_factor_est==0){
        L_0 <- matrix(0, N, 1)
        R_0 <- matrix(0, T, 1)
    }
    
    fe_fit_data = MLE_poisson(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe_data <- fe_fit_data$beta_est
    L_fe_data = fe_fit_data$L_est
    R_fe_data = fe_fit_data$R_est
    
    # beta_fe analytical bias correction with number of factor (data_driven)
    bias_corr_data <- Compute_bias_corr_poisson(Y, X, beta_fe_data, L_fe_data, R_fe_data, 0)
    beta_corr_data <- bias_corr_data$beta_corr
    if(is.na(beta_corr_data[1])){
        beta_corr_data <- beta_fe_data
    }
    std_beta_corr_data <- bias_corr_data$std_est
    
    # beta_fe sample splitting bias correction with number of factor (data_driven)
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_poisson(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left_data <- fe_fit_left$beta_est
    fe_fit_right = MLE_poisson(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right_data <- fe_fit_right$beta_est
    fe_fit_upper = MLE_poisson(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper_data <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_poisson(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom_data <- fe_fit_bottom$beta_est
    
    beta_corr_sp_data <- 3* beta_fe_data - (beta_fe_left_data + beta_fe_right_data + beta_fe_upper_data + beta_fe_bottom_data)/2
    
    
    
    result <- list(
        beta_pool = beta_pool,
        beta_nnr = beta_nnr,
        num_factor_est = num_factor_est,
        beta_fe = beta_fe,
        beta_corr = beta_corr,
        beta_corr_sp = beta_corr_sp,
        std_beta_corr = std_beta_corr,
        beta_fe_data = beta_fe_data, 
        beta_corr_data = beta_corr_data, 
        beta_corr_sp_data = beta_corr_sp_data,
        std_beta_corr_data = std_beta_corr_data)
    
}
stopCluster(myCluster)  
summary_pe_100_100 <- rearrange_result(raw_data, beta, num_X)


seed <- sample.int(100000, sim, replace = FALSE)

N = 200
T = 200

myCluster <- makeCluster(detectCores()-1, type = "FORK")
registerDoParallel(myCluster)
raw_data <- foreach (k=1:sim) %dopar% {
    set.seed(seed[k])
    data <- gen_data_poisson(N, T, num_factor, num_Z, beta)
    
    Yvec <- data$Y %>% c() 
    X1vec <- data$X[[1]]%>% c()
    id <- rep(1:N, T) %>% matrix(N, T) %>% c()
    time <- rep(1:T, N) %>% matrix(T, N) %>% t() %>% c()
    df <- data.frame(Yvec, X1vec, id, time)
    
    X <- data$X
    Y <- data$Y
    
    pool <- glm(Yvec ~ X1vec, data = df, family = poisson(link = "log"))
    
    beta_pool <- pool$coefficients[-1]
    beta_pool <- c(unlist(beta_pool))
    print(beta_pool)
    
    # Compute tuning parameter
    
    TW_model <- fepois(Yvec ~ X1vec | id + time,  data = df)
    Y_fit <- fitted(TW_model) %>% matrix(N, T)
    phi <- SVD(Y -Y_fit)$sigma[1] * (1 + delta)
    print(sqrt(N) * mean(Y))
    print(phi)
    
    nnr_fit = fit_poisson(X, Y,  phi, 1, 1000, 1e-8)
    beta_nnr <-  nnr_fit$beta_est
    
    sigma <-  SVD(nnr_fit$Theta_est)$sigma
    num_factor_est <- determine_num_factor(sigma, R_max)
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    
    fe_fit = MLE_poisson(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe <- fe_fit$beta_est
    L_fe = fe_fit$L_est
    R_fe = fe_fit$R_est
    
    
    index <- Compute_index_with_LR(X, beta_fe, L_fe, R_fe)
    index <- exp(index)
    phi <- SVD(Y -index)$sigma[1] * (1 + delta)
    
    
    
    
    
    # NNR 
    
    nnr_fit = fit_poisson(X, Y,  phi, 1, 1000, 1e-8)
    beta_nnr <-  nnr_fit$beta_est
    print(beta_pool)
    
    # determine the number of factor
    
    sigma <-  SVD(nnr_fit$Theta_est)$sigma
    num_factor_est <- determine_num_factor(sigma, R_max)
    print(num_factor_est)
    
    # beta_fe with correct number of factor
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    
    fe_fit = MLE_poisson(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe <- fe_fit$beta_est
    L_fe = fe_fit$L_est
    R_fe = fe_fit$R_est
    
    # beta_fe analytical bias correction with correct number of factor
    
    
    bias_corr <- Compute_bias_corr_poisson(Y, X, beta_fe, L_fe, R_fe, 0)
    beta_corr <- bias_corr$beta_corr
    if(is.na(beta_corr[1])){
        beta_corr <- beta_fe
    }
    std_beta_corr <- bias_corr$std_est
    
    # beta_fe sample splitting bias correction with correct number of factor
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_poisson(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left <- fe_fit_left$beta_est
    fe_fit_right = MLE_poisson(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right <- fe_fit_right$beta_est
    fe_fit_upper = MLE_poisson(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_poisson(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom <- fe_fit_bottom$beta_est
    
    beta_corr_sp <- 3* beta_fe - (beta_fe_left + beta_fe_right + beta_fe_upper + beta_fe_bottom)/2
    
    
    # beta_fe with number of factor (data_driven)
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor_est);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    if (num_factor_est==0){
        L_0 <- matrix(0, N, 1)
        R_0 <- matrix(0, T, 1)
    }
    
    fe_fit_data = MLE_poisson(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe_data <- fe_fit_data$beta_est
    L_fe_data = fe_fit_data$L_est
    R_fe_data = fe_fit_data$R_est
    
    # beta_fe analytical bias correction with number of factor (data_driven)
    bias_corr_data <- Compute_bias_corr_poisson(Y, X, beta_fe_data, L_fe_data, R_fe_data, 0)
    beta_corr_data <- bias_corr_data$beta_corr
    if(is.na(beta_corr_data[1])){
        beta_corr_data <- beta_fe_data
    }
    std_beta_corr_data <- bias_corr_data$std_est
    
    # beta_fe sample splitting bias correction with number of factor (data_driven)
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_poisson(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left_data <- fe_fit_left$beta_est
    fe_fit_right = MLE_poisson(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right_data <- fe_fit_right$beta_est
    fe_fit_upper = MLE_poisson(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper_data <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_poisson(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom_data <- fe_fit_bottom$beta_est
    
    beta_corr_sp_data <- 3* beta_fe_data - (beta_fe_left_data + beta_fe_right_data + beta_fe_upper_data + beta_fe_bottom_data)/2
    
    
    
    result <- list(
        beta_pool = beta_pool,
        beta_nnr = beta_nnr,
        num_factor_est = num_factor_est,
        beta_fe = beta_fe,
        beta_corr = beta_corr,
        beta_corr_sp = beta_corr_sp,
        std_beta_corr = std_beta_corr,
        beta_fe_data = beta_fe_data, 
        beta_corr_data = beta_corr_data, 
        beta_corr_sp_data = beta_corr_sp_data,
        std_beta_corr_data = std_beta_corr_data)
    
}
stopCluster(myCluster)   
summary_pe_200_200 <- rearrange_result(raw_data, beta, num_X)



#=======================================================



library(NNRPanel)
library(dplyr)
library(fixest)
library(doParallel)

set.seed(199505)

num_Z = 1
num_W = 1
num_X = num_Z + num_W
num_factor = 2
beta_W = c(0.5)
beta_Z = c(0.2)
beta = c(0.5, 0.2)
R_max = 5

delta = 0.05


sim = 1000


seed <- sample.int(100000, sim, replace = FALSE)

N = 50
T = 40

myCluster <- makeCluster(detectCores()-1, type = "FORK")
registerDoParallel(myCluster)
raw_data <- foreach (k=1:sim)%dopar%{
    set.seed(seed[k])
    
    data <- gen_data_dynamic_logit(N, T, num_factor, num_Z, num_W, beta_W, beta_Z)
    
    # Pool
    Yvec <- data$Y %>% c() 
    Xvec <- data$X[[1]] %>% c()
    Zvec <- data$X[[2]]%>% c()
    fvec <- data$interact_f %>% c()
    id <- rep(1:N, T) %>% matrix(N, T) %>% c()
    time <- rep(1:T, N) %>% matrix(T, N) %>% t() %>% c()
    df <- data.frame(Yvec, Xvec, Zvec, id, time)
    
    X <- data$X
    Y <- data$Y
    print(mean(Y))
    print(mean(X[[1]]))
    print(mean(X[[2]]))
    
    pool <- glm(Yvec ~ Xvec + Zvec, data = df, family = binomial(link = "logit"))
    summary(pool)
    beta_pool <- pool$coefficients[-1]
    beta_pool <- c(unlist(beta_pool))
    
    
    # Compute tuning parameter
    
    TW_model <- feglm(Yvec ~ Xvec + Zvec | id + time, family = binomial(link = "logit"), data = df)
    Y_fit <- fitted(TW_model) %>% matrix(N, T)
    phi <- SVD(Y -Y_fit)$sigma[1] * (1 + delta)
    print(phi)
    
    nnr_fit = fit_logit(X, Y,  phi, 10, 300, 1e-8)
    beta_nnr <-  nnr_fit$beta_est
    
    sigma <-  SVD(nnr_fit$Theta_est)$sigma
    num_factor_est <- determine_num_factor(sigma, R_max)
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    
    fe_fit = MLE_logit(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe <- fe_fit$beta_est
    L_fe = fe_fit$L_est
    R_fe = fe_fit$R_est
    
    
    index <- Compute_index_with_LR(X, beta_fe, L_fe, R_fe)
    P <-  1 - 1/(exp(index) + 1) 
    phi <- SVD(Y -P)$sigma[1] * (1 + delta)
    
    
    
    
    # Compute tuning parameter
    
    TW_model <- feglm(Yvec ~ Xvec + Zvec | id + time, family = binomial(link = "logit"), data = df)
    Y_fit <- fitted(TW_model) %>% matrix(N, T)
    phi <- SVD(Y -Y_fit)$sigma[1] * (1 + delta)
    print(phi)
    
    nnr_fit = fit_logit(X, Y,  phi, 10, 300, 1e-8)
    beta_nnr <-  nnr_fit$beta_est
    
    sigma <-  SVD(nnr_fit$Theta_est)$sigma
    num_factor_est <- determine_num_factor(sigma, R_max)
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    
    fe_fit = MLE_logit(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe <- fe_fit$beta_est
    L_fe = fe_fit$L_est
    R_fe = fe_fit$R_est
    
    
    index <- Compute_index_with_LR(X, beta_fe, L_fe, R_fe)
    P <-  1 - 1/(exp(index) + 1) 
    phi <- SVD(Y -P)$sigma[1] * (1 + delta)
    
    
    # NNR 
    
    nnr_fit = fit_logit(X, Y,  mean(Y)*sqrt(N), 10, 1000, 1e-8)
    beta_nnr <-  nnr_fit$beta_est
    print(beta_pool)
    
    # determine the number of factor
    
    sigma <-  SVD(nnr_fit$Theta_est)$sigma
    num_factor_est <- determine_num_factor(sigma, R_max)
    print(num_factor_est)
    
    # beta_fe with correct number of factor
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    
    fe_fit = MLE_logit(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe <- fe_fit$beta_est
    L_fe = fe_fit$L_est
    R_fe = fe_fit$R_est
    
    # beta_fe analytical bias correction with correct number of factor
    
    
    bias_corr <- Compute_bias_corr_logit(Y, X, beta_fe, L_fe, R_fe, 1)
    beta_corr <- bias_corr$beta_corr
    if(is.na(beta_corr[1])){
        beta_corr <- beta_fe
    }
    std_beta_corr <- bias_corr$std_est
    
    # beta_fe sample splitting bias correction with correct number of factor
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_logit(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left <- fe_fit_left$beta_est
    fe_fit_right = MLE_logit(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right <- fe_fit_right$beta_est
    fe_fit_upper = MLE_logit(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_logit(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom <- fe_fit_bottom$beta_est
    
    beta_corr_sp <- 3* beta_fe - (beta_fe_left + beta_fe_right + beta_fe_upper + beta_fe_bottom)/2
    
    
    # beta_fe with number of factor (data_driven)
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor_est);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    if (num_factor_est==0){
        L_0 <- matrix(0, N, 1)
        R_0 <- matrix(0, T, 1)
    }
    
    fe_fit_data = MLE_logit(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe_data <- fe_fit_data$beta_est
    L_fe_data = fe_fit_data$L_est
    R_fe_data = fe_fit_data$R_est
    
    # beta_fe analytical bias correction with number of factor (data_driven)
    bias_corr_data <- Compute_bias_corr_logit(Y, X, beta_fe_data, L_fe_data, R_fe_data, 1)
    beta_corr_data <- bias_corr_data$beta_corr
    if(is.na(beta_corr_data[1])){
        beta_corr_data <- beta_fe_data
    }
    std_beta_corr_data <- bias_corr_data$std_est
    
    # beta_fe sample splitting bias correction with number of factor (data_driven)
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_logit(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left_data <- fe_fit_left$beta_est
    fe_fit_right = MLE_logit(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right_data <- fe_fit_right$beta_est
    fe_fit_upper = MLE_logit(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper_data <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_logit(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom_data <- fe_fit_bottom$beta_est
    
    beta_corr_sp_data <- 3* beta_fe_data - (beta_fe_left_data + beta_fe_right_data + beta_fe_upper_data + beta_fe_bottom_data)/2
    
    
    
    result <- list(
        beta_pool = beta_pool,
        beta_nnr = beta_nnr,
        num_factor_est = num_factor_est,
        beta_fe = beta_fe,
        beta_corr = beta_corr,
        beta_corr_sp = beta_corr_sp,
        std_beta_corr = std_beta_corr,
        beta_fe_data = beta_fe_data, 
        beta_corr_data = beta_corr_data, 
        beta_corr_sp_data = beta_corr_sp_data,
        std_beta_corr_data = std_beta_corr_data)
    
}
stopCluster(myCluster)
summary_ld_50_40 <- rearrange_result(raw_data, beta, num_X)


seed <- sample.int(100000, sim, replace = FALSE)

N = 100
T = 40

myCluster <- makeCluster(detectCores()-1, type = "FORK")
registerDoParallel(myCluster)
raw_data <- foreach (k=1:sim)%dopar%{
    set.seed(seed[k])
    
    data <- gen_data_dynamic_logit(N, T, num_factor, num_Z, num_W, beta_W, beta_Z)
    
    # Pool
    Yvec <- data$Y %>% c() 
    Xvec <- data$X[[1]] %>% c()
    Zvec <- data$X[[2]]%>% c()
    fvec <- data$interact_f %>% c()
    id <- rep(1:N, T) %>% matrix(N, T) %>% c()
    time <- rep(1:T, N) %>% matrix(T, N) %>% t() %>% c()
    df <- data.frame(Yvec, Xvec, Zvec, id, time)
    
    X <- data$X
    Y <- data$Y
    print(mean(Y))
    print(mean(X[[1]]))
    print(mean(X[[2]]))
    
    pool <- glm(Yvec ~ Xvec + Zvec, data = df, family = binomial(link = "logit"))
    summary(pool)
    beta_pool <- pool$coefficients[-1]
    beta_pool <- c(unlist(beta_pool))
    
    
    # Compute tuning parameter
    
    TW_model <- feglm(Yvec ~ Xvec + Zvec | id + time, family = binomial(link = "logit"), data = df)
    Y_fit <- fitted(TW_model) %>% matrix(N, T)
    phi <- SVD(Y -Y_fit)$sigma[1] * (1 + delta)
    print(phi)
    
    nnr_fit = fit_logit(X, Y,  phi, 10, 300, 1e-8)
    beta_nnr <-  nnr_fit$beta_est
    
    sigma <-  SVD(nnr_fit$Theta_est)$sigma
    num_factor_est <- determine_num_factor(sigma, R_max)
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    
    fe_fit = MLE_logit(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe <- fe_fit$beta_est
    L_fe = fe_fit$L_est
    R_fe = fe_fit$R_est
    
    
    index <- Compute_index_with_LR(X, beta_fe, L_fe, R_fe)
    P <-  1 - 1/(exp(index) + 1) 
    phi <- SVD(Y -P)$sigma[1] * (1 + delta)
    
    
    # NNR 
    
    nnr_fit = fit_logit(X, Y,  mean(Y)*sqrt(N), 10, 1000, 1e-8)
    beta_nnr <-  nnr_fit$beta_est
    print(beta_pool)
    
    # determine the number of factor
    
    sigma <-  SVD(nnr_fit$Theta_est)$sigma
    num_factor_est <- determine_num_factor(sigma, R_max)
    print(num_factor_est)
    
    # beta_fe with correct number of factor
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    
    fe_fit = MLE_logit(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe <- fe_fit$beta_est
    L_fe = fe_fit$L_est
    R_fe = fe_fit$R_est
    
    # beta_fe analytical bias correction with correct number of factor
    
    
    bias_corr <- Compute_bias_corr_logit(Y, X, beta_fe, L_fe, R_fe, 0)
    beta_corr <- bias_corr$beta_corr
    std_beta_corr <- bias_corr$std_est
    
    # beta_fe sample splitting bias correction with correct number of factor
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_logit(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left <- fe_fit_left$beta_est
    fe_fit_right = MLE_logit(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right <- fe_fit_right$beta_est
    fe_fit_upper = MLE_logit(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_logit(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom <- fe_fit_bottom$beta_est
    
    beta_corr_sp <- 3* beta_fe - (beta_fe_left + beta_fe_right + beta_fe_upper + beta_fe_bottom)/2
    
    
    # beta_fe with number of factor (data_driven)
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor_est);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    if (num_factor_est==0){
        L_0 <- matrix(0, N, 1)
        R_0 <- matrix(0, T, 1)
    }
    
    fe_fit_data = MLE_logit(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe_data <- fe_fit_data$beta_est
    L_fe_data = fe_fit_data$L_est
    R_fe_data = fe_fit_data$R_est
    
    # beta_fe analytical bias correction with number of factor (data_driven)
    bias_corr_data <- Compute_bias_corr_logit(Y, X, beta_fe_data, L_fe_data, R_fe_data, 0)
    beta_corr_data <- bias_corr_data$beta_corr
    std_beta_corr_data <- bias_corr_data$std_est
    
    # beta_fe sample splitting bias correction with number of factor (data_driven)
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_logit(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left_data <- fe_fit_left$beta_est
    fe_fit_right = MLE_logit(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right_data <- fe_fit_right$beta_est
    fe_fit_upper = MLE_logit(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper_data <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_logit(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom_data <- fe_fit_bottom$beta_est
    
    beta_corr_sp_data <- 3* beta_fe_data - (beta_fe_left_data + beta_fe_right_data + beta_fe_upper_data + beta_fe_bottom_data)/2
    
    
    
    result <- list(
        beta_pool = beta_pool,
        beta_nnr = beta_nnr,
        num_factor_est = num_factor_est,
        beta_fe = beta_fe,
        beta_corr = beta_corr,
        beta_corr_sp = beta_corr_sp,
        std_beta_corr = std_beta_corr,
        beta_fe_data = beta_fe_data, 
        beta_corr_data = beta_corr_data, 
        beta_corr_sp_data = beta_corr_sp_data,
        std_beta_corr_data = std_beta_corr_data)
    
}
stopCluster(myCluster)
summary_ld_100_40 <- rearrange_result(raw_data, beta, num_X)



seed <- sample.int(100000, sim, replace = FALSE)

N = 200
T = 40

myCluster <- makeCluster(detectCores()-1, type = "FORK")
registerDoParallel(myCluster)
raw_data <- foreach (k=1:sim)%dopar%{
    set.seed(seed[k])
    
    data <- gen_data_dynamic_logit(N, T, num_factor, num_Z, num_W, beta_W, beta_Z)
    
    # Pool
    Yvec <- data$Y %>% c() 
    Xvec <- data$X[[1]] %>% c()
    Zvec <- data$X[[2]]%>% c()
    fvec <- data$interact_f %>% c()
    id <- rep(1:N, T) %>% matrix(N, T) %>% c()
    time <- rep(1:T, N) %>% matrix(T, N) %>% t() %>% c()
    df <- data.frame(Yvec, Xvec, Zvec, id, time)
    
    X <- data$X
    Y <- data$Y
    print(mean(Y))
    print(mean(X[[1]]))
    print(mean(X[[2]]))
    
    pool <- glm(Yvec ~ Xvec + Zvec, data = df, family = binomial(link = "logit"))
    summary(pool)
    beta_pool <- pool$coefficients[-1]
    beta_pool <- c(unlist(beta_pool))
    
    
    # Compute tuning parameter
    
    TW_model <- feglm(Yvec ~ Xvec + Zvec | id + time, family = binomial(link = "logit"), data = df)
    Y_fit <- fitted(TW_model) %>% matrix(N, T)
    phi <- SVD(Y -Y_fit)$sigma[1] * (1 + delta)
    print(phi)
    
    nnr_fit = fit_logit(X, Y,  phi, 10, 300, 1e-8)
    beta_nnr <-  nnr_fit$beta_est
    
    sigma <-  SVD(nnr_fit$Theta_est)$sigma
    num_factor_est <- determine_num_factor(sigma, R_max)
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    
    fe_fit = MLE_logit(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe <- fe_fit$beta_est
    L_fe = fe_fit$L_est
    R_fe = fe_fit$R_est
    
    
    index <- Compute_index_with_LR(X, beta_fe, L_fe, R_fe)
    P <-  1 - 1/(exp(index) + 1) 
    phi <- SVD(Y -P)$sigma[1] * (1 + delta)
    
    
    # NNR 
    
    nnr_fit = fit_logit(X, Y,  mean(Y)*sqrt(N), 10, 1000, 1e-8)
    beta_nnr <-  nnr_fit$beta_est
    print(beta_pool)
    
    # determine the number of factor
    
    sigma <-  SVD(nnr_fit$Theta_est)$sigma
    num_factor_est <- determine_num_factor(sigma, R_max)
    print(num_factor_est)
    
    # beta_fe with correct number of factor
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    
    fe_fit = MLE_logit(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe <- fe_fit$beta_est
    L_fe = fe_fit$L_est
    R_fe = fe_fit$R_est
    
    # beta_fe analytical bias correction with correct number of factor
    
    
    bias_corr <- Compute_bias_corr_logit(Y, X, beta_fe, L_fe, R_fe, 0)
    beta_corr <- bias_corr$beta_corr
    std_beta_corr <- bias_corr$std_est
    
    # beta_fe sample splitting bias correction with correct number of factor
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_logit(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left <- fe_fit_left$beta_est
    fe_fit_right = MLE_logit(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right <- fe_fit_right$beta_est
    fe_fit_upper = MLE_logit(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_logit(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom <- fe_fit_bottom$beta_est
    
    beta_corr_sp <- 3* beta_fe - (beta_fe_left + beta_fe_right + beta_fe_upper + beta_fe_bottom)/2
    
    
    # beta_fe with number of factor (data_driven)
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor_est);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    if (num_factor_est==0){
        L_0 <- matrix(0, N, 1)
        R_0 <- matrix(0, T, 1)
    }
    
    fe_fit_data = MLE_logit(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe_data <- fe_fit_data$beta_est
    L_fe_data = fe_fit_data$L_est
    R_fe_data = fe_fit_data$R_est
    
    # beta_fe analytical bias correction with number of factor (data_driven)
    bias_corr_data <- Compute_bias_corr_logit(Y, X, beta_fe_data, L_fe_data, R_fe_data, 0)
    beta_corr_data <- bias_corr_data$beta_corr
    std_beta_corr_data <- bias_corr_data$std_est
    
    # beta_fe sample splitting bias correction with number of factor (data_driven)
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_logit(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left_data <- fe_fit_left$beta_est
    fe_fit_right = MLE_logit(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right_data <- fe_fit_right$beta_est
    fe_fit_upper = MLE_logit(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper_data <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_logit(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom_data <- fe_fit_bottom$beta_est
    
    beta_corr_sp_data <- 3* beta_fe_data - (beta_fe_left_data + beta_fe_right_data + beta_fe_upper_data + beta_fe_bottom_data)/2
    
    
    
    result <- list(
        beta_pool = beta_pool,
        beta_nnr = beta_nnr,
        num_factor_est = num_factor_est,
        beta_fe = beta_fe,
        beta_corr = beta_corr,
        beta_corr_sp = beta_corr_sp,
        std_beta_corr = std_beta_corr,
        beta_fe_data = beta_fe_data, 
        beta_corr_data = beta_corr_data, 
        beta_corr_sp_data = beta_corr_sp_data,
        std_beta_corr_data = std_beta_corr_data)
    
}
stopCluster(myCluster)
summary_ld_200_40 <- rearrange_result(raw_data, beta, num_X)


seed <- sample.int(100000, sim, replace = FALSE)

N = 100
T = 100

myCluster <- makeCluster(detectCores()-1, type = "FORK")
registerDoParallel(myCluster)
raw_data <- foreach (k=1:sim)%dopar%{
    set.seed(seed[k])
    
    data <- gen_data_dynamic_logit(N, T, num_factor, num_Z, num_W, beta_W, beta_Z)
    
    # Pool
    Yvec <- data$Y %>% c() 
    Xvec <- data$X[[1]] %>% c()
    Zvec <- data$X[[2]]%>% c()
    fvec <- data$interact_f %>% c()
    id <- rep(1:N, T) %>% matrix(N, T) %>% c()
    time <- rep(1:T, N) %>% matrix(T, N) %>% t() %>% c()
    df <- data.frame(Yvec, Xvec, Zvec, id, time)
    
    X <- data$X
    Y <- data$Y
    print(mean(Y))
    print(mean(X[[1]]))
    print(mean(X[[2]]))
    
    pool <- glm(Yvec ~ Xvec + Zvec, data = df, family = binomial(link = "logit"))
    summary(pool)
    beta_pool <- pool$coefficients[-1]
    beta_pool <- c(unlist(beta_pool))
    
    
    # Compute tuning parameter
    
    TW_model <- feglm(Yvec ~ Xvec + Zvec | id + time, family = binomial(link = "logit"), data = df)
    Y_fit <- fitted(TW_model) %>% matrix(N, T)
    phi <- SVD(Y -Y_fit)$sigma[1] * (1 + delta)
    print(phi)
    
    nnr_fit = fit_logit(X, Y,  phi, 10, 300, 1e-8)
    beta_nnr <-  nnr_fit$beta_est
    
    sigma <-  SVD(nnr_fit$Theta_est)$sigma
    num_factor_est <- determine_num_factor(sigma, R_max)
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    
    fe_fit = MLE_logit(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe <- fe_fit$beta_est
    L_fe = fe_fit$L_est
    R_fe = fe_fit$R_est
    
    
    index <- Compute_index_with_LR(X, beta_fe, L_fe, R_fe)
    P <-  1 - 1/(exp(index) + 1) 
    phi <- SVD(Y -P)$sigma[1] * (1 + delta)
    
    
    # NNR 
    
    nnr_fit = fit_logit(X, Y,  mean(Y)*sqrt(N), 10, 1000, 1e-8)
    beta_nnr <-  nnr_fit$beta_est
    print(beta_pool)
    
    # determine the number of factor
    
    sigma <-  SVD(nnr_fit$Theta_est)$sigma
    num_factor_est <- determine_num_factor(sigma, R_max)
    print(num_factor_est)
    
    # beta_fe with correct number of factor
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    
    fe_fit = MLE_logit(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe <- fe_fit$beta_est
    L_fe = fe_fit$L_est
    R_fe = fe_fit$R_est
    
    # beta_fe analytical bias correction with correct number of factor
    
    
    bias_corr <- Compute_bias_corr_logit(Y, X, beta_fe, L_fe, R_fe, 0)
    beta_corr <- bias_corr$beta_corr
    std_beta_corr <- bias_corr$std_est
    
    # beta_fe sample splitting bias correction with correct number of factor
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_logit(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left <- fe_fit_left$beta_est
    fe_fit_right = MLE_logit(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right <- fe_fit_right$beta_est
    fe_fit_upper = MLE_logit(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_logit(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom <- fe_fit_bottom$beta_est
    
    beta_corr_sp <- 3* beta_fe - (beta_fe_left + beta_fe_right + beta_fe_upper + beta_fe_bottom)/2
    
    
    # beta_fe with number of factor (data_driven)
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor_est);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    if (num_factor_est==0){
        L_0 <- matrix(0, N, 1)
        R_0 <- matrix(0, T, 1)
    }
    
    fe_fit_data = MLE_logit(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe_data <- fe_fit_data$beta_est
    L_fe_data = fe_fit_data$L_est
    R_fe_data = fe_fit_data$R_est
    
    # beta_fe analytical bias correction with number of factor (data_driven)
    bias_corr_data <- Compute_bias_corr_logit(Y, X, beta_fe_data, L_fe_data, R_fe_data, 0)
    beta_corr_data <- bias_corr_data$beta_corr
    std_beta_corr_data <- bias_corr_data$std_est
    
    # beta_fe sample splitting bias correction with number of factor (data_driven)
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_logit(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left_data <- fe_fit_left$beta_est
    fe_fit_right = MLE_logit(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right_data <- fe_fit_right$beta_est
    fe_fit_upper = MLE_logit(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper_data <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_logit(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom_data <- fe_fit_bottom$beta_est
    
    beta_corr_sp_data <- 3* beta_fe_data - (beta_fe_left_data + beta_fe_right_data + beta_fe_upper_data + beta_fe_bottom_data)/2
    
    
    
    result <- list(
        beta_pool = beta_pool,
        beta_nnr = beta_nnr,
        num_factor_est = num_factor_est,
        beta_fe = beta_fe,
        beta_corr = beta_corr,
        beta_corr_sp = beta_corr_sp,
        std_beta_corr = std_beta_corr,
        beta_fe_data = beta_fe_data, 
        beta_corr_data = beta_corr_data, 
        beta_corr_sp_data = beta_corr_sp_data,
        std_beta_corr_data = std_beta_corr_data)
    
}
stopCluster(myCluster)
summary_ld_100_100 <- rearrange_result(raw_data, beta, num_X)



seed <- sample.int(100000, sim, replace = FALSE)

N = 200
T = 200

myCluster <- makeCluster(detectCores()-1, type = "FORK")
registerDoParallel(myCluster)
raw_data <- foreach (k=1:sim)%dopar%{
    set.seed(seed[k])
    
    data <- gen_data_dynamic_logit(N, T, num_factor, num_Z, num_W, beta_W, beta_Z)
    
    # Pool
    Yvec <- data$Y %>% c() 
    Xvec <- data$X[[1]] %>% c()
    Zvec <- data$X[[2]]%>% c()
    fvec <- data$interact_f %>% c()
    id <- rep(1:N, T) %>% matrix(N, T) %>% c()
    time <- rep(1:T, N) %>% matrix(T, N) %>% t() %>% c()
    df <- data.frame(Yvec, Xvec, Zvec, id, time)
    
    X <- data$X
    Y <- data$Y
    print(mean(Y))
    print(mean(X[[1]]))
    print(mean(X[[2]]))
    
    pool <- glm(Yvec ~ Xvec + Zvec, data = df, family = binomial(link = "logit"))
    summary(pool)
    beta_pool <- pool$coefficients[-1]
    beta_pool <- c(unlist(beta_pool))
    
    
    # Compute tuning parameter
    
    TW_model <- feglm(Yvec ~ Xvec + Zvec | id + time, family = binomial(link = "logit"), data = df)
    Y_fit <- fitted(TW_model) %>% matrix(N, T)
    phi <- SVD(Y -Y_fit)$sigma[1] * (1 + delta)
    print(phi)
    
    nnr_fit = fit_logit(X, Y,  phi, 10, 300, 1e-8)
    beta_nnr <-  nnr_fit$beta_est
    
    sigma <-  SVD(nnr_fit$Theta_est)$sigma
    num_factor_est <- determine_num_factor(sigma, R_max)
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    
    fe_fit = MLE_logit(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe <- fe_fit$beta_est
    L_fe = fe_fit$L_est
    R_fe = fe_fit$R_est
    
    
    index <- Compute_index_with_LR(X, beta_fe, L_fe, R_fe)
    P <-  1 - 1/(exp(index) + 1) 
    phi <- SVD(Y -P)$sigma[1] * (1 + delta)
    
    
    # NNR 
    
    nnr_fit = fit_logit(X, Y,  mean(Y)*sqrt(N), 10, 1000, 1e-8)
    beta_nnr <-  nnr_fit$beta_est
    print(beta_pool)
    
    # determine the number of factor
    
    sigma <-  SVD(nnr_fit$Theta_est)$sigma
    num_factor_est <- determine_num_factor(sigma, R_max)
    print(num_factor_est)
    
    # beta_fe with correct number of factor
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    
    fe_fit = MLE_logit(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe <- fe_fit$beta_est
    L_fe = fe_fit$L_est
    R_fe = fe_fit$R_est
    
    # beta_fe analytical bias correction with correct number of factor
    
    
    bias_corr <- Compute_bias_corr_logit(Y, X, beta_fe, L_fe, R_fe, 0)
    beta_corr <- bias_corr$beta_corr
    std_beta_corr <- bias_corr$std_est
    
    # beta_fe sample splitting bias correction with correct number of factor
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_logit(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left <- fe_fit_left$beta_est
    fe_fit_right = MLE_logit(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right <- fe_fit_right$beta_est
    fe_fit_upper = MLE_logit(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_logit(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom <- fe_fit_bottom$beta_est
    
    beta_corr_sp <- 3* beta_fe - (beta_fe_left + beta_fe_right + beta_fe_upper + beta_fe_bottom)/2
    
    
    # beta_fe with number of factor (data_driven)
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor_est);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    if (num_factor_est==0){
        L_0 <- matrix(0, N, 1)
        R_0 <- matrix(0, T, 1)
    }
    
    fe_fit_data = MLE_logit(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe_data <- fe_fit_data$beta_est
    L_fe_data = fe_fit_data$L_est
    R_fe_data = fe_fit_data$R_est
    
    # beta_fe analytical bias correction with number of factor (data_driven)
    bias_corr_data <- Compute_bias_corr_logit(Y, X, beta_fe_data, L_fe_data, R_fe_data, 0)
    beta_corr_data <- bias_corr_data$beta_corr
    std_beta_corr_data <- bias_corr_data$std_est
    
    # beta_fe sample splitting bias correction with number of factor (data_driven)
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_logit(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left_data <- fe_fit_left$beta_est
    fe_fit_right = MLE_logit(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right_data <- fe_fit_right$beta_est
    fe_fit_upper = MLE_logit(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper_data <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_logit(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom_data <- fe_fit_bottom$beta_est
    
    beta_corr_sp_data <- 3* beta_fe_data - (beta_fe_left_data + beta_fe_right_data + beta_fe_upper_data + beta_fe_bottom_data)/2
    
    
    
    result <- list(
        beta_pool = beta_pool,
        beta_nnr = beta_nnr,
        num_factor_est = num_factor_est,
        beta_fe = beta_fe,
        beta_corr = beta_corr,
        beta_corr_sp = beta_corr_sp,
        std_beta_corr = std_beta_corr,
        beta_fe_data = beta_fe_data, 
        beta_corr_data = beta_corr_data, 
        beta_corr_sp_data = beta_corr_sp_data,
        std_beta_corr_data = std_beta_corr_data)
    
}
stopCluster(myCluster)
summary_ld_200_200 <- rearrange_result(raw_data, beta, num_X)




#================================================


num_Z = 1
num_W = 1
num_X = num_Z + num_W
num_factor = 2
beta_W = c(0.25)
beta_Z = c(0)
beta = c(0.25, 0)
R_max = 5



seed <- sample.int(100000, sim, replace = FALSE)

N = 50
T = 40

myCluster <- makeCluster(detectCores()-1, type = "FORK")
registerDoParallel(myCluster)
raw_data <- foreach (k=1:sim)%dopar%{
    set.seed(seed[k])
    
    data <- gen_data_dynamic_poisson(N, T, num_factor, num_Z, num_W, beta_W, beta_Z)
    
    # Pool
    Yvec <- data$Y %>% c() 
    Xvec <- data$X[[1]] %>% c()
    Zvec <- data$X[[2]]%>% c()
    fvec <- data$interact_f %>% c()
    df <- data.frame(Yvec, Xvec, Xvec)
    
    X <- data$X
    Y <- data$Y
    
    pool <- glm(Yvec ~ Xvec + Zvec, data = df, family = poisson(link = "log"))
    summary(pool)
    beta_pool <- pool$coefficients[-1]
    beta_pool <- c(unlist(beta_pool))
    
    # NNR 
    
    nnr_fit = fit_poisson(X, Y,  mean(Y)*sqrt(N), 10, 300, 1e-8)
    beta_nnr <-  nnr_fit$beta_est
    print(beta_pool)
    
    # determine the number of factor
    
    sigma <-  SVD(nnr_fit$Theta_est)$sigma
    num_factor_est <- determine_num_factor(sigma, R_max)
    print(num_factor_est)
    
    # beta_fe with correct number of factor
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    
    fe_fit = MLE_poisson(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe <- fe_fit$beta_est
    L_fe = fe_fit$L_est
    R_fe = fe_fit$R_est
    
    # beta_fe analytical bias correction with correct number of factor
    
    
    bias_corr <- Compute_bias_corr_poisson(Y, X, beta_fe, L_fe, R_fe, 0)
    beta_corr <- bias_corr$beta_corr
    std_beta_corr <- bias_corr$std_est
    
    # beta_fe sample splitting bias correction with correct number of factor
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_poisson(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left <- fe_fit_left$beta_est
    fe_fit_right = MLE_poisson(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right <- fe_fit_right$beta_est
    fe_fit_upper = MLE_poisson(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_poisson(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom <- fe_fit_bottom$beta_est
    
    beta_corr_sp <- 3* beta_fe - (beta_fe_left + beta_fe_right + beta_fe_upper + beta_fe_bottom)/2
    
    
    # beta_fe with number of factor (data_driven)
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor_est);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    if (num_factor_est==0){
        L_0 <- matrix(0, N, 1)
        R_0 <- matrix(0, T, 1)
    }
    
    fe_fit_data = MLE_poisson(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe_data <- fe_fit_data$beta_est
    L_fe_data = fe_fit_data$L_est
    R_fe_data = fe_fit_data$R_est
    
    # beta_fe analytical bias correction with number of factor (data_driven)
    bias_corr_data <- Compute_bias_corr_poisson(Y, X, beta_fe_data, L_fe_data, R_fe_data, 0)
    beta_corr_data <- bias_corr_data$beta_corr
    std_beta_corr_data <- bias_corr_data$std_est
    
    # beta_fe sample splitting bias correction with number of factor (data_driven)
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_poisson(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left_data <- fe_fit_left$beta_est
    fe_fit_right = MLE_poisson(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right_data <- fe_fit_right$beta_est
    fe_fit_upper = MLE_poisson(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper_data <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_poisson(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom_data <- fe_fit_bottom$beta_est
    
    beta_corr_sp_data <- 3* beta_fe_data - (beta_fe_left_data + beta_fe_right_data + beta_fe_upper_data + beta_fe_bottom_data)/2
    
    
    
    result <- list(
        beta_pool = beta_pool,
        beta_nnr = beta_nnr,
        num_factor_est = num_factor_est,
        beta_fe = beta_fe,
        beta_corr = beta_corr,
        beta_corr_sp = beta_corr_sp,
        std_beta_corr = std_beta_corr,
        beta_fe_data = beta_fe_data, 
        beta_corr_data = beta_corr_data, 
        beta_corr_sp_data = beta_corr_sp_data,
        std_beta_corr_data = std_beta_corr_data)
}
summary_pd_50_40 <- rearrange_result(raw_data, beta, num_X)




seed <- sample.int(100000, sim, replace = FALSE)

N = 100
T = 40

myCluster <- makeCluster(detectCores()-1, type = "FORK")
registerDoParallel(myCluster)
raw_data <- foreach (k=1:sim)%dopar%{
    set.seed(seed[k])
    
    data <- gen_data_dynamic_poisson(N, T, num_factor, num_Z, num_W, beta_W, beta_Z)
    
    # Pool
    Yvec <- data$Y %>% c() 
    Xvec <- data$X[[1]] %>% c()
    Zvec <- data$X[[2]]%>% c()
    fvec <- data$interact_f %>% c()
    df <- data.frame(Yvec, Xvec, Xvec)
    
    X <- data$X
    Y <- data$Y
    
    pool <- glm(Yvec ~ Xvec + Zvec, data = df, family = poisson(link = "log"))
    summary(pool)
    beta_pool <- pool$coefficients[-1]
    beta_pool <- c(unlist(beta_pool))
    
    # NNR 
    
    nnr_fit = fit_poisson(X, Y,  mean(Y)*sqrt(N), 10, 300, 1e-8)
    beta_nnr <-  nnr_fit$beta_est
    print(beta_pool)
    
    # determine the number of factor
    
    sigma <-  SVD(nnr_fit$Theta_est)$sigma
    num_factor_est <- determine_num_factor(sigma, R_max)
    print(num_factor_est)
    
    # beta_fe with correct number of factor
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    
    fe_fit = MLE_poisson(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe <- fe_fit$beta_est
    L_fe = fe_fit$L_est
    R_fe = fe_fit$R_est
    
    # beta_fe analytical bias correction with correct number of factor
    
    
    bias_corr <- Compute_bias_corr_poisson(Y, X, beta_fe, L_fe, R_fe, 0)
    beta_corr <- bias_corr$beta_corr
    std_beta_corr <- bias_corr$std_est
    
    # beta_fe sample splitting bias correction with correct number of factor
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_poisson(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left <- fe_fit_left$beta_est
    fe_fit_right = MLE_poisson(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right <- fe_fit_right$beta_est
    fe_fit_upper = MLE_poisson(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_poisson(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom <- fe_fit_bottom$beta_est
    
    beta_corr_sp <- 3* beta_fe - (beta_fe_left + beta_fe_right + beta_fe_upper + beta_fe_bottom)/2
    
    
    # beta_fe with number of factor (data_driven)
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor_est);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    if (num_factor_est==0){
        L_0 <- matrix(0, N, 1)
        R_0 <- matrix(0, T, 1)
    }
    
    fe_fit_data = MLE_poisson(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe_data <- fe_fit_data$beta_est
    L_fe_data = fe_fit_data$L_est
    R_fe_data = fe_fit_data$R_est
    
    # beta_fe analytical bias correction with number of factor (data_driven)
    bias_corr_data <- Compute_bias_corr_poisson(Y, X, beta_fe_data, L_fe_data, R_fe_data, 0)
    beta_corr_data <- bias_corr_data$beta_corr
    std_beta_corr_data <- bias_corr_data$std_est
    
    # beta_fe sample splitting bias correction with number of factor (data_driven)
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_poisson(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left_data <- fe_fit_left$beta_est
    fe_fit_right = MLE_poisson(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right_data <- fe_fit_right$beta_est
    fe_fit_upper = MLE_poisson(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper_data <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_poisson(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom_data <- fe_fit_bottom$beta_est
    
    beta_corr_sp_data <- 3* beta_fe_data - (beta_fe_left_data + beta_fe_right_data + beta_fe_upper_data + beta_fe_bottom_data)/2
    
    
    
    result <- list(
        beta_pool = beta_pool,
        beta_nnr = beta_nnr,
        num_factor_est = num_factor_est,
        beta_fe = beta_fe,
        beta_corr = beta_corr,
        beta_corr_sp = beta_corr_sp,
        std_beta_corr = std_beta_corr,
        beta_fe_data = beta_fe_data, 
        beta_corr_data = beta_corr_data, 
        beta_corr_sp_data = beta_corr_sp_data,
        std_beta_corr_data = std_beta_corr_data)
}
summary_pd_100_40 <- rearrange_result(raw_data, beta, num_X)



seed <- sample.int(100000, sim, replace = FALSE)

N = 200
T = 40

myCluster <- makeCluster(detectCores()-1, type = "FORK")
registerDoParallel(myCluster)
raw_data <- foreach (k=1:sim)%dopar%{
    set.seed(seed[k])
    
    data <- gen_data_dynamic_poisson(N, T, num_factor, num_Z, num_W, beta_W, beta_Z)
    
    # Pool
    Yvec <- data$Y %>% c() 
    Xvec <- data$X[[1]] %>% c()
    Zvec <- data$X[[2]]%>% c()
    fvec <- data$interact_f %>% c()
    df <- data.frame(Yvec, Xvec, Xvec)
    
    X <- data$X
    Y <- data$Y
    
    pool <- glm(Yvec ~ Xvec + Zvec, data = df, family = poisson(link = "log"))
    summary(pool)
    beta_pool <- pool$coefficients[-1]
    beta_pool <- c(unlist(beta_pool))
    
    # NNR 
    
    nnr_fit = fit_poisson(X, Y,  mean(Y)*sqrt(N), 10, 300, 1e-8)
    beta_nnr <-  nnr_fit$beta_est
    print(beta_pool)
    
    # determine the number of factor
    
    sigma <-  SVD(nnr_fit$Theta_est)$sigma
    num_factor_est <- determine_num_factor(sigma, R_max)
    print(num_factor_est)
    
    # beta_fe with correct number of factor
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    
    fe_fit = MLE_poisson(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe <- fe_fit$beta_est
    L_fe = fe_fit$L_est
    R_fe = fe_fit$R_est
    
    # beta_fe analytical bias correction with correct number of factor
    
    
    bias_corr <- Compute_bias_corr_poisson(Y, X, beta_fe, L_fe, R_fe, 0)
    beta_corr <- bias_corr$beta_corr
    std_beta_corr <- bias_corr$std_est
    
    # beta_fe sample splitting bias correction with correct number of factor
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_poisson(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left <- fe_fit_left$beta_est
    fe_fit_right = MLE_poisson(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right <- fe_fit_right$beta_est
    fe_fit_upper = MLE_poisson(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_poisson(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom <- fe_fit_bottom$beta_est
    
    beta_corr_sp <- 3* beta_fe - (beta_fe_left + beta_fe_right + beta_fe_upper + beta_fe_bottom)/2
    
    
    # beta_fe with number of factor (data_driven)
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor_est);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    if (num_factor_est==0){
        L_0 <- matrix(0, N, 1)
        R_0 <- matrix(0, T, 1)
    }
    
    fe_fit_data = MLE_poisson(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe_data <- fe_fit_data$beta_est
    L_fe_data = fe_fit_data$L_est
    R_fe_data = fe_fit_data$R_est
    
    # beta_fe analytical bias correction with number of factor (data_driven)
    bias_corr_data <- Compute_bias_corr_poisson(Y, X, beta_fe_data, L_fe_data, R_fe_data, 0)
    beta_corr_data <- bias_corr_data$beta_corr
    std_beta_corr_data <- bias_corr_data$std_est
    
    # beta_fe sample splitting bias correction with number of factor (data_driven)
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_poisson(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left_data <- fe_fit_left$beta_est
    fe_fit_right = MLE_poisson(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right_data <- fe_fit_right$beta_est
    fe_fit_upper = MLE_poisson(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper_data <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_poisson(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom_data <- fe_fit_bottom$beta_est
    
    beta_corr_sp_data <- 3* beta_fe_data - (beta_fe_left_data + beta_fe_right_data + beta_fe_upper_data + beta_fe_bottom_data)/2
    
    
    
    result <- list(
        beta_pool = beta_pool,
        beta_nnr = beta_nnr,
        num_factor_est = num_factor_est,
        beta_fe = beta_fe,
        beta_corr = beta_corr,
        beta_corr_sp = beta_corr_sp,
        std_beta_corr = std_beta_corr,
        beta_fe_data = beta_fe_data, 
        beta_corr_data = beta_corr_data, 
        beta_corr_sp_data = beta_corr_sp_data,
        std_beta_corr_data = std_beta_corr_data)
}
summary_pd_200_40 <- rearrange_result(raw_data, beta, num_X)




seed <- sample.int(100000, sim, replace = FALSE)

N = 100
T = 100

myCluster <- makeCluster(detectCores()-1, type = "FORK")
registerDoParallel(myCluster)
raw_data <- foreach (k=1:sim)%dopar%{
    set.seed(seed[k])
    
    data <- gen_data_dynamic_poisson(N, T, num_factor, num_Z, num_W, beta_W, beta_Z)
    
    # Pool
    Yvec <- data$Y %>% c() 
    Xvec <- data$X[[1]] %>% c()
    Zvec <- data$X[[2]]%>% c()
    fvec <- data$interact_f %>% c()
    df <- data.frame(Yvec, Xvec, Xvec)
    
    X <- data$X
    Y <- data$Y
    
    pool <- glm(Yvec ~ Xvec + Zvec, data = df, family = poisson(link = "log"))
    summary(pool)
    beta_pool <- pool$coefficients[-1]
    beta_pool <- c(unlist(beta_pool))
    
    # NNR 
    
    nnr_fit = fit_poisson(X, Y,  mean(Y)*sqrt(N), 10, 300, 1e-8)
    beta_nnr <-  nnr_fit$beta_est
    print(beta_pool)
    
    # determine the number of factor
    
    sigma <-  SVD(nnr_fit$Theta_est)$sigma
    num_factor_est <- determine_num_factor(sigma, R_max)
    print(num_factor_est)
    
    # beta_fe with correct number of factor
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    
    fe_fit = MLE_poisson(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe <- fe_fit$beta_est
    L_fe = fe_fit$L_est
    R_fe = fe_fit$R_est
    
    # beta_fe analytical bias correction with correct number of factor
    
    
    bias_corr <- Compute_bias_corr_poisson(Y, X, beta_fe, L_fe, R_fe, 0)
    beta_corr <- bias_corr$beta_corr
    std_beta_corr <- bias_corr$std_est
    
    # beta_fe sample splitting bias correction with correct number of factor
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_poisson(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left <- fe_fit_left$beta_est
    fe_fit_right = MLE_poisson(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right <- fe_fit_right$beta_est
    fe_fit_upper = MLE_poisson(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_poisson(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom <- fe_fit_bottom$beta_est
    
    beta_corr_sp <- 3* beta_fe - (beta_fe_left + beta_fe_right + beta_fe_upper + beta_fe_bottom)/2
    
    
    # beta_fe with number of factor (data_driven)
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor_est);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    if (num_factor_est==0){
        L_0 <- matrix(0, N, 1)
        R_0 <- matrix(0, T, 1)
    }
    
    fe_fit_data = MLE_poisson(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe_data <- fe_fit_data$beta_est
    L_fe_data = fe_fit_data$L_est
    R_fe_data = fe_fit_data$R_est
    
    # beta_fe analytical bias correction with number of factor (data_driven)
    bias_corr_data <- Compute_bias_corr_poisson(Y, X, beta_fe_data, L_fe_data, R_fe_data, 0)
    beta_corr_data <- bias_corr_data$beta_corr
    std_beta_corr_data <- bias_corr_data$std_est
    
    # beta_fe sample splitting bias correction with number of factor (data_driven)
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_poisson(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left_data <- fe_fit_left$beta_est
    fe_fit_right = MLE_poisson(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right_data <- fe_fit_right$beta_est
    fe_fit_upper = MLE_poisson(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper_data <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_poisson(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom_data <- fe_fit_bottom$beta_est
    
    beta_corr_sp_data <- 3* beta_fe_data - (beta_fe_left_data + beta_fe_right_data + beta_fe_upper_data + beta_fe_bottom_data)/2
    
    
    
    result <- list(
        beta_pool = beta_pool,
        beta_nnr = beta_nnr,
        num_factor_est = num_factor_est,
        beta_fe = beta_fe,
        beta_corr = beta_corr,
        beta_corr_sp = beta_corr_sp,
        std_beta_corr = std_beta_corr,
        beta_fe_data = beta_fe_data, 
        beta_corr_data = beta_corr_data, 
        beta_corr_sp_data = beta_corr_sp_data,
        std_beta_corr_data = std_beta_corr_data)
}
summary_pd_100_100 <- rearrange_result(raw_data, beta, num_X)


seed <- sample.int(100000, sim, replace = FALSE)

N = 200
T = 200

myCluster <- makeCluster(detectCores()-1, type = "FORK")
registerDoParallel(myCluster)
raw_data <- foreach (k=1:sim)%dopar%{
    set.seed(seed[k])
    
    data <- gen_data_dynamic_poisson(N, T, num_factor, num_Z, num_W, beta_W, beta_Z)
    
    # Pool
    Yvec <- data$Y %>% c() 
    Xvec <- data$X[[1]] %>% c()
    Zvec <- data$X[[2]]%>% c()
    fvec <- data$interact_f %>% c()
    df <- data.frame(Yvec, Xvec, Xvec)
    
    X <- data$X
    Y <- data$Y
    
    pool <- glm(Yvec ~ Xvec + Zvec, data = df, family = poisson(link = "log"))
    summary(pool)
    beta_pool <- pool$coefficients[-1]
    beta_pool <- c(unlist(beta_pool))
    
    # NNR 
    
    nnr_fit = fit_poisson(X, Y,  mean(Y)*sqrt(N), 10, 300, 1e-8)
    beta_nnr <-  nnr_fit$beta_est
    print(beta_pool)
    
    # determine the number of factor
    
    sigma <-  SVD(nnr_fit$Theta_est)$sigma
    num_factor_est <- determine_num_factor(sigma, R_max)
    print(num_factor_est)
    
    # beta_fe with correct number of factor
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    
    fe_fit = MLE_poisson(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe <- fe_fit$beta_est
    L_fe = fe_fit$L_est
    R_fe = fe_fit$R_est
    
    # beta_fe analytical bias correction with correct number of factor
    
    
    bias_corr <- Compute_bias_corr_poisson(Y, X, beta_fe, L_fe, R_fe, 0)
    beta_corr <- bias_corr$beta_corr
    std_beta_corr <- bias_corr$std_est
    
    # beta_fe sample splitting bias correction with correct number of factor
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_poisson(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left <- fe_fit_left$beta_est
    fe_fit_right = MLE_poisson(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right <- fe_fit_right$beta_est
    fe_fit_upper = MLE_poisson(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_poisson(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom <- fe_fit_bottom$beta_est
    
    beta_corr_sp <- 3* beta_fe - (beta_fe_left + beta_fe_right + beta_fe_upper + beta_fe_bottom)/2
    
    
    # beta_fe with number of factor (data_driven)
    
    list_LR = Low_rank_appro(nnr_fit$Theta_est, num_factor_est);
    L_0 = list_LR$L
    R_0 = list_LR$R
    beta_0 = beta_nnr
    if (num_factor_est==0){
        L_0 <- matrix(0, N, 1)
        R_0 <- matrix(0, T, 1)
    }
    
    fe_fit_data = MLE_poisson(X, Y, beta_0, L_0, R_0, 10, 10000, 1e-8)
    beta_fe_data <- fe_fit_data$beta_est
    L_fe_data = fe_fit_data$L_est
    R_fe_data = fe_fit_data$R_est
    
    # beta_fe analytical bias correction with number of factor (data_driven)
    bias_corr_data <- Compute_bias_corr_poisson(Y, X, beta_fe_data, L_fe_data, R_fe_data, 0)
    beta_corr_data <- bias_corr_data$beta_corr
    std_beta_corr_data <- bias_corr_data$std_est
    
    # beta_fe sample splitting bias correction with number of factor (data_driven)
    
    data_sample_split <-  sample_split(data, L_0, R_0)
    data_left <- data_sample_split$data_left
    data_right <- data_sample_split$data_right
    data_upper <- data_sample_split$data_upper
    data_bottom <- data_sample_split$data_bottom
    
    fe_fit_left = MLE_poisson(data_left$X, data_left$Y, beta_0, data_left$L, data_left$R, 10, 10000, 1e-8)
    beta_fe_left_data <- fe_fit_left$beta_est
    fe_fit_right = MLE_poisson(data_right$X, data_right$Y, beta_0, data_right$L, data_right$R, 10, 10000, 1e-8)
    beta_fe_right_data <- fe_fit_right$beta_est
    fe_fit_upper = MLE_poisson(data_upper$X, data_upper$Y, beta_0, data_upper$L, data_upper$R, 10, 10000, 1e-8)
    beta_fe_upper_data <- fe_fit_upper$beta_est
    fe_fit_bottom = MLE_poisson(data_bottom$X, data_bottom$Y, beta_0, data_bottom$L, data_bottom$R, 10, 10000, 1e-8)
    beta_fe_bottom_data <- fe_fit_bottom$beta_est
    
    beta_corr_sp_data <- 3* beta_fe_data - (beta_fe_left_data + beta_fe_right_data + beta_fe_upper_data + beta_fe_bottom_data)/2
    
    
    
    result <- list(
        beta_pool = beta_pool,
        beta_nnr = beta_nnr,
        num_factor_est = num_factor_est,
        beta_fe = beta_fe,
        beta_corr = beta_corr,
        beta_corr_sp = beta_corr_sp,
        std_beta_corr = std_beta_corr,
        beta_fe_data = beta_fe_data, 
        beta_corr_data = beta_corr_data, 
        beta_corr_sp_data = beta_corr_sp_data,
        std_beta_corr_data = std_beta_corr_data)
}
summary_pd_200_200 <- rearrange_result(raw_data, beta, num_X)









