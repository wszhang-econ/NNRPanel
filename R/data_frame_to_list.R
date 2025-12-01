data_frame_to_list <-function(data_frame){
    
    # delete missing rows 
    data_frame <- data_frame[complete.cases(data_frame), , drop = FALSE]
    
    # sort by id (the first column) and year (the second column)
    data_frame <-data_frame[order(data_frame[, 1], data_frame[, 2]), , drop = FALSE]
    
    # check unbalanced panel
    N <- length(unique(data_frame[, 1]))
    T <- length(unique(data_frame[, 2]))
    
    if (nrow(data_frame) != N * T) {
        stop("Unbalanced panel")
    }
    
    # transform the data_frame to list 
    
    num_X <- ncol(data_frame) -3 
        
    Y <- data_frame[, 3, drop = TRUE] 
    Y <- matrix(Y, T, N)
    Y <- t(Y)
    
    X <- lapply(4:ncol(data_frame), function(j) data_frame[, j, drop = TRUE])
    for (i in num_X){
        X[[i]] <- matrix(X[[i]], T, N)
        X[[i]] <- t(X[[i]])
    }
    
    data_list <- list(Y = Y, 
                 X = X)
    
    return(data_list) 
}




list_to_data_frame <-function(data_list){
    
    N <- nrow(data_list$Y)
    T <- ncol(data_list$Y) 
    
    num_X <- length(data_list$X) 
    
    
    data_frame <- matrix(0, N * T, 3 + num_X) 
    data_frame[, 1] <- rep(1:N, each = T)
    data_frame[, 2] <- rep(1:T, times = N)
    data_frame[, 3] <- t(data_list$Y)
    for(i in num_X){
        data_frame[, i + 3] <- t(data_list$X[[i]])
    }
    
    
    
    return(data_frame)  
}

