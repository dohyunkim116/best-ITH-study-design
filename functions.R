suppressMessages(library(doParallel))
suppressMessages(library(foreach))
Rcpp::sourceCpp("recursive_search.cpp")

simulate_tumors <- function(nTumors,nSamp=NA,nSamp_vec=NA,nSeg,mutationRates=NA,theta_lb=NA,theta_ub=NA,thetas=NA,
                               sim=NA){
    if(!is.na(sim)) set.seed(sim)
    if(is.na(nSamp)) {
        if(is.na(nSamp_vec[1]))
            nSamp_vec <- sample(2:10,nTumors,replace = T)
    }
    else {
        nSamp_vec <- rep(nSamp,nTumors)
    }
    if(is.na(mutationRates[1])){
        P <- runif(nSeg)
    }
    else {
        P <- mutationRates
    }
    
    Deltas <- lapply(1:nTumors,function(i){
        rbernoulli(length(P),P) %>% as.numeric()
    })
    
    if(!is.na(theta_lb) && !is.na(theta_ub))
        thetas <- runif(nTumors,theta_lb,theta_ub)
    
    tumor_mat_list <- lapply(1:nTumors,function(i){
        delta <- Deltas[[i]]
        theta <- thetas[i]
        nSamp <- nSamp_vec[i]
        temp <- lapply(delta,function(s){
            if(s == 1){
                rbernoulli(nSamp,theta)
            }
            else {
                rep(0,nSamp)
            }
        }) %>% lapply(as.numeric) %>% do.call(rbind,.)
        colnames(temp) <- paste0("s",1:nSamp)
        temp
    })
    attr(tumor_mat_list,"mutation_rate_vec") <- P
    attr(tumor_mat_list,"nSamp_vec") <- nSamp_vec
    tumor_mat_list
}

funWeightV3_mc <- function(tumor_mat_list,num_cores){
    distMatList <- lapply(tumor_mat_list,function(tumor_mat){
        dist(t(tumor_mat),method = "manhattan",diag = T,upper = T) %>% as.matrix()
    })
    
    tempFun <- function(distMat){
        k <- nrow(distMat)
        if(k <= 1)stop("Something is wrong!")
        distVec <- distMat[lower.tri(distMat)]
        k2 <- cbind(distVec, distVec)
        k3 <- NULL
        if(k >= 3){
            tempCombn <- combn(k, 3)
            k3 <- mclapply(1:ncol(tempCombn),function(i){
                tempIdx <- tempCombn[, i]
                tempMat <- cbind(c(distMat[tempIdx[1], tempIdx[2]], 
                                   distMat[tempIdx[1], tempIdx[2]], 
                                   distMat[tempIdx[1], tempIdx[3]]), 
                                 c(distMat[tempIdx[1], tempIdx[3]], 
                                   distMat[tempIdx[2], tempIdx[3]], 
                                   distMat[tempIdx[2], tempIdx[3]]))
                tempMat
            },mc.cores = num_cores) %>% do.call(rbind,.)
        }
        k4 <- NULL
        if(k >= 4){
            tempCombn <- combn(k, 4)
            k4 <- mclapply(1:ncol(tempCombn),function(i){
                tempIdx <- tempCombn[, i]
                tempMat <- cbind(c(distMat[tempIdx[1], tempIdx[2]], distMat[tempIdx[1], tempIdx[3]], distMat[tempIdx[1], tempIdx[4]]), c(distMat[tempIdx[3], tempIdx[4]], distMat[tempIdx[2], tempIdx[4]], distMat[tempIdx[2], tempIdx[3]]))
                tempMat
            },mc.cores = num_cores) %>% do.call(rbind,.)
        }
        return(list(k2 = k2, k3 = k3, k4 = k4))
    }
    tempList <- lapply(distMatList, tempFun)
    #	str(tempList)
    
    k2 <- k3 <- k4 <- NULL
    for(i in 1:length(tempList)){
        k2 <- rbind(k2, tempList[[i]]$k2)
        k3 <- rbind(k3, tempList[[i]]$k3)
        k4 <- rbind(k4, tempList[[i]]$k4)
    }
    #	str(k2)
    #	str(k3)
    #	str(k4)
    
    tempFun2 <- function(distMatList,i,j){
        tempDistMat1 <- distMatList[[i]]
        tempDistMat2 <- distMatList[[j]]
        tempDistVec1 <- tempDistMat1[lower.tri(tempDistMat1)]
        tempDistVec2 <- tempDistMat2[lower.tri(tempDistMat2)]
        tempMat <- cbind(rep(tempDistVec1, each = length(tempDistVec2)), 
                         rep(tempDistVec2, length(tempDistVec1)))
        tempMat        
    }
    registerDoParallel(cores = num_cores)
    kk <- foreach(i=1:(length(distMatList) - 1),.combine = rbind) %:%
        foreach(j=(i+1):length(distMatList),.combine = rbind) %dopar% {
            tempFun2(distMatList,i,j)
        }
    #	str(kk)
    
    a1 <- mean((k4[, 1] - k4[, 2])^2)
    a2 <- mean((k3[, 1] - k3[, 2])^2)
    a3 <- mean((kk[, 1] - kk[, 2])^2)
    b1 <- a1 / 2
    b2 <- (a1 - a2) / 2
    
    m <- unlist(lapply(distMatList, nrow))
    sigmaVec <- 2 / m / (m - 1) * b1 + 4 * (m - 2) / m / (m - 1) * b2
    sigmaVec[sigmaVec < 0] <- 0 #remove negative values
    xVec <- unlist(lapply(distMatList, function(x)mean(x[lower.tri(x)])))
    funLogL <- function(tau2, xVec, sigmaVec){
        wVec <- 1 / (sigmaVec^2 + tau2)
        mu <- sum(wVec * xVec) / sum(wVec)
        return(-0.5 * sum(dnorm(xVec, mu, sqrt(sigmaVec^2 + tau2), log = TRUE)))
    }
    tempOptim <- optim(0, funLogL, xVec = xVec, sigmaVec = sigmaVec, method = "L-BFGS-B", lower = 1e-7, upper = Inf)
    b3 <- tempOptim$par
    
    w <- 1 / (b3 + sigmaVec)
    attr(w, "par") <- c("sigma_square"=b1, "rho"=b2, "tau_square"=b3)
    return(w)
}

estimate_parameters <- function(tumor_mat_list){
    weight <- funWeightV3_mc(tumor_mat_list,1)
    attr(weight,"par")
}

phi1 <- function(tau_sq, sigma_sq, rho, nSamp, allowed_total_nSamp){
    between_var <- tau_sq
    avg_within_var <- sigma_sq
    avg_within_covar <- rho
    var_est_ith <- 2*avg_within_var/(nSamp^2-nSamp) + (4*nSamp-8)*avg_within_covar/(nSamp^2-nSamp)
    nSubjects <- allowed_total_nSamp/nSamp
    res <- nSubjects/(between_var + var_est_ith)
    attr(res,"sigma_n_sq") <- var_est_ith
    res
}

power1 <- function(sigma_sq,rho,tau_sq,beta,sigma_y,nSamp,allowed_total_nSamp){
    between_var <- tau_sq
    avg_within_var <- sigma_sq
    avg_within_covar <- rho
    var_est_ith <- 2*avg_within_var/(nSamp^2-nSamp) + (4*nSamp-8)*avg_within_covar/(nSamp^2-nSamp)
    nSubjects <- allowed_total_nSamp/nSamp
    phi1 <- nSubjects/(between_var + var_est_ith)
    xi <- beta/sigma_y*tau_sq*sqrt(phi1)
    res <- 2 - pnorm(1.96 - xi) - pnorm(1.96 + xi)
    attr(res,"sigma_n_sq") <- var_est_ith
    attr(res,"phi1") <- phi1
    attr(res,"xi") <- xi
    res
}

























