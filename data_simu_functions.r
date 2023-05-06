#########################################################################################################
######## Project 2 Hierarchical SAE with CRC estimator for areal population size estimation #############
#########################################################################################################


# library("copula")
if (!require("mvtnorm")) install.packages("mvtnorm") else (require("mvtnorm", quietly = TRUE)) 
if (!require("matrixsampling")) install.packages("matrixsampling") else (require("matrixsampling", quietly = TRUE)) 
if (!require("LaplacesDemon")) install.packages("LaplacesDemon") else (require("LaplacesDemon", quietly = TRUE)) 
if (!require("spam")) install.packages("spam") else (require("spam", quietly = TRUE)) 


##############################
### Create location matrix ###
##############################
### Create and returns an adjacent matrix of 2D grid graph
prepGridAdj <- function(n){
  ret <- diag(rep(1,n-1))
  ret <- rbind(rep(0,n-1),ret)
  ret <- cbind(ret,rep(0,n))
  return(ret + t(ret))
}
genGridAdj <- function(m1,m2){
  n1 <- length(m1[,1])
  n2 <- length(m2[,1])
  n <- n1*n2
  N <- as.matrix(expand.grid(1:n1,1:n2))
  M <- matrix(FALSE,n,n)
  for(i in 1:n){
    for(j in 1:n){
      if((N[i,1] == N[j,1]) & m2[N[i,2],N[j,2]]){
        M[i,j] <- TRUE
      }else if((N[i,2] == N[j,2]) & m1[N[i,1],N[j,1]]){
        M[i,j] <- TRUE
      }
    }
  }
  return(M)
} # return an adjacency matrix with FALSE = not boundary connected, TRUE = boundary connected
getAdjMatrix <- function(size1,size2=size1){
  # size1, size2 are integers indicating grid size, no need to be size1=size2, but size1=size2 makes a squared grid
  # e.g. size1 = 2, size2 = 8/size1 to make a 2*4 dimension grid
  N <- prepGridAdj(size1)
  M <- prepGridAdj(size2)
  AdjMatrix <- genGridAdj(N,M)
  W <- ifelse(AdjMatrix == TRUE, 1, 0)
  return(W)
} # return an adjacency matrix with 0 = not boundary connected, 1 = boundary connected

### Determine the log-normal distribution parameter by two quantiles
getParams_LogNorm <- function(x1, x2, p1, p2) {
  # x1 = first location, corresponding to the percentile p1 (p1 is a cdf)
  # x2 = second location, corresponding to the percentile p2 (p2 is a cdf)
  logx1 <- log(x1)
  logx2 <- log(x2)
  sigma <- (logx2 - logx1)/(qnorm(p2) - qnorm(p1))
  mu <- (logx1*qnorm(p2) - logx2*qnorm(p1))/(qnorm(p2) - qnorm(p1))
  res <- c(mu, sigma)
  return(res)
}

### Compute inverse of spatial covariance matrix (CAR)
getSpatialCov <- function(W, tau2 = 0.1, rho = 0.95){
  # W = neighbor matrix, main diag = n of neighbors
  # tau2 = the variance of spatial random variables
  # rho = control overall level of spatial autocorrelation
  
  ### Sources for computing Q.W
  # Paper Mitzi Morris(2019) T3.2S3N31
  # Rstan: B <- solve(D_w) %*% W; Q.W <- D_w %*% (I - rho*B) ## Here B is scaled adj matrix, the scale is inverse matrix of the number of neighbors
  # CAR.simGLM: https://rdrr.io/cran/mclcar/man/CAR.simLM.html
  
  ### Compute Precision Matrix
  D_w <- diag(apply(W, 2, sum)) # Adjacency matrix with number of neighbors
  I <- diag(nrow(W)) # identify matrix
  ### If just regular CAR, get precision matrix
  rho <- rho 
  prec_scaler <- 1/tau2
  Q.W <- D_w %*% (I - rho * solve(D_w) %*% W) 
  Q.W_spam <- as.spam(prec_scaler * Q.W) # act like as.matrix, returns a vector of non-zero values (column-wise listed) of the result from previous step of Q.W
  ### Create the spatial covariance matrix and simulate a set of spatial random effects
  Sigma <- solve(Q.W) # For plotting the summary stats
  res <- list(Sigma = Sigma,
              Q.W = Q.W,
              Q.W_spam = Q.W_spam)
  return(res)
} # cov.inv_mat and spam version of cov.inv_mat

##############################################################
### Model 1: Sample prevalence probabilities for locations ###
##############################################################
getPrev <- function(cov.inv_mat, beta0, eps_sd, sample_size = 1){ 
  # cov.inv_mat = inverse spatial covariance matrix (Q.W)
  # beta0 = overall baseline prevalence
  # eps_sd = unmeasured dispersion (unstructed iid random effect)on prevalence
  # sample_size = number of simulation (number of datasets)
  
  k <- nrow(cov.inv_mat)
  n <- nrow(cov.inv_mat)
  cholQ.W <- chol.spam(cov.inv_mat, pivot = "MMD", memory = list(nnzcolindices = 6.25 * n))
  # we can therefore calculate the Cholesky decomposition, making the simulation remarkably fast and efficient.
  phi <- backsolve(cholQ.W, rnorm(k)) # returns random draws for spatial random effects
  # alternatively,   
  # cov_mat <- solve(cov.inv_mat)
  # phi <- rmvnorm(n = sample_size, mean = rep(0, k), sigma = (cov_mat)) # phi is spatial random variable
  
  if(eps_sd > 0){
    eps_sd2 <- eps_sd^2
    eps <- rmvnorm(n = sample_size, mean = rep(0, k), sigma = eps_sd2 * diag(k)) # eps is error term that explain other heterogeneity
    # alternatively,
    # eps <- rnorm(n = k*sample_size, mean = 0, sd = eps_sd)
    # eps <- matrix(eps, nrow = sample_size, ncol = k)
  }
  else{eps <- rep(0,k)}
  
  logit_p <- beta0 + phi + eps # prevalence is a logit function of linear combination of overall prev, spatial random effect, and IID RE term
  p <- exp(logit_p) / (1 + exp(logit_p)) # equivalent to plogis(logit_p) or invlogit(logit_p) under LaplacesDemon package
  
  res <- data.frame(Prev = p, SpatialRE = phi)
  return(res)
} # return regional prevalence, length of K

### Get prevalence from copula
# getPrevCopula <- function(cov_mat, eps_sd = 0.1, tau2 = 0.01, sample_size = 1){
# cov_mat = spatial covariance matrix
# mean_phi = prev centered at mean_phi level
# eps_sd = unmeasured dispersion on prevalence
# tau2 = variance for spatial structure
# sample_size = number of simulation (number of datasets)

#   k <- nrow(cov_mat)
#   I <- diag(k)
#   eps_mat <- eps_sd^2*I
#   corr_mat <- cov2cor(tau2*cov_mat + eps_mat)
#   param_vec <- P2p(corr_mat) 
#   cop <- normalCopula(param = param_vec, dim = k, dispstr = "un")
#   p <- rCopula(sample_size, cop) # pnorm(phi)
#   # alternatively,
#   # phi <- rmvnorm(n = sample_size, mean = rep(0, k), sigma = (tau2 * cov_mat))
#   # p <- pnorm(phi, sd = sqrt(tau2))
#   return(p)
# }

###########################################
###/// Variation of Prevalence Model ///###
###########################################
### Ignore spatial term and error term, only apply universal mean prevalence, no location difference 
getPrev_1 <- function(eps0, eps0_sd, sample_size = 1){
  # eps0 = mean prevalence
  # eps0_sd = sd of prevalence
  logit_p <- rnorm(n = sample_size, mean = eps0, sd = eps0_sd)
  p <- exp(logit_p) / (1 + exp(logit_p)) 
  return(p)
}
# ### Ignore spatial term, only apply random effect of location with a common mean and variance
# getPrev_1.3.1 <- function(eps_mean, eps_sd, sample_size = 1){
#   # eps_mean = mean prevalence
#   # eps_sd = sd of prevalence
#   logit_p <- rnorm(n = sample_size, mean = eps_mean, sd = eps_sd)
#   p <- exp(logit_p) / (1 + exp(logit_p)) 
#   return(p)
# }
# ### Ignore spatial term, only apply fixed effect of location (coefficients), i.e. assign fixed value of the coefficients
# getPrev_1.3.2 <- function(beta_mean, beta_sd, sample_size = 1){
#   # beta_mean = mean prevalence
#   # beta_sd = sd of prevalence
#   logit_p <- rnorm(n = sample_size, mean = beta_mean, sd = beta_sd)
#   p <- exp(logit_p) / (1 + exp(logit_p)) 
#   return(p)
# }


######################################
### Model 2: detection probability ###
######################################
### Generate list effects mu_j
genListEffects <- function(p_j){
  # p_j = wanted detection probabilities
  
  logit_p_j <- log(p_j/(1-p_j))
  return(logit_p_j)
}

### Generate alpha - List interaction with Individual Heterogeneity
genAlpha <- function(J, alpha_sd = 100, no_interact = TRUE){ 
  # J = number of Datasets
  # alpha_sd = sd of loading factors
  # no_interact = no behavioral response
  
  if(no_interact){
    mat <- diag(J)
  }
  else{
    mat <- diag(J)
    n <- J*(J-1)/2
    alpha_jrs <- rnorm(n, 0, alpha_sd)
    mat[lower.tri(mat, diag = FALSE)] <- alpha_jrs 
  }
  return(mat) # loading matrix of alpha
}

### Generate individual heterogeneity vector 
genTheta <- function(D = n_latent_factor, sample_size, theta_mean = theta_val, theta_sd = theta_sd_val, theta_corr = theta_corr_val, no_interact = FALSE){
  # p = number of latent factors 
  # sample_size = number of individuals need to draw values
  # theta_mean = mean of individual heterogeneity latent factors, dim = 1*D
  # theta_sd = sd of individual heterogeneity latent factors, dim = 1*D
  # theta_corr = correlation of latent factors, dim = 1*D
  # no_interact = logical operation to have interaction between list and individual heterogeneity
  if(D == 1){
    if(!no_interact){
      theta_vec <- rnorm(sample_size, theta_mean, sd = theta_sd) # Generate theta for all individuals
      theta_mat <- matrix(rep(theta_vec, each = D), ncol = D, byrow = TRUE) # Repeat for J lists
    }
  }
  if(D > 1){
    if(no_interact){ 
      # Model Mt+h
      theta_vec <- rnorm(sample_size, theta_mean, sd = theta_sd) # Generate theta for all individuals
      theta_mat <- matrix(rep(theta_vec, each = D), ncol = D, byrow = TRUE) # Repeat for J lists
    }
    else{
      # Generalized Mth
      # theta_Sigma = Covariance matrix of the latent factors
      # option 1: generate theta Sigma and use inverse wishart distribution
      if(is.null(theta_sd) & is.null(theta_corr)){
        theta_Sigma_array <- rinvwishart(n = 1, nu = D+1, Omega = diag(D))
        theta_Sigma <- apply(theta_Sigma_array, 2, c)  
      }
      if(!is.null(theta_sd) & !is.null(theta_corr)){
        # option 2: import fixed theta Sigma and use multivariate normal distribution
        theta_cov <- theta_sd*theta_sd*theta_corr # compute covariance between the latent factors
        theta_vcov <- diag(D)
        theta_vcov <- ifelse(theta_vcov == 0, theta_cov, theta_vcov) # covariance between latent factors 
        diag(theta_vcov) <-  theta_sd^2 # variance of each latent factor
        theta_Sigma <- theta_vcov
      }
      theta_mean <- theta_mean
      theta_mat <- rmvnorm(sample_size, theta_mean, sigma = theta_Sigma)
    } 
  }  
  return(theta_mat)
} # theta_val matrix for each k

### Get detection probability
getDetectProb <- function(n_target, list_effects, alpha, theta){ 
  # n_target = N_targets[i]
  # list_effects = mu_j
  # alpha = Alpha matrix for factor loading, dim = J*D (n lists * n latent factors)
  # theta = Individual heterogeneity for each k, dim = P_val_k * D  
  
  if(!is.null(alpha)){
    ## With individual heterogeneity, whether interaction or not is controlled inside of getAlpha
    theta_target <- theta[sample(nrow(theta), n_target),] # only keep the target population, dim = n_target*D
    theta_target <- matrix(theta_target, nrow = n_target, ncol = ncol(theta)) # make sure the dim of the theta for target population
    theta_target <- t(theta_target) # dim = D*n_target
    alpha <- as.matrix(alpha)
    loading_mat <- t(alpha%*%theta_target) # dim = n_target*J
    list_effects_mat <- matrix(rep(list_effects, each = n_target), nrow = n_target, byrow = FALSE) # dim = n_target * J
    p <- plogis(list_effects_mat + loading_mat) # dim = number of OUD at single K (n_target) * number of lists (J)
  }
  if(is.null(alpha) & is.null(theta)){
    ## Only lists effect
    list_effects_mat <- matrix(rep(list_effects, each = n_target), nrow = n_target, byrow = FALSE) # dim = n_target * J
    p <- plogis(list_effects_mat) # dim = number of OUD at single K (n_target) * number of lists (J)  
  }
  
  return(p) # individual detection probability across lists at single location k
} # p_detect



###############################
###/// For model variation ///###
###############################

### Generate individual effect ###
## Case 1.2.0
genTheta_1.2.0 <- function(p, sample_size, theta_sd = 0.1){
  # p = number of latent factors 
  # sample_size = number of individuals need to draw values
  # theta_sd = sd of individual heterogeneity
  
  theta_vec <- rnorm(sample_size, 0, sd = theta_sd) # Generate theta for all individuals
  theta_mat <- matrix(rep(theta_vec, each = p), ncol = p, byrow = TRUE) # Repeat for J lists
  return(theta_mat)
}
## Case 1.2.1
genTheta_1.2.1 <- function(p, sample_size, theta_mean = logit(0.01), theta_sd = 0){
  # p = number of latent factors 
  # sample_size = number of individuals need to draw values
  # theta_sd = sd of individual heterogeneity
  theta_vec <- rnorm(sample_size, mean = theta_mean, sd = theta_sd) # Generate theta for all individuals
  theta_mat <- matrix(rep(theta_vec, each = p), ncol = p, byrow = TRUE) # Repeat for J lists
  return(theta_mat)
}
## Case 1.2.2
genTheta_1.2.2 <- function(p, sample_size, betaX1_mean, betaX1_sd, betaX1, fix_beta = TRUE){
  # p = number of latent factors 
  # sample_size = number of individuals need to draw values
  # X1 = a vector of covariate values (fixed effect)
  X1 <- rep(0, sample_size)
  X1[1:sample_size/2] <- 0; X1[(sample_size/2+1):sample_size] <- 1 # generate 50% for X1 = 0 and 50% for X1 = 1
  X1.cat <- length(unique(X1)) # number of groups in covariate X
  if(fix_beta){
    if(X1.cat != length(betaX1)){stop("Error in genFixedEff_Theta: length of betaX1 is not the same as X1.cat")}
    else{
      theta_vec <- model.matrix(~-1 + as.factor(X1)) %*% betaX1 # fixed effect X1^T * betaX1
    }
  }
  if(fix_beta == FALSE){
    if(X1.cat != length(betaX1_mean)| X1.cat != length(betaX1_sd)){stop("Error in genFixedEff_Theta: length of betaX1_mean is not the same as X1.cat")}
    else{
      theta_vec <- rnorm(X1.cat, mean = betaX1_mean, sd = betaX1_sd) # Generate theta for all individuals
    }
  }
  theta_mat <- matrix(rep(theta_vec, each = p), ncol = p, byrow = TRUE) # Repeat for J lists
  return(theta_mat)
}

### Ignore individual heterogeneity, only apply list effects
getDetectProb_1.1 <- function(n_target, list_effects){
  # n_target = N_target[i]
  # list_effects = mu_j
  list_effects_mat <- matrix(rep(list_effects, each = n_target), nrow = n_target, byrow = FALSE) # dim = n_target * J
  p <- plogis(list_effects_mat) # dim = number of OUD at single K (n_target) * number of lists (J)
  return(p) # individual detection probability across lists at single location k
} # p_detect

### Both individual heterogeneity and list effects, no interaction
getDetectProb_1.2 <- function(n_target, list_effects, alpha, theta){
  # n_target = N_targets[i]
  # list_effects = mu_j
  # alpha = Alpha matrix for factor loading
  # theta = Individual heterogeneity 
  theta_target <- theta[sample(nrow(theta), n_target), ] # only keep the target population, dim = n_target * J
  list_effects_mat <- matrix(rep(list_effects, each = n_target), nrow = n_target, byrow = FALSE) # dim = n_target * J
  p <- plogis(list_effects_mat + theta_target%*%alpha) # dim = number of OUD at single K (n_target) * number of lists (J)
  return(p) # individual detection probability across lists at single location k
} # p_detect






###########################################################
### Draw full of all targets across lists and locations ###
###########################################################
genYfull <- function(p){
  # p = p_detect, for one location, given p matrix sample 0/1 across lists
  
  p_vec <- as.vector(p) # append by columns
  yfull_vec <- rbinom(n = length(p_vec), size = 1, prob = p_vec) # Across all detection prob, draw 0/1 for each
  yfull_mat <- matrix(yfull_vec, nrow = nrow(p), ncol = ncol(p), byrow = FALSE) # convert back to p matrix
  return(yfull_mat) # detection history in one location across lists
}

### Manipulate full of all target data shape 
makeYfullTable <- function(yfull){ 
  # yfull = yfull list of yfull_mat, list of all target people, some observed, some not
  
  datalist <- list()
  for(i in 1:length(yfull)) {
    yfull_i <- yfull[[i]]
    obs_label <- apply(yfull_i, 1, max) # Mark if observed or not
    dt_i <- data.frame(yfull_i, target_label = 1, obs_label = obs_label, loc_label = i) 
    datalist[[i]] <- dt_i  
  }
  YfullTable <- do.call(rbind,datalist)
  return(YfullTable)
}
