#########################################################################################################
######## Project 2 Hierarchical SAE with CRC estimator for areal population size estimation #############
#########################################################################################################

## Create initial values for loading factor matrix ##
genLoading_mat_init <- function(seed, J, D){
  set.seed(seed = seed)
  loading_mat_init <- matrix(rnorm(n = J*D, mean = 0, sd = 1), nrow = J, ncol = D)
  loading_mat_init <- round(loading_mat_init, digits = 3)
  loading_mat_init[upper.tri(loading_mat_init)] <- 0  
  if(D == 1){
    loading_mat_init[1,1] <- 1
  }
  return(loading_mat_init)
}

## Create initial values for individual latent factors ##
genTheta_init <- function(seed, aug_size = aug_size, D){
  if(D == 1){
    set.seed(seed)
    theta <- rnorm(n = aug_size, mean = 0, sd = 0.1)
  }
  if(D > 1){
    set.seed(seed)
    theta <-  matrix(rnorm(n = aug_size*D, mean = 0, sd = 0.1), nrow = aug_size, ncol = D)
  }
  return(theta)
}