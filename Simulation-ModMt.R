# Jianing Wang
# Boston University Dissertation
# Chapter 2


########################################################
##################### Fit Model Mt #####################
########################################################

if (!require("dplyr")) install.packages("dplyr") else (require("dplyr", quietly = TRUE)) 
if (!require("reshape2")) install.packages("reshape2") else (require("reshape2", quietly = TRUE)) 
if (!require("stringr")) install.packages("stringr") else (require("stringr", quietly = TRUE)) 
if (!require("nimble")) install.packages("nimble") else (require("nimble", quietly = TRUE)) 
if (!require("coda")) install.packages("coda") else (require("coda", quietly = TRUE)) 
# if (!require("CARBayesdata")) install.packages("CARBayesdata") else (require("CARBayesdata", quietly = TRUE)) 
# if (!require("sp")) install.packages("sp") else (require("sp", quietly = TRUE)) 
# if (!require("spdep")) install.packages("spdep") else (require("spdep", quietly = TRUE)) 
# if (!require("bayesplot")) install.packages("bayesplot") else (require("bayesplot", quietly = TRUE)) 
if (!require("LaplacesDemon")) install.packages("LaplacesDemon") else (require("LaplacesDemon", quietly = TRUE)) 
if (!require("ggmcmc")) install.packages("ggmcmc") else (require("ggmcmc", quietly = TRUE)) 
if (!require("jagsUI")) install.packages("jagsUI") else (require("jagsUI", quietly = TRUE)) # parallel computing
if (!require("parallel")) install.packages("parallel") else (require("parallel", quietly = TRUE)) # parallel computing

#### Choose the Scenario ####
scenario <- 1.3 # Mt model starts with 1.
nDatasets_val <- 3
## If varying or not varying the detection probabilities
DetectProb_scenario <- "Vary" # "Same_large" # "Same_small" # "Vary"
## If including or excluding the extremely small data source
Vary_Extreme_scenario <- "woExtreme" # "woExtreme" # "wExtreme" 

if(scenario < 1.5){
  case <- paste("Simu", scenario, "_J", nDatasets_val, DetectProb_scenario, "PDetect", sep = "")
}else{case <- paste("Simu", scenario, "_J", nDatasets_val, DetectProb_scenario, Vary_Extreme_scenario, "PDetect", sep = "")}

#### Run with R-intel ####
run_with_intel <- "yes" # yes

#### Parallel the Chains ####
parallel_run <- "yes"

#### Remove the list with lowest capturability ####
Remove_Smallest_List <- "yes"

#### Set Working Directory and Path to Output ####
if(scenario < 1.5){
  path_to_input <-  paste(".../Simulation/GenDt_Mt/Simu",scenario,"_J", nDatasets_val, "_", DetectProb_scenario, "PDetect/", sep = "")
}else{path_to_input <-  paste(".../Simulation/GenDt_Mt/Simu",scenario,"_J", nDatasets_val, "_", DetectProb_scenario, Vary_Extreme_scenario, "PDetect/", sep = "")}

if(Remove_Smallest_List == "no"){
  path_to_output <- paste(".../Simulation/Simu_MCMC_Out/Mt/Simu",scenario,"_J", nDatasets_val, "_", DetectProb_scenario, Vary_Extreme_scenario, "PDetect/", sep = "")
}
if(Remove_Smallest_List == "yes"){
  path_to_output <- paste(".../Simulation/Simu_MCMC_Out/Mt/Simu",scenario,"_J", nDatasets_val, "_", DetectProb_scenario, Vary_Extreme_scenario, "PDetect", "woSmallestJ/", sep = "")
}

path_to_funcs <- ".../Simulation/"

#### Load Functions ####
source(paste0(path_to_funcs, "data_simu_functions.r"))
source(paste0(path_to_funcs, "postprocess_functions.r"))
source(paste0(path_to_funcs, "rmSmallestList_functions.r"))

#### Borrow the job ID to be seed ####
task_id <- as.integer(Sys.getenv("SGE_TASK_ID"))
# task_id <- 1

#### Import the Datasets ####
dtYfullTable_one_simu <- readRDS(paste0(path_to_input, "dtYfullTable_simu", task_id, ".RData"))
dtNtarget_one_simu <- readRDS(paste0(path_to_input, "dtNtarget_simu", task_id, ".RData"))
dtyObs_by_loc_one_simu <- readRDS(paste0(path_to_input, "dtyObs_by_loc_simu", task_id, ".RData"))
dtyObs_by_list_loc_simu <- readRDS(paste0(path_to_input, "dtyObs_by_list_loc_simu", task_id, ".RData"))
dtSpatialRE_one_simu <- readRDS(paste0(path_to_input, "dtSpatialRE_simu", task_id, ".RData"))
dtBasicSetUp_one_simu <- readRDS(paste0(path_to_input, "BasicSetUp.RData"))
W_val <- dtBasicSetUp_one_simu$W
D_w_val <- dtBasicSetUp_one_simu$D_w
PDetect_val <- readRDS(paste0(path_to_input, "dtPDetect_simu", task_id, ".RData"))
names_of_k <- dtBasicSetUp_one_simu$names_of_k 

#### Remove the List with the Smallest Catchability ####
if(Remove_Smallest_List == "yes"){
dtYfullTable_one_simu <- remove_smallest_list(dtYfullTable = dtYfullTable_one_simu, nDatasets_val = nDatasets_val)
dtyObs_by_loc_one_simu <- aggregate(dtYfullTable_one_simu$obs_label, by=list(Location=dtYfullTable_one_simu$loc_label), FUN=sum)
PDetect_val <- PDetect_val[-length(PDetect_val)]
}

#### Data Manipulation ####
### Rename ###
single_dt <- dtYfullTable_one_simu # temporarily use one simulated dataset
K <- dtBasicSetUp_one_simu$K # number of subregions
if(Remove_Smallest_List == "yes"){
  J <- dtBasicSetUp_one_simu$J - 1 # number of lists  
}else{
  J <- dtBasicSetUp_one_simu$J # number of lists
  }
D <- dtBasicSetUp_one_simu$D # number of latent factors
P <- dtBasicSetUp_one_simu$P # number of normal people by region
N_target <- dtNtarget_one_simu # number of target people by region
N_obs <- dtyObs_by_loc_one_simu$x # number of observed people by region
TotalN_target <- sum(N_target)
TotalN_obs <- sum(N_obs)
Prev <- N_target/P
### Compute theoretical/empirical list effects in real generated data for posterior check ###
empirical_listeffect <- logit(PDetect_val)
if(Remove_Smallest_List == "yes"){
  theoretical_listeffect <- dtBasicSetUp_one_simu$list_effects[-length(dtBasicSetUp_one_simu$list_effects)]
}else{
  theoretical_listeffect <- dtBasicSetUp_one_simu$list_effects
}

### Create the observed dataset for one single simulated dataset ###
dtYObs <- subset(dtYfullTable_one_simu, obs_label == 1)

#### Data Augmentation ####
prev_ub <- 0.3
M_aug <- round(P*prev_ub,digits = 0)
aug_size <- sum(M_aug) # The augmented size
pseudo_prev <- N_target/M_aug # Effectively what is the true value of prevalence after DA
logit_pseudo_prev <- logit(pseudo_prev) # Effectively what is the true value of lambda after DA

all_zeros_size <- M_aug - N_obs # N of additional zero rows by location
dt_all_zeros <- list()
for(i in 1:K){
  dt_all_zeros_loc_i <- array(0, dim = c(all_zeros_size[i],ncol(dtYObs)-1))
  dt_all_zeros_loc_i <- cbind(dt_all_zeros_loc_i, rep(i, nrow(dt_all_zeros_loc_i)))
  dt_all_zeros_loc_i <- as.data.frame(dt_all_zeros_loc_i)
  colnames(dt_all_zeros_loc_i) <- colnames(dtYObs)
  dt_all_zeros[[i]] <- dt_all_zeros_loc_i
}
dt_all_zeros <- do.call(rbind,dt_all_zeros)
## Augmentation to make all locations have untarget + target people
dtYObs_aug <- rbind(dtYObs, dt_all_zeros)
## Order by observed first, and then unobserved, nrow = sum(P)
dtYObs_aug <- arrange(dtYObs_aug, desc(obs_label))
dtYObs_yaug <- dtYObs_aug[,seq(J)] # this is the dataset to be passed into the model
########## END WITH DATA AUGEMENTATION ##########


#### Prepare Intrinsic CAR (ICAR) model specification ####
## Reshape the adjacency matrix to get adjacency list 
## (i.e. A vector of indices indicating which regions are neighbors of which.)
## Change the character name of Adj matrix to numeric name 
W_numeric_name <- data.frame(name = names_of_k, id = seq(1:length(names_of_k)))
row.names(W_val) <- ifelse(row.names(W_val) %in% W_numeric_name$name, W_numeric_name$id,NA)
colnames(W_val) <- ifelse(colnames(W_val) %in% W_numeric_name$name, W_numeric_name$id,NA)
W <- W_val
W_ls <- reshape2::melt(W)
Adj_vec <- W_ls$Var1[W_ls$value==1]

## A vector of weights. In this case, all weights are 1.
Weights_vec <- W_ls$value[W_ls$value==1]

## A vector of length N. num[n] indicates how many neighbors region n contains.
## This helps map the adj vector to the starting region.
D_w <- diag(apply(W, 2, sum))
D_w_ls <- reshape2::melt(D_w)
D_w_ls <- D_w_ls$value[which(D_w_ls$value>0)]


#### Set up parallel computation ####
## First you create a cluster using the total amount of cores you have but one to make sure your computer can go on working:
nbcores <- detectCores() - 1 
my_cluster <- makeCluster(nbcores)

## Then you wrap your workflow in a function to be run in parallel:
workflow <- function(seed, config_param){
  
  if (!require("nimble")) install.packages("nimble") else (require("nimble", quietly = TRUE)) 
  if (!require("matrixsampling")) install.packages("matrixsampling") else (require("matrixsampling", quietly = TRUE)) 
  
  ### Create a NIMBLE model ###
  ## Define the model code, all nodes, variables and relationship, and all constants
  Mod_Mt_ICAR_Code <- nimbleCode( {
    ############
    ## Priors ##
    ############
    
    # Spatial RE [ICAR] #
    alpha ~ dflat() # vague uniform prior due to centering constraint at 0
    
    # variance parameter of the phi #
    ## Option 1: Variance ~ Inverse Gamma, Precision tau2 ~ Gamma
    # tau2 ~ dgamma(shape = 0.5, scale = 1/0.5) # inverse variance of spatial component (precision of the ICAR Component, equivalent variance ~ IG(0.5,0.5))
    ## Option 2: Scale ~ Half-Cauchy (t distribution with df = 1)
    inv_tau  ~ T(dt(mu = 0, tau = sd_spRE_prior, df = 1), 0, ) # half-t truncated at 0 and set no upper limit
    tau2 <- 1/(inv_tau^2) # compute precision  
    
    # phi from ICAR #
    phi[1:K] ~ dcar_normal(adj[1:L], weights[1:L], num[1:K], tau2, zero_mean = 1) # L = length of Adj_vec, K = n of regions, zero_mean = 1 put constraint centers at 0
    
    # List_Effects FE #
    for(j in 1:J){ 
      logit_list_effects[j] ~ dnorm(0, sd = sigma_mu) 
    }
    
    ################
    ## Likelihood ##
    ################
    
    # Prevalence Level #
    for (i in 1:M){ # M = augmented large population
      
      logit(lambda[i]) <- alpha + phi[prev_re[i]] # pseudo logit prev via nested indicators, dim = length of M
      
      z[i] ~ dbern(lambda[i]) # inclusion indicator to be a target subject, same location's subjects share the same pseudo prev
      
      for(j in 1:J){
        # detection prob model #
        logit(p.detect[i,j]) <- logit_list_effects[j] # detection probability
        
        # compute effective p #
        p.eff[i,j] <- z[i] * p.detect[i,j]  # being a potential target * being observed
        y[i,j] ~ dbern(p.eff[i,j]) # individual capture histories  
      } # j
    } # i
    
    #########################
    ## Derive the quantity ##
    #########################
    # <1> Number of target people
    N <- sum(z[1:M]) 
    # <2> Overall baseline value of lambda (to compute baseline prevalence)
    beta0 <- alpha # compute overall baseline value of lambda
    # <3> Variance of spatial component
    # sigma2.phi <- 1/tau2 # variance of spatial component using inverse gamma for variance parameter
    sigma2.phi <- inv_tau^2 # variance of spatial component using half-cauchy for variance parameter
    # <4> Area-specific residual
    resPhi[1:K] <- phi[1:K]
    # <5> Overall residuals at individual level
    # res[1:M] <- (z[1:M]-lambda[1:M])/sqrt(lambda[1:M])          
    # <6> Lambda
    for(k in 1:K){
      pesudo_logit_prev[k] <- alpha + phi[k]
    }
    
  })
  ########## FINISH BAYESIAN MODEL WRITTING ##########
  
  run_stepwise_mcmc <- "yes" # "yes", if run stepwise, then the initial values cannot be specified with nested list.
  
  
  ## Use only when parallel computation ##
  J <- config_param$Constants$J
  K <- config_param$Constants$K
  aug_size <- config_param$Constants$M
  
  ### Configuration of MCMC ###
  ## Initials
  set.seed(123) # ensure reproducibility while drawing randomly initial values.
  
  if(run_stepwise_mcmc == "yes"){
    Mod_Mt_ICAR_Inits <- list(alpha = runif(1,0.1,0.5), 
                              tau2 = rgamma(1, shape = 0.5, scale = 1/0.5),  # if using inverse gamma prior distribution for spRE var parameter
                              inv_tau = runif(1,1,10), # if using half-cauchy prior distribution for spRE var parameter
                              logit_list_effects = rep(runif(1,-2,2),J), 
                              phi = runif(K,-0.2,0.5),
                              z = rep(1, aug_size))
  }
  
  #############################
  ### Run the MCMC sampling ###
  #############################
  ## Run MCMC step by step ##
  ## Build model ##
  Mod_Mt_ICAR <- nimbleModel(code = Mod_Mt_ICAR_Code, name = "Mod_Mt_ICAR",
                             constants = config_param$Constants,
                             data = config_param$Data, inits = Mod_GenMt_ICAR_Inits,
                             calculate = FALSE # disable the calculation of all deterministic nodes and log-likelihood
  )
  ## Compile the nimble model ##
  cMod_Mt_ICAR <- compileNimble(Mod_Mt_ICAR)
  ## (1) Configuration
  Mod_Mt_ICARConf <- configureMCMC(Mod_Mt_ICAR,
                                   useConjugacy = FALSE, # because the model is anyway not conjugacy, no need to check, speed up
                                   monitors = c("N","pesudo_logit_prev", "sigma2.phi", "resPhi", "beta0","logit_list_effects"),
                                   enableWAIC = nimbleOptions('enableWAIC' = TRUE),
                                   print = TRUE)
  ## (2) Build MCMC algorithm for the model
  Mod_Mt_ICARMCMC <- buildMCMC(Mod_Mt_ICARConf)
  ## (3) Compile C++ after configuration
  C_mcmc_Mod_Mt_ICAR <- compileNimble(Mod_Mt_ICARMCMC, project = Mod_Mt_ICAR)
  ## (4) Run MCMC sampling
  samples_Mod_Mt_ICAR <- runMCMC(C_mcmc_Mod_Mt_ICAR, 
                                 niter = config_param$niter, 
                                 nburnin = config_param$nburnin, 
                                 nchains = config_param$nchains,  # comment this line because we set three chains outside
                                 thin = config_param$nthin, 
                                 setSeed = seed,  # control the state of R the random number generator for reproducibility
                                 # progressBar = TRUE,
                                 summary = TRUE,
                                 WAIC = TRUE,
                                 samplesAsCodaMCMC = TRUE)
  ########## END OF RUNNING MODELS ##########
  
  
  return(samples_Mod_Mt_ICAR)
  
} ######## END OF PARALLEL COMPUTATION WRAP ##########

## Now we run the code using parLapply(), which uses cluster nodes to execute our workflow:

## Configuration of MCMC sampling 
nchains <- 3  # comment this line because we set three chains outside
niter <- 100000 
nburnin <- 50000 
nthin <- 100

## Constants 
Mod_Mt_ICAR_Consts <- list(M = aug_size,
                           J = J,
                           K = K,
                           L = length(Adj_vec),
                           adj = Adj_vec,
                           weights = Weights_vec,
                           num = D_w_ls, # number of neighbors
                           sd_spRE_prior = 10, # standard deviation of half-cauchy prior distribution for spRE variance parameter
                           sigma_mu = 10, # list effects sigma
                           prev_re =  dtYObs_aug$loc_label # location label for each individual
)


## Data, pass in y 
Mod_Mt_ICAR_Data <- list(y = dtYObs_yaug)

## Put together into a list 
Mod_Mt_ICAR_config <- list(nchains = nchains,
                               niter = niter,
                               nburnin = nburnin,
                               nthin = nthin,
                               Constants = Mod_Mt_ICAR_Consts,
                               Data = Mod_Mt_ICAR_Data)

## Run parallel versions of lapply
output <- parLapply(cl = my_cluster, 
                    X = 9, 
                    fun = workflow, 
                    config_param = Mod_Mt_ICAR_config ## Configuration, Data and Constants ##
)

## Close the cluster with stopCluster() so that processes do not continue to run in the background and slow down other processes:
stopCluster(my_cluster)
########## END WITH PARALLEL COMPUTATION ###########

#### Post-Process MCMC Samples ####
# length(output) = 1
samples_Mod_Mt_ICAR <- output[[1]] # length = 1



#############################################
### Save the outputs and plot the results ###
#############################################

chain1 <- samples_Mod_Mt_ICAR$samples$chain1
chain2 <- samples_Mod_Mt_ICAR$samples$chain2
chain3 <- samples_Mod_Mt_ICAR$samples$chain3

names_chain1 <- colnames(chain1)
names_chain2 <- colnames(chain2)
names_chain3 <- colnames(chain3)

nsamples <- (niter - nburnin)/nthin

###################
## Target size N ##
###################
samp1_N <- chain1[,'N']
samp2_N <- chain2[,'N']
samp3_N <- chain3[,'N']

## Combine chains ##
samp_combo_N <- c(samp1_N, samp2_N, samp3_N)

################
## Prevalence ##
################
# Keep only lambda parameters' samp
samp1_lambda <- save_lambda(samples = chain1, names_in_chain = names_chain1)
samp2_lambda <- save_lambda(samples = chain2, names_in_chain = names_chain2)
samp3_lambda <- save_lambda(samples = chain3, names_in_chain = names_chain3)

ls_samp_lambda <- list(chain1 = samp1_lambda,
                       chain2 = samp2_lambda,
                       chain3 = samp3_lambda)

## Effective prevalence ##
samp1_effprev <- invlogit(samp1_lambda) * prev_ub
samp2_effprev <- invlogit(samp2_lambda) * prev_ub
samp3_effprev <- invlogit(samp3_lambda) * prev_ub

## Combine chains
samp_combo_effprev <- as.data.frame(rbind(samp1_effprev, samp2_effprev, samp3_effprev))

################################
## Variance parameter of spRE ##
################################
samp1_sigma2.phi <- as.vector(chain1[,'sigma2.phi'])
samp2_sigma2.phi <- as.vector(chain2[,'sigma2.phi'])
samp3_sigma2.phi <- as.vector(chain3[,'sigma2.phi'])

## Combine chains
samp_combo_Sigma2.phi <- c(samp1_sigma2.phi, samp2_sigma2.phi, samp3_sigma2.phi)

###########################################
## Overall effective baseline prevalence ##
###########################################
samp1_beta0 <- as.vector(chain1[,'beta0'])
samp2_beta0 <- as.vector(chain2[,'beta0'])
samp3_beta0 <- as.vector(chain3[,'beta0'])

## Effective overall prevalence 
samp1_eff_baselineprev <- invlogit(samp1_beta0) * prev_ub
samp2_eff_baselineprev <- invlogit(samp2_beta0) * prev_ub
samp3_eff_baselineprev <- invlogit(samp3_beta0) * prev_ub

## Combine chains
samp_combo_Baseline_prev <- c(samp1_eff_baselineprev, samp2_eff_baselineprev, samp3_eff_baselineprev)

#################################
## Spatial residual phi values ##
#################################
samp1_resPhi <- save_resPhi(samp = chain1, names_in_chain = names_chain1) # nsamp * K
samp2_resPhi <- save_resPhi(samp = chain2, names_in_chain = names_chain2)
samp3_resPhi <- save_resPhi(samp = chain3, names_in_chain = names_chain3)

## Combine chains
samp_combo_resPhi <- data.frame(rbind(samp1_resPhi, samp2_resPhi, samp3_resPhi)) # (nchains * nsamp) * K

######################
## List Effect (mu) ##
######################
samp1_logitlisteffect <- as.matrix(save_logitlisteffect(samples = chain1, names_in_chain = names_chain1))
samp2_logitlisteffect <- as.matrix(save_logitlisteffect(samples = chain2, names_in_chain = names_chain2))
samp3_logitlisteffect <- as.matrix(save_logitlisteffect(samples = chain3, names_in_chain = names_chain3))

## Combine chains
samp_combo_logitlisteffect <- as.data.frame(rbind(samp1_logitlisteffect, samp2_logitlisteffect, samp3_logitlisteffect))


##################################
### Save the samples from MCMC ###
##################################
## Total size N 
chain_index <- rep(1:nchains,each = nsamples)
combodt_samp_N <- data.frame(Chains = chain_index,
                             Combo_samp_N = samp_combo_N)
colnames(combodt_samp_N) <- c("Chains","Ntarget")

## Computed prevalence 
chain_index <- rep(1:nchains,each = nsamples)
combodt_samp_prev <- data.frame(Chains = chain_index,
                                Combo_samp_prev = samp_combo_effprev)
colnames(combodt_samp_prev) <- c("Chains",paste("effprev_K",seq(1,K),sep=""))

## Sigma2 in Phi 
chain_index <- rep(1:nchains,each = nsamples)
combodt_samp_sigma2.phi <- data.frame(Chains = chain_index,
                                      Combo_samp_sigma2.phi = samp_combo_Sigma2.phi)
colnames(combodt_samp_sigma2.phi) <- c("Chains","sigma2.phi")

## Baseline overall prevalence 
chain_index <- rep(1:nchains,each = nsamples)
combodt_samp_effbaselineprev <- data.frame(Chains = chain_index,
                                           Combo_samp_eff_baselineprev = samp_combo_Baseline_prev)
colnames(combodt_samp_effbaselineprev) <- c("Chains","effbaselineprev")

## Residual Phi 
chain_index <- rep(1:nchains,each = nsamples)
combodt_samp_resPhi <- data.frame(Chains = chain_index,
                                  Combo_samp_resPhi = samp_combo_resPhi)
colnames(combodt_samp_resPhi) <- c("Chains",names_of_k) # Order resPhi1:K

## Logit list effect
chain_index <- rep(1:nchains,each = nsamples)
combodt_samp_logitlisteffect <- data.frame(Chains = chain_index,
                                           Combo_samp_logitlisteffect = samp_combo_logitlisteffect)
colnames(combodt_samp_logitlisteffect) <- c("Chains",paste("logitlisteffect_J",seq(1,J),sep = ""))

## Extract the WAIC
mcmc.out.WAIC <- samples_Mod_Mt_ICAR$WAIC

## Extract all chains summary
mcmc.out.summary <- samples_Mod_Mt_ICAR$summary$all.chains

## Save the true values
TrueVals <- list(Parallel_Run = parallel_run,
                 Intel_Run = run_with_intel,
                 TotalN_target = TotalN_target,
                 TotalN_obs = TotalN_obs,
                 N_target = N_target,
                 Prev = Prev,
                 Theoretical_list_effects = theoretical_listeffect, # theoretical list effect designed in simulation
                 Empirical_list_effects = empirical_listeffect, # empirical log of pdetect in simulated data
                 Basline_prev = unique(invlogit(dtBasicSetUp_one_simu$beta0)),
                 Sigma2.phi = dtBasicSetUp_one_simu$tau2, # true spatial RE variance
                 SpatialRE = dtSpatialRE_one_simu
                 )

## Save MCMC Setup 
NIMBLE_Setup <- list(case = case,
                     Remove_Smallest_List = Remove_Smallest_List,
                     nchains = nchains,
                     niter = niter,
                     nburnin = nburnin,
                     nthin = nthin,
                     nsamples = nsamples,
                     sd_spRE_prior = Mod_Mt_ICAR_Consts$sd_spRE_prior, # prior value for spatial RE half-cauchy variance
                     sigma_mu = Mod_Mt_ICAR_Consts$sigma_mu # prior value for list effect variance,
                     )

## Save data
saveRDS(samples_Mod_Mt_ICAR$samples, file = paste0(path_to_output, "mcmc.out_ModMtICAR_samples", "_simu", task_id, ".Rdata"))
saveRDS(combodt_samp_N, file = paste0(path_to_output, "mcmc.out_ModMtICAR_N", "_simu", task_id, ".Rdata"))
saveRDS(combodt_samp_prev, file = paste0(path_to_output, "mcmc.out_ModMtICAR_prev", "_simu", task_id, ".Rdata"))
saveRDS(combodt_samp_logitlisteffect, file = paste0(path_to_output, "mcmc.out_ModMtICAR_logitlisteffect", "_simu", task_id, ".Rdata"))
saveRDS(combodt_samp_sigma2.phi, file = paste0(path_to_output, "mcmc.out_ModMtICAR_sigma2.phi", "_simu", task_id, ".Rdata"))
saveRDS(combodt_samp_effbaselineprev, file = paste0(path_to_output, "mcmc.out_ModMtICAR_effbaselineprev", "_simu", task_id, ".Rdata"))
saveRDS(combodt_samp_resPhi, file = paste0(path_to_output, "mcmc.out_ModMtICAR_resPhi", "_simu", task_id, ".Rdata"))
saveRDS(mcmc.out.WAIC, file = paste0(path_to_output, "mcmc.out_ModMtICAR_WAIC", "_simu", task_id, ".Rdata"))
saveRDS(mcmc.out.summary, file = paste0(path_to_output, "mcmc.out_ModMtICAR_summary", "_simu", task_id, ".Rdata"))
saveRDS(TrueVals, file = paste0(path_to_output, "TrueVals", "_simu", task_id, ".Rdata"))
saveRDS(NIMBLE_Setup, file = paste0(path_to_output, "NIMBLE_Setup", "_simu", task_id, ".Rdata"))
