# Jianing Wang
# Boston University Dissertation
# Chapter 2

######################################################
##### Run Candidate Model: Log-linear CRC model ######
######################################################

if (!require("MASS")) install.packages("MASS") else (require("MASS", quietly = TRUE)) 
if (!require("dplyr")) install.packages("dplyr") else (require("dplyr", quietly = TRUE)) 
if (!require("tidyr")) install.packages("tidyr") else (require("tidyr", quietly = TRUE)) 
if (!require("data.table")) install.packages("data.table") else (require("data.table", quietly = TRUE)) 
if (!require("reshape2")) install.packages("reshape2") else (require("reshape2", quietly = TRUE)) 
if (!require("ggplot2")) install.packages("ggplot2") else (require("ggplot2", quietly = TRUE)) 
if (!require("stringr")) install.packages("stringr") else (require("stringr", quietly = TRUE)) 

## Path to functions
path_to_funcs <- ".../Simulation/"
source(paste(path_to_funcs, "data_simu_functions.r", sep = ""))
source(paste(path_to_funcs, "postprocess_functions.r", sep = ""))
source(paste0(path_to_funcs, "modelfit_functions.r"))

## Save data
path_to_data <- ".../Simulation/MapData/"

# Note
# To show the superiority of the proposed model, we compare with the independent log-linear model that assume independent fit to each location.

## Import the counties files ##
AdjMat_Ohio <- read.csv(paste0(path_to_data, "OHData/OHcounties_adjacency2010.csv", sep = ""))
### Adjacency matrix ###
W_val <- as.matrix(AdjMat_Ohio)
rownames(W_val) <- colnames(AdjMat_Ohio)
colnames(W_val) <- colnames(AdjMat_Ohio)
### Compute the Neighbors matrix and output for model fitting use ###
D_w_val <- read.csv(paste0(path_to_data, "OHData/nNeighbor_mat_OHcounties.csv", sep = "")) # Ajacency matrix with number of neighbors
D_w_val <- as.matrix(D_w_val)
rownames(D_w_val) <- colnames(D_w_val)
names_all_areas <- colnames(W_val) # names of all counties from adjacency matrix


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


#### Set Working Directory ####
if(scenario < 1.8){
  if(scenario < 1.5){
    path_to_input <-  paste(".../Simulation/GenDt_Mt/Simu",scenario,"_J", nDatasets_val, "_", DetectProb_scenario, "PDetect/", sep = "")
    path_to_BHM_output <- paste(".../Simulation/Simu_MCMC_Out/Mt/Simu",scenario,"_J", nDatasets_val, "_", DetectProb_scenario, "PDetect/", sep = "")
    path_to_LogLinear_output <- ".../Simulation/Simu_Loglinear_Out/Mt/"
  }else{
    path_to_input <-  paste(".../Simulation/GenDt_Mt/Simu",scenario,"_J", nDatasets_val, "_", DetectProb_scenario, Vary_Extreme_scenario, "PDetect/", sep = "")
    path_to_BHM_output <- paste(".../Simulation/Simu_MCMC_Out/Mt/Simu",scenario,"_J", nDatasets_val, "_", DetectProb_scenario, Vary_Extreme_scenario, "PDetect/", sep = "")
    path_to_LogLinear_output <- ".../Simulation/Simu_Loglinear_Out/Mt/"
    }
}
if(scenario > 1.8){
  if(scenario < 1.5){
    path_to_input <-  paste(".../Simulation/GenDt_Mth/Simu",scenario,"_J", nDatasets_val, "_", DetectProb_scenario, "PDetect/", sep = "")
    path_to_BHM_output <- paste(".../Simulation/Simu_MCMC_Out/Mth/Simu",scenario,"_J", nDatasets_val, "_", DetectProb_scenario, "PDetect/", sep = "")
    path_to_LogLinear_output <- ".../Simulation/Simu_Loglinear_Out/Mth/"
  }else{
    path_to_input <-  paste(".../Simulation/GenDt_Mth/Simu",scenario,"_J", nDatasets_val, "_", DetectProb_scenario, Vary_Extreme_scenario, "PDetect/", sep = "")
    path_to_BHM_output <- paste(".../Simulation/Simu_MCMC_Out/Mth/Simu",scenario,"_J", nDatasets_val, "_", DetectProb_scenario, Vary_Extreme_scenario, "PDetect/", sep = "")
    path_to_LogLinear_output <- ".../Simulation/Simu_Loglinear_Out/Mth/"
  }
}


### Find the number of simulations that BHM ran, use the same data here
TrueVals <- list.files(path = path_to_BHM_output, pattern = 'TrueVals_(simu\\d+).Rdata$')
# Find number of simulations
nsim <- length(TrueVals)
# Find which simulation runs through
simu_ran_through <- c()
for(task_id in 1:nsim){
  which_ran_through <- str_extract_all(TrueVals[task_id],"\\(?[0-9]+\\)?")
  simu_ran_through <- c(simu_ran_through,which_ran_through)
}
simu_ran_through <- sort(as.numeric(do.call(c,simu_ran_through)))



############################
### Fit Log-Linear Model ###
############################

### For each generated data, fit log-linear model with Negative Binomial Distribution assumption ###
rel.bias <- c()
AbsoluteError <- c()
SquaredError <- c()
rel.bias_prev <- matrix(0, nrow = length(names_all_areas), ncol = nsim)
AbsoluteError_prev <- matrix(0, nrow = length(names_all_areas), ncol = nsim)
SquaredError_prev <- matrix(0, nrow = length(names_all_areas), ncol = nsim)

for(i in 1:nsim){
  task_id <- simu_ran_through[i]
  ### Read the data  
  dtYfullTable_one_simu <- readRDS(paste0(path_to_input, "dtYfullTable_simu", task_id, ".RData"))
  dtNtarget_one_simu <- readRDS(paste0(path_to_input, "dtNtarget_simu", task_id, ".RData"))
  dtBasicSetUp_one_simu <- readRDS(paste0(path_to_input, "BasicSetUp.RData"))
  P <- dtBasicSetUp_one_simu$P
  J <- dtBasicSetUp_one_simu$J
  K <- dtBasicSetUp_one_simu$K
  
  ### Create the contingency table 
  if(J == 3){
    empty_grid_list <- expand.grid(List1 = c(0,1), List2 = c(0,1), List3 = c(0,1), K = seq(1:K)) 
    empty_grid_list <- empty_grid_list[-which(empty_grid_list$List1 == 0 & empty_grid_list$List2 == 0 & empty_grid_list$List3 == 0),]
    yObs_by_list_loc <- aggregate(dtYfullTable_one_simu$obs_label, by=list(List1 = dtYfullTable_one_simu$X1,
                                                                           List2 = dtYfullTable_one_simu$X2,
                                                                           List3 = dtYfullTable_one_simu$X3,
                                                                           K = dtYfullTable_one_simu$loc_label), FUN=sum)
    # Get complete yObs table by list and location
    yObs_contingency_final <- merge(x = empty_grid_list, y = yObs_by_list_loc, by.x = c("List1","List2","List3", "K"),
                                    by.y = c("List1","List2","List3", "K"), 
                                    all.x = TRUE)
    yObs_contingency_final$x <- ifelse(is.na(yObs_contingency_final$x),0, yObs_contingency_final$x)
    colnames(yObs_contingency_final) <- c("List1","List2","List3","K","n_obs")
    yObs_contingency_final <- with(yObs_contingency_final, yObs_contingency_final[order(K, List1, List2, List3),])
  }
  if(J == 5){
    empty_grid_list <- expand.grid(List1 = c(0,1), List2 = c(0,1), List3 = c(0,1), List4 = c(0,1), List5 = c(0,1), K = seq(1:K)) 
    empty_grid_list <- empty_grid_list[-which(empty_grid_list$List1 == 0 & empty_grid_list$List2 == 0 & empty_grid_list$List3 == 0 & empty_grid_list$List4 == 0 & empty_grid_list$List5 == 0),]
    yObs_by_list_loc <- aggregate(dtYfullTable_one_simu$obs_label, by=list(List1 = dtYfullTable_one_simu$X1,
                                                                           List2 = dtYfullTable_one_simu$X2,
                                                                           List3 = dtYfullTable_one_simu$X3,
                                                                           List4 = dtYfullTable_one_simu$X4,
                                                                           List5 = dtYfullTable_one_simu$X5,
                                                                           K = dtYfullTable_one_simu$loc_label), FUN=sum)
    # Get complete yObs table by list and location
    yObs_contingency_final <- merge(x = empty_grid_list, y = yObs_by_list_loc, by.x = c("List1","List2","List3","List4","List5", "K"),
                                    by.y = c("List1","List2","List3","List4","List5", "K"), 
                                    all.x = TRUE)
    yObs_contingency_final$x <- ifelse(is.na(yObs_contingency_final$x),0, yObs_contingency_final$x)
    colnames(yObs_contingency_final) <- c("List1","List2","List3","List4","List5","K","n_obs")
    yObs_contingency_final <- with(yObs_contingency_final, yObs_contingency_final[order(K, List1, List2, List3, List4, List5),])
  }
  if(J == 7){
    empty_grid_list <- expand.grid(List1 = c(0,1), List2 = c(0,1), List3 = c(0,1), List4 = c(0,1), List5 = c(0,1), List6 = c(0,1), List7 = c(0,1), K = seq(1:K)) 
    empty_grid_list <- empty_grid_list[-which(empty_grid_list$List1 == 0 & empty_grid_list$List2 == 0 & empty_grid_list$List3 == 0 & empty_grid_list$List4 == 0 & empty_grid_list$List5 == 0 & empty_grid_list$List6 == 0 & empty_grid_list$List7 == 0),]
    yObs_by_list_loc <- aggregate(dtYfullTable_one_simu$obs_label, by=list(List1 = dtYfullTable_one_simu$X1,
                                                                           List2 = dtYfullTable_one_simu$X2,
                                                                           List3 = dtYfullTable_one_simu$X3,
                                                                           List4 = dtYfullTable_one_simu$X4,
                                                                           List5 = dtYfullTable_one_simu$X5,
                                                                           List6 = dtYfullTable_one_simu$X6,
                                                                           List7 = dtYfullTable_one_simu$X7,
                                                                           K = dtYfullTable_one_simu$loc_label), FUN=sum)
    # Get complete yObs table by list and location
    yObs_contingency_final <- merge(x = empty_grid_list, y = yObs_by_list_loc, by.x = c("List1","List2","List3","List4","List5","List6","List7","K"),
                                    by.y = c("List1","List2","List3","List4","List5","List6","List7", "K"), 
                                    all.x = TRUE)
    yObs_contingency_final$x <- ifelse(is.na(yObs_contingency_final$x),0, yObs_contingency_final$x)
    colnames(yObs_contingency_final) <- c("List1","List2","List3","List4","List5","List6","List7","K","n_obs")
    yObs_contingency_final <- with(yObs_contingency_final, yObs_contingency_final[order(K, List1, List2, List3, List4, List5, List6, List7),])
  }
  
  
  ### Separate into K location-wise contingency table ###
  yObs_contingency_by_loc_ls <- split(yObs_contingency_final, f = yObs_contingency_final$K)    
  
  ## NB Models
  tryCatch(
    expr = {if(J == 3){
      est_total_by_loc <- sapply(1:K, FUN = function(k){
        dt <- as.data.frame( yObs_contingency_by_loc_ls[k])
        colnames(dt) <- c("List1","List2","List3","K","n_obs")
        dt$List1 <- as.character(dt$List1)
        dt$List2 <- as.character(dt$List2)
        dt$List3 <- as.character(dt$List3)
        if(sum(dt$n_obs)>0){
          summary_nb_mod_k <- summary(nb_mod <- glm.nb(n_obs ~ List1 + List2 + List3, data = dt))
          intercept_k <- summary_nb_mod_k$coefficients[1,1]
          est_unknown_k <- exp(intercept_k)
        }
        if(sum(dt$n_obs) == 0){
          est_unknown_k <- 0
        }
        known_k <- sum(dt$n_obs)
        est_total_k <- round(est_unknown_k + known_k, digits = 0)
        return(est_total_k)
      })
    }
    },
    error = function(e){
      print(paste("Error in model fit", i, sep = " "))
    }
  )
  
  tryCatch(
    expr = {if(J == 5){
      est_total_by_loc <- sapply(1:K, FUN = function(k){
        dt <- as.data.frame( yObs_contingency_by_loc_ls[k])
        colnames(dt) <- c("List1","List2","List3","List4","List5","K","n_obs")
        dt$List1 <- as.character(dt$List1)
        dt$List2 <- as.character(dt$List2)
        dt$List3 <- as.character(dt$List3)
        dt$List4 <- as.character(dt$List4)
        dt$List5 <- as.character(dt$List5)
        if(sum(dt$n_obs)>0){
          summary_nb_mod_k <- summary(nb_mod <- glm.nb(n_obs ~ List1 + List2 + List3 + List4 + List5, data = dt))
          intercept_k <- summary_nb_mod_k$coefficients[1,1]
          est_unknown_k <- exp(intercept_k)
        }
        if(sum(dt$n_obs) == 0){
          est_unknown_k <- 0
        }
        known_k <- sum(dt$n_obs)
        est_total_k <- round(est_unknown_k + known_k, digits = 0)
        return(est_total_k)
      })
    }
    },
    error = function(e){
      print(paste("Error in model fit", i, sep = " "))
    }
  )
  
  tryCatch(
    expr = {if(J == 7){
      est_total_by_loc <- sapply(1:K, FUN = function(k){
        dt <- as.data.frame( yObs_contingency_by_loc_ls[k])
        colnames(dt) <- c("List1","List2","List3","List4","List5","List6","List7","K","n_obs")
        dt$List1 <- as.character(dt$List1)
        dt$List2 <- as.character(dt$List2)
        dt$List3 <- as.character(dt$List3)
        dt$List4 <- as.character(dt$List4)
        dt$List5 <- as.character(dt$List5)
        dt$List6 <- as.character(dt$List6)
        dt$List7 <- as.character(dt$List7)
        if(sum(dt$n_obs)>0){
          summary_nb_mod_k <- summary(nb_mod <- glm.nb(n_obs ~ List1 + List2 + List3 + List4 + List5 + List6 + List7, data = dt))
          intercept_k <- summary_nb_mod_k$coefficients[1,1]
          est_unknown_k <- exp(intercept_k)
        }
        if(sum(dt$n_obs) == 0){
          est_unknown_k <- 0
        }
        known_k <- sum(dt$n_obs)
        est_total_k <- round(est_unknown_k + known_k, digits = 0)
        return(est_total_k)
      })
    }
    },
    error = function(e){
      print(paste("Error in model fit", i, sep = " "))
    }
  )
  
  
  ### Compute relative bias on total Ntarget
  rel.bias_task_id <- sum(est_total_by_loc)/sum(dtNtarget_one_simu)
  rel.bias <- c(rel.bias, rel.bias_task_id)
  
  ### Compute RMSE/MAE on total Ntarget
  Error_task_id <- sum(est_total_by_loc) - sum(dtNtarget_one_simu)
  SquaredError_task_id <- (Error_task_id)^2
  AbsoluteError <- c(AbsoluteError, abs(Error_task_id))
  SquaredError <- c(SquaredError, SquaredError_task_id)
  
  ### Compute relative bias on area-specific prevalence
  est_prev_by_loc <- est_total_by_loc/P
  true_prev_by_loc <- dtNtarget_one_simu/P
  rel.bias_prev_task_id <- est_prev_by_loc/true_prev_by_loc
  rel.bias_prev[,i] <- rel.bias_prev_task_id
  
  ### Compute RMSE/MAE on area-specific prevalence
  Error_prev_task_id <- est_prev_by_loc - true_prev_by_loc
  AbsoluteError_prev[,i] <- abs(Error_prev_task_id)
  SquaredError_prev_task_id <- (Error_prev_task_id)^2
  SquaredError_prev[,i] <- SquaredError_prev_task_id
  
  print(i)
} # end of the loop for task_id

### Compute Mean of the relative bias for total Ntarget
rel.bias_final <- rel.bias[which(rel.bias < 100000)]
Mean.rel.bias <- mean(rel.bias_final)
### Compute RMSE for total Ntarget
SquaredError_final <- SquaredError[which(SquaredError < 100000)]
RMSE <- sqrt(mean(SquaredError_final))
### Compute MAE for total Ntarget
AbsoluteError_final <- AbsoluteError[which(AbsoluteError < 100000)]
MAE <- median(AbsoluteError_final)


### Summary for area-specific prevalence
rel.bias_prev_final <- apply(rel.bias_prev, c(1,2), FUN = function(x) {
  if(x > 5) {x <- NA}
  else{x <- x}
  return(x)})
### Number of times across nsim for each area to produce unstable estimates
freq_infinite_est_by_loc <- apply(rel.bias_prev_final, 1, FUN = function(x) sum(is.na(x)))
freq_infinite_est_by_loc_df <- data.frame(K = names_all_areas, freq_infinite_est_by_loc = freq_infinite_est_by_loc)
### Number of areas that never have unstable estimates
length(freq_infinite_est_by_loc_df$K[freq_infinite_est_by_loc_df$freq_infinite_est_by_loc== 0])
### Compute Mean of the relative bias for area-specific prevalence
Mean.rel.bias_prev <- apply(rel.bias_prev_final, 1, FUN = function(x){mean(x, na.rm = TRUE)})

AbsoluteError_prev_final <- apply(AbsoluteError_prev, c(1,2), FUN = function(x) {
  if(x > 1) {x <- NA}
  else{x <- x}
  return(x)})
### Compute Mean of the relative bias for area-specific prevalence
MAE_prev <- apply(AbsoluteError_prev_final, 1, FUN = function(x){median(x, na.rm = TRUE)})

SquaredError_prev_final <- apply(AbsoluteError_prev, c(1,2), FUN = function(x) {
  if(x > 1) {x <- NA}
  else{x <- x}
  return(x)})
### Compute Mean of the relative bias for area-specific prevalence
RMSE_prev <- apply(SquaredError_prev_final, 1, FUN = function(x){sqrt(mean(x, na.rm = TRUE))})


### Print the summary of total Ntarget ###
print(paste("The case is ", case, "; the mean relative bias is ", round(Mean.rel.bias, digits = 3), 
            "; the RMSE is ", round(RMSE, digits = 3), 
            "; the MAE is ", round(MAE, digits = 3), sep = ""))

length(which(rel.bias > 100))

summary_Ntarget <- data.frame(case = case,
                              Mean.Rel.Bias_Ntarget = round(Mean.rel.bias, digits = 3),
                              MAE_Ntarget = round(MAE, digits = 3),
                              RMSE_Ntarget = round(RMSE, digits = 3),
                              n_infinite_est = length(which(rel.bias > 100))
)

### Print the summary of prevalence ###
summary_prev <- data.frame(Case = rep(case, length(names_all_areas)),
                           K = names_all_areas, 
                           Mean.Rel.Bias_prev = Mean.rel.bias_prev,
                           MAE_prev = MAE_prev,
                           RMSE_prev = RMSE_prev, 
                           freq_infinite_est_by_loc = freq_infinite_est_by_loc)


write.csv(summary_Ntarget, file = paste(path_to_LogLinear_output, "Summary_TotalNtarget.csv", sep = ""), row.names = FALSE)
write.csv(summary_prev, file = paste(path_to_LogLinear_output, "Summary_Prev.csv", sep = ""), row.names = FALSE)
