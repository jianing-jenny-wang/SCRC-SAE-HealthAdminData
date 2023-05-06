#####################################
######### Post Processing ###########
#####################################
library(ggplot2)
library(tidyr)


color1 <- rgb(255,192,203, max = 255, alpha = 90, names = "lt.pink")
color2 <- rgb(158,162,163,max = 255, alpha = 90, names = "lt.gray")

##################################
### Summary of True Prevalence ###
##################################
summary_truePrev <- function(true_samples){
 nsim <- dim(true_samples)[1] 
 K <- dim(TruePrev_all_simu)[2]
 summary_true_prev_across_simu <- apply(true_samples, 2, FUN = function(x){
   mean_x <- mean(x)
   sd_x <- sd(x)
   median_x <- median(x)
   lb95_x <- quantile(x, probs = 0.025)
   ub95_x <- quantile(x, probs = 0.975)
   min_x <- min(x)
   max_x <- max(x)
   res_vec <- c(mean_x, sd_x, median_x, lb95_x, ub95_x, min_x, max_x)
   return(res_vec)
 })
 summary_true_prev_across_simu <- t(summary_true_prev_across_simu)
 colnames(summary_true_prev_across_simu) <- c("Mean_TruePrev_k", "SD_TruePrev_k","Median_TruePrev_k",
                                              "LB95_TruePrev_k", "UB95_TruePrev_k", "Min_TruePrev_k", "Max_TruePrev_k")
 summary_true_prev_across_simu <- data.frame(K = seq(1:K), summary_true_prev_across_simu)
 summary_true_prev <- colMeans(summary_true_prev_across_simu[,c(2:ncol(summary_true_prev_across_simu))])
 names(summary_true_prev) <- c("Ave_Mean_TruePrev_k", "Ave_SD_TruePrev_k","Ave_Median_TruePrev_k",
                                  "Ave_LB95_TruePrev_k", "Ave_UB95_TruePrev_k", "Ave_Min_TruePrev_k", "Ave_Max_TruePrev_k")
 summary_true_prev <- round(summary_true_prev, digits = 4)
 return(summary_true_prev)
 }



#####################
### Target size N ###
#####################

## Draw trace plot of total N target ##
trace_postNtarget <- function(samples, title = "(C1) Number of population"){
  plot(jitter(samples), xlab = "iteration", ylab = "N", main = title, type = "l")
}

## Draw posterior disribution of total N target ##
plot_postNtarget <- function(samples, xlimit, TotalN_obs, TotalN_target, names_of_chain = "(C1)"){
  hist(samples, nclass = 50, col = "gray95", main = paste(names_of_chain, "Posterior of Total N target", sep = " "), xlab = "Posterior of total target N",
       las = 1, xlim = xlimit)
  abline(v = quantile(samples, prob = 0.25), col="black", lwd = 3) # Lower bound of samples
  abline(v = quantile(samples, prob = 0.75), col="black", lwd = 3) # Upper bound of samples
  abline(v = quantile(samples, prob = 0.025), col="gray50", lwd = 3) # Lower bound of samples
  abline(v = quantile(samples, prob = 0.975), col="gray50", lwd = 3) # Upper bound of samples
  abline(v = mean(samples), col="darkorange", lwd = 3) # Mean of samples
  abline(v = median(samples), col="red", lwd = 3) # Median of samples
  abline(v = TotalN_obs, col="black", lwd = 4) # Observed (Total)
  abline(v = TotalN_target, col = "blue", lwd = 4) # Truth N target (Total)
  legend("topright", c("95% CrI", "50% CrI", "Median N", "Mean N", "True Total N"), fill=c("gray75", "gray40", "red", "darkorange", "blue"), 
         xpd=TRUE, cex=0.6, bty='n')
}

## Compute summary statistics of post samples N ##
stat_fun <- function(x) {
  mean <- mean(x)
  med <- median(x)
  Qt2.5 <- quantile(x, prob = 0.025)
  Qt97.5 <- quantile(x, prob = 0.975)
  res <- c(mean, med, Qt2.5, Qt97.5)
  names(res) <- c("Mean", "Median", "Qt 2.5%", "Qt 97.5%")
  return(res)
}

#################################
### Pseudo prevalence: lambda ###
#################################

## Save samples from lambda when lambda was generated individually ##
save_lambda <- function(samples, names_in_chain) {
  index_lambda <- grep("^.*pesudo_logit_prev*.*",names_in_chain)
  chain_c_lambda <- samples[,index_lambda] # K columns
  return(chain_c_lambda)
}
plot_prev <- function(K, samples, truevalue, names_of_chain, names_of_k){
  # if((K %% 4) == 0){
  # par(mfrow=c(4,4))  
  # }
  # if((K %% 4) != 0){
  #   par(mfrow=c(4,4))  
  # }
  par(mfrow=c(4,4)) 
  for(k in 1:K){
    hist(samples[,k], col = "gray", main = paste(names_of_chain, " Post Prev. at ", names_of_k[k], sep = ""), xlab = "Prevalence", las = 1, cex=1)
    abline(v = quantile(samples[,k], prob = 0.25), col="black", lwd = 3) # Lower bound of samples
    abline(v = quantile(samples[,k], prob = 0.75), col="black", lwd = 3) # Upper bound of samples
    abline(v = quantile(samples[,k], prob = 0.025), col="gray50", lwd = 3) # Lower bound of samples
    abline(v = quantile(samples[,k], prob = 0.975), col="gray50", lwd = 3) # Upper bound of samples
    abline(v = mean(samples[,k]), col="darkorange", lwd = 3) # Mean of samples
    abline(v = median(samples[,k]), col="red", lwd = 3) # Median of samples
    abline(v = truevalue[k], col = "blue", lwd = 4) # Truth prevalence = N/P
  }
  # legend("topright", c("(95%)50% CrI", "Med-Prev","Mean-Prev", "True-Prev"), fill=c("gray75", "red", "darkorange", "blue"), 
  #        xpd=TRUE, cex=0.7, bty='n')
}

##################
### Prevalence ###
##################

## Draw posterior distribution of effective prevalence ##
## If location is stratified ##
plot_postprev_k <- function(K, samples, title = "(C1) Post Prevalence at location = "){
  if((K %% 2) == 0){
    par(mfrow=c(2,K/2))  
  }
  if((K %% 2) != 0){
    par(mfrow=c(2,(K+1)/2))  
  }
  for(k in 1:K){
    hist(samples, col = "gray85", main = paste(title, k, sep = ""), xlab = "Prevalence", las = 1)
    abline(v = quantile(samples, prob = 0.25), col="black", lwd = 3) # Lower bound of samples
    abline(v = quantile(samples, prob = 0.75), col="black", lwd = 3) # Upper bound of samples
    abline(v = quantile(samples, prob = 0.025), col="gray50", lwd = 3) # Lower bound of samples
    abline(v = quantile(samples, prob = 0.975), col="gray50", lwd = 3) # Upper bound of samples
    abline(v = mean(samples), col="darkorange", lwd = 3) # Mean of samples
    abline(v = median(samples), col="red", lwd = 3) # Median of samples
    abline(v = Prev[k], col = "blue", lwd = 4) # Truth prevalence = N/P
    legend("topright", c("(95%)50% CrI", "Median Prev","Mean Prev", "True Prev"), fill=c("gray75", "red", "darkorange", "blue"), 
           xpd=TRUE, cex=0.7, bty='n')
  }
  par(mfrow=c(1,1))
}
## If location is not stratified
plot_postprev_nosp <- function(samples, title = "(C1) Post Prevalence across all location"){
  hist(samples, col = "gray85", main = title, xlab = "Prevalence", las = 1)
  abline(v = quantile(samples, prob = 0.25), col="black", lwd = 3) # Lower bound of samples
  abline(v = quantile(samples, prob = 0.75), col="black", lwd = 3) # Upper bound of samples
  abline(v = quantile(samples, prob = 0.025), col="gray50", lwd = 3) # Lower bound of samples
  abline(v = quantile(samples, prob = 0.975), col="gray50", lwd = 3) # Upper bound of samples
  abline(v = median(samples), col="red", lwd = 3) # Median of samples
  abline(v = mean(samples), col="darkorange", lwd = 3) # Mean of samples
  abline(v = mean(Prev), col = "blue", lwd = 4) # Truth prevalence = N/P
  legend("topright", c("(95%)50% CrI", "Median Prev","Mean Prev", "True Prev"), fill=c("gray75", "red", "darkorange", "blue"), 
         xpd=TRUE, cex=0.7, bty='n')
}

## Compute summary statistics of post samples prevalence ##
stat_fun_prev <- function(samples, stat_fun = stat_fun){
  stat_effprev <- apply(samples, 2, FUN = stat_fun)
  stat_effprev <- round(stat_effprev, digits = 3)
  # colnames(stat_effprev) <- c("Prev")
  stat_effprev <- t(stat_effprev)
  return(stat_effprev)
}

##################################
### Prevalence terms - epsilon/beta ###
##################################

## Plot the comparison between posterior and prior distribution of epsilon
plot_eps <- function(K, post_eps, prior_eps, truevalue, names_of_chain){
  if((K %% 2) == 0){
    par(mfrow=c(2,K/2))  
  }
  if((K %% 2) != 0){
    par(mfrow=c(2,(K+1)/2))  
  }
  for(k in 1:K){
    post <- hist(post_eps[,k], plot = FALSE)
    prior <- hist(prior_eps[,k], plot = FALSE)
    xlim <- range(c(post$breaks, prior$breaks))
    ylim <- max(c(post$count, prior$count)) 
    plot(post, col = color2, xlim = xlim, ylim = c(0,ylim), main = paste(names_of_chain, " eps at k = ", k, sep = ""), xlab = expression(epsilon), las = 1, cex = 0.8)
    plot(prior, add = TRUE, col = color1)
    abline(v = truevalue[k], col = "blue", lwd = 3) #  pseudo logit of truth prevalence = N/P
    # legend("topright", expression("Post","Prior",paste("True ",epsilon, sep = " ")), fill=c(color2,color1,"blue"), 
    #        xpd=TRUE, cex=0.7, bty='n')
  }
}

plot_beta <- function(K, post_beta, prior_beta, truevalue, names_of_chain){
  if((K %% 2) == 0){
    par(mfrow=c(2,K/2))  
  }
  if((K %% 2) != 0){
    par(mfrow=c(2,(K+1)/2))  
  }
  for(k in 1:K){
    post <- hist(post_beta[,k], plot = FALSE)
    prior <- hist(prior_beta[,k], plot = FALSE)
    xlim <- range(c(post$breaks, prior$breaks))
    ylim <- max(c(post$count, prior$count)) 
    plot(post, col = color2, xlim = xlim, ylim = c(0,ylim), main = paste(names_of_chain, " beta at k = ", k, sep = ""), xlab = expression(beta), las = 1, cex = 0.8)
    plot(prior, add = TRUE, col = color1)
    abline(v = truevalue[k], col = "blue", lwd = 3) #  pseudo logit of truth prevalence = N/P
    # legend("topright", expression("Post","Prior",paste("True ",betailon, sep = " ")), fill=c(color2,color1,"blue"), 
    #        xpd=TRUE, cex=0.7, bty='n')
  }
}



###############################
### Detection Probabilities ###
###############################
## Compute summary statistics of post samples p.detect ##

stat_fun_p.detect <- function(samples, name_in_chain){
  # Keep only pdetect parameters' samples
  index_pdetect <- grep("^.*p.detect*.*",name_in_chain)
  chain_c_pdetect <- samples[,index_pdetect]
  # Remove the duplication because mcmc gives the pdetect for each individual
  mean_chain_c_pdetect <- colMeans(chain_c_pdetect) # vector
  mean_chain_c_pdetect_1 <- mean_chain_c_pdetect[1:aug_size]
  mean_chain_c_pdetect_2 <- mean_chain_c_pdetect[(aug_size + 1):(2*aug_size)]
  mean_chain_c_pdetect_3 <- mean_chain_c_pdetect[(2*aug_size + 1):length(mean_chain_c_pdetect)]
  # Combine list detection probabilities together and compute summary statistics
  mean_chain_c_pdetect_ls <- list(List1 = mean_chain_c_pdetect_1,
                                  List2 = mean_chain_c_pdetect_2,
                                  List3 = mean_chain_c_pdetect_3)
  stat_pdetect_ls <- lapply(mean_chain_c_pdetect_ls, FUN = stat_fun)
  stat_pdetect <- do.call(rbind, stat_pdetect_ls)
  stat_pdetect <- round(stat_pdetect, digits = 3)
  row.names(stat_pdetect) <- c("List1","List2","List3")
  return(stat_pdetect)
}

###########################
### logit list effect ###
###########################
## Obtain post samples of logit_list_effect ##
save_logitlisteffect <- function(samples, names_in_chain) {
  index_logitlisteff <- grep("^.*logit_list_effects*.*",names_in_chain)
  chain_c_logitlisteff <- samples[,index_logitlisteff] # K columns
  return(chain_c_logitlisteff)
}
## Compute summary statistics of post samples logit_list_effect ##
stat_fun_logitlisteffect <- function(J = J, chain = 1, samples, true_logitlisteff = dtListEffects_one_simu){
  mean <- colMeans(samples)
  median <- apply(samples,2,median)
  Qt2.5 <- apply(samples, 2, FUN = function(x) quantile(x, prob = 0.025))
  Qt97.5 <- apply(samples, 2, FUN = function(x) quantile(x, prob = 0.975))
  res <- data.frame(Chain = rep(chain, length(mean)), List = seq(1:J), Mean = mean, Median = median, Qt2.5 = Qt2.5, Qt97.5 = Qt97.5, true_logitlisteff = true_logitlisteff)
  names(res) <- c("Chain", "List", "Mean", "Median", "Qt 2.5%", "Qt 97.5%","True logit lists effect")
  row.names(res) <- NULL
  return(res)     
}

## Plot the comparison between posterior and prior distributions
plot_listeffect <- function(J, post_listeffect, prior_listeffect, truevalue, names_of_chain){
  if(J %% 2 == 0){
  par(mfrow=c(J/2,2))  
  }
  if(J %% 2 != 0){
  par(mfrow=c((J+1)/2,2))   
  }

  for(j in 1:J){
    post <- hist(post_listeffect[,j], plot = FALSE)
    prior <- hist(prior_listeffect[,j], plot = FALSE)
    xlim <- range(c(post$breaks, prior$breaks))
    ylim <- max(c(post$count, prior$count)) 
    plot(post, col = color2, xlim = xlim, ylim = c(0,ylim), main = paste(names_of_chain, " logit list effect j =", j, sep = ""), xlab = "logit list effect", las = 1)
    plot(prior, add = TRUE, col = color1)
    abline(v = truevalue[j], col = "blue", lwd = 3) # true value of logit list effect
    legend("topright", c("Post","Prior"), fill=c(color2,color1), 
           xpd=TRUE, cex=1.2, bty='n')
  }
  
}


#############################################
### Spatial variance parameter sigma2.phi ###
#############################################
plot_Sigma2.phi <- function(post_sigma2.phi, prior_sigma2.phi, truevalue, names_of_chain){
post <- hist(post_sigma2.phi, plot = FALSE)
if(!is.null(prior_sigma2.phi)){
  prior <- hist(prior_sigma2.phi, plot = FALSE) 
}
if(is.null(prior_sigma2.phi)){
  prior <- NULL
}
xlim <- if(is.null(prior_sigma2.phi)){
range(c(post$breaks))
}else{range(c(post$breaks, prior$breaks))}
ylim <- if(is.null(prior_sigma2.phi)){
  max(post$count)
  }else{max(c(post$count, prior$count))}
plot(post, col = color2, xlim = xlim, ylim = c(0,ylim), main = paste(names_of_chain, " Sigma2.phi", sep = ""), xlab = "Sigma2 for spRE", las = 1)
if(!is.null(prior_sigma2.phi)){
plot(prior, add = TRUE, col = color1)
}
abline(v = quantile(post_sigma2.phi, prob = 0.25), col="black", lwd = 3) # Lower bound of samples
abline(v = quantile(post_sigma2.phi, prob = 0.75), col="black", lwd = 3) # Upper bound of samples
abline(v = quantile(post_sigma2.phi, prob = 0.025), col="gray50", lwd = 3) # Lower bound of samples
abline(v = quantile(post_sigma2.phi, prob = 0.975), col="gray50", lwd = 3) # Upper bound of samples
abline(v = median(post_sigma2.phi), col="red", lwd = 3) # Median of samples
abline(v = mean(post_sigma2.phi), col="darkorange", lwd = 3) # Mean of samples
abline(v = truevalue, col = "blue", lwd = 3) # true value 
legend("topright", c("(95%)50% CrI", "Median Prev","Mean Prev", "True Prev"), fill=c("gray75", "red", "darkorange", "blue"), 
       xpd=TRUE, cex=0.7, bty='n')
}


################################################
### Spatial component residuals (phi) values ###
################################################
## Save samples from phi when phi was generated by location ##
save_resPhi <- function(samples, names_in_chain) {
  index_resPhi <- grep("^.*resPhi*.*",names_in_chain)
  chain_c_resPhi <- samples[,index_resPhi] # Dim = K columns
  return(chain_c_resPhi)
}

plot_resPhi <- function(K, nsamples, names_of_k, post, true_resPhi, names_of_chain, nchains){
  truedt_resPhi <- data.frame(true_res_Phi = true_resPhi,
                              names_of_k = names_of_k)
  post <- as.data.frame(post)
  if(nchains == 1){
    post$chains <- rep(nchains, nsamples)
  }
  if(nchains > 1){
    post$chains <- rep(seq(1:nchains), each = nsamples)  
  }
  colnames(post) <- c(paste("Area", seq(1:K), sep = ""), "chains")
  post <- post[,c("chains", paste("Area", seq(1:K), sep = ""))]
  post_resPhi <- gather(post, Areas, resPhi, Area1:paste("Area", K, sep = ""), factor_key=TRUE)
  loc_label <- rep(names_of_k, each = nchains*nsamples)
  postdt_resPhi <- data.frame(post_resPhi = post_resPhi,
                              names_of_k = loc_label)
  colnames(postdt_resPhi) <- c("chains", "Areas", "post_resPhi", "names_of_k")
  p <- ggplot(data = postdt_resPhi, aes(x = names_of_k, y = post_resPhi)) +
    geom_violin() +
    geom_point(data = truedt_resPhi, aes(x = names_of_k, y = true_res_Phi), color = "black") +
    theme_bw() +
    theme(
      legend.position="none",
      plot.title = element_text(size=13),
      axis.text.x = element_text(size=13, angle = 60, vjust =0.5),
      axis.text.y = element_text(size=13),
      axis.title = element_text(size = 13)) +
    xlab("Areas") +
    ylab(paste("Posterior of Spatial RE (Phi)", names_of_chain, sep = ""))
  return(p)
}


##################################
## Baselines overall prevalence ##
##################################
# 
# 
plot_baseline_prev <- function(post, true_baseline_prev, names_of_chain){
   post_baseline_prev <- as.vector(post)
   post <- hist(post_baseline_prev, plot = FALSE)
   xlim <- range(c(true_baseline_prev,post$breaks))
   ylim <- max(c(post$count))
   plot(post, col = color2, xlim = xlim, ylim = c(0,ylim), main = paste(names_of_chain, "Posterior of Baseline Prev", sep = ""), xlab = "Baseline Prev", las = 1)
   abline(v = quantile(post_baseline_prev, prob = 0.25), col="black", lwd = 3) # Lower bound of samples
   abline(v = quantile(post_baseline_prev, prob = 0.75), col="black", lwd = 3) # Upper bound of samples
   abline(v = quantile(post_baseline_prev, prob = 0.025), col="gray50", lwd = 3) # Lower bound of samples
   abline(v = quantile(post_baseline_prev, prob = 0.975), col="gray50", lwd = 3) # Upper bound of samples
   abline(v = median(post_baseline_prev), col="red", lwd = 3) # Median of samples
   abline(v = mean(post_baseline_prev), col="darkorange", lwd = 3) # Mean of samples
   abline(v = true_baseline_prev, col = "blue", lwd = 3) # true value 
   legend("topright", c("(95%)50% CrI", "Median Prev","Mean Prev", "True Prev"), fill=c("gray75", "red", "darkorange", "blue"), 
          xpd=TRUE, cex=0.7, bty='n')
 }
 

###########################
## Factor Loading Matrix ##
###########################
save_loading_mat <- function(samples, names_in_chain) {
  index_loading <- grep("^.*loading_mat_upper_tri*.*",names_in_chain)
  chain_c_loading <- samples[,index_loading] # K columns
  return(chain_c_loading)
}

## Plot the comparison between posterior and prior distributions
plot_loading_mat <- function(J, D, post_loading_mat, truevalue, names_of_chain){
  if(J*D %% 2 == 0){
    par(mfrow=c(J*D/2,2))  
  }
  if(J*D %% 2 != 0){
    par(mfrow=c((J*D+1)/2,2))   
  }
  # post_loading_mat has dim = nsamples * D*J
  # at column-wise, it will go through 1:J for D=1, then go through 1:J for D=2, etc
  # Therefore, when D = 1, post_loading_mat has J columns
  # when D > 1, post_loading_mat has J*D columns
  if(D == 1){
    for(j in 1:J){
      post <- hist(post_loading_mat[,j], plot = FALSE)
      xlim <- range(c(post$breaks))
      ylim <- max(c(post$count)) 
      plot(post, col = color2, xlim = xlim, ylim = c(0,ylim), main = paste(names_of_chain, "loading matrix row j =", j, sep = ""), xlab = "loading matrix", las = 1)
      abline(v = truevalue[j,], col = "blue", lwd = 3) # true value of loading matrix
      legend("topright", c("Post"), fill=c(color2), 
             xpd=TRUE, cex=1.2, bty='n')
    }
  }
  if(D > 1){
    truevalue_vec <- as.vector(truevalue)
    element_index_j <- seq(1:J)
    element_index_d <- seq(1:D)
    permutation_jd <- expand.grid(element_index_j,element_index_d)
    element_index <- paste("[J = ", permutation_jd[,1], ", D = ", permutation_jd[,2], "]", sep = "")
    
    for(jd in 1:(J*D)){
      post <- hist(post_loading_mat[,jd], plot = FALSE)
      xlim <- range(c(post$breaks))
      ylim <- max(c(post$count)) 
      plot(post, col = color2, xlim = xlim, ylim = c(0,ylim), 
           main = paste(names_of_chain, " loading matrix element ", element_index[jd], sep = ""), 
           xlab = "loading matrix element value", las = 1)
      abline(v = truevalue_vec[jd], col = "blue", lwd = 3) # true value of loading matrix
      legend("topright", c("Post"), fill=c(color2), 
             xpd=TRUE, cex=1.2, bty='n')
      print(jd)
    }
  }
  par(mfrow=c(1,1)) 
}





########################################
### Post process to compute measures ###
########################################
### Compute the Mean Relative Bias, RMSE, and MAE across nsim ###
Comp.Bias.RMSE.MAE.vec <- function(est,truevalue, nsim){
  if(length(est) == length(truevalue)){
    MSE <- mean(est - truevalue)^2  
    RMSE <- sqrt(MSE)
    MAE <- median(abs(est - truevalue))
    var.rel.bias <- sd(est/truevalue)
    mean.rel.bias <- mean(est/truevalue)
    Qt.2.5.rel.bias <- quantile(est/truevalue, probs = c(0.025, 0.975))[1]
    Qt.97.5.rel.bias <- quantile(est/truevalue, probs = c(0.025, 0.975))[2]
    T.Test.Diff.bias <- t.test(est-truevalue)
    T.Test.P.Value.Diff.bias <- T.Test.Diff.bias$p.value
    comp.bias <- c(MSE, RMSE, MAE, mean.rel.bias, var.rel.bias, 
                    Qt.2.5.rel.bias, Qt.97.5.rel.bias,
                    T.Test.P.Value.Diff.bias)
    names(comp.bias) <- c("MSE","RMSE","MAE", "Mean.Rel.Bias", "SD.Rel.Bias",
                          "Qt2.5.Rel.Bias", "Qt97.5.Rel.Bias", "T.Test.P.Value.Diff.Bias")
  }
  if(length(est) != length(truevalue)){
    message("estimated values must be the same length of number of true values")
  }
  return(comp.bias)
}


###############################################################
### Compute the overall capturability given independent p_j ###
###############################################################
theory_capturability_Mt_J3 <- function(N = 5000, p1, p2, p3){
  n_obs <- N*p1+N*p2+N*p3-N*(p1*p2+p1*p3+p2*p3+p1*p2*p3)
  prop <- n_obs/N
  obs1 <- sample(1:N,  N*p1, replace = FALSE)
  obs2 <- sample(1:N,  N*p2, replace = FALSE)
  obs3 <- sample(1:N,  N*p3, replace = FALSE)
  obs <- unique(c(obs1, obs2, obs3))
  empirical_prop <- length(obs)/N
  var_p <- var(c(p1,p2,p3))
  return(empirical_prop)
}
## Example 
# p1 <- 0.45
# p2 <- 0.145
# p3 <- 0.005
# N <- 5000
# obs1 <- sample(1:N,  N*p1, replace = FALSE)
# obs2 <- sample(1:N,  N*p2, replace = FALSE)
# obs3 <- sample(1:N,  N*p3, replace = FALSE)
# obs <- unique(c(obs1, obs2, obs3))
# # length(obs)/N
