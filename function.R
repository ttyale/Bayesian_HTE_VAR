######################################applying HTE Bayesian approaches to the KDPP data example######################################
#function clean the kdpp data set
cleandataset <- function(kdpp1){
  #Input:
  #kdpp1:original kdpp data set
  
  y_temp<-kdpp1$ddiabp
  
  cluster_id<-kdpp1$cluster
  mean(sort(unique(cluster_id))  == c(1:max(cluster_id)))
  n_c<-length(unique(cluster_id))
  n<-rep(NA,
         times = n_c)
  z<-rep(NA,
         times = n_c)
  v1<-rep(NA,
          times = n_c)
  v2<-rep(NA,
          times = n_c)
  v3<-rep(NA,
          times = n_c)
  v4<-rep(NA,
          times = n_c)
  for(j in 1:n_c){
    
    n[j]<-sum(cluster_id == sort(unique(cluster_id))[j]) #cluster size
    z[j]<-as.numeric(unique(kdpp1$arms0[cluster_id == sort(unique(cluster_id))[j]]) == "Intervention")
    v1[j]<-kdpp1$avschyr[cluster_id == sort(unique(cluster_id))[j]][1]
    v2[j]<-kdpp1$avage[cluster_id == sort(unique(cluster_id))[j]][1]
    v3[j]<-mean(kdpp1$religion[cluster_id == sort(unique(cluster_id))[j]] == "Christian")  #Categories:  Hindu, Muslim, Christian
    v4[j]<-mean(kdpp1$religion[cluster_id == sort(unique(cluster_id))[j]] == "Muslim")  #Categories:  Hindu, Muslim, Christian
    
  }
  
  v<-cbind(1,
           z,
           n,
           v1,
           v2,
           v3,
           v4)
  
  unique(kdpp1$religion)
  x<-cbind(1,
           scale(kdpp1$schoolyear),
           scale(kdpp1$age),
           as.numeric(kdpp1$religion == "Hindu"),   #Reference:  Christian
           as.numeric(kdpp1$religion == "Muslim"))  #Reference:  Christian 
  
  #Ignore Individual-Level Predictors
  x<-matrix(1,
            nrow = nrow(x),
            ncol = 1)
  
  y<-matrix(NA,
            nrow = n_c,
            ncol = max(n))
  for(j in 1:n_c){
    y[j, 1:n[j]]<-y_temp[cluster_id == sort(unique(cluster_id))[j]]
  }
  
  return(list(n_c,n,z,v,x,y))
}
#function producing posterior samples using 4 proposed models
model_mcmc <- function(n_c,n,z,v,x,y,seed=75425,niter,burnin,thin){
  #Input:
  #n_c: number of clusters
  #n: subject id
  #z: treatment assignment
  #v: vector of covariates
  #x: fixed effect intercept
  #y: observed outcome
  #seed: seeds for reproducibility
  #niter: total number of iterations
  #burnin: number of burn-in period
  #thin: thinning parameter
  
  require("rjags")
  
  #################################################
  #JAGS Code
  #################################################
  model_string<-"
model{
     
     #Likelihood
     for(i in 1:n_c){

        for(j in 1:n[i]){

           y[i,j] ~ dnorm(mu_y[i,j], 
                          sigma2_inv[i])

           mu_y[i,j] <- x[i,]%*%beta +
                        z[i]*eta + 
                        theta[i]
 
           }

        #Residual Variances
        sigma2_inv[i] <- 1.00/sigma2[i]
        sigma2[i] <- exp(2.00*log_sigma[i])
        log_sigma[i] ~ dnorm(mu_log_sigma[i],
                             phi2_inv)
        mu_log_sigma[i] <- v[i,]%*%gamma
        
        #Random Effects
        theta[i] ~ dnorm(0.00,
                         tau2_inv[i])
        tau2_inv[i] <- 1.00/tau2[i]
        tau2[i] <- tau[i]*tau[i]
        tau[i] <- omega0*(1.00 - z[i]) +
                  omega1*z[i]
        
        #ICCs 
        delta[i] <- tau2[i]/(tau2[i] + sigma2[i])

        }

     #Priors
     for(j in 1:n_col_x){
        beta[j] ~ dnorm(0.00,
                        0.0001)
        }

     eta ~ dnorm(0.00,
                 0.0001)

     for(j in 1:n_col_v){
        gamma[j] ~ dnorm(0.00,
                         0.0001)
        }

     omega0 ~ dunif(0.00, 
                    100.00)
     omega1 ~ dunif(0.00, 
                    100.00)
     omega_diff <- omega1 - 
                   omega0
     phi2_inv <- 1.00/(phi*phi)
     phi2 <- phi^2
     phi ~ dunif(0.00, 
                 100.00)

     }
"
  
  #######################################################
  #Model Organization
  #######################################################
  model_jags<-jags.model(textConnection(model_string),
                         data=list('n_c' = n_c,
                                   'n' = n,
                                   'y' = y,
                                   'x' = x,
                                   'n_col_x' = ncol(x),
                                   'z' = z,
                                   'v' = v,
                                   'n_col_v' = ncol(v)),
                         n.chains = 3,
                         inits = list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed))  #Number of Chains  
  
  ##############################################################
  #Posterior Sampling
  ##############################################################
  update(model_jags,
         n.iter = burnin)  #Burnin
  posterior_samples<-coda.samples(model_jags, 
                                  variable.names=c("beta",
                                                   "eta", 
                                                   "theta",
                                                   "sigma2",
                                                   "tau2",
                                                   "phi",
                                                   "phi2",
                                                   "gamma",
                                                   "omega_diff",
                                                   "delta"),
                                  thin = thin,        #Thinning
                                  n.iter = niter,
                                  seed=seed)  #Samples
  
  ###############################################################
  model_string1<-"
model{
     
     #Likelihood
     for(i in 1:n_c){

        for(j in 1:n[i]){

           y[i,j] ~ dnorm(mu_y[i,j], 
                          sigma2_inv)

           mu_y[i,j] <- beta0 +
                        eta*z[i] + 
                        theta[i]
 
           }

        #Random Effects
        theta[i] ~ dnorm(0.00,
                         tau2_inv)

        }
    
     #Priors
     beta0 ~ dnorm(0.00,
                   0.0001)
     eta ~ dnorm(0.00,
                 0.0001)

     sigma2_inv <- 1.00/(sigma*sigma)
     sigma ~ dunif(0.00, 
                   100.00)

     tau2_inv <- 1.00/(tau*tau)
     tau ~ dunif(0.00, 
                 100.00)

     #ICC
     delta <- (tau*tau)/(tau*tau + sigma*sigma)

     }
"
  
  ############################################################
  #2.2 JAGS code for invariable outcome variance but variable tau model##
  ############################################################
  ############################################################
  #2.2 JAGS code for invariable outcome variance but variable tau model##
  ############################################################
  model_string2<-"
model{
     
     #Likelihood
     for(i in 1:n_c){

        for(j in 1:n[i]){

           y[i,j] ~ dnorm(mu_y[i,j], 
                          sigma2_inv)

           mu_y[i,j] <- beta0 +
                        eta*z[i] + 
                        theta[i]
 
           }
  
        #Random Effects
        theta[i] ~ dnorm(0.00,
                         tau2_inv[i])
        tau2_inv[i] <- 1.00/tau2[i]
        tau2[i] <- tau[i]*tau[i]
        tau[i] <- omega0*(1.00 - z[i]) +
                  omega1*z[i]
        
        #ICCs 
        delta[i] <- tau2[i]/(tau2[i] + sigma2)

        }

     #Priors
     beta0 ~ dnorm(0.00,
                   0.0001)
     eta ~ dnorm(0.00,
                 0.0001)
    
     sigma2_inv <- 1.00/(sigma*sigma)
     sigma2 <- sigma*sigma
     sigma ~ dunif(0.00, 
                   100.00)
    
     omega0 ~ dunif(0.00, 
                    100.00)
     omega1 ~ dunif(0.00, 
                    100.00)
     omega_diff <- omega1 - 
                   omega0
    
     }
"
  ##################################################################################
  #2.3 JAGS code for variable out variance but invariant re variance model #########
  ##################################################################################
  model_string3<-"
model{
     
     #Likelihood
     for(i in 1:n_c){

        for(j in 1:n[i]){

           y[i,j] ~ dnorm(mu_y[i,j], 
                          sigma2_inv[i])

           mu_y[i,j] <- beta0 +
                        eta*z[i] + 
                        theta[i]
 
           }

        #Residual Variances
        sigma2_inv[i] <- 1.00/sigma2[i]
        sigma2[i] <- exp(2.00*log_sigma[i])
        log_sigma[i] ~ dnorm(mu_log_sigma[i],
                             phi2_inv)
        mu_log_sigma[i] <- v[i,]%*%gamma
        
        #Random Effects
        theta[i] ~ dnorm(0.00,
                         tau2_inv)
       
        #ICCs 
        delta[i] <- tau*tau/(tau*tau + sigma2[i])

        }

     #Priors
     beta0 ~ dnorm(0.00,
                   0.0001)
     eta ~ dnorm(0.00,
                 0.0001)
     mu ~ dnorm(0.00,
                0.0001)
     for(j in 1:n_col_v){
        gamma[j] ~ dnorm(0.00,
                         0.0001)
        }
     tau2_inv <- 1.00/(tau*tau)
     tau ~ dunif(0.00, 
                 100.00)
        
     phi2_inv <- 1.00/(phi*phi)
     phi2 <- phi^2
     phi ~ dunif(0.00, 
                 100.00)

     }
"
  
  
  #######################################################
  #Model 3 Variable ICC Fixed RE Var
  #######################################################
  model_jags3<-jags.model(textConnection(model_string3),
                          data=list('n_c' = n_c,
                                    'n' = n,
                                    'y' = y,
                                    'z' = z,
                                    'v' = v,
                                    'n_col_v' = ncol(v)),
                          n.chains = 3,
                          inits = list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed))  #Number of Chains  
  ##############################################################
  #Posterior Sampling Variable ICC
  ##############################################################
  update(model_jags3,
         n.iter = burnin)  #Burnin
  
  posterior_sample3<-coda.samples(model_jags3, 
                                  variable.names=c("beta0",
                                                   "eta", 
                                                   "theta",
                                                   "sigma2",
                                                   "tau",
                                                   "phi",
                                                   "phi2",
                                                   "gamma",
                                                   "delta"),
                                  thin = thin,        #Thinning
                                  n.iter = niter)  #Samples
  
  
  #######################################################
  #Model 2 Fixed residual variance Variable RE Var
  #######################################################
  model_jags2<-jags.model(textConnection(model_string2),
                          data=list('n_c' = n_c,
                                    'n' = n,
                                    'y' = y,
                                    'z' = z),
                          n.chains = 3,
                          inits = list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed))  #Number of Chains  
  
  ##############################################################
  #Posterior Sampling Variable ICC
  ##############################################################
  update(model_jags2,
         n.iter = burnin)  #Burnin
  
  posterior_sample2<-coda.samples(model_jags2, 
                                  variable.names=c("beta0",
                                                   "eta", 
                                                   "theta",
                                                   "sigma",
                                                   "tau2",
                                                   "omega_diff",
                                                   "delta"),
                                  thin = thin,      #Thinning
                                  n.iter = niter)  #Samples
  
  
  
  
  
  #######################################################
  #Model Organization Traditional
  #######################################################
  model_jags1<-jags.model(textConnection(model_string1),
                          data=list('n_c' = n_c,
                                    'n' = n,
                                    'y' = y,
                                    'z' = z),
                          n.chains = 3,
                          inits = list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed))  #Number of Chains  
  
  ############################################################
  #Posterior Sampling traditional
  ############################################################
  update(model_jags1,
         n.iter = burnin)  #Burnin
  
  posterior_sample1<-coda.samples(model_jags1, 
                                  variable.names=c("beta0",
                                                   "eta",
                                                   "sigma", 
                                                   "theta",
                                                   "tau", 
                                                   "delta"),
                                  thin = thin,        #Thinning
                                  n.iter = niter)  #Samples
  
  return(list(posterior_samples,posterior_sample3,posterior_sample2,posterior_sample1))
}

#function calculating WAIC for 4 proposed models using posterior samples
WAIC_Cal <- function(posterior,n_c,n,z,v,x,y){
  #Input:
  #posterior: posterior samples 
  #n_c: number of clusters
  #n: subject id
  #z: treatment assignment
  #v: vector of covariates
  #x: fixed effect intercept
  #y: observed outcome
  
  posterior_samples <- posterior[[1]]
  posterior_sample3 <- posterior[[2]]
  posterior_sample2 <- posterior[[3]]
  posterior_sample1 <- posterior[[4]]
  
  if(sum((substr(colnames(posterior_samples[[1]]),1,4) == "beta")) > 1){
    
    beta<-posterior_samples[[1]][,(substr(colnames(posterior_samples[[1]]),1,4) == "beta")]
    for(j in 2:length(posterior_samples)){
      beta<-rbind(beta,
                  posterior_samples[[2]][,(substr(colnames(posterior_samples[[2]]),1,4) == "beta")])
    }
    
  }
  
  if(sum((substr(colnames(posterior_samples[[1]]),1,4) == "beta")) == 1){
    
    beta<-posterior_samples[[1]][,(substr(colnames(posterior_samples[[1]]),1,4) == "beta")]
    for(j in 2:length(posterior_samples)){
      beta<-c(beta,
              posterior_samples[[2]][,(substr(colnames(posterior_samples[[2]]),1,4) == "beta")])
    }
    beta<-matrix(beta,
                 ncol = 1)
    
  }
  
  eta<-posterior_samples[[1]][,(substr(colnames(posterior_samples[[1]]),1,3) == "eta")]
  for(j in 2:length(posterior_samples)){
    eta<-c(eta,
           posterior_samples[[2]][,(substr(colnames(posterior_samples[[2]]),1,3) == "eta")])
  }
  
  theta<-posterior_samples[[1]][,(substr(colnames(posterior_samples[[1]]),1,5) == "theta")]
  for(j in 2:length(posterior_samples)){
    theta<-rbind(theta,
                 posterior_samples[[2]][,(substr(colnames(posterior_samples[[2]]),1,5) == "theta")])
  }
  
  tau2<-posterior_samples[[1]][,(substr(colnames(posterior_samples[[1]]),1,4) == "tau2")]
  for(j in 2:length(posterior_samples)){
    tau2<-rbind(tau2,
                posterior_samples[[2]][,(substr(colnames(posterior_samples[[2]]),1,4) == "tau2")])
  }
  
  sigma2<-posterior_samples[[1]][,(substr(colnames(posterior_samples[[1]]),1,6) == "sigma2")]
  for(j in 2:length(posterior_samples)){
    sigma2<-rbind(sigma2,
                  posterior_samples[[2]][,(substr(colnames(posterior_samples[[2]]),1,6) == "sigma2")])
  }
  
  #####################################################################################################
  #WAIC Calculation
  #####################################################################################################
  piece<-matrix(0, 
                nrow = nrow(beta), 
                ncol = (nrow(y)*ncol(y)))
  for(j in 1:nrow(beta)){
    
    mean<-z*eta[j] +
      theta[j,]
    
    sd<-sqrt(sigma2[j,])
    
    mu_y<-rep(NA,
              times = (nrow(y)*ncol(y)))
    sd_y<-rep(NA,
              times = (nrow(y)*ncol(y)))
    
    counter<-0
    for(k in 1:nrow(y)){
      
      #Following Code Gives and Error But It Is OK
      if(sum((substr(colnames(posterior_samples[[1]]),1,4) == "beta")) > 1){
        mu_y[(1 + (k-1)*ncol(y)):(k*ncol(y))]<-mean[k] +
          x[c((1 + counter):(n[k] + counter)),]%*%beta[j,]
      }
      
      if(sum((substr(colnames(posterior_samples[[1]]),1,4) == "beta")) == 1){
        mu_y[(1 + (k-1)*ncol(y)):(k*ncol(y))]<-mean[k] +
          x[c((1 + counter):(n[k] + counter))]*beta[j]
      }
      
      sd_y[(1 + (k-1)*ncol(y)):(k*ncol(y))]<-sd[k]
      counter<-counter +
        n[k]
      
    }
    
    piece[j,]<-dnorm(x = c(t(y)),
                     mean = mu_y, 
                     sd = sd_y, 
                     log = TRUE)
    
   # print(j/nrow(beta))
    
  }
  
  llpd<-sum(log(colMeans(exp(piece), na.rm = TRUE)), na.rm = TRUE)
  
  PWAIC_1<-2*sum(log(colMeans(exp(piece), na.rm = TRUE)) - colMeans(piece, na.rm = TRUE), na.rm = TRUE)
  WAIC_1<- -2*(llpd - PWAIC_1)
  WAIC_1
  
  temp<-rep(0, 
            times = (nrow(y)*ncol(y)))
  for(j in 1:(nrow(y)*ncol(y))){
    temp[j]<-var(piece[,j], na.rm = TRUE)
  }
  PWAIC_2<-sum(temp, na.rm = TRUE)
  WAIC_2<- -2*(llpd - PWAIC_2)
  
  
  #WAIC and WAIC2
  #Model 3
  beta0<-posterior_sample3[[1]][,(substr(colnames(posterior_sample3[[1]]),1,5) == "beta0")]
  eta<-posterior_sample3[[1]][,(substr(colnames(posterior_sample3[[1]]),1,3) == "eta")]
  theta<-posterior_sample3[[1]][,(substr(colnames(posterior_sample3[[1]]),1,5) == "theta")]
  sigma2<-posterior_sample3[[1]][,(substr(colnames(posterior_sample3[[1]]),1,6) == "sigma2")]
  #WAIC Calculation
  piece<-matrix(0, 
                nrow = length(beta0), 
                ncol = (nrow(y)*ncol(y)))
  
  for(j in 1:length(beta0)){
    
    mean<-beta0[j] +
      eta[j]*z +
      theta[j,] 
    sd<-sqrt(sigma2[j,])
    
    mu_y<-rep(NA,
              times = (nrow(y)*ncol(y)))
    sd_y<-rep(NA,
              times = (nrow(y)*ncol(y)))
    
    counter<-0
    for(k in 1:nrow(y)){
      
      mu_y[(1 + (k-1)*ncol(y)):(k*ncol(y))]<-mean[k]
      sd_y[(1 + (k-1)*ncol(y)):(k*ncol(y))]<-sd[k]
      
    }
    
    piece[j,]<-dnorm(x = c(t(y)),
                     mean = mu_y, 
                     sd = sd_y, 
                     log = TRUE)
    
    #print(j/length(beta0))
  }
  
  llpd<-sum(log(colMeans(exp(piece), na.rm = TRUE)), na.rm = TRUE)
  
  PWAIC_1<-2*sum(log(colMeans(exp(piece), na.rm = TRUE)) - colMeans(piece, na.rm = TRUE), na.rm = TRUE)
  WAIC31<- -2*(llpd - PWAIC_1)
  
  temp<-rep(0, 
            times = (nrow(y)*ncol(y)))
  for(j in 1:(nrow(y)*ncol(y))){
    temp[j]<-var(piece[,j], na.rm = TRUE)
  }
  PWAIC_2<-sum(temp, na.rm = TRUE)
  WAIC32<- -2*(llpd - PWAIC_2) 
  
  
  
  #Model 2
  beta0<-posterior_sample2[[1]][,(substr(colnames(posterior_sample2[[1]]),1,5) == "beta0")]
  eta<-posterior_sample2[[1]][,(substr(colnames(posterior_sample2[[1]]),1,3) == "eta")]
  theta<-posterior_sample2[[1]][,(substr(colnames(posterior_sample2[[1]]),1,5) == "theta")]
  sigma<-posterior_sample2[[1]][,(substr(colnames(posterior_sample2[[1]]),1,5) == "sigma")]
  #WAIC Calculation
  piece<-matrix(0, 
                nrow = length(beta0), 
                ncol = (nrow(y)*ncol(y)))
  
  for(j in 1:length(beta0)){
    
    mean<-beta0[j] +
      eta[j]*z +
      theta[j,] 
    
    mu_y<-rep(NA,
              times = (nrow(y)*ncol(y)))
    
    counter<-0
    for(k in 1:nrow(y)){
      mu_y[(1 + (k-1)*ncol(y)):(k*ncol(y))]<-mean[k]
    }
    
    piece[j,]<-dnorm(x = c(t(y)),
                     mean = mu_y, 
                     sd = sigma[j], 
                     log = TRUE)
  }
  
  
  llpd<-sum(log(colMeans(exp(piece), na.rm = TRUE)), na.rm = TRUE)
  
  PWAIC_1<-2*sum(log(colMeans(exp(piece), na.rm = TRUE)) - colMeans(piece, na.rm = TRUE), na.rm = TRUE)
  WAIC21<- -2*(llpd - PWAIC_1)
  
  temp<-rep(0, 
            times = (nrow(y)*ncol(y)))
  for(j in 1:(nrow(y)*ncol(y))){
    temp[j]<-var(piece[,j], na.rm = TRUE)
  }
  PWAIC_2<-sum(temp, na.rm = TRUE)
  WAIC22<- -2*(llpd - PWAIC_2) 
  
  
  #Model 1
  beta0<-posterior_sample1[[1]][,(substr(colnames(posterior_sample1[[1]]),1,5) == "beta0")]
  eta<-posterior_sample1[[1]][,(substr(colnames(posterior_sample1[[1]]),1,3) == "eta")]
  theta<-posterior_sample1[[1]][,(substr(colnames(posterior_sample1[[1]]),1,5) == "theta")]
  sigma<-posterior_sample1[[1]][,(substr(colnames(posterior_sample1[[1]]),1,5) == "sigma")]
  #################
  #WAIC Calculation
  piece<-matrix(0, 
                nrow = length(beta0), 
                ncol = (nrow(y)*ncol(y)))
  for(j in 1:length(beta0)){
    
    mean<-beta0[j] +
      eta[j]*z +
      theta[j,] 
    
    mu_y<-rep(NA,
              times = (nrow(y)*ncol(y)))
    counter<-0
    for(k in 1:nrow(y)){
      mu_y[(1 + (k-1)*ncol(y)):(k*ncol(y))]<-mean[k]
    }
    
    piece[j,]<-dnorm(x = c(t(y)),
                     mean = mu_y, 
                     sd = sigma[j], 
                     log = TRUE)
    
    #print(j/length(beta0))
  }
  
  llpd<-sum(log(colMeans(exp(piece), na.rm = TRUE)), na.rm = TRUE)
  
  PWAIC_1<-2*sum(log(colMeans(exp(piece), na.rm = TRUE)) - colMeans(piece, na.rm = TRUE), na.rm = TRUE)
  WAIC11<- -2*(llpd - PWAIC_1)
  
  
  temp<-rep(0, 
            times = (nrow(y)*ncol(y)))
  for(j in 1:(nrow(y)*ncol(y))){
    temp[j]<-var(piece[,j], na.rm = TRUE)
  }
  PWAIC_2<-sum(temp, na.rm = TRUE)
  WAIC12<- -2*(llpd - PWAIC_2)

  
  return(c(WAIC_2,WAIC32,WAIC22,WAIC12))
}
#function calculating means, medians and credible intervals for the intervention effect and regression coefficients
summary_statistics <- function(posterior,n_c,n,z,v,x,y){
  #Input:
  #posterior: posterior samples 
  #n_c: number of clusters
  #n: subject id
  #z: treatment assignment
  #v: vector of covariates
  #x: fixed effect intercept
  #y: observed outcome
  
  posterior_samples <- posterior[[1]]
  posterior_sample3 <- posterior[[2]]
  posterior_sample2 <- posterior[[3]]
  posterior_sample1 <- posterior[[4]]
  
  WAIC <- WAIC_Cal(posterior,n_c,n,z,v,x,y)
  
  #ATE estimates by model 1
  atemean1 <- mean(posterior_sample1[[1]][,(substr(colnames(posterior_sample1[[1]]),1,3) == "eta")])
  atemedian1 <- median(posterior_sample1[[1]][,(substr(colnames(posterior_sample1[[1]]),1,3) == "eta")])
  #credible interval
  atelq1 <- summary(posterior_sample1[[1]])$quantiles[(substr(rownames(summary(posterior_sample1[[1]])$quantiles),1,3) == "eta"),][1]
  ateuq1 <- summary(posterior_sample1[[1]])$quantiles[(substr(rownames(summary(posterior_sample1[[1]])$quantiles),1,3) == "eta"),][5]
  
  #ATE estimates by model 2
  atemean2 <- mean(posterior_sample2[[1]][,(substr(colnames(posterior_sample2[[1]]),1,3) == "eta")])
  atemedian2 <- median(posterior_sample2[[1]][,(substr(colnames(posterior_sample2[[1]]),1,3) == "eta")])
  #credible interval
  atelq2 <- summary(posterior_sample2[[1]])$quantiles[(substr(rownames(summary(posterior_sample2[[1]])$quantiles),1,3) == "eta"),][1]
  ateuq2 <- summary(posterior_sample2[[1]])$quantiles[(substr(rownames(summary(posterior_sample2[[1]])$quantiles),1,3) == "eta"),][5]
  
  #ATE estimates by model 3
  atemean3 <- mean(posterior_sample3[[1]][,(substr(colnames(posterior_sample3[[1]]),1,3) == "eta")])
  atemedian3 <- atemedian2 <- median(posterior_sample3[[1]][,(substr(colnames(posterior_sample3[[1]]),1,3) == "eta")])
  #credible interval
  atelq3 <- summary(posterior_sample3[[1]])$quantiles[(substr(rownames(summary(posterior_sample3[[1]])$quantiles),1,3) == "eta"),][1]
  ateuq3 <- summary(posterior_sample3[[1]])$quantiles[(substr(rownames(summary(posterior_sample3[[1]])$quantiles),1,3) == "eta"),][5]
  
  #ATE estimates by model 4
  atemean4 <- mean(posterior_samples[[1]][,(substr(colnames(posterior_samples[[1]]),1,3) == "eta")])
  atemedian4 <- median(posterior_samples[[1]][,(substr(colnames(posterior_samples[[1]]),1,3) == "eta")])
  #credible interval
  atelq4 <- summary(posterior_samples[[1]])$quantiles[(substr(rownames(summary(posterior_samples[[1]])$quantiles),1,3) == "eta"),][1]
  ateuq4 <- summary(posterior_samples[[1]])$quantiles[(substr(rownames(summary(posterior_samples[[1]])$quantiles),1,3) == "eta"),][5]
  
  
  #coefficients of model 1
  coefmean1 <- rep(NA,7)
  coefmedian1 <- rep(NA,7)
  coeflq1 <- rep(NA,7)
  coefuq1 <- rep(NA,7)
  
  #coefficients of model 2
  coefmean2 <- rep(NA,7)
  coefmedian2 <- rep(NA,7)
  coeflq2 <- rep(NA,7)
  coefuq2 <- rep(NA,7)
  
  #coefficients of model 3
  coefmean3 <- round(colMeans(posterior_sample3[[1]][,(substr(colnames(posterior_sample3[[1]]),1,5) == "gamma")]),2)
  coefmedian3 <- round(apply(posterior_sample3[[1]][,(substr(colnames(posterior_sample3[[1]]),1,5) == "gamma")],2,median),2)
  #credible interval
  coefci3 <- round(summary(posterior_sample3[[1]])$quantiles[(substr(rownames(summary(posterior_sample3[[1]])$quantiles),1,5) == "gamma"),],2)
  coeflq3 <- coefci3[,1]
  coefuq3 <- coefci3[,5]
  
  #coefficients of model 4
  coefmean4 <- round(colMeans(posterior_samples[[1]][,(substr(colnames(posterior_samples[[1]]),1,5) == "gamma")]),2)
  coefmedian4 <- round(apply(posterior_samples[[1]][,(substr(colnames(posterior_samples[[1]]),1,5) == "gamma")],2,median),2)
  #credible interval
  coefci4 <- round(summary(posterior_samples[[1]])$quantiles[(substr(rownames(summary(posterior_samples[[1]])$quantiles),1,5) == "gamma"),],2)
  coeflq4 <- coefci4[,1]
  coefuq4 <- coefci4[,5]
  
  r4 <- c(atemean4,atemedian4,atelq4,ateuq4,coefmean4,coefmedian4,coeflq4,coefuq4,WAIC[1])
  r3 <- c(atemean3,atemedian3,atelq3,ateuq3,coefmean3,coefmedian3,coeflq3,coefuq3,WAIC[2])
  r2 <- c(atemean2,atemedian2,atelq2,ateuq2,coefmean2,coefmedian2,coeflq2,coefuq2,WAIC[3])
  r1 <- c(atemean1,atemedian1,atelq1,ateuq1,coefmean1,coefmedian1,coeflq1,coefuq1,WAIC[4])
  
  df <- rbind(r1,r2,r3,r4)
  colnames(df) <- c("ate_mean","ate_median","ate_cil","ate_ciu",
                    paste("coefficients_mean", 1:7, sep=""),
                    paste("coefficients_median", 1:7, sep=""),
                    paste("coefficients_cil", 1:7, sep=""),
                    paste("coefficients_ciu", 1:7, sep=""),
                    "WAIC")
  
  return(df)
}
