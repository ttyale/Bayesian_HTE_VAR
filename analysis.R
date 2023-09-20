######################################applying HTE Bayesian approaches to the KDPP data example######################################
source("./function.R")
kdpp1 <- read.csv("kdpp_cleaned.csv")
clearned_kdpp <- cleandataset(kdpp1)
n_c <- clearned_kdpp[[1]]
n <- clearned_kdpp[[2]]
z <- clearned_kdpp[[3]]
v <- clearned_kdpp[[4]]
x <- clearned_kdpp[[5]]
y <- clearned_kdpp[[6]]
burnin <- 30000
niter <- 300000
thin <- 1
posterior <- model_mcmc(n_c,n,z,v,x,y,seed=75425,niter,burnin,thin)
result <- summary_statistics(posterior,n_c,n,z,v,x,y)
result <- as.data.frame(result)
write.csv(result,file="table4.csv")
