##==============================================================================
## time series sampling
##
## aewsma@rit.edu
##==============================================================================

# par_time_center and par_time_stdev are the guys

f_var_red <- matrix( c(0,0.82,0.80,0.84,0.20,0.55,0.55,0.60,0,0.75,0.65,0.55), 12, 1)
rownames(f_var_red) <- colnames(par_time_stdev)

# pick degrees of freedom for each parameter so the variance matches the
# reduced variance (improved model failure rate) experiments of Royer et al 2014
df <- matrix(rep(NA,12),12,1)
rownames(df) <- colnames(par_time_stdev)
for (pp in 1:12) {
  df[pp] <- n_time + 3 + ceil(2*max(par_time_stdev[,pp]^2)/(1-f_var_red[pp]))
}

# STOP -- think about how large this matrix might be for a large number of Latin
# hypercube samples. Throw out anythign with %outbound > .50 immediately?
n_sample <- 1000
samples <- array(dim=c(n_time,n_time,n_sample))
for (pp in 1:12) {
  samples <- rInvWishart(n_sample, df[pp], (df[pp]-(n_time+1))*diag(par_time_stdev^2))
}



if(FALSE){
# how it's done, accounting for the degrees of freedom and keeping the mean
# of the sampled matrices equal to the diagonal covariance
p <- 4
df1 <- p+2; A1 <- rInvWishart(1000, df1, (df1-(p+1))*diag(rep(16,p)))
df2 <- 2*p; A2 <- rInvWishart(1000, df2, (df2-(p+1))*diag(rep(16,p)))
df3 <- 6*p; A3 <- rInvWishart(1000, df3, (df3-(p+1))*diag(rep(16,p)))
print(mean(A1[1,1,])); hist(A1[1,1,], breaks=50, xlim=c(0,200))
print(mean(A2[1,1,])); hist(A2[1,1,])
print(mean(A3[1,1,])); hist(A3[1,1,])
}


##==============================================================================
## End
##==============================================================================
