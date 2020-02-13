##==============================================================================
## lhs.R
##
## todo...
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================
## Copyright 2019 Tony Wong
## This file is part of GEOCARB-calibration.
## GEOCARB-calibration is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## GEOCARB-calibration is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with GEOCARB-calibration.  If not, see <http://www.gnu.org/licenses/>.
##==============================================================================

## Clear workspace
rm(list=ls())

# needed packages
library(foreach)
library(doParallel)
library(lhs)
library(Hmisc)
library(CholWishart)
library(abind)
library(MASS)

## Set testing number of samples and file name appendix here
## if there aren't enough samples on the MCMC output file, will break.
n_sample <- 1e5 # note this will be doubled if using LHS precalibration
appen <- 'test'

n_node <- 6 # use parallel evaluation of ensembles in Sobol' integration?

param_choice <- 'all'   # Calibrate all 68 parameters? ("all") or only the 6 from Park and Royer 2011 ("PR2011")
data_choice <- 'F2017'    # Which data set?  PR2011 = Park and Royer (2011), or F2017 = Foster et al (2017)
fSR_choice <- 'DT2019'     # Which fSR time series? ("PR2011", "LENTON", "DT2019")
plot.dir <- '../figures/'
prcout_threshold <- 0.5  # only save simulations with percent_outbound < this

# calibration parameters and input data
filename.data <- '../input_data/CO2_Proxy_Foster2017_calib_SN-co2_25Sep2018.csv'
filename.calibinput <- '../input_data/GEOCARB_input_summaries_calib_all.csv'

if(Sys.info()['user']=='tony') {
  # Tony's local machine (if you aren't me, you almost certainly need to change this...)
  machine <- 'local'
  setwd('/Users/tony/codes/GEOCARB-calibration/R')
} else {
  # assume on a cluster of some kind...
  machine <- 'remote'
  setwd('~/work/codes/GEOCARB-calibration/R')
}

source("model_setup.R")
n_parameters <- length(parnames_calib)

## preliminary simulation to get the length
model_out <- model_forMCMC( par_calib=par_calib0,
                            par_time=par_time_center,
                            par_fixed=par_fixed0,
                            parnames_calib=parnames_calib,
                            parnames_fixed=parnames_fixed,
                            parnames_time=parnames_time,
                            age=age,
                            ageN=ageN,
                            ind_const_calib=ind_const_calib,
                            ind_time_calib=ind_time_calib,
                            ind_const_fixed=ind_const_fixed,
                            ind_time_fixed=ind_time_fixed,
                            ind_expected_time=ind_expected_time,
                            ind_expected_const=ind_expected_const,
                            iteration_threshold=iteration_threshold)
n_time <- length(model_out[,"co2"])
time <- model_out[,1]
source('percent_outbound.R')
source('constraints.R')

##==============================================================================



##==============================================================================
## Below here will be in a separate sampling routine
## vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv


n_sample <- 1000
n_sample_per_chunk <- 400 # maximum number of time series samples to consider at once
n_node <- 1

if (n_sample <= n_sample_per_chunk) {
  # business as usual

  set.seed(2020)
  parameters_lhs <- randomLHS(n_sample, n_parameters)
  par_calib <- parameters_lhs  # initialize

#TODO
print("here for some reason?")
#do simulations... only return the results with at most 50% %outbound

} else {
  # divide the samples into chunks  <-- TODO: could try to distribute the load more evenly if parallelized
  n_chunk <- ceiling(n_sample/n_sample_per_chunk)
  n_sample_this_chunk <- rep(n_sample_per_chunk, n_chunk)
  n_sample_this_chunk[n_chunk] <- n_sample - (n_chunk-1)*n_sample_per_chunk

  # sample the constant parameters from one large LHS and divide
  set.seed(2020)
  parameters_lhs <- randomLHS(n_sample, n_parameters)
  # constant parameter sampling
  n_const_calib <- length(ind_const_calib)
  par_calib_tmp <- parameters_lhs  # initialize
  colnames(par_calib_tmp) <- parnames_calib
  for (i in 1:n_const_calib) {
    row_num <- match(parnames_calib[i],input$parameter)
    if(input[row_num, 'distribution_type']=='gaussian') {
      par_calib_tmp[,i] <- qnorm(p=parameters_lhs[,ind_const_calib[i]], mean=input[row_num,"mean"], sd=(0.5*input[row_num,"two_sigma"]))
    } else if(input[row_num, 'distribution_type']=='lognormal') {
      par_calib_tmp[,i] <- qlnorm(p=parameters_lhs[,ind_const_calib[i]], meanlog=log(input[row_num,"mean"]), sdlog=log(0.5*input[row_num,"two_sigma"]))
    } else {
      print('ERROR - unknown prior distribution type')
    }
  }

  # divide the scalar parameters across the chunks
  par_calib <- vector('list', n_chunk)
  for (cc in 1:(n_chunk-1)) {
    par_calib[[cc]] <- par_calib_tmp[((cc-1)*n_sample_per_chunk+1):(cc*n_sample_per_chunk),]
  }
  par_calib[[n_chunk]] <- par_calib_tmp[(cc*n_sample_per_chunk+1):n_sample,]
  rm(list=c("par_calib_tmp"))

  # set up arrays to hold the parameters and covariances with low enough %outbound
  par_calib_save <- NULL
  par_time_save <- NULL
  par_covar_save <- NULL

  for (cc in 1:n_chunk) {
    # do simulations... only return the results with at most 50% %outbound

    # within simulation loop, do the time series sampling
    covariance_samples <- array(dim=c(n_time,n_time,n_sample_this_chunk[cc],n_parameters_time))
    time_series_samples <- array(dim=c(n_time,n_sample_this_chunk[cc],n_parameters_time))
    for (pp in 1:n_parameters_time) {
      # these degrees of freedom give variances in line with the reduced variances
      # of Royer et al (2014, AJS), and set the sampled mean covariance on the
      # diagonal matrix of variances (df-(n_time+1))
      covariance_samples[,,,pp] <- rInvWishart(n_sample_this_chunk[cc], df[pp], (df[pp]-(n_time+1))*diag(par_time_stdev[,pp]^2))
      # now draw the actual time series
      for (ii in 1:n_sample_this_chunk[cc]) {
        time_series_samples[,ii,pp] <- mvrnorm(n=1, mu=par_time_center[,pp], Sigma=covariance_samples[,,ii,pp])
      }
    }

    # run the simulations
    model_co2_this_chunk <- model_temp_this_chunk <- mat.or.vec(nr=n_time, nc=n_sample_this_chunk[cc])
    tbeg <- proc.time()
    for (ii in 1:n_sample_this_chunk[cc]) {
      model_out <- model_forMCMC(par_calib=par_calib[[cc]][ii,],      ## <--------------------- UPDATE HERE NOW!!  combien this loop with the one above, and check for %outbound after each
                                 par_time=time_series_samples[,ii,],
                                 par_fixed=par_fixed0,
                                 parnames_calib=parnames_calib,
                                 parnames_fixed=parnames_fixed,
                                 parnames_time=parnames_time,
                                 age=age,
                                 ageN=ageN,
                                 ind_const_calib=ind_const_calib,
                                 ind_time_calib=ind_time_calib,
                                 ind_const_fixed=ind_const_fixed,
                                 ind_time_fixed=ind_time_fixed,
                                 ind_expected_time=ind_expected_time,
                                 ind_expected_const=ind_expected_const,
                                 iteration_threshold=iteration_threshold)
      model_co2_this_chunk[,ii] <- model_out[,"co2"]
      model_temp_this_chunk[,ii] <- model_out[,"temp"] + 15 # Note: Berner (2004; Eq 2.8 and 2.28) assumes present CO2 is 280 ppmv and T is 15 deg C
      ##    # normalize relative to "present" (t=0 Mya)
      ##    for (ii in 1:ncol(model_temp)) {model_temp[,ii] <- model_temp[,ii] - model_temp[58,ii]}
    }
    # use the percent-outbound approach of Mill et al 2019 (Gondwana Research)
    prcout_co2 <- sapply(1:n_sample_this_chunk[cc], function(ss) {percout(model_co2_this_chunk[,ss], windows$co2)})
    prcout_temp <- sapply(1:n_sample_this_chunk[cc], function(ss) {percout(model_temp_this_chunk[,ss], windows$temp)})
    idx_save <- which((prcout_co2 <= prcout_threshold) & (prcout_temp <= prcout_threshold))
      #### ^---- or, do this within the chunk's sampling loop?

    # combine with previous good estimates
    par_calib_save <- rbind(par_calib_save, par_calib[[cc]][idx_save,])
    par_time_save <- abind(par_time_save, time_series_samples[,idx_save,], along=2)
    par_covar_save <- abind(par_covar_save, covariance_samples[,,idx_save,], along=3)
  }
}
tend <- proc.time()
print(paste('precalibration took ',(tend[3]-tbeg[3])/60,' minutes', sep=''))


## ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
##==============================================================================



##==============================================================================
## parameter precalibration
##=========================

## preliminary simulation to get the length
model_out <- model_forMCMC( par_calib=par_calib0,
                            par_time=par_time_center,
                            par_fixed=par_fixed0,
                            parnames_calib=parnames_calib,
                            parnames_fixed=parnames_fixed,
                            parnames_time=parnames_time,
                            age=age,
                            ageN=ageN,
                            ind_const_calib=ind_const_calib,
                            ind_time_calib=ind_time_calib,
                            ind_const_fixed=ind_const_fixed,
                            ind_time_fixed=ind_time_fixed,
                            ind_expected_time=ind_expected_time,
                            ind_expected_const=ind_expected_const,
                            iteration_threshold=iteration_threshold)
n_time <- length(model_out[,"co2"])
time <- model_out[,1]

## run the ensemble

## TODO -- should parallelize this
model_co2 <- model_temp <- mat.or.vec(nr=n_time, nc=n_sample)
tbeg <- proc.time()
for (ii in 1:n_sample) {
  model_out <- model_forMCMC(par_calib=par_calib[ii,],
                             par_fixed=par_fixed0,
                             parnames_calib=parnames_calib,
                             parnames_fixed=parnames_fixed,
                             parnames_time=parnames_time,
                             age=age,
                             ageN=ageN,
                             ind_const_calib=ind_const_calib,
                             ind_time_calib=ind_time_calib,
                             ind_const_fixed=ind_const_fixed,
                             ind_time_fixed=ind_time_fixed,
                             ind_expected_time=ind_expected_time,
                             ind_expected_const=ind_expected_const,
                             iteration_threshold=iteration_threshold,
                             do_sample_tvq=DO_SAMPLE_TVQ,
                             par_time_center=par_time_center,
                             par_time_stdev=par_time_stdev)
  model_co2[,ii] <- model_out[,"co2"]
  model_temp[,ii] <- model_out[,"temp"]
}
tend <- proc.time()

print(paste('precalibration took ',(tend[3]-tbeg[3])/60,' minutes', sep=''))

# Note: Berner (2004; Eq 2.8 and 2.28) assumes present CO2 is 280 ppmv and T is 15 deg C

# normalize relative to "present" (t=0 Mya)
for (ii in 1:ncol(model_temp)) {model_temp[,ii] <- model_temp[,ii] - model_temp[58,ii]}

# use the percent-outbound approach of Mill et al 2019 (Gondwana Research)
source('percent_outbound.R')
source('constraints.R')
prcout_co2 <- sapply(1:n_sample, function(ss) {percout(model_co2[,ss], windows$co2)})
prcout_temp <- sapply(1:n_sample, function(ss) {percout(model_temp[,ss], windows$temp)})





##==============================================================================
## End
##==============================================================================
