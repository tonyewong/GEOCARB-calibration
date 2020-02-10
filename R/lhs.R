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
library(invgamma)

## Set testing number of samples and file name appendix here
## if there aren't enough samples on the MCMC output file, will break.
n_sample <- 1e5 # note this will be doubled if using LHS precalibration
appen <- 'test'

n_node <- 6 # use parallel evaluation of ensembles in Sobol' integration?

param_choice <- 'all'   # Calibrate all 68 parameters? ("all") or only the 6 from Park and Royer 2011 ("PR2011")
data_choice <- 'F2017'    # Which data set?  PR2011 = Park and Royer (2011), or F2017 = Foster et al (2017)
fSR_choice <- 'DT2019'     # Which fSR time series? ("PR2011", "LENTON", "DT2019")
plot.dir <- '../figures/'

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

##==============================================================================



##==============================================================================
## Draw parameter samples
##=======================

## scale up to the actual parameter distributions

## draw parameters by Latin Hypercube (Sample)
set.seed(2020)
parameters_lhs <- randomLHS(n_sample, n_parameters)
par_calib <- parameters_lhs  # initialize

n_const_calib <- length(ind_const_calib)
par_calib <- parameters_lhs  # initialize
colnames(par_calib) <- parnames_calib
for (i in 1:n_const_calib) {
  row_num <- match(parnames_calib[i],input$parameter)
  if(input[row_num, 'distribution_type']=='gaussian') {
    par_calib[,i] <- qnorm(p=parameters_lhs[,ind_const_calib[i]], mean=input[row_num,"mean"], sd=(0.5*input[row_num,"two_sigma"]))
  } else if(input[row_num, 'distribution_type']=='lognormal') {
    par_calib[,i] <- qlnorm(p=parameters_lhs[,ind_const_calib[i]], meanlog=log(input[row_num,"mean"]), sdlog=log(0.5*input[row_num,"two_sigma"]))
  } else {
    print('ERROR - unknown prior distribution type')
  }
}
for (i in (n_const_calib+1):length(parnames_calib)) {
  par_calib[,i] <- qbeta(p=parameters_lhs[,i], shape1=5, shape2=5)
}
##==============================================================================



##==============================================================================
## parameter precalibration
##=========================

## preliminary simulation to get the length
model_out <- model_forMCMC( par_calib=par_calib0,
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

# try the percent-outbound approach of Mill et al 2019 (Gondwana Research)
source('percent_outbound.R')
source('constraints.R')
source('temperature_constraints.R')
prcout_co2 <- sapply(1:n_sample, function(ss) {percout(model_co2[,ss], windows)})
prcout_temp <- sapply(1:n_sample, function(ss) {percout(model_temp[,ss], windows_temp)})





##==============================================================================
## End
##==============================================================================
