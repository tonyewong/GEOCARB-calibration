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
n_sample <- 500000
n_sample_per_chunk <- 15000 # maximum number of time series samples to consider at once
n_sample_min <- 40000 # minimum number of samples we'd be happy with; stop after this to avoid overrunning RAM
n_node <- 1 # distribute the chunks across multiple nodes?

param_choice <- 'all'   # Calibrate all 68 parameters? ("all") or only the 6 from Park and Royer 2011 ("PR2011")
data_choice <- 'F2017'    # Which data set?  PR2011 = Park and Royer (2011), or F2017 = Foster et al (2017)
fSR_choice <- 'DT2019'     # Which fSR time series? ("PR2011", "LENTON", "DT2019")
plot.dir <- '../figures/'
prcout_threshold <- 0.4  # only save simulations with percent_outbound < this

# calibration parameters and input data
filename.data <- '../input_data/CO2_Proxy_Foster2017_calib_SN-co2_25Sep2018.csv'
filename.calibinput <- '../input_data/GEOCARB_input_summaries_calib_all.csv'
filename.covarout <- paste('../output/lhs_covar_out',prcout_threshold*100,'.RData',sep='')
filename.paramout <- paste('../output/lhs_param_out',prcout_threshold*100,'.RData',sep='')
                      
if(Sys.info()['user']=='tony') {
  # Tony's local machine (if you aren't me, you almost certainly need to change this...)
  machine <- 'local'
  setwd('/Users/tony/codes/GEOCARB-calibration/R')
} else if(Sys.info()['user']=='aewsma') {
  machine <- 'office'
  setwd('/Users/aewsma/codes/GEOCARB-calibration/R')
} else {
  # assume on another cluster of some kind...
  machine <- 'remote'
  setwd('~/work/codes/GEOCARB-calibration/R')
}

source("model_setup.R")
n_parameters <- length(parnames_const_calib)
n_parameters_time <- length(parnames_time_calib)

## preliminary simulation to get the length
model_out <- model_forMCMC( par_calib=par_calib0,
                            par_time=par_time_center,
                            par_fixed=par_fixed0,
                            parnames_calib=parnames_const_calib,
                            parnames_fixed=parnames_const_fixed,
                            parnames_time=parnames_time_calib,
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
## Do the actual sampling
##========================

source("lhs_sampling.R")

##==============================================================================



##==============================================================================
## Save
##========================

save(list=c("par_covar_save"), file=filename.covarout)
save(list=c("par_calib_save","par_time_save"), file=filename.paramout)

##==============================================================================



##==============================================================================
## End
##==============================================================================
