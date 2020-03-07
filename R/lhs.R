##==============================================================================
## lhs.R
##
## todo...
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================
## Copyright 2020 Tony Wong
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

# calibration parameters and input data
filename.data <- '../input_data/CO2_Proxy_Foster2017_calib_SN-co2_25Sep2018.csv'
filename.calibinput <- '../input_data/GEOCARB_input_summaries_calib_all.csv'

# create output file names
data_appen <- ""
if (use_co2) {data_appen <- paste(data_appen,"c",sep="")}
if (use_temperature) {data_appen <- paste(data_appen,"t",sep="")}
filename.covarout <- paste('../output/lhs_covar_',data_appen,'_out',prcout_threshold*100,'.RData',sep='')
filename.paramout <- paste('../output/lhs_param_',data_appen,'_out',prcout_threshold*100,'.RData',sep='')

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
