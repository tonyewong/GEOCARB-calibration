##==============================================================================
## model_forMCMC_supp.R
##
## Input:
##  par_calib        vector of scalar input parameters
##  par_time         array of time-varying parameters
##  par_fixed        vector structured like par_calib, but for fixed arrays
##                   (not calibrated)
##  age              age, in millions of years
##  ageN             number of time steps
##  ind_const_calib  number of calibration parameters constant in time
##  ind_time_calib   number of calibration parameters varying in time (with ageN
##                   different values)
##  ind_const_fixed  number of fixed parameters constant in time
##  ind_time_fixed   number of fixed parameters varying in time (with ageN
##                   different values)
##  parnames_calib   calibration parameter names
##  parnames_fixed   fixed parameter names
##  parnames_time    time-varying parameter names
##
## Output:
##  age              age, in millinos of years
##  co2              CO2 concentration, ppmv
##  o2               O2 concentration, ppmv
##  temp             temperature, deg C
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

model_forMCMC_supp <- function(par_calib, par_time, par_fixed, parnames_calib,
                          parnames_fixed, parnames_time, age, ageN,
                          ind_const_calib, ind_time_calib,
                          ind_const_fixed, ind_time_fixed,
                          ind_expected_time, ind_expected_const,
                          iteration_threshold) {

  # here, we are assuming that we will calibrate all of the parameters
  # par_calib = a vector of all of the constant parameters
  # par_time = a matrix of all of the time-varying parameters
  # par_fixed will be NULL in this case
  # parnames_calib = names of all of the constant parameters
  # parnames_time = names of all of the time-varying parameters

  N_const_total <- length(c(ind_const_calib, ind_const_fixed))
  N_time_total <- length(c(ind_time_calib, ind_time_fixed))#/ageN

  # set up the time-constant parameter matrix (vector, only 1 simulation set)
  # first length(ind_const_calib) values are the calibration parameters
  # then the fixed values come at the end
  Matrix_56_unordered <- matrix(c(par_calib[ind_const_calib], par_fixed[ind_const_fixed]), nrow=N_const_total, ncol=1)
  rownames(Matrix_56_unordered) <- c( parnames_calib[ind_const_calib],
                                      parnames_fixed[ind_const_fixed] )

  # set up the time-varying parameter matrix
  Matrix_12_unordered <- par_time

  geoRes <- run_geocarb_suppF(Matrix_56=Matrix_56_unordered,
                         Matrix_12=Matrix_12_unordered,
                         age=age,
                         ageN=ageN,
                         iteration_threshold=iteration_threshold,
                         ind_expected_time=ind_expected_time,
                         ind_expected_const=ind_expected_const)

  GEOCARB_output <- cbind(geoRes$age, geoRes$CO2_out, geoRes$O2_out, geoRes$temp_out)
  colnames(GEOCARB_output) <- c('age','co2','o2','temp')

  return(GEOCARB_output)
}

##==============================================================================
## End
##==============================================================================
