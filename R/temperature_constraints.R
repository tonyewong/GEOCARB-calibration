##==============================================================================
## temperature_constraints.R
##
## void function. Returns an object `windows` that is ntime x 2, where the two
## columns are a lower and upper bound for windows at each time slice. The
## precalibration will only permit simulations that pass through all of these
## windows (or some minimum number of the windows).
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


filename.temperature <- "../input_data/Mills_GR_2019_temp_co2.csv"
data_temps <- read.csv(filename.temperature, col.names=c('time','co2_max','co2_min','T_max','T_min','T_avg'))

# comparisons of temperatures - all relative to "present" (t=0)
# For the Mills et al, upper and lower are the average +/- 1 sigma.
# We want nsig (from `constraints.R`), so take as bounds:
#    upper |--> average + nsig*(upper - average)
#    lower |--> average - nsig*(average - lower)

temp_upper <- data_temps[,"T_avg"] + nsig*(data_temps[,"T_max"]-data_temps[,"T_avg"])
temp_lower <- data_temps[,"T_avg"] - nsig*(data_temps[,"T_avg"]-data_temps[,"T_min"])

# normalize relative to "present"
subtract <- data_temps[1,"Tavg"]
temp_upper <- temp_upper - subtract
temp_lower <- temp_lower - subtract


##==============================================================================
## End
##==============================================================================
