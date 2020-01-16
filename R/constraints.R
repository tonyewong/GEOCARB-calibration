##==============================================================================
## constraints.R
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

if (TRUE) {

## for each time period, pool the observations from that 10 Myr slice

nsig <- 4
windows <- mat.or.vec(n_time, 2)
colnames(windows) <- c("low","high")

for (tt in 1:n_time) {
    # start each window at 0 to 50000 ppmv range
    windows[tt,"low"] <- 0
    windows[tt,"high"] <- 50000
    idx <- which(time[tt]-data_calib$age < 10 & time[tt]-data_calib$age >= 0)
    if (length(idx) > 0) {
        sigmas <- log(data_calib$co2_high[idx]/data_calib$co2[idx])
        centers <- log(data_calib$co2[idx])
        tops <- centers + nsig*sigmas
        bots <- centers - nsig*sigmas
        tops_max <- exp(max(tops))
        bots_min <- exp(min(bots))
        windows[tt,"low"] <- max(c(bots_min,windows[tt,"low"]))
        windows[tt,"high"] <- min(c(tops_max,windows[tt,"high"]))
    }
}
# make sure the last window is wider
windows[58,"low"] <- 180
windows[58,"high"] <- 400

}


if (FALSE) {

ibad <- NULL
for (ss in 1:n_sample) {
  if( any(model_out[,ss] < .lower_bound_co2) |
      any(model_out[,ss] > .upper_bound_co2) |
      model_out[58,ss] < 180 | model_out[58,ss] > 400) {
      ##model_out[58,ss] < 280 | model_out[58,ss] > 400) {
    ibad <- c(ibad,ss)
  }
}

n_time <- nrow(model_out)
parameters_good <- par_calib[-ibad,]
model_good <- model_out[,-ibad]
tend <- proc.time()
}


##==============================================================================
## End
##==============================================================================
