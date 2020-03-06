##==============================================================================
## percent_outbound.R
##
## Compute the proportion of model simulation points from `time_series`
## that lie outside the given `windows`.
##
## In GEOCARB, many model simulations using naive parameters (sampled from
## the prior distributions, for example) will lead to model failure in the
## sense that some of the timesteps have infinities or NANs. In those cases
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

percout <- function(time_series, windows) {
  if (any(is.infinite(time_series))) {return(1)} else {
    p_out <- sum(time_series < windows[,1] | time_series > windows[,2]) / length(time_series)
    return(p_out)
  }
}

##==============================================================================
## End
##==============================================================================
