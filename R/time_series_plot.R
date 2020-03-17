##==============================================================================
## time_series_plot.R
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

time_series_plot <- function() {

  bb <- "30"
  dd <- "ct"

  plot(-time, time_series_from_priors_quantiles[,"0.5",pp], type='l', xlim=c(-450,0), ylim=ylims, xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
  grid()
  polygon(-c(time,rev(time)), c(time_series_from_priors_quantiles[,"0.025",pp],rev(time_series_from_priors_quantiles[,"0.975",pp])), col=rgb(.5,.5,.5,.25), border=NA)
  lines(-time, time_series_from_priors_quantiles[,"0.5",pp], lwd=1, lty=1, col=rgb(.5,.5,.5))
  polygon(-c(time,rev(time)), c(par_time_quantiles[[bb]][[dd]][,"0.025",pp],rev(par_time_quantiles[[bb]][[dd]][,"0.975",pp])), col=rgb(.6,0,0,.25), border=NA)
  lines(-time, par_time_quantiles[[bb]][[dd]][,"0.5",pp], lwd=1, lty=1, col=rgb(.6,0,0))
  mtext('Time [Myr ago]', side=1, line=2.8)
  #mtext(expression("        Global average\nsurface temperature ["*degree*"C]"), side=2, line=2.2)
  mtext(paste(parnames_time[pp],units), side=2, line=3.2)
  axis(1, at=seq(-400,0,100), labels=c('400','300','200','100','0'))
  ticks=seq(from=ylims[1], to=ylims[2], by=dy)
  axis(2, at=ticks, labels=ticks, las=1)

}

##==============================================================================
## End
##==============================================================================
