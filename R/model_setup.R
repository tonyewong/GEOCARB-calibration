##==============================================================================
## model_setup.R
##
## Short-cut routine for configuring a model simulation.
##
## Tony Wong (aewsma@rit.edu)
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

# If using the Foster et al 2017 data set, which proxy sets to assimilate?
#   (set what you want to "TRUE", others to "FALSE")
data_to_assim <- cbind( c("paleosols" , TRUE),
                        c("alkenones" , TRUE),
                        c("stomata"   , TRUE),
                        c("boron"     , TRUE),
                        c("liverworts", TRUE) )

if (fSR_choice=="PR2011") {
  USE_LENTON_FSR <- FALSE
  USE_DT2019_FSR <- FALSE
} else if (fSR_choice=="LENTON") {
  USE_LENTON_FSR <- TRUE
  USE_DT2019_FSR <- FALSE
} else if (fSR_choice=="DT2019") {
  USE_LENTON_FSR <- FALSE
  USE_DT2019_FSR <- TRUE
} else {
  print("ERROR: unknown fSR_choice")
}

# create appendix tag and file name for this simulation set
appen <- ""
if (data_choice=="PR2011") {appen <- paste(appen,"dP", sep="")} else if (data_choice=="F2017") {appen <- paste(appen,"dF", sep="")}
if (substr(param_choice,1,6)=="PR2011") {appen <- paste(appen,"pP", sep="")} else if (substr(param_choice,1,3)=="all") {appen <- paste(appen,"pA", sep="")}
if (fSR_choice=="PR2011") {appen <- paste(appen,"sO", sep="")} else if (fSR_choice=="LENTON") {appen <- paste(appen,"sL", sep="")} else if (fSR_choice=="DT2019") {appen <- paste(appen,"sR", sep="")}
today <- Sys.Date(); today <- format(today,format="%d%b%Y")
filename.lhs = paste('../output/geocarb_mcmcoutput_',appen,'_',today,'.RData',sep="")

# upper bound from Royer et al 2014 (should be yielding a failed run anyhow)
# lower bound relaxed in light of additional proxy data
.upper_bound_co2 <- 50000
.lower_bound_co2 <- 0

# sampling the time-varying arrays
if (substr(param_choice,1,6)=="PR2011") {
  DO_SAMPLE_TVQ <- FALSE
} else {
  DO_SAMPLE_TVQ <- TRUE
}

# Now, set up the parameters as desired. Possibly the same.
filename.calibinput <- paste('../input_data/GEOCARB_input_summaries_calib_',param_choice,'.csv', sep='')
source('parameterSetup_tvq.R')
source('model_forMCMC_tvq.R')
source('run_geocarbF.R')
##==============================================================================



##==============================================================================
## Data
##=====

source("getData.R")
##==============================================================================


##==============================================================================
## Calibration parameter prior distributions
##==========================================

# Get model parameter prior distributions
names <- as.character(input$parameter)
bound_lower <- rep(NA, length(names))
bound_upper <- rep(NA, length(names))

ind_neg_inf <- which(input[,'lower_limit']=='_inf')
bound_lower[ind_neg_inf] <- -Inf
bound_lower[setdiff(1:length(names), ind_neg_inf)] <- as.numeric(as.character(input$lower_limit[setdiff(1:length(names), ind_neg_inf)]))
bound_upper <- input$upper_limit

bounds <- cbind(bound_lower, bound_upper)
rownames(bounds) <- as.character(input$parameter)

# only actually need the calibration parameters' bounds, so reformat the bounds
# array to match the vector of calibration parameters
bounds_calib <- mat.or.vec(nr=length(parnames_calib), nc=2)
colnames(bounds_calib) <- c('lower','upper')
rownames(bounds_calib) <- parnames_calib
for (i in 1:length(parnames_calib)) {
  bounds_calib[i,'lower'] <- bounds[parnames_calib[i],'bound_lower']
  bounds_calib[i,'upper'] <- bounds[parnames_calib[i],'bound_upper']
}

rm(list=c('bound_lower','bound_upper','bounds'))

##==============================================================================
## End
##==============================================================================
