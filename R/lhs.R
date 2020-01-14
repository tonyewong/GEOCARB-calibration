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
library(sensitivity)
library(sn)
library(foreach)
library(doParallel)
library(lhs)
library(Hmisc)

## Set testing number of samples and file name appendix here
## if there aren't enough samples on the MCMC output file, will break.
n_sample <- 1e5 # note this will be doubled if using LHS precalibration
appen <- 'final'

n_node <- 6 # use parallel evaluation of ensembles in Sobol' integration?

param_choice <- 'all_stdev'   # Calibrate all 69 parameters? ("all") or only the 6 from Park and Royer 2011 ("PR2011")
data_choice <- 'F2017'    # Which data set?  PR2011 = Park and Royer (2011), or F2017 = Foster et al (2017)
lhood_choice <- 'mixture'  # Mixture model ("mixture") or unimodal ("unimodal")?
fSR_choice <- 'DT2019'     # Which fSR time series? ("PR2011", "LENTON", "DT2019")
dist <- 'sn'               # kernel choice for each data point (sn (skew-normal), ln (log-normal), nm (normal))
plot.dir <- '../figures/'

# latin hypercube precalibration
sens='NS' # valid values:  L1, L2, NS, NSL, pres

# calibration parameters and input data
filename.data <- '../input_data/CO2_Proxy_Foster2017_calib_SN-co2_25Sep2018.csv'
filename.calibinput <- '../input_data/GEOCARB_input_summaries_calib_all_stdev.csv'

# upper bound from Royer et al 2014 (should be yielding a failed run anyhow)
# lower bound relaxed in light of additional proxy data
.upper_bound_co2 <- 50000
.lower_bound_co2 <- 0

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

# Other characteristics are in 'input' and 'time_arrays', which are fed into
# the sensitivity analysis call below.
##==============================================================================



##==============================================================================
## Draw parameter samples
##=======================

## scale up to the actual parameter distributions

## draw parameters by Latin Hypercube (Sample)
set.seed(2019)
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
tbeg <- proc.time()
model_out <- sapply(1:n_sample, function(ss) {
                    model_forMCMC(par_calib=par_calib[ss,],
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
                                  par_time_stdev=par_time_stdev)[,'co2']})
ibad <- NULL
for (ss in 1:n_sample) {
  if( any(model_out[,ss] < windows[,"low"]) | any(model_out[,ss] > windows[,"high"]) ) {
    ibad <- c(ibad,ss)
  }
}
n_time <- nrow(model_out)
parameters_good <- par_calib[-ibad,]
model_good <- model_out[,-ibad]
tend <- proc.time()

# report success rate
print(paste('precalibration took ',(tend[3]-tbeg[3])/60,' minutes', sep=''))
print(paste('with success rate of ',nrow(parameters_good),'/',n_sample, sep=''))

# show model hindcasts against the observational data

# need a fix for the NAN that results when co2_low for data points is 0
idx_low <- which(data_calib$co2_low == 0)
data_calib$co2_low[idx_low] <- 1

# get 5-95% range and median  are cols 1-3; max-post will be 4
quantiles_i_want <- c(0,0.005,.025,.05,.5,.95,.975,0.995,1)
model_quantiles <- mat.or.vec(nr=n_time, nc=(length(quantiles_i_want)+1))
colnames(model_quantiles) <- c('q000','q005','q025','q05','q50','q95','q975','q995','q100','maxpost')
for (t in 1:n_time) {
    model_quantiles[t,1:length(quantiles_i_want)] <- quantile(model_good[t,], quantiles_i_want)
}

## Log scale (model and points, with CO2 and age error bars, no likelihood surface, but with Royer et al 2014 too)
## MS FIGURE 2

pdf(paste(plot.dir,'precal_out.pdf',sep=''),width=4,height=3,colormodel='cmyk', pointsize=11)
par(mfrow=c(1,1), mai=c(.65,.9,.15,.15))
plot(-time, log10(model_quantiles[,'q50']), type='l', xlim=c(-450,0), ylim=c(0.7,log10(30000)), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
#polygon(-c(time,rev(time)), log10(c(model_quantiles[,'q025'],rev(model_quantiles[,'q975']))), col=rgb(.2,.6,.6,.5), border=NA)
polygon(-c(time,rev(time)), log10(c(model_quantiles[,'q000'],rev(model_quantiles[,'q100']))), col="gray60", border=NA)
polygon(-c(time,rev(time)), log10(c(model_quantiles[,'q025'],rev(model_quantiles[,'q975']))), col="gray80", border=NA)
lines(-time, log10(model_quantiles[,'q50']), lwd=2, col='black')
#lines(-time, log10(model_quantiles[,'q50']), lwd=2, col=rgb(.2,.6,.6))
points(-data_calib$age, log10(data_calib$co2), pch=16, cex=0.4, lwd=.4)
for (ii in 1:nrow(data_calib)) {
    arrows(-data_calib$age[ii], log10(data_calib$co2_low[ii]), -data_calib$age[ii], log10(data_calib$co2_high[ii]), length=0.02, angle=90, code=3, lwd=0.25)
    arrows(-data_calib$age_old[ii], log10(data_calib$co2[ii]), -data_calib$age_young[ii], log10(data_calib$co2[ii]), length=0.02, angle=90, code=3, lwd=0.25)
}
mtext('Time [Myr ago]', side=1, line=2.1)
mtext(expression('CO'[2]*' concentration [ppmv]'), side=2, line=3.2)
axis(1, at=seq(-400,0,100), labels=c('400','300','200','100','0'))
ticks=log10(c(seq(10,100,10),seq(200,1000,100),seq(2000,10000,1000),seq(20000,30000,10000)))
axis(2, at=ticks, labels=rep('',length(ticks)))
axis(2, at=log10(c(10,30,100,300,1000,3000,10000,30000)), labels=c('10','30','100','300','1000','3000','10000','30000'), las=1)
legend(-452, log10(85), c('data','median'), pch=c(16,NA), lty=c(NA,1), col=c('black','black'), pt.cex=0.8, cex=.6, bty='n')
legend(-452, log10(35), c('min-max range','95% range'), pch=c(15), col=c('gray60','gray80'), pt.cex=1.2, cex=.6, bty='n', y.intersp=1.6)
minor.tick(nx=5, ny=0, tick.ratio=0.5)
dev.off()




# save precalibration ensemble
today=Sys.Date(); today=format(today,format="%d%b%Y")
saveRDS(parameters_good, paste('../output/precal_parameters_',today,'.rds', sep=''))


if(FALSE) {


set.seed(9102)
colnames(parameters_good) <- parnames_calib
indAvailable <- 1:nrow(parameters_good)
indA <- sample(indAvailable, size=floor(length(indAvailable)/2), replace=FALSE)
indAvailable <- indAvailable[-indA]
indB <- sample(indAvailable, size=length(indA), replace=FALSE)
parameters_sampleA <- parameters_good[indA,]
parameters_sampleB <- parameters_good[indB,]
colnames(parameters_sampleA) <- colnames(parameters_sampleB) <- parnames_calib
##==============================================================================



##==============================================================================
## Sobol analysis
##===============

## Get a reference simulation for integrated sensitivity measure (if using L1, e.g.)
model_out <- model_forMCMC(par_calib=par_calib0,
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
model_ref <- model_out[,'co2']
model_age <- model_out[,'age']

l_scaled <- TRUE
export_names <- c('model_forMCMC', 'run_geocarbF', 'age', 'ageN',
                  'par_fixed0', 'parnames_calib', 'parnames_fixed', 'parnames_time',
                  'ind_const_calib', 'ind_time_calib', 'ind_const_fixed',
                  'ind_time_fixed', 'input', 'ind_expected_time',
                  'ind_expected_const', 'iteration_threshold', 'l_scaled', 'sens',
                  'model_ref', 'data_calib', 'ind_mod2obs',
                  'DO_SAMPLE_TVQ', 'par_time_center', 'par_time_stdev')

n_boot <- .Nboot
conf <- .confidence
if (n_node > 1) {DO_PARALLEL <- TRUE} else {DO_PARALLEL <- FALSE}

## Actually run the Sobol'

source('sobolTony.R')

if (!DO_PARALLEL) {
  t.out <- system.time(s.out <- sobolTony(parameters_sampleA, parameters_sampleB, sens,
                     par_fixed=par_fixed0, parnames_calib=parnames_calib,
                     parnames_fixed=parnames_fixed, parnames_time=parnames_time,
                     age=age, ageN=ageN,
                     ind_const_calib=ind_const_calib, ind_time_calib=ind_time_calib,
                     ind_const_fixed=ind_const_fixed, ind_time_fixed=ind_time_fixed,
                     input=input, ind_expected_time=ind_expected_time,
                     ind_expected_const=ind_expected_const, model_ref=model_ref,
                     iteration_threshold=iteration_threshold, data_calib=data_calib,
                     do_sample_tvq=DO_SAMPLE_TVQ, par_time_center=par_time_center,
                     par_time_stdev=par_time_stdev, n_boot=n_boot, conf=conf, second=.second))
} else {
  t.out <- system.time(s.out <- sobolTony(parameters_sampleA, parameters_sampleB, sens,
                     par_fixed=par_fixed0, parnames_calib=parnames_calib,
                     parnames_fixed=parnames_fixed, parnames_time=parnames_time,
                     age=age, ageN=ageN,
                     ind_const_calib=ind_const_calib, ind_time_calib=ind_time_calib,
                     ind_const_fixed=ind_const_fixed, ind_time_fixed=ind_time_fixed,
                     input=input, ind_expected_time=ind_expected_time,
                     ind_expected_const=ind_expected_const, model_ref=model_ref,
                     iteration_threshold=iteration_threshold, data_calib=data_calib,
                     do_sample_tvq=DO_SAMPLE_TVQ, par_time_center=par_time_center,
                     par_time_stdev=par_time_stdev, parallel=TRUE, n_core=n_node,
                     export_names=export_names, n_boot=n_boot, conf=conf, second=.second))
}

print(paste('Sobol simulations took ',t.out[3],' seconds', sep=''))

today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.sobol <- paste('../output/sobol_sens',sens,'_',appen,'_',today,'.rds', sep='')
saveRDS(s.out, filename.sobol)

## Check convergence - is maximum confidence interval width < 10% of the highest
## total order index?
max_sens_ind <- 0.1*max(s.out$T)
max_conf_int <- max(max(s.out$S[,3]-s.out$S[,2]), max(s.out$T[,3]-s.out$T[,2]))
print(paste('max. conf int=',max_conf_int, ' want less than:  0.1 * max. sensitivity index=',max_sens_ind, sep=''))
##==============================================================================



##==============================================================================
## Write indices, results we'd need to file

today=Sys.Date(); today=format(today,format="%d%b%Y")
file.sobolout1 <- paste('../output/geocarb_sobol-1-tot_sens',sens,'_',appen,'_',today,'.txt',sep='')
file.sobolout2 <- paste('../output/geocarb_sobol-2_sens',sens,'_',appen,'_',today,'.txt',sep='')

headers.1st.tot <- matrix(c('Parameter', 'S1', 'S1_conf_low', 'S1_conf_high',
                            'ST', 'ST_conf_low', 'ST_conf_high'), nrow=1)
output.1st.tot  <- data.frame(cbind( parnames_calib, s.out$S, s.out$T))
write.table(headers.1st.tot, file=file.sobolout1, append=FALSE, sep = " ",
            quote=FALSE    , row.names = FALSE , col.names=FALSE)
write.table(output.1st.tot , file=file.sobolout1, append=TRUE , sep = " ",
            quote=FALSE    , row.names = FALSE , col.names=FALSE)

if (.second) {
  headers.2nd <- matrix(c('Parameter_1', 'Parameter_2', 'S2', 'S2_conf_low',
                          'S2_conf_high'), nrow=1)

  output.2nd <- data.frame(cbind( s.out$S2.names,s.out$S2 ))
  write.table(headers.2nd    , file=file.sobolout2, append=FALSE , sep = " ",
              quote=FALSE    , row.names = FALSE , col.names=FALSE)
  write.table(output.2nd     , file=file.sobolout2, append=TRUE , sep = " ",
              quote=FALSE    , row.names = FALSE , col.names=FALSE)
}
##==============================================================================



}



##==============================================================================
## End
##==============================================================================
