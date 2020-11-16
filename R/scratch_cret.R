
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

idx_cret <- c(48:50)
experiments <- c("ctrl","cret"); n_experiments <- length(experiments)
kdes <- model_temp <- model_co2 <- par_calib <- par_time <- vector("list", length=n_experiments)
names(kdes) <- names(model_temp) <- names(model_co2) <- names(par_calib) <- names(par_time) <- experiments

## load both sets of model outputs
load("../output/lhs_cret_param_ct_out50_ctrl.RData")
par_calib$ctrl <- par_calib_save
par_time$ctrl <- par_time_save

load("../output/lhs_cret_param_ct_out50_cret.RData")
par_calib$cret <- par_calib_save
par_time$cret <- par_time_save

## set up stuff for running model
param_choice <- 'all'   # Calibrate all 68 parameters? ("all") or only the 6 from Park and Royer 2011 ("PR2011")
data_choice <- 'F2017'    # Which data set?  PR2011 = Park and Royer (2011), or F2017 = Foster et al (2017)
fSR_choice <- 'DT2019'     # Which fSR time series? ("PR2011", "LENTON", "DT2019")
plot.dir <- '../figures/'
threshold_choices <- c(.5)
use_co2 <- use_temperature <- TRUE

# calibration parameters and input data
filename.data <- '../input_data/CO2_Proxy_Foster2017_calib_SN-co2_25Sep2018.csv'
filename.calibinput <- '../input_data/GEOCARB_input_summaries_calib_all.csv'

# create output file names
data_appen <- ""
if (use_co2) {data_appen <- paste(data_appen,"c",sep="")}
if (use_temperature) {data_appen <- paste(data_appen,"t",sep="")}
if (!exists("appen")) {appen <- ""}

source("model_setup.R")
source('run_geocarbF.R')
source('model_forMCMC.R')
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
## Establish upper and lower bounds for each parameter
## (NA if none (inf))
##====================================================

bounds <- mat.or.vec(nr=length(ind_const_calib), nc=2)
colnames(bounds) <- c("lower","upper")
rownames(bounds) <- parnames_const_calib
for (ii in 1:nrow(bounds)) {
  row_num <- match(parnames_const_calib[ii],input$parameter)
  if (as.vector(input[row_num,"lower_limit"]) != "_inf") {
    # there is a lower bound
    bounds[ii,"lower"] <- as.numeric(as.vector(input[row_num,"lower_limit"]))
  } else {
    # no lower bound
    bounds[ii,"lower"] <- NA
  }
  if (is.finite(input[row_num,"upper_limit"])) {
    # there is an upper bound
    bounds[ii,"upper"] <- input[row_num,"upper_limit"]
  } else {
    # no upper bound
    bounds[ii,"upper"] <- NA
  }
}
##==============================================================================

## run the ensembles using the control/cretaceous experiment results
n_const_calib <- length(ind_const_calib)
for (ee in experiments) {
  model_co2[[ee]] <- model_temp[[ee]] <- mat.or.vec(nr=n_time, nc=nrow(par_calib[[ee]]))
  prcout_co2 <- prcout_temp <- rep(NA, nrow(par_calib[[ee]]))
  for (ii in 1:nrow(par_calib[[ee]])) {
      model_out <- model_forMCMC(par_calib=par_calib[[ee]][ii,],
                                 par_time=par_time[[ee]][,ii,],
                                 par_fixed=par_const_fixed0,
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
      model_co2[[ee]][,ii] <- model_out[,"co2"]
      model_temp[[ee]][,ii] <- model_out[,"temp"] + 15 # Note: Berner (2004; Eq 2.8 and 2.28) assumes present CO2 is 280 ppmv and T is 15 deg C
      ##    # normalize relative to "present" (t=0 Mya)
      ##    for (ii in 1:ncol(model_temp)) {model_temp[,ii] <- model_temp[,ii] - model_temp[58,ii]}
      # use the percent-outbound approach of Mill et al 2019 (Gondwana Research)
      prcout_co2[ii] <- percout(model_co2[[ee]][,ii], windows$co2)
      prcout_temp[ii] <- percout(model_temp[[ee]][,ii], windows$temp)
  }
}

##==============================================================================
## PLOTTING

# preliminary plot of time series for temperature
par(mfrow=c(2,1))
ee <- "ctrl"
plot(-age, model_temp[[ee]][,1], type='l', ylim=c(0,50)); for (ii in 1:nrow(par_calib[[ee]])) {lines(-age, model_temp[[ee]][,ii], col='gray80')}
lines(-age, windows$temp[,"low"],col="red"); lines(-age, windows$temp[,"high"],col="red")
ee <- "cret"
plot(-age, model_temp[[ee]][,1], type='l', ylim=c(0,50)); for (ii in 1:nrow(par_calib[[ee]])) {lines(-age, model_temp[[ee]][,ii], col='gray80')}
lines(-age, windows$temp[,"low"],col="red"); lines(-age, windows$temp[,"high"],col="red")

##==============================================================================
## NUMBERS

write.csv(cor(par_calib$ctrl), file="corr_ctrl.csv")
write.csv(cor(par_calib$cret), file="corr_cret.csv")

##==============================================================================
## parameter distributions

for (ee in experiments) {
  kdes[[ee]] <- vector("list", ncol(par_calib[[ee]]))
  names(kdes[[ee]]) <- parnames_const
  for (pp in parnames_const) {
    kdes[[ee]][[pp]] <- density(par_calib[[ee]][,match(pp,parnames_const)])
  }
}

# boxplots
pdf('../figures/par_calib_comparison_cret.pdf',width=7.5,height=9, colormodel='cmyk', pointsize=11)
par(mfrow=c(10,6), mai=c(.4,.35,.05,.05))
for (pp in 1:length(parnames_const)) {
    boxplot(par_calib$ctrl[,pp], par_calib$cret[,pp], names=c("Ctrl", "Cret"),
            horizontal=TRUE, las=1, boxwex=0.5, pch=4, lty=1, lwd=1)
    mtext(side=1, text=parnames_const[pp], line=2.1, cex=0.75)
}
dev.off()

# deltaT2X density plots
pdf('../figures/deltaT2X_cret.pdf',width=4,height=3, colormodel='cmyk', pointsize=11)
par(mfrow=c(1,1), mai=c(.7,.35,.13,.15))
plot(kdes$ctrl$deltaT2X$x, kdes$ctrl$deltaT2X$y, type='l', lwd=1.5, col='black',
     xlab='', ylab='', yaxt='n', xlim=c(1.65,6.65), ylim=c(0,0.65), yaxs='i')
lines(kdes$cret$deltaT2X$x, kdes$cret$deltaT2X$y, lwd=1.5, lty=5, col='firebrick')
mtext(side=1, text=expression(Delta*"T2X [deg C]"), line=2.3)
mtext(side=2, text="Probability density", line=0.5)
dev.off()

##==============================================================================
## End
##==============================================================================
