
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

idx_cret <- c(49)
experiments <- c("ctrl","cret"); n_experiments <- length(experiments)
kdes <- model_temp <- model_co2 <- par_calib <- par_time <- vector("list", length=n_experiments)
names(kdes) <- names(model_temp) <- names(model_co2) <- names(par_calib) <- names(par_time) <- experiments

## load both sets of model outputs
#load("../output/lhs_cret_param_ct_out50_ctrl.RData")
load("../output/lhs_param_ct_out50.RData")
par_calib$ctrl <- par_calib_save
par_time$ctrl <- par_time_save

load("../output/lhs_cret_param_ct_out50_cret.RData")
par_calib$cret <- par_calib_save
par_time$cret <- par_time_save

## time-varying parameter means
par_time_mean <- vector("list", n_experiments); names(par_time_mean) <- experiments
for (ee in experiments) {
  par_time_mean[[ee]] <- mat.or.vec(nr=dim(par_time[[ee]])[1], nc=dim(par_time[[ee]])[3])
  for (pp in 1:dim(par_time[[ee]])[3]) {
    par_time_mean[[ee]][,pp] <- apply(X=par_time[[ee]][,,pp], MARGIN=1, FUN=mean)
  }
}
write.csv(par_time_mean$cret, file="par_time_mean_cret.csv")

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
prcout_co2 <- prcout_temp <- vector('list', n_experiments)
names(prcout_co2) <- names(prcout_temp) <- experiments
n_const_calib <- length(ind_const_calib)
for (ee in experiments) {
  model_co2[[ee]] <- model_temp[[ee]] <- mat.or.vec(nr=n_time, nc=nrow(par_calib[[ee]]))
  prcout_co2[[ee]] <- prcout_temp[[ee]] <- rep(NA, nrow(par_calib[[ee]]))
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
      prcout_co2[[ee]][ii] <- percout(model_co2[[ee]][,ii], windows$co2)
      prcout_temp[[ee]][ii] <- percout(model_temp[[ee]][,ii], windows$temp)
  }
}

##==============================================================================
## PLOTTING

# preliminary plot of time series for temperature
pdf('../figures/temp_cret.pdf',width=4,height=6, colormodel='cmyk', pointsize=11)
par(mfrow=c(2,1), mai=c(.7,.65,.3,.15))
ee <- "ctrl"
plot(-age, model_temp[[ee]][,1], type='l', xlim=c(-420,0), ylim=c(0,50), xaxt='n', yaxt='n', xlab='', ylab='')
for (ii in 1:nrow(par_calib[[ee]])) {lines(-age, model_temp[[ee]][,ii], col='gray80')}
lines(-age, windows$temp[,"low"],col="red", lwd=1.5); lines(-age, windows$temp[,"high"],col="red", lwd=1.5)
mtext(side=1, text="Time [Myr ago]", line=2.2)
mtext(side=2, text="Temperature [deg C]", line=2.2)
mtext(side=3, text="Control", line=0.15)
axis(side=1, at=seq(from=-400, to=0, by=100), labels=seq(from=400, to=0, by=-100), padj=-0.4)
axis(side=1, at=seq(from=-350, to=0, by=100), labels=rep('', 4), lwd=0.5)
axis(side=2, at=seq(from=0, to= 50, by=10), las=1, hadj=0.8)
ee <- "cret"
plot(-age, model_temp[[ee]][,1], type='l', xlim=c(-420,0), ylim=c(0,50), xaxt='n', yaxt='n', xlab='', ylab='')
for (ii in 1:nrow(par_calib[[ee]])) {lines(-age, model_temp[[ee]][,ii], col='gray80')}
lines(-age, windows$temp[,"low"],col="red", lwd=1.5); lines(-age, windows$temp[,"high"],col="red", lwd=1.5)
mtext(side=1, text="Time [Myr ago]", line=2.2)
mtext(side=2, text="Temperature [deg C]", line=2.2)
mtext(side=3, text="Cretaceous-matching", line=0.15)
axis(side=1, at=seq(from=-400, to=0, by=100), labels=seq(from=400, to=0, by=-100), padj=-0.4)
axis(side=1, at=seq(from=-350, to=0, by=100), labels=rep('', 4), lwd=0.5)
axis(side=2, at=seq(from=0, to= 50, by=10), las=1, hadj=0.8)
dev.off()

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
legend(4, 0.65, c("Control", "Cretaceous-matching"), lty=c(1,5), lwd=1.5, col=c("black", "firebrick"), bty="n")
dev.off()

##==============================================================================
## changes in deltaT2X with experiment, and prcout

idx30 <- idx40 <- idx50 <- vector('list', length=n_experiments)
names(idx30) <- names(idx40) <- names(idx50) <- experiments
for (ee in experiments) {
    idx30[[ee]] <- which(prcout_temp[[ee]] <= 0.3 & prcout_co2[[ee]] <= 0.3)
    idx40[[ee]] <- which(prcout_temp[[ee]] <= 0.4 & prcout_co2[[ee]] <= 0.4)
    idx50[[ee]] <- which(prcout_temp[[ee]] <= 0.5 & prcout_co2[[ee]] <= 0.5)
}
print(quantile(par_calib$cret[,10], c(.05,.5,.95)))
print(quantile(par_calib$cret[idx50$cret,10], c(.05,.5,.95)))
print(quantile(par_calib$cret[idx40$cret,10], c(.05,.5,.95)))
print(quantile(par_calib$cret[idx30$cret,10], c(.05,.5,.95)))

bp <- vector('list', length=6); names(bp) <- c("Cret30", "Ctrl30", "Cret40", "Ctrl40", "Cret50", "Ctrl50")
bp$Ctrl50 <- par_calib$ctrl[idx50$ctrl,10]
bp$Cret50 <- par_calib$cret[idx50$cret,10]
bp$Ctrl40 <- par_calib$ctrl[idx40$ctrl,10]
bp$Cret40 <- par_calib$cret[idx40$cret,10]
bp$Ctrl30 <- par_calib$ctrl[idx30$ctrl,10]
bp$Cret30 <- c(par_calib$cret[idx30$cret,10], par_calib_save[,10])

pdf('../figures/deltaT2X_boxplots_cret.pdf',width=4,height=3, colormodel='cmyk', pointsize=11)
par(mfrow=c(1,1), mai=c(.7,.35,.13,.15))
boxplot(bp, lty=1, pch=4, horizontal=TRUE, las=1, col=c("lightblue","gray70","lightblue","gray70","lightblue","gray70"), yaxt='n', cex.axis=1.2)
mtext(side=1, expression(Delta * "T2X ["^o *"C]"), line=2.6, cex=0.85)
axis(side=2, at=1:6, las=1, labels=c("Cret-30%","Ctrl-30%","Cret-40%","Ctrl-40%","Cret-50%","Ctrl-50%"), cex.axis=1.1)
dev.off()


##==============================================================================
## time series plots

c(7,8,9,10)

pp <- 7; plot(-age, par_time_mean$ctrl[,pp], col="black", lwd=1.5, type='l', ylab=parnames_time[pp]); lines(-age, par_time_mean$cret[,pp], lwd=1.5, col="firebrick")


##==============================================================================
## NUMBERS

Tcret <- vector('list', n_experiments); names(Tcret) <- experiments
for (ee in experiments) {Tcret[[ee]] <- model_temp[[ee]][idx_cret,]; par_calib[[ee]] <- cbind(par_calib[[ee]], Tcret[[ee]])}

idx_to_use <- idx40
corr <- corr_p <- vector('list', n_experiments); names(corr) <- names(corr_p) <- experiments
idx_sample <- sample(1:length(idx_to_use$ctrl), size=length(idx_to_use$cret), replace=FALSE) # subsample to avoid sample size effects
for (ee in experiments) {
    if (ee=="ctrl") {idx <- idx_sample} else {idx <- 1:length(idx_to_use$cret)}
    corr[[ee]] <- corr_p[[ee]] <- mat.or.vec(nr=ncol(par_calib[[ee]]), nc=ncol(par_calib[[ee]]))
    colnames(corr[[ee]]) <- colnames(corr_p[[ee]]) <- rownames(corr[[ee]]) <- rownames(corr_p[[ee]]) <- colnames(par_calib[[ee]])
    for (p1 in 1:(ncol(par_calib[[ee]])-1)) {
        for (p2 in (p1+1):ncol(par_calib[[ee]])) {
            corr_result <- cor.test(par_calib[[ee]][idx_to_use[[ee]][idx], p1], par_calib[[ee]][idx_to_use[[ee]][idx], p2], method="pearson")
            corr[[ee]][p1,p2] <- corr_result$estimate
            corr_p[[ee]][p1,p2] <- corr_result$p.value
        }
    }
}
write.csv(corr$ctrl, file="corr_ctrl40.csv")
write.csv(corr_p$ctrl, file="corr_p_ctrl40.csv")
write.csv(corr$cret, file="corr_cret40.csv")
write.csv(corr_p$cret, file="corr_p_cret40.csv")

##==============================================================================
## End
##==============================================================================
