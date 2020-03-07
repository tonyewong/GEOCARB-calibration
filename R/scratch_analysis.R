##==============================================================================
## analysis.R
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

rm(list=ls())

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

##
## read and save quantiles of all the parameters for the experiments
## when computing the quantiles, sample down to only 10000 from each experiment
## (so number of samples does not bias results)
##

num_samples <- 10000
quantiles_i_want <- c(0, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 1)
prc_outbound <- as.character(c(30,35,40,45,50))
data_sets <- c("c","t","ct")
par_calib <- par_quantiles <- par_time <- idx_sample <- vector("list", length(prc_outbound))
names(par_calib) <- names(par_quantiles) <- names(par_time) <- names(idx_sample) <- prc_outbound

for (bb in prc_outbound) {
    par_calib[[bb]] <- par_quantiles[[bb]] <- par_time[[bb]] <- vector("list", length(data_sets))
    names(par_calib[[bb]]) <- names(par_quantiles[[bb]]) <- names(par_time[[bb]]) <- data_sets
    for (dd in data_sets) {
        load(paste("../output/lhs_param_",dd,"_out",bb,".RData", sep=""))
        par_calib[[bb]][[dd]] <- cbind(par_calib_save, rep(0,nrow(par_calib_save)))
        colnames(par_calib[[bb]][[dd]]) <- c(colnames(par_calib_save), "deltaT2Xglac")
        par_time[[bb]][[dd]] <- par_time_save
        par_quantiles[[bb]][[dd]] <- mat.or.vec(nr=ncol(par_calib[[bb]][[dd]])+1, nc=length(quantiles_i_want))
        rownames(par_quantiles[[bb]][[dd]]) <- c(colnames(par_calib[[bb]][[dd]]), "deltaT2Xglac")
        colnames(par_quantiles[[bb]][[dd]]) <- as.character(quantiles_i_want)
        idx_sample[[bb]][[dd]] <- sample(1:nrow(par_calib[[bb]][[dd]]), size=num_samples, replace=FALSE)
        for (pp in 1:ncol(par_calib[[bb]][[dd]])) {
            par_quantiles[[bb]][[dd]][pp,] <- quantile(par_calib[[bb]][[dd]][idx_sample[[bb]][[dd]],pp], quantiles_i_want)
        }
        par_calib[[bb]][[dd]][,"deltaT2Xglac"] <- par_calib[[bb]][[dd]][,"deltaT2X"]*par_calib[[bb]][[dd]][,"GLAC"]
        par_quantiles[[bb]][[dd]]["deltaT2Xglac",] <- quantile(par_calib[[bb]][[dd]][idx_sample[[bb]][[dd]],"deltaT2Xglac"], quantiles_i_want)
    }
}

##
## run ensemble with the 30-both parameters
## --> will need to read the time series parameters from lhs_covar_ct_out30.RData
##

#load("../output/lhs_covar_ct_out30.RData")  #> dim(par_covar_save)  [1]    58    58 10000    12
## make sure the samples are concomitant with the constant parameters
#par_covar <- par_covar_save[,,idx_sample$`30`$`ct`,]
#rm(list=c("par_covar_save")) # to save ram

## model setup

param_choice <- 'all'   # Calibrate all 68 parameters? ("all") or only the 6 from Park and Royer 2011 ("PR2011")
data_choice <- 'F2017'    # Which data set?  PR2011 = Park and Royer (2011), or F2017 = Foster et al (2017)
fSR_choice <- 'DT2019'     # Which fSR time series? ("PR2011", "LENTON", "DT2019")
filename.data <- '../input_data/CO2_Proxy_Foster2017_calib_SN-co2_25Sep2018.csv'
source("model_setup.R")

## run hindcasts

model_hindcast <- vector("list", length(prc_outbound))
names(model_hindcast) <- prc_outbound

for (bb in prc_outbound) {
    model_hindcast[[bb]] <- vector("list", length(data_sets))
    names(model_hindcast[[bb]]) <- data_sets
    for (dd in data_sets) {
        model_hindcast[[bb]][[dd]] <- vector('list', 3)
        names(model_hindcast[[bb]][[dd]]) <- c("co2","temp")
        model_hindcast[[bb]][[dd]]$co2 <- model_hindcast[[bb]][[dd]]$temp <- mat.or.vec(nr=58, nc=num_samples)
        for (ii in 1:num_samples) {
            model_out<- model_forMCMC(par_calib=par_calib[[bb]][[dd]][idx_sample[[bb]][[dd]][ii],],
                                       par_time=par_time[[bb]][[dd]][,idx_sample[[bb]][[dd]][ii],],
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
            model_hindcast[[bb]][[dd]]$co2[,ii] <- model_out[,"co2"]
            model_hindcast[[bb]][[dd]]$temp[,ii] <- model_out[,"temp"] + 15 # as in Berner 2004
        }
    }
}
n_time <- length(model_out[,"co2"])
time <- model_out[,1]

## compute model quantiles for plotting against proxy data, for each experiment

model_quantiles <- vector("list", length(prc_outbound))
names(model_quantiles) <- prc_outbound
for (bb in prc_outbound) {
		model_quantiles[[bb]] <- vector("list", length(data_sets))
    names(model_quantiles[[bb]]) <- data_sets
    for (dd in data_sets) {
    		model_quantiles[[bb]][[dd]] <- vector("list", 2)
        names(model_quantiles[[bb]][[dd]]) <- c("co2","temp")
        for (oo in c("co2", "temp")) {
        		model_quantiles[[bb]][[dd]][[oo]] <- mat.or.vec(nr=58, nc=length(quantiles_i_want))
            for (tt in 1:58) {

## FOR NOW GET RID OF FAILED RUNS - LATER GO BACK AND RE-RUN, DISCARDING ANY THAT ARE INFINITY
idx <- which(is.finite(model_hindcast[[bb]][[dd]][[oo]][tt,]))
model_quantiles[[bb]][[dd]][[oo]][tt,] <- quantile(model_hindcast[[bb]][[dd]][[oo]][tt,idx], quantiles_i_want)

            		#model_quantiles[[bb]][[dd]][[oo]][tt,] <- quantile(model_hindcast[[bb]][[dd]][[oo]][tt,], quantiles_i_want)
            }
            colnames(model_quantiles[[bb]][[dd]][[oo]]) <- as.character(quantiles_i_want)
        }
    }
}

## make a preliminary plot of the hindcasts relative to proxy data

todo

source("getData.R")
source("constraints.R")

# pick which simulation set to display
bb <- "30"
dd <- "ct"

#pdf(paste(plot.dir,'model_ensemble_vs_obspts_logscale_2errbars+Royer14.pdf',sep=''),width=4,height=3,colormodel='cmyk', pointsize=11)

par(mfrow=c(2,1), mai=c(1,1,.15,.15))
plot(-time, log10(model_quantiles$`50`$ct$co2[,"0.5"]), type='l', xlim=c(-450,0), ylim=c(0,log10(40000)), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
polygon(-c(time,rev(time)), log10(c(model_quantiles[[bb]][[dd]]$co2[,"0.025"],rev(model_quantiles[[bb]][[dd]]$co2[,"0.975"]))), col=rgb(.6,.2,.6,.25), border=NA)
polygon(-c(time,rev(time)), log10(c(model_quantiles[[bb]][[dd]]$co2[,"0.05"],rev(model_quantiles[[bb]][[dd]]$co2[,"0.95"]))), col=rgb(.6,.2,.6,.45), border=NA)
polygon(-c(time,rev(time)), log10(c(model_quantiles[[bb]][[dd]]$co2[,"0.25"],rev(model_quantiles[[bb]][[dd]]$co2[,"0.75"]))), col=rgb(.6,.2,.6,.65), border=NA)
polygon(-c(time,rev(time)), log10(c(windows$co2[,"high"],rev(windows$co2[,"low"]))), col=rgb(.5,.5,.5,.5), border=NA)
lines(-time, log10(model_quantiles$`50`$ct$co2[,"0.5"]), lwd=1, lty=5, col=rgb(.6,.2,.6))
mtext('Time [Myr ago]', side=1, line=2.1)
mtext(expression('CO'[2]*' concentration [ppmv]'), side=2, line=3.2)
axis(1, at=seq(-400,0,100), labels=c('400','300','200','100','0'))
ticks=log10(c(seq(1,10,1),seq(10,100,10),seq(200,1000,100),seq(2000,10000,1000),seq(20000,100000,10000)))
axis(2, at=ticks, labels=rep('',length(ticks)))
axis(2, at=log10(c(1,10,100,1000,10000)), labels=c('1','10','100','1000','10000'), las=1)
#axis(2, at=log10(c(3,10,30,100,300,1000,3000,10000)), labels=c('3','10','30','100','300','1000','3000','10000'), las=1)
legend(-452, log10(7), c('This work','Foster et al [2017]'), pch=c(15,15), col=c(rgb(.6,.2,.6,.5),rgb(.5,.5,.5,.5)), pt.cex=1.2, cex=.6, bty='n', y.intersp=1.6)
legend(-452, log10(7), c('This work','Foster et al [2017]'), pch=c('-',''), col=c(rgb(.6,.2,.6,.5),rgb(.5,.5,.5,.5)), cex=.6, bty='n', y.intersp=1.6)

plot(-time, model_quantiles$`50`$ct$temp[,"0.5"], type='l', xlim=c(-450,0), ylim=c(0,50), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
polygon(-c(time,rev(time)), c(model_quantiles[[bb]][[dd]]$temp[,"0.025"],rev(model_quantiles[[bb]][[dd]]$temp[,"0.975"])), col=rgb(.6,.2,.6,.25), border=NA)
polygon(-c(time,rev(time)), c(model_quantiles[[bb]][[dd]]$temp[,"0.05"],rev(model_quantiles[[bb]][[dd]]$temp[,"0.95"])), col=rgb(.6,.2,.6,.45), border=NA)
polygon(-c(time,rev(time)), c(model_quantiles[[bb]][[dd]]$temp[,"0.25"],rev(model_quantiles[[bb]][[dd]]$temp[,"0.75"])), col=rgb(.6,.2,.6,.65), border=NA)
polygon(-c(time,rev(time)), c(windows$temp[,"high"],rev(windows$temp[,"low"])), col=rgb(.5,.5,.5,.5), border=NA)
lines(-time, model_quantiles$`50`$ct$temp[,"0.5"], lwd=1, lty=5, col=rgb(.6,.2,.6))
mtext('Time [Myr ago]', side=1, line=2.1)
mtext('Global average surface\n temperature [deg C]', side=2, line=2.2)
axis(1, at=seq(-400,0,100), labels=c('400','300','200','100','0'))
ticks=seq(from=0, to=60, by=5)
axis(2, at=ticks, labels=rep('',length(ticks)))
axis(2, at=seq(from=0, to=60, by=10), las=1)
legend(-452, 50, c('This work','Mills et al [2019]'), pch=c(15,15), col=c(rgb(.6,.2,.6,.5),rgb(.5,.5,.5,.5)), pt.cex=1.2, cex=.6, bty='n', y.intersp=1.6)
legend(-452, 50, c('This work','Mills et al [2019]'), pch=c('-',''), col=c(rgb(.6,.2,.6,.5),rgb(.5,.5,.5,.5)), cex=.6, bty='n', y.intersp=1.6)

#dev.off()




## compute uncertainty ranges for ESS for each experiment

todo



##==============================================================================
## End
##==============================================================================
