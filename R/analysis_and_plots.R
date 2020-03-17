##==============================================================================
## analysis_and_plots.R
##
## Read results from the Latin Hypercube Sampling experiments.
## Generate figures for manuscript and supplemental material.
## Calculate numbers for analysis and comparison among experiments and previous
## works.
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

library(abind)
library(Hmisc)

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

## model setup
param_choice <- 'all'   # Calibrate all 68 parameters? ("all") or only the 6 from Park and Royer 2011 ("PR2011")
data_choice <- 'F2017'    # Which data set?  PR2011 = Park and Royer (2011), or F2017 = Foster et al (2017)
fSR_choice <- 'DT2019'     # Which fSR time series? ("PR2011", "LENTON", "DT2019")
filename.data <- '../input_data/CO2_Proxy_Foster2017_calib_SN-co2_25Sep2018.csv'
source("model_setup.R")
n_time <- nrow(par_time_center)

num_samples <- 10000
quantiles_i_want <- c(0, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 1)
prc_outbound <- as.character(c(30,35,40,45,50))
data_sets <- c("c","t","ct")
par_calib <- par_quantiles <- par_time <- par_time_quantiles <- idx_sample <- vector("list", length(prc_outbound))
names(par_calib) <- names(par_quantiles) <- names(par_time) <- names(par_time_quantiles) <- names(idx_sample) <- prc_outbound

for (bb in prc_outbound) {
    par_calib[[bb]] <- par_quantiles[[bb]] <- par_time[[bb]] <- par_time_quantiles[[bb]] <- vector("list", length(data_sets))
    names(par_calib[[bb]]) <- names(par_quantiles[[bb]]) <- names(par_time[[bb]]) <- names(par_time_quantiles[[bb]]) <- data_sets
    for (dd in data_sets) {
        if ((bb!="30") | (dd!="ct")) {
            load(paste("../output/lhs_param_",dd,"_out",bb,".RData", sep=""))
            par_calib[[bb]][[dd]] <- cbind(par_calib_save, rep(0,nrow(par_calib_save)))
            colnames(par_calib[[bb]][[dd]]) <- c(colnames(par_calib_save), "deltaT2Xglac")
            par_time[[bb]][[dd]] <- par_time_save
        } else {
            ## need to read the parameters for the 30% outbound and carbon+temperature data
            ## experiment separately because the success rate is quite low
            # first simulation set, with seed=1234*ii=13574 (ii=11 is this simulation set)
            load(paste("../output/lhs_param_ct_out30_seed13574.RData", sep=""))
            par_calib[[bb]][[dd]] <- cbind(par_calib_save, rep(0,nrow(par_calib_save)))
            colnames(par_calib[[bb]][[dd]]) <- c(colnames(par_calib_save), "deltaT2Xglac")
            par_time[[bb]][[dd]] <- par_time_save
            # second simulation set, with seed=ii=11
            load(paste("../output/lhs_param_ct_out30_seed11.RData", sep=""))
            par_calib[[bb]][[dd]] <- rbind(par_calib[[bb]][[dd]], cbind(par_calib_save, rep(0,nrow(par_calib_save))))
            par_time[[bb]][[dd]] <- abind(par_time[[bb]][[dd]], par_time_save, along=2)
            # third simulation set, with seed=2020
            #TODO...
            # trim down to 10,000 (or whatever num_samples is above) to match the other simulations
        }
        ## quantiles for the constant parameters
        par_quantiles[[bb]][[dd]] <- mat.or.vec(nr=ncol(par_calib[[bb]][[dd]])+1, nc=length(quantiles_i_want))
        rownames(par_quantiles[[bb]][[dd]]) <- c(colnames(par_calib[[bb]][[dd]]), "deltaT2Xglac")
        colnames(par_quantiles[[bb]][[dd]]) <- as.character(quantiles_i_want)
        if(nrow(par_calib[[bb]][[dd]]) >= num_samples) {
            idx_sample[[bb]][[dd]] <- sample(1:nrow(par_calib[[bb]][[dd]]), size=num_samples, replace=FALSE)
        } else {
            idx_sample[[bb]][[dd]] <- sample(1:nrow(par_calib[[bb]][[dd]]), size=num_samples, replace=TRUE)
        }
        for (pp in 1:ncol(par_calib[[bb]][[dd]])) {
            par_quantiles[[bb]][[dd]][pp,] <- quantile(par_calib[[bb]][[dd]][idx_sample[[bb]][[dd]],pp], quantiles_i_want)
        }
        par_calib[[bb]][[dd]][,"deltaT2Xglac"] <- par_calib[[bb]][[dd]][,"deltaT2X"]*par_calib[[bb]][[dd]][,"GLAC"]
        par_quantiles[[bb]][[dd]]["deltaT2Xglac",] <- quantile(par_calib[[bb]][[dd]][idx_sample[[bb]][[dd]],"deltaT2Xglac"], quantiles_i_want)
        ## quantiles for the time-varying parameters
        par_time_quantiles[[bb]][[dd]] <- array(NA, dim=c(dim(par_time[[bb]][[dd]])[1],length(quantiles_i_want),dim(par_time[[bb]][[dd]])[3]), dimnames=list(1:n_time, quantiles_i_want, parnames_time))
        for (pp in 1:dim(par_time[[bb]][[dd]])[3]) {
            for (tt in 1:dim(par_time[[bb]][[dd]])[1]) {
                par_time_quantiles[[bb]][[dd]][tt,,pp] <- quantile(par_time[[bb]][[dd]][tt,,pp], quantiles_i_want)
            }
        }
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
            model_hindcast[[bb]][[dd]]$temp[,ii] <- model_out[,"temp"] + par_calib[[bb]][[dd]][idx_sample[[bb]][[dd]][ii],"Ws"]*model_out[,1]/570 + 15 # as in Berner 2004
            # ^-- adding the Ws*t/570 solar luminosity contribution back in, so we have actual temperatures
        }
    }
}
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
            		model_quantiles[[bb]][[dd]][[oo]][tt,] <- quantile(model_hindcast[[bb]][[dd]][[oo]][tt,], quantiles_i_want)
            }
            colnames(model_quantiles[[bb]][[dd]][[oo]]) <- as.character(quantiles_i_want)
        }
    }
}



##==============================================================================
##
## make a plot of the hindcasts relative to proxy data
##

source("getData.R")
source("constraints.R")
windows$temp_sol <- windows$temp + array(rep(7.4*time/570,2), dim=c(58,2))

# pick which simulation set to display
bb <- "30"
dd <- "ct"

pdf('../figures/model_vs_proxy.pdf',width=4,height=6,colormodel='cmyk', pointsize=11)

par(mfrow=c(2,1), mai=c(0.6,.9,.15,.15))
plot(-time, log10(model_quantiles[[bb]][[dd]]$co2[,"0.5"]), type='l', xlim=c(-450,0), ylim=c(0,log10(40000)), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
grid()
polygon(-c(time,rev(time)), log10(c(model_quantiles[[bb]][[dd]]$co2[,"0.025"],rev(model_quantiles[[bb]][[dd]]$co2[,"0.975"]))), col=rgb(0,0,.6,.25), border=NA)
polygon(-c(time,rev(time)), log10(c(model_quantiles[[bb]][[dd]]$co2[,"0.05"],rev(model_quantiles[[bb]][[dd]]$co2[,"0.95"]))), col=rgb(0,0,.6,.45), border=NA)
polygon(-c(time,rev(time)), log10(c(model_quantiles[[bb]][[dd]]$co2[,"0.25"],rev(model_quantiles[[bb]][[dd]]$co2[,"0.75"]))), col=rgb(0,0,.6,.65), border=NA)
polygon(-c(time,rev(time)), log10(c(windows$co2[,"high"],rev(windows$co2[,"low"]))), col=rgb(.5,.5,.5,.5), border=NA)
lines(-time, log10(model_quantiles[[bb]][[dd]]$co2[,"0.5"]), lwd=1, lty=1, col=rgb(0,0,.6))
mtext('Time [Myr ago]', side=1, line=2.1)
mtext(expression('CO'[2]*' concentration [ppmv]'), side=2, line=3.5)
axis(1, at=seq(-400,0,100), labels=c('400','300','200','100','0'))
ticks=log10(c(seq(1,10,1),seq(10,100,10),seq(200,1000,100),seq(2000,10000,1000),seq(20000,100000,10000)))
axis(2, at=ticks, labels=rep('',length(ticks)))
axis(2, at=log10(c(1,10,100,1000,10000)), labels=c('1','10','100','1000','10000'), las=1)
#axis(2, at=log10(c(3,10,30,100,300,1000,3000,10000)), labels=c('3','10','30','100','300','1000','3000','10000'), las=1)
legend(-450, log10(8), c('This work','Foster et al [2017]'), pch=c(15,15), col=c(rgb(0,0,.6,.5),rgb(.5,.5,.5,.5)), pt.cex=1.2, cex=.6, bty='n', y.intersp=1.6)
legend(-450, log10(8), c('This work','Foster et al [2017]'), pch=c('-',''), col=c(rgb(0,0,.6,.5),rgb(.5,.5,.5,.5)), cex=.6, bty='n', y.intersp=1.6)

plot(-time, model_quantiles[[bb]][[dd]]$temp[,"0.5"], type='l', xlim=c(-450,0), ylim=c(0,50), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
grid()
polygon(-c(time,rev(time)), c(model_quantiles[[bb]][[dd]]$temp[,"0.025"],rev(model_quantiles[[bb]][[dd]]$temp[,"0.975"])), col=rgb(.6,0,0,.25), border=NA)
polygon(-c(time,rev(time)), c(model_quantiles[[bb]][[dd]]$temp[,"0.05"],rev(model_quantiles[[bb]][[dd]]$temp[,"0.95"])), col=rgb(.6,0,0,.45), border=NA)
polygon(-c(time,rev(time)), c(model_quantiles[[bb]][[dd]]$temp[,"0.25"],rev(model_quantiles[[bb]][[dd]]$temp[,"0.75"])), col=rgb(.6,0,0,.65), border=NA)
polygon(-c(time,rev(time)), c(windows$temp_sol[,"high"],rev(windows$temp_sol[,"low"])), col=rgb(.5,.5,.5,.5), border=NA)
lines(-time, model_quantiles[[bb]][[dd]]$temp[,"0.5"], lwd=1, lty=1, col=rgb(.6,0,0))
mtext('Time [Myr ago]', side=1, line=2.1)
mtext(expression("        Global average\nsurface temperature ["*degree*"C]"), side=2, line=2.2)
axis(1, at=seq(-400,0,100), labels=c('400','300','200','100','0'))
ticks=seq(from=0, to=60, by=5)
axis(2, at=ticks, labels=rep('',length(ticks)))
axis(2, at=seq(from=0, to=60, by=10), las=1)
legend(-450, 50, c('This work','Mills et al [2019]'), pch=c(15,15), col=c(rgb(.6,0,0,.5),rgb(.5,.5,.5,.5)), pt.cex=1.2, cex=.6, bty='n', y.intersp=1.6)
legend(-450, 50, c('This work','Mills et al [2019]'), pch=c('-',''), col=c(rgb(.6,0,0,.5),rgb(.5,.5,.5,.5)), cex=.6, bty='n', y.intersp=1.6)

dev.off()
##==============================================================================



##==============================================================================
##
## compute uncertainty ranges for ESS for each experiment
##

# density estimates for the distributions of deltaT2X

# experiment you want
bb <- "30"
dd <- "ct"

# density of deltaT2X from this work
parnames_calib <- colnames(par_calib[[bb]][[dd]])
ics <- match('deltaT2X', parnames_calib)
deltaT2X_density <- density(par_calib[[bb]][[dd]][,ics], from=0, to=10)

iglac <- match('GLAC', parnames_calib)
glac_density <- density(par_calib[[bb]][[dd]][,iglac], from=1, to=5)

icsg <- match('deltaT2Xglac', parnames_calib)
deltaT2Xglac_density <- density(par_calib[[bb]][[dd]][,icsg], from=0, to=25)


pr2011_dat <- read.csv('../input_data/ParkRoyer2011_Fig3_85varred.csv')
pr2011_cdf <- approxfun(pr2011_dat[,1], pr2011_dat[,4])
pr2011_icdf <- approxfun(pr2011_dat[,4], pr2011_dat[,1])
pr2011_pdf <- approxfun(pr2011_dat[,1], pr2011_dat[,3])

deltaT2X_density_pr2011 <- vector('list', 2); names(deltaT2X_density_pr2011) <- c('x','y')
deltaT2X_density_pr2011$x <- deltaT2X_density$x
deltaT2X_density_pr2011$y <- pr2011_pdf(deltaT2X_density_pr2011$x)

# Park and Royer 2011 have about 16% probability above deltaT2X = 6 deg C
print(1-pr2011_cdf(6))
print(1-pr2011_cdf(7))

# get priors too
row_num <- match('deltaT2X',input$parameter)
x_cs <- seq(from=0, to=10, by=0.1)
f_cs <- dlnorm(x=x_cs, meanlog=log(input[row_num,"mean"]), sdlog=log(0.5*input[row_num,"two_sigma"]))
row_num <- match('GLAC',input$parameter)
x_gl <- seq(from=0, to=10, by=0.1)
f_gl <- dnorm(x=x_gl, mean=input[row_num,"mean"], sd=(0.5*input[row_num,"two_sigma"]))

# Royer et al 2007:  1.5 and 6.2 deg C (5–95% range), 2.8 best fit
x_5_95_royer2007 <- c(1.6, 2.8, 5.5)
# Park and Royer 2011 from CSV/Excel table
x_5_95_pr2011 <- pr2011_icdf(c(.05,.5,.95))
x_5_95_thisstudy <- quantile(par_calib[[bb]][[dd]][,ics], c(.05,.5,.95))  #
x_5_95_glac <- quantile(par_calib[[bb]][[dd]][,ics]*par_calib[[bb]][[dd]][,iglac], c(.05,.5,.95))  #
x_ktc2017 <- c(3.7, 5.6, 7.5)

##
## figure comparing the 30%-outbound ct experiment
##

offset <- 0.1

pdf('../figures/deltaT2X_distributions.pdf',width=4,height=3, colormodel='cmyk', pointsize=11)
par(mfrow=c(1,1), mai=c(.7,.3,.13,.15))
plot(deltaT2X_density$x, deltaT2X_density$y + offset, type='l', lwd=1.7, xlim=c(0.9,10.1), ylim=c(0,.9+offset),
     xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n', axes=FALSE, col="steelblue")
#polygon(-c(time,rev(time)), c(model_quantiles[[bb]][[dd]][,'q025'],rev(model_quantiles[,'q975'])), col='aquamarine1', border=NA)
#polygon(-c(time,rev(time)), c(model_quantiles[,'q05'],rev(model_quantiles[,'q95'])), col='aquamarine3', border=NA)
#lines(deltaT2X_density_nm$x, deltaT2X_density_nm$y + offset, lwd=2, lty=3)
lines(deltaT2X_density_pr2011$x, deltaT2X_density_pr2011$y + offset, lwd=1.7, lty=2)
lines(x_cs, f_cs + offset, lwd=1.7, lty=3, col="steelblue")
mtext(expression(Delta*"T(2x) ["*degree*"C]"), side=1, line=2.3)
mtext('Density', side=2, line=0.3)
arrows(1, 0, 1, .85+offset, length=0.08, angle=30, code=2)
axis(1, at=seq(0,10))
minor.tick(nx=4, ny=0, tick.ratio=0.5)
y0 <- 0.7*offset; arrows(x_5_95_thisstudy[1], y0, x_5_95_thisstudy[3], y0, lwd=1.5, length=0.04, angle=90, code=3, col="steelblue"); points(x_5_95_thisstudy[2], y0, pch=16, col="steelblue")
#y1 <- 0.35*offset; arrows(x_5_95_royer2007[1], y1, x_5_95_royer2007[3], y1, lwd=1.5, length=0.04, angle=90, code=3); points(x_5_95_royer2007[2], y1, pch=15)
y1 <- 0.3*offset; arrows(x_5_95_pr2011[1], y1, x_5_95_pr2011[3], y1, lwd=1.5, length=0.04, angle=90, code=3); points(x_5_95_pr2011[2], y1, pch=15)
#y2 <- 0.08; arrows(x_5_95_ktc2017[1], y2, x_5_95_ktc2017[3], y2, length=0.04, angle=90, code=3); points(x_5_95_ktc2017[2], y2, pch=17)
legend(5.1,1.02, c('5-95% range, PR2011','PR2011','5-95% range, this study','a posteriori, this study','a priori, both studies'),
       pch=c(15,NA,16,NA,NA), lty=c(1,2,1,1,3), col=c("black","black","steelblue","steelblue","steelblue"), bty='n', lwd=1.7, cex=0.9)
dev.off()
##==============================================================================



##==============================================================================
##
## figure comparing the quantiles across the different experiments
##

y0 <- 0
ylims <- c(0,10)
width <- 0.5
offset <- 0.15
colors <- vector("list", 3); names(colors) <- data_sets
alphas <- c(0.45, 0.65)
colors$c <- c(rgb(0,0,0.5,alphas[1]), rgb(0,0,0.9,alphas[2]))
colors$t <- c(rgb(0.5,0,0,alphas[1]), rgb(0.8,0,0,alphas[2]))
colors$ct <- c(rgb(0.5,0,0.5,alphas[1]), rgb(0.85,0,0.85,alphas[2]))

bb <- "30"
pdf('../figures/deltaT2X_experiments.pdf',width=3,height=4, colormodel='cmyk', pointsize=11)
par(mfrow=c(1,1), mai=c(.7,.3,.13,.15))
plot(-10,-10, xlim=c(1.5,5.5), ylim=c(0+0.5*offset,10+0.5*offset), yaxt='n', ylab='', xlab='')
for (xx in seq(0,10)) {lines(rep(xx,2), 2*ylims, col=rgb(.65,.65,.65,1), lty=3, lwd=0.5)}
for (dd in rev(data_sets)) {
    for (bb in prc_outbound) {
        tmp <- par_quantiles[[bb]][[dd]]["deltaT2X",c("0.25","0.75","0.05","0.95","0.5")]
        polygon(c(tmp[3:4], rev(tmp[3:4])), c(rep(y0,2),rep(y0+width,2)), lwd=1, col=colors[[dd]][1])
        polygon(c(tmp[1:2], rev(tmp[1:2])), c(rep(y0,2),rep(y0+width,2)), lwd=1, col=colors[[dd]][2])
        lines(c(tmp[5],tmp[5]), c(y0,y0+width), type='l', lwd=2)
        if(dd=="ct") {text(5.3,y0+0.5*width,paste(bb,"%",sep=""))}
        y0 <- y0+width+offset
    }
    y0 <- y0+2*offset
}
mtext(expression(Delta*"T(2x) ["*degree*"C]"), side=1, line=2.3)
minor.tick(ny=0)
mtext(expression('CO'[2]), side=2, line=.3, adj=0.87)
mtext("T", side=2, line=.5, adj=0.5)
mtext(expression("CO"[2]*" & T"), side=2, line=.3, adj=0.1)
dev.off()
##==============================================================================



##==============================================================================
##
## figure comparing the quantiles as sample size increases
##

y0 <- 0
ylims <- c(0,8.5)
width <- 0.7
offset <- 0.2
alphas <- c(0.45, 0.65)

bb <- "30"
dd <- "ct"
pdf('../figures/deltaT2X_samplesize.pdf',width=3,height=3, colormodel='cmyk', pointsize=11)
par(mfrow=c(1,1), mai=c(.7, .85, .15, .15))
plot(-10,-10, xlim=c(1.5,5.5), ylim=c(ylims[1]+1.0*offset,ylims[2]+0.5*offset), yaxt='n', ylab='', xlab='')
for (xx in seq(0,10)) {lines(rep(xx,2), 2*ylims, col=rgb(.65,.65,.65,1), lty=3, lwd=0.5)}
nn_test <- seq(1000,10000,1000)
yvals <- c()
for (nn in nn_test) {
    tmp_sample <- par_calib[[bb]][[dd]][1:nn,"deltaT2X"]
    tmp <- quantile(tmp_sample, c(0.25,0.75,0.05,0.95,0.5))
    polygon(c(tmp[3:4], rev(tmp[3:4])), c(rep(y0,2),rep(y0+width,2)), lwd=1, col=rgb(.5,.5,.5,alphas[1]))
    polygon(c(tmp[1:2], rev(tmp[1:2])), c(rep(y0,2),rep(y0+width,2)), lwd=1, col=rgb(.5,.5,.5,alphas[2]))
    lines(c(tmp[5],tmp[5]), c(y0,y0+width), type='l', lwd=2)
    yvals <- c(yvals, y0+0.5*width)
    y0 <- y0+width+offset
}
mtext(expression(Delta*"T(2x) ["*degree*"C]"), side=1, line=2.3)
mtext("Sample size", side=2, line=3.3)
axis(side=2, at=yvals, labels=nn_test[1:length(yvals)], las=1)
minor.tick(ny=0)
dev.off()
##==============================================================================



##==============================================================================
##
## time series plots - priors vs posteriors
##

library(CholWishart)
library(MASS)

# get quantiles for the time series parameters
# generate samples from multivariate normal sampling with diagonal covariance

n_parameters_time <- dim(par_time[[bb]][[dd]])[3]
time_series_from_priors <- array(dim=c(n_time,num_samples,n_parameters_time), dimnames=list(1:n_time, 1:num_samples, parnames_time))
time_series_from_priors_quantiles <- array(dim=c(n_time,length(quantiles_i_want),n_parameters_time), dimnames=list(1:n_time, quantiles_i_want, parnames_time))
#covariance_from_priors <- array(dim=c(n_time,n_time,num_samples,n_parameters_time))
covariance_from_priors <- array(dim=c(n_time, n_time))
source("time_series_df.R")

for (pp in 1:n_parameters_time) {
  # these degrees of freedom give variances in line with the reduced variances
  # of Royer et al (2014, AJS), and set the sampled mean covariance on the
  # diagonal matrix of variances (df-(n_time+1))
  #covariance_from_priors[,,,pp] <- rInvWishart(num_samples, df[pp], (df[pp]-(n_time+1))*diag(par_time_stdev[,pp]^2))
  # now draw the actual time series
  for (ii in 1:num_samples) {
    # draw the covariance matrices JIT
    covariance_from_priors <- rInvWishart(1, df[pp], (df[pp]-(n_time+1))*diag(par_time_stdev[,pp]^2))
    #time_series_from_priors[,ii,pp] <- mvrnorm(n=1, mu=par_time_center[,pp], Sigma=covariance_from_priors[,,ii,pp])
    time_series_from_priors[,ii,pp] <- mvrnorm(n=1, mu=par_time_center[,pp], Sigma=covariance_from_priors[,,1])
    # normalize the ones that must be normalized (first batch with 1 as present, second with 0 as present)
    if (pp %in% c("fR", "fL", "fA", "fAw_fA", "fC")) {
      time_series_from_priors[,ii,pp] <- time_series_from_priors[,ii,pp]/time_series_from_priors[n_time,ii,pp]
    } else if (pp %in% c("GEOG")) {
      time_series_from_priors[,ii,pp] <- time_series_from_priors[,ii,pp] - time_series_from_priors[n_time,ii,pp]
    }
  }
  for (tt in 1:n_time) {
    time_series_from_priors_quantiles[tt,,pp] <- quantile(time_series_from_priors[tt,,pp], quantiles_i_want)
  }
}

source("time_series_plot.R")

pdf('../figures/time_series_parameters.pdf',width=6,height=8,colormodel='cmyk', pointsize=11)
par(mfrow=c(4,3), mai=c(0.6,.65,.15,.15))

pp <- 1; ylims <- c(60,105); dy <- 5; units <- "[units]"; time_series_plot()
legend(-450, 106, c("a priori","a posteriori"), pch=c(15,15), col=c(rgb(.5,.5,.5,.5), rgb(0,0,.6,.5)), pt.cex=1.2, cex=1.1, bty='n')
pp <- 2; ylims <- c(-2,6); dy <- 1; units <- "[units]"; time_series_plot()
pp <- 3; ylims <- c(10,35); dy <- 5; units <- "[units]"; time_series_plot()
pp <- 4; ylims <- c(0,1.5); dy <- 0.25; units <- "[units]"; time_series_plot()
pp <- 5; ylims <- c(0,3); dy <- 0.5; units <- "[units]"; time_series_plot()
pp <- 6; ylims <- c(0,2); dy <- 0.5; units <- "[units]"; time_series_plot()
pp <- 7; ylims <- c(0,2); dy <- 0.5; units <- "[units]"; time_series_plot()
pp <- 8; ylims <- c(0,2); dy <- 0.5; units <- "[units]"; time_series_plot()
pp <- 9; ylims <- c(0,.12); dy <- 0.02; units <- "[units]"; time_series_plot()
pp <- 10; ylims <- c(-6,8); dy <- 2; units <- "[units]"; time_series_plot()
pp <- 11; ylims <- c(0,4); dy <- 1; units <- "[units]"; time_series_plot()
pp <- 12; ylims <- c(0,2); dy <- 0.5; units <- "[units]"; time_series_plot()

dev.off()

##==============================================================================



##==============================================================================
##
## what variable is correlated to temperature at age = 100 Myr (or thereabouts)
##

age_of_interest <- c(90,100,110)
idx_of_interest <- match(age_of_interest, time)

bb <- "40"
dd <- "t"
mat <- cbind(par_calib[[bb]][[dd]][idx_sample[[bb]][[dd]],], t(model_hindcast[[bb]][[dd]]$temp[idx_of_interest,]))

res2 <-cor.test(mat,  method = "spearman")

ii <- 10
jj <- 1+jj
res1 <- cor.test(mat[,ii], mat[,jj],  method = "spearman")
plot(mat[,ii], mat[,jj], main=paste(round(res1$estimate,4),"   ",round(res1$p.value,4)), xlab=ii, ylab=jj)


todo

##==============================================================================


##==============================================================================
## End
##==============================================================================
