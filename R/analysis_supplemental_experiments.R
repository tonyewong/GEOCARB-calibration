##==============================================================================
## analysis_supplemental_experiments.R
##
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================

source("model_forMCMC_supp.R")
source("run_geocarb_suppF.R")

##
## read the parameters from the set of supplemental experiments
## -- compare against par_calib$`50`$ct, the control
##

supp_names <- c("control", "gym", "timing", "gym+timing")
par_calib_supp <- par_quantiles_supp <- par_time_supp <- par_time_quantiles_supp <- prcout_supp <- vector("list", length(supp_names))
names(par_calib_supp) <- names(par_quantiles_supp) <- names(par_time_supp) <- names(par_time_quantiles_supp) <- names(prcout_supp) <- supp_names

# control
appen <- supp_names[1]
# -- get parameters
load(paste("../output/lhs_param_ct_out50.RData", sep=""))
par_calib_supp[[appen]] <- cbind(par_calib_save, rep(0,nrow(par_calib_save)))
colnames(par_calib_supp[[appen]]) <- c(colnames(par_calib_save), "deltaT2Xglac")
par_time_supp[[appen]] <- par_time_save
# -- quantiles for the constant parameters
par_quantiles_supp[[appen]] <- mat.or.vec(nr=ncol(par_calib_supp[[appen]])+1, nc=length(quantiles_i_want))
rownames(par_quantiles_supp[[appen]]) <- c(colnames(par_calib_supp[[appen]]), "deltaT2Xglac")
colnames(par_quantiles_supp[[appen]]) <- as.character(quantiles_i_want)
for (pp in 1:ncol(par_calib_supp[[appen]])) {
    par_quantiles_supp[[appen]][pp,] <- quantile(par_calib_supp[[appen]][,pp], quantiles_i_want)
}
par_calib_supp[[appen]][,"deltaT2Xglac"] <- par_calib_supp[[appen]][,"deltaT2X"]*par_calib_supp[[appen]][,"GLAC"]
par_quantiles_supp[[appen]]["deltaT2Xglac",] <- quantile(par_calib_supp[[appen]][,"deltaT2Xglac"], quantiles_i_want)
# -- quantiles for the time-varying parameters
par_time_quantiles_supp[[appen]] <- array(NA, dim=c(dim(par_time_supp[[appen]])[1],length(quantiles_i_want),dim(par_time_supp[[appen]])[3]), dimnames=list(1:n_time, quantiles_i_want, parnames_time))
for (pp in 1:dim(par_time_supp[[appen]])[3]) {
    for (tt in 1:dim(par_time_supp[[appen]])[1]) {
        par_time_quantiles_supp[[appen]][tt,,pp] <- quantile(par_time_supp[[appen]][tt,,pp], quantiles_i_want)
    }
}

# GYM-only
appen <- supp_names[2]
# -- get parameters
load(paste("../output/lhs_param_ct_out50_",appen,".RData", sep=""))
par_calib_supp[[appen]] <- cbind(par_calib_save, rep(0,nrow(par_calib_save)))
colnames(par_calib_supp[[appen]]) <- c(colnames(par_calib_save), "deltaT2Xglac")
par_time_supp[[appen]] <- par_time_save
# -- quantiles for the constant parameters
par_quantiles_supp[[appen]] <- mat.or.vec(nr=ncol(par_calib_supp[[appen]])+1, nc=length(quantiles_i_want))
rownames(par_quantiles_supp[[appen]]) <- c(colnames(par_calib_supp[[appen]]), "deltaT2Xglac")
colnames(par_quantiles_supp[[appen]]) <- as.character(quantiles_i_want)
for (pp in 1:ncol(par_calib_supp[[appen]])) {
    par_quantiles_supp[[appen]][pp,] <- quantile(par_calib_supp[[appen]][,pp], quantiles_i_want)
}
par_calib_supp[[appen]][,"deltaT2Xglac"] <- par_calib_supp[[appen]][,"deltaT2X"]*par_calib_supp[[appen]][,"GLAC"]
par_quantiles_supp[[appen]]["deltaT2Xglac",] <- quantile(par_calib_supp[[appen]][,"deltaT2Xglac"], quantiles_i_want)
# -- quantiles for the time-varying parameters
par_time_quantiles_supp[[appen]] <- array(NA, dim=c(dim(par_time_supp[[appen]])[1],length(quantiles_i_want),dim(par_time_supp[[appen]])[3]), dimnames=list(1:n_time, quantiles_i_want, parnames_time))
for (pp in 1:dim(par_time_supp[[appen]])[3]) {
    for (tt in 1:dim(par_time_supp[[appen]])[1]) {
        par_time_quantiles_supp[[appen]][tt,,pp] <- quantile(par_time_supp[[appen]][tt,,pp], quantiles_i_want)
    }
}

# timing-only
appen <- supp_names[3]
# -- get parameters
load(paste("../output/lhs_param_ct_out50_",appen,".RData", sep=""))
par_calib_supp[[appen]] <- cbind(par_calib_save, rep(0,nrow(par_calib_save)))
colnames(par_calib_supp[[appen]]) <- c(colnames(par_calib_save), "deltaT2Xglac")
par_time_supp[[appen]] <- par_time_save
# -- quantiles for the constant parameters
par_quantiles_supp[[appen]] <- mat.or.vec(nr=ncol(par_calib_supp[[appen]])+1, nc=length(quantiles_i_want))
rownames(par_quantiles_supp[[appen]]) <- c(colnames(par_calib_supp[[appen]]), "deltaT2Xglac")
colnames(par_quantiles_supp[[appen]]) <- as.character(quantiles_i_want)
for (pp in 1:ncol(par_calib_supp[[appen]])) {
    par_quantiles_supp[[appen]][pp,] <- quantile(par_calib_supp[[appen]][,pp], quantiles_i_want)
}
par_calib_supp[[appen]][,"deltaT2Xglac"] <- par_calib_supp[[appen]][,"deltaT2X"]*par_calib_supp[[appen]][,"GLAC"]
par_quantiles_supp[[appen]]["deltaT2Xglac",] <- quantile(par_calib_supp[[appen]][,"deltaT2Xglac"], quantiles_i_want)
# -- quantiles for the time-varying parameters
par_time_quantiles_supp[[appen]] <- array(NA, dim=c(dim(par_time_supp[[appen]])[1],length(quantiles_i_want),dim(par_time_supp[[appen]])[3]), dimnames=list(1:n_time, quantiles_i_want, parnames_time))
for (pp in 1:dim(par_time_supp[[appen]])[3]) {
    for (tt in 1:dim(par_time_supp[[appen]])[1]) {
        par_time_quantiles_supp[[appen]][tt,,pp] <- quantile(par_time_supp[[appen]][tt,,pp], quantiles_i_want)
    }
}

# both
appen <- supp_names[4]
# -- get parameters
load(paste("../output/lhs_param_ct_out50_",appen,".RData", sep=""))
par_calib_supp[[appen]] <- cbind(par_calib_save, rep(0,nrow(par_calib_save)))
colnames(par_calib_supp[[appen]]) <- c(colnames(par_calib_save), "deltaT2Xglac")
par_time_supp[[appen]] <- par_time_save
# -- quantiles for the constant parameters
par_quantiles_supp[[appen]] <- mat.or.vec(nr=ncol(par_calib_supp[[appen]])+1, nc=length(quantiles_i_want))
rownames(par_quantiles_supp[[appen]]) <- c(colnames(par_calib_supp[[appen]]), "deltaT2Xglac")
colnames(par_quantiles_supp[[appen]]) <- as.character(quantiles_i_want)
for (pp in 1:ncol(par_calib_supp[[appen]])) {
    par_quantiles_supp[[appen]][pp,] <- quantile(par_calib_supp[[appen]][,pp], quantiles_i_want)
}
par_calib_supp[[appen]][,"deltaT2Xglac"] <- par_calib_supp[[appen]][,"deltaT2X"]*par_calib_supp[[appen]][,"GLAC"]
par_quantiles_supp[[appen]]["deltaT2Xglac",] <- quantile(par_calib_supp[[appen]][,"deltaT2Xglac"], quantiles_i_want)
# -- quantiles for the time-varying parameters
par_time_quantiles_supp[[appen]] <- array(NA, dim=c(dim(par_time_supp[[appen]])[1],length(quantiles_i_want),dim(par_time_supp[[appen]])[3]), dimnames=list(1:n_time, quantiles_i_want, parnames_time))
for (pp in 1:dim(par_time_supp[[appen]])[3]) {
    for (tt in 1:dim(par_time_supp[[appen]])[1]) {
        par_time_quantiles_supp[[appen]][tt,,pp] <- quantile(par_time_supp[[appen]][tt,,pp], quantiles_i_want)
    }
}


##
## run hindcasts
##

model_hindcast_supp <- vector('list', length(supp_names))
names(model_hindcast_supp) <- supp_names

for (ee in supp_names) {
  model_hindcast_supp[[ee]] <- vector('list', 2)
  names(model_hindcast_supp[[ee]]) <- c("co2","temp")
  model_hindcast_supp[[ee]]$co2 <- model_hindcast_supp[[ee]]$temp <- mat.or.vec(nr=58, nc=num_samples)
  if (ee=="control") {supp_experiment_parameters=c(130,80,1)
  } else if (ee=="gym") {supp_experiment_parameters=c(130,80,0.25)
  } else if (ee=="timing") {supp_experiment_parameters=c(110,60,1)
  } else if (ee=="gym+timing") {supp_experiment_parameters=c(110,60,0.25)
  } else {print("ERROR")}
  prcout_supp[[ee]] <- mat.or.vec(nr=num_samples, nc=2)
  colnames(prcout_supp[[ee]]) <- c("co2","temp")
  for (ii in 1:num_samples) {
      model_out <- model_forMCMC_supp(par_calib=par_calib_supp[[ee]][ii,],
                                 par_time=par_time_supp[[ee]][,ii,],
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
                                 iteration_threshold=iteration_threshold,
                                 supp_experiment_parameters=supp_experiment_parameters)
      # compute %outbound
      prcout_supp[[ee]][ii,"co2"] <- percout(model_out[,"co2"], windows$co2)
      prcout_supp[[ee]][ii,"temp"] <- percout(model_out[,"temp"] + 15, windows$temp)
      # save hindcasts
      model_hindcast_supp[[ee]]$co2[,ii] <- model_out[,"co2"]
      model_hindcast_supp[[ee]]$temp[,ii] <- model_out[,"temp"] + par_calib_supp[[ee]][ii,"Ws"]*model_out[,"age"]/570 + 15 # as in Berner 2004
      # ^-- adding the Ws*t/570 solar luminosity contribution back in, so we have actual temperatures
      #     it is subtracted out in run_geocarb.f90 so it is easier to compare with the Mills temperatures
  }
}


##
## compute model quantiles for plotting against proxy data, for each experiment
##

model_quantiles_supp <- vector("list", length(supp_names))
names(model_quantiles_supp) <- supp_names
for (ee in supp_names) {
  model_quantiles_supp[[ee]] <- vector("list", 2)
  names(model_quantiles_supp[[ee]]) <- c("co2","temp")
  for (oo in c("co2", "temp")) {
  		model_quantiles_supp[[ee]][[oo]] <- mat.or.vec(nr=58, nc=length(quantiles_i_want))
      for (tt in 1:58) {
      		model_quantiles_supp[[ee]][[oo]][tt,] <- quantile(model_hindcast_supp[[ee]][[oo]][tt,], quantiles_i_want)
      }
      colnames(model_quantiles_supp[[ee]][[oo]]) <- as.character(quantiles_i_want)
  }
}


##
## compute empirical survival function for model %outbound
##

prcout_supp_sorted <- prcout_supp
for (ee in supp_names) {
  prcout_supp_sorted[[ee]][,"co2"] <- sort(prcout_supp[[ee]][,"co2"])
  prcout_supp_sorted[[ee]][,"temp"] <- sort(prcout_supp[[ee]][,"temp"])
}
ecdf_values <- seq(from=1, to=num_samples, by=1)/(num_samples+1)

# plot of survival functions (not very compelling - not included in manuscript)
plot(prcout_supp_sorted$control[,"temp"], log10(1-ecdf_values), type="l", ylim=c(-2,0))
lines(prcout_supp_sorted$`gym`[,"temp"], log10(1-ecdf_values), col='blue')
lines(prcout_supp_sorted$`timing`[,"temp"], log10(1-ecdf_values), col='red')
lines(prcout_supp_sorted$`gym+timing`[,"temp"], log10(1-ecdf_values), col='purple')


##
## plot showing hindcast comparison
##

TODO

TODO - make usre you use windows$temp_sol






##
##
##
##
##

##
## make a plot of the hindcasts relative to proxy data
##

# don't want to plot proxy windows for time periods without data
idx_no_data <- which(windows$co2[,"low"]==0 | windows$co2[,"high"]==50000)
idx_data <- setdiff(1:n_time, idx_no_data)
ifirst <- idx_data[1]
idx_data <- c(ifirst-1, idx_data) # start time series 1 earlier for continuity in the figure

pdf('../figures/model_vs_proxy_supp.pdf',width=6,height=8,colormodel='cmyk', pointsize=11)
par(mfrow=c(4,2), mai=c(0.6,.9,.18,.15))

## Control (ct prcout=50)
ee <- 1
# -- co2
plot(-time, log10(model_quantiles_supp[[ee]]$co2[,"0.5"]), type='l', xlim=c(-425,0), ylim=c(-0.3,log10(40000)), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
grid()
polygon(-c(time[idx_data],rev(time[idx_data])), log10(c(windows$co2[idx_data,"high"],rev(windows$co2[idx_data,"low"]))), col=rgb(.5,.5,.5,.5), border=1, lty=1)
#igood <- which(is.finite(model_quantiles[[bb]]$t$co2[,"0.025"])); lines(-time[igood], log10(model_quantiles[[bb]]$t$co2[igood,"0.025"]), lwd=1, lty=5)
#igood <- which(is.finite(model_quantiles[[bb]]$t$co2[,"0.975"])); lines(-time[igood], log10(model_quantiles[[bb]]$t$co2[igood,"0.975"]), lwd=1, lty=5)
polygon(-c(time,rev(time)), log10(c(model_quantiles_supp[[ee]]$co2[,"0.05"],rev(model_quantiles_supp[[ee]]$co2[,"0.95"]))), col=rgb(0,.15,.7,.45), border=NA)
polygon(-c(time,rev(time)), log10(c(model_quantiles_supp[[ee]]$co2[,"0.25"],rev(model_quantiles_supp[[ee]]$co2[,"0.75"]))), col=rgb(0,.15,.7,.65), border=NA)
#polygon(-c(time,rev(time)), log10(c(windows$co2[,"high"],rev(windows$co2[,"low"]))), col=rgb(.5,.5,.5,.5), border=NA)
lines(-time, log10(model_quantiles_supp[[ee]]$co2[,"0.5"]), lwd=1, lty=1, col=rgb(0,.15,.7))
mtext('Time [Myr ago]', side=1, line=2.1)
mtext(expression('CO'[2]*' concentration [ppmv]'), side=2, line=3.5)
axis(1, at=seq(-400,0,100), labels=c('400','300','200','100','0'))
ticks=log10(c(seq(1,10,1),seq(10,100,10),seq(200,1000,100),seq(2000,10000,1000),seq(20000,100000,10000)))
axis(2, at=ticks, labels=rep('',length(ticks)))
axis(2, at=log10(c(1,10,100,1000,10000)), labels=c('1','10','100','1000','10000'), las=1)
#axis(2, at=log10(c(3,10,30,100,300,1000,3000,10000)), labels=c('3','10','30','100','300','1000','3000','10000'), las=1)
legend(-420, 0.6, c('This work (median and 50%, 90% and 95% ranges)',expression('95% range without CO'[2]*' data'),'Data from Foster et al [2017]'), pch=c(15,NA,15), lty=c(NA,5,NA), col=c(rgb(0,0,.6,.5),'black',rgb(.5,.5,.5,.5)), pt.cex=1.2, cex=.6, bty='n', y.intersp=1)
legend(-420, 0.6, c('This work (median and 50%, 90% and 95% ranges)',expression('95% range without CO'[2]*' data'),'Data from Foster et al [2017]'), pch=c('-','',''), lty=c(NA,5,NA), col=c(rgb(0,0,.6,.5),'black',rgb(.5,.5,.5,.5)), cex=.6, bty='n', y.intersp=1)
mtext(side=3, text=expression(bold('   a')), line=0, cex=1, adj=-0.24);
# -- temperature
plot(-time, model_quantiles_supp[[ee]]$temp[,"0.5"], type='l', xlim=c(-425,0), ylim=c(0,51), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
grid()
#igood <- which(is.finite(model_quantiles[[bb]]$c$temp[,"0.025"])); lines(-time[igood], model_quantiles[[bb]]$c$temp[igood,"0.025"], lwd=1, lty=5)
#igood <- which(is.finite(model_quantiles[[bb]]$c$temp[,"0.975"])); lines(-time[igood], model_quantiles[[bb]]$c$temp[igood,"0.975"], lwd=1, lty=5)
polygon(-c(time,rev(time)), c(windows$temp_sol[,"high"],rev(windows$temp_sol[,"low"])), col=rgb(.5,.5,.5,.5), border=1, lty=1)
polygon(-c(time,rev(time)), c(model_quantiles_supp[[ee]]$temp[,"0.025"],rev(model_quantiles_supp[[ee]]$temp[,"0.975"])), col=rgb(.6,0,0,.25), border=NA)
polygon(-c(time,rev(time)), c(model_quantiles_supp[[ee]]$temp[,"0.05"],rev(model_quantiles_supp[[ee]]$temp[,"0.95"])), col=rgb(.6,0,0,.45), border=NA)
polygon(-c(time,rev(time)), c(model_quantiles_supp[[ee]]$temp[,"0.25"],rev(model_quantiles_supp[[ee]]$temp[,"0.75"])), col=rgb(.6,0,0,.65), border=NA)
lines(-time, model_quantiles_supp$temp[,"0.5"], lwd=1, lty=1, col=rgb(.6,0,0))
mtext('Time [Myr ago]', side=1, line=2.1)
mtext(expression("        Global average\nsurface temperature ["*degree*"C]"), side=2, line=2.2)
axis(1, at=seq(-400,0,100), labels=c('400','300','200','100','0'))
ticks=seq(from=0, to=60, by=5)
axis(2, at=ticks, labels=rep('',length(ticks)))
axis(2, at=seq(from=0, to=60, by=10), las=1)
legend(-420, 51, c('This work (median and 50%, 90% and 95% ranges)','95% range without temperature data','Data from Mills et al [2019]'), pch=c(15,NA,15), lty=c(NA,5,NA), col=c(rgb(.6,0,0,.5),'black',rgb(.5,.5,.5,.5)), pt.cex=1.2, cex=.6, bty='n', y.intersp=1)
legend(-420, 51, c('This work (median and 50%, 90% and 95% ranges)','95% range without temperature data','Data from Mills et al [2019]'), pch=c('-','',''), lty=c(NA,5,NA), col=c(rgb(.6,0,0,.5),'black',rgb(.5,.5,.5,.5)), cex=.6, bty='n', y.intersp=1)
mtext(side=3, text=expression(bold('   b')), line=0, cex=1, adj=-0.24);

## GYM only
ee <- 2
# -- co2
plot(-time, log10(model_quantiles_supp[[ee]]$co2[,"0.5"]), type='l', xlim=c(-425,0), ylim=c(-0.3,log10(40000)), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
grid()
polygon(-c(time[idx_data],rev(time[idx_data])), log10(c(windows$co2[idx_data,"high"],rev(windows$co2[idx_data,"low"]))), col=rgb(.5,.5,.5,.5), border=1, lty=1)
#igood <- which(is.finite(model_quantiles[[bb]]$t$co2[,"0.025"])); lines(-time[igood], log10(model_quantiles[[bb]]$t$co2[igood,"0.025"]), lwd=1, lty=5)
#igood <- which(is.finite(model_quantiles[[bb]]$t$co2[,"0.975"])); lines(-time[igood], log10(model_quantiles[[bb]]$t$co2[igood,"0.975"]), lwd=1, lty=5)
polygon(-c(time,rev(time)), log10(c(model_quantiles_supp[[ee]]$co2[,"0.05"],rev(model_quantiles_supp[[ee]]$co2[,"0.95"]))), col=rgb(0,.15,.7,.45), border=NA)
polygon(-c(time,rev(time)), log10(c(model_quantiles_supp[[ee]]$co2[,"0.25"],rev(model_quantiles_supp[[ee]]$co2[,"0.75"]))), col=rgb(0,.15,.7,.65), border=NA)
#polygon(-c(time,rev(time)), log10(c(windows$co2[,"high"],rev(windows$co2[,"low"]))), col=rgb(.5,.5,.5,.5), border=NA)
lines(-time, log10(model_quantiles_supp[[ee]]$co2[,"0.5"]), lwd=1, lty=1, col=rgb(0,.15,.7))
mtext('Time [Myr ago]', side=1, line=2.1)
mtext(expression('CO'[2]*' concentration [ppmv]'), side=2, line=3.5)
axis(1, at=seq(-400,0,100), labels=c('400','300','200','100','0'))
ticks=log10(c(seq(1,10,1),seq(10,100,10),seq(200,1000,100),seq(2000,10000,1000),seq(20000,100000,10000)))
axis(2, at=ticks, labels=rep('',length(ticks)))
axis(2, at=log10(c(1,10,100,1000,10000)), labels=c('1','10','100','1000','10000'), las=1)
#axis(2, at=log10(c(3,10,30,100,300,1000,3000,10000)), labels=c('3','10','30','100','300','1000','3000','10000'), las=1)
mtext(side=3, text=expression(bold('   c')), line=0, cex=1, adj=-0.24);
# -- temperature
plot(-time, model_quantiles_supp[[ee]]$temp[,"0.5"], type='l', xlim=c(-425,0), ylim=c(0,51), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
grid()
#igood <- which(is.finite(model_quantiles[[bb]]$c$temp[,"0.025"])); lines(-time[igood], model_quantiles[[bb]]$c$temp[igood,"0.025"], lwd=1, lty=5)
#igood <- which(is.finite(model_quantiles[[bb]]$c$temp[,"0.975"])); lines(-time[igood], model_quantiles[[bb]]$c$temp[igood,"0.975"], lwd=1, lty=5)
polygon(-c(time,rev(time)), c(windows$temp_sol[,"high"],rev(windows$temp_sol[,"low"])), col=rgb(.5,.5,.5,.5), border=1, lty=1)
polygon(-c(time,rev(time)), c(model_quantiles_supp[[ee]]$temp[,"0.025"],rev(model_quantiles_supp[[ee]]$temp[,"0.975"])), col=rgb(.6,0,0,.25), border=NA)
polygon(-c(time,rev(time)), c(model_quantiles_supp[[ee]]$temp[,"0.05"],rev(model_quantiles_supp[[ee]]$temp[,"0.95"])), col=rgb(.6,0,0,.45), border=NA)
polygon(-c(time,rev(time)), c(model_quantiles_supp[[ee]]$temp[,"0.25"],rev(model_quantiles_supp[[ee]]$temp[,"0.75"])), col=rgb(.6,0,0,.65), border=NA)
lines(-time, model_quantiles_supp$temp[,"0.5"], lwd=1, lty=1, col=rgb(.6,0,0))
mtext('Time [Myr ago]', side=1, line=2.1)
mtext(expression("        Global average\nsurface temperature ["*degree*"C]"), side=2, line=2.2)
axis(1, at=seq(-400,0,100), labels=c('400','300','200','100','0'))
ticks=seq(from=0, to=60, by=5)
axis(2, at=ticks, labels=rep('',length(ticks)))
axis(2, at=seq(from=0, to=60, by=10), las=1)
mtext(side=3, text=expression(bold('   d')), line=0, cex=1, adj=-0.24);

## timing only
ee <- 3
# -- co2
plot(-time, log10(model_quantiles_supp[[ee]]$co2[,"0.5"]), type='l', xlim=c(-425,0), ylim=c(-0.3,log10(40000)), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
grid()
polygon(-c(time[idx_data],rev(time[idx_data])), log10(c(windows$co2[idx_data,"high"],rev(windows$co2[idx_data,"low"]))), col=rgb(.5,.5,.5,.5), border=1, lty=1)
#igood <- which(is.finite(model_quantiles[[bb]]$t$co2[,"0.025"])); lines(-time[igood], log10(model_quantiles[[bb]]$t$co2[igood,"0.025"]), lwd=1, lty=5)
#igood <- which(is.finite(model_quantiles[[bb]]$t$co2[,"0.975"])); lines(-time[igood], log10(model_quantiles[[bb]]$t$co2[igood,"0.975"]), lwd=1, lty=5)
polygon(-c(time,rev(time)), log10(c(model_quantiles_supp[[ee]]$co2[,"0.05"],rev(model_quantiles_supp[[ee]]$co2[,"0.95"]))), col=rgb(0,.15,.7,.45), border=NA)
polygon(-c(time,rev(time)), log10(c(model_quantiles_supp[[ee]]$co2[,"0.25"],rev(model_quantiles_supp[[ee]]$co2[,"0.75"]))), col=rgb(0,.15,.7,.65), border=NA)
#polygon(-c(time,rev(time)), log10(c(windows$co2[,"high"],rev(windows$co2[,"low"]))), col=rgb(.5,.5,.5,.5), border=NA)
lines(-time, log10(model_quantiles_supp[[ee]]$co2[,"0.5"]), lwd=1, lty=1, col=rgb(0,.15,.7))
mtext('Time [Myr ago]', side=1, line=2.1)
mtext(expression('CO'[2]*' concentration [ppmv]'), side=2, line=3.5)
axis(1, at=seq(-400,0,100), labels=c('400','300','200','100','0'))
ticks=log10(c(seq(1,10,1),seq(10,100,10),seq(200,1000,100),seq(2000,10000,1000),seq(20000,100000,10000)))
axis(2, at=ticks, labels=rep('',length(ticks)))
axis(2, at=log10(c(1,10,100,1000,10000)), labels=c('1','10','100','1000','10000'), las=1)
#axis(2, at=log10(c(3,10,30,100,300,1000,3000,10000)), labels=c('3','10','30','100','300','1000','3000','10000'), las=1)
mtext(side=3, text=expression(bold('   e')), line=0, cex=1, adj=-0.24);
# -- temperature
plot(-time, model_quantiles_supp[[ee]]$temp[,"0.5"], type='l', xlim=c(-425,0), ylim=c(0,51), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
grid()
#igood <- which(is.finite(model_quantiles[[bb]]$c$temp[,"0.025"])); lines(-time[igood], model_quantiles[[bb]]$c$temp[igood,"0.025"], lwd=1, lty=5)
#igood <- which(is.finite(model_quantiles[[bb]]$c$temp[,"0.975"])); lines(-time[igood], model_quantiles[[bb]]$c$temp[igood,"0.975"], lwd=1, lty=5)
polygon(-c(time,rev(time)), c(windows$temp_sol[,"high"],rev(windows$temp_sol[,"low"])), col=rgb(.5,.5,.5,.5), border=1, lty=1)
polygon(-c(time,rev(time)), c(model_quantiles_supp[[ee]]$temp[,"0.025"],rev(model_quantiles_supp[[ee]]$temp[,"0.975"])), col=rgb(.6,0,0,.25), border=NA)
polygon(-c(time,rev(time)), c(model_quantiles_supp[[ee]]$temp[,"0.05"],rev(model_quantiles_supp[[ee]]$temp[,"0.95"])), col=rgb(.6,0,0,.45), border=NA)
polygon(-c(time,rev(time)), c(model_quantiles_supp[[ee]]$temp[,"0.25"],rev(model_quantiles_supp[[ee]]$temp[,"0.75"])), col=rgb(.6,0,0,.65), border=NA)
lines(-time, model_quantiles_supp$temp[,"0.5"], lwd=1, lty=1, col=rgb(.6,0,0))
mtext('Time [Myr ago]', side=1, line=2.1)
mtext(expression("        Global average\nsurface temperature ["*degree*"C]"), side=2, line=2.2)
axis(1, at=seq(-400,0,100), labels=c('400','300','200','100','0'))
ticks=seq(from=0, to=60, by=5)
axis(2, at=ticks, labels=rep('',length(ticks)))
axis(2, at=seq(from=0, to=60, by=10), las=1)
mtext(side=3, text=expression(bold('   f')), line=0, cex=1, adj=-0.24);

## GYM and timing
ee <- 4
# -- co2
plot(-time, log10(model_quantiles_supp[[ee]]$co2[,"0.5"]), type='l', xlim=c(-425,0), ylim=c(-0.3,log10(40000)), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
grid()
polygon(-c(time[idx_data],rev(time[idx_data])), log10(c(windows$co2[idx_data,"high"],rev(windows$co2[idx_data,"low"]))), col=rgb(.5,.5,.5,.5), border=1, lty=1)
#igood <- which(is.finite(model_quantiles[[bb]]$t$co2[,"0.025"])); lines(-time[igood], log10(model_quantiles[[bb]]$t$co2[igood,"0.025"]), lwd=1, lty=5)
#igood <- which(is.finite(model_quantiles[[bb]]$t$co2[,"0.975"])); lines(-time[igood], log10(model_quantiles[[bb]]$t$co2[igood,"0.975"]), lwd=1, lty=5)
polygon(-c(time,rev(time)), log10(c(model_quantiles_supp[[ee]]$co2[,"0.05"],rev(model_quantiles_supp[[ee]]$co2[,"0.95"]))), col=rgb(0,.15,.7,.45), border=NA)
polygon(-c(time,rev(time)), log10(c(model_quantiles_supp[[ee]]$co2[,"0.25"],rev(model_quantiles_supp[[ee]]$co2[,"0.75"]))), col=rgb(0,.15,.7,.65), border=NA)
#polygon(-c(time,rev(time)), log10(c(windows$co2[,"high"],rev(windows$co2[,"low"]))), col=rgb(.5,.5,.5,.5), border=NA)
lines(-time, log10(model_quantiles_supp[[ee]]$co2[,"0.5"]), lwd=1, lty=1, col=rgb(0,.15,.7))
mtext('Time [Myr ago]', side=1, line=2.1)
mtext(expression('CO'[2]*' concentration [ppmv]'), side=2, line=3.5)
axis(1, at=seq(-400,0,100), labels=c('400','300','200','100','0'))
ticks=log10(c(seq(1,10,1),seq(10,100,10),seq(200,1000,100),seq(2000,10000,1000),seq(20000,100000,10000)))
axis(2, at=ticks, labels=rep('',length(ticks)))
axis(2, at=log10(c(1,10,100,1000,10000)), labels=c('1','10','100','1000','10000'), las=1)
#axis(2, at=log10(c(3,10,30,100,300,1000,3000,10000)), labels=c('3','10','30','100','300','1000','3000','10000'), las=1)
mtext(side=3, text=expression(bold('   g')), line=0, cex=1, adj=-0.24);
# -- temperature
plot(-time, model_quantiles_supp[[ee]]$temp[,"0.5"], type='l', xlim=c(-425,0), ylim=c(0,51), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
grid()
#igood <- which(is.finite(model_quantiles[[bb]]$c$temp[,"0.025"])); lines(-time[igood], model_quantiles[[bb]]$c$temp[igood,"0.025"], lwd=1, lty=5)
#igood <- which(is.finite(model_quantiles[[bb]]$c$temp[,"0.975"])); lines(-time[igood], model_quantiles[[bb]]$c$temp[igood,"0.975"], lwd=1, lty=5)
polygon(-c(time,rev(time)), c(windows$temp_sol[,"high"],rev(windows$temp_sol[,"low"])), col=rgb(.5,.5,.5,.5), border=1, lty=1)
polygon(-c(time,rev(time)), c(model_quantiles_supp[[ee]]$temp[,"0.025"],rev(model_quantiles_supp[[ee]]$temp[,"0.975"])), col=rgb(.6,0,0,.25), border=NA)
polygon(-c(time,rev(time)), c(model_quantiles_supp[[ee]]$temp[,"0.05"],rev(model_quantiles_supp[[ee]]$temp[,"0.95"])), col=rgb(.6,0,0,.45), border=NA)
polygon(-c(time,rev(time)), c(model_quantiles_supp[[ee]]$temp[,"0.25"],rev(model_quantiles_supp[[ee]]$temp[,"0.75"])), col=rgb(.6,0,0,.65), border=NA)
lines(-time, model_quantiles_supp$temp[,"0.5"], lwd=1, lty=1, col=rgb(.6,0,0))
mtext('Time [Myr ago]', side=1, line=2.1)
mtext(expression("        Global average\nsurface temperature ["*degree*"C]"), side=2, line=2.2)
axis(1, at=seq(-400,0,100), labels=c('400','300','200','100','0'))
ticks=seq(from=0, to=60, by=5)
axis(2, at=ticks, labels=rep('',length(ticks)))
axis(2, at=seq(from=0, to=60, by=10), las=1)
mtext(side=3, text=expression(bold('   h')), line=0, cex=1, adj=-0.24);

dev.off()

##
## End
##
