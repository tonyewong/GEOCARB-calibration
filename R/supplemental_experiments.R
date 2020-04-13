##
##
##

source("model_forMCMC_supp.R")
source("run_geocarb_suppF.R")


## run hindcasts

bb <- "30"
dd <- "ct"

model_hindcast_supp <- vector('list', 3)
names(model_hindcast_supp) <- c("co2","temp")
model_hindcast_supp$co2 <- model_hindcast_supp$temp <- mat.or.vec(nr=58, nc=num_samples)
for (ii in 1:num_samples) {
    model_out<- model_forMCMC_supp(par_calib=par_calib[[bb]][[dd]][idx_sample[[bb]][[dd]][ii],],
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
    model_hindcast_supp$co2[,ii] <- model_out[,"co2"]
    model_hindcast_supp$temp[,ii] <- model_out[,"temp"] + par_calib[[bb]][[dd]][idx_sample[[bb]][[dd]][ii],"Ws"]*model_out[,1]/570 + 15 # as in Berner 2004
    # ^-- adding the Ws*t/570 solar luminosity contribution back in, so we have actual temperatures
}

## compute model quantiles for plotting against proxy data, for each experiment

model_quantiles_supp <- vector("list", 2)
names(model_quantiles_supp) <- c("co2","temp")
for (oo in c("co2", "temp")) {
		model_quantiles_supp[[oo]] <- mat.or.vec(nr=58, nc=length(quantiles_i_want))
    for (tt in 1:58) {
    		model_quantiles_supp[[oo]][tt,] <- quantile(model_hindcast_supp[[oo]][tt,], quantiles_i_want)
    }
    colnames(model_quantiles_supp[[oo]]) <- as.character(quantiles_i_want)
}

#cbind(age, model_quantiles_supp$temp[,c("0.025","0.975")], model_quantiles$`30`$ct$temp[,c("0.025","0.975")], model_quantiles_supp$temp[,c("0.025","0.975")]-model_quantiles$`30`$ct$temp[,c("0.025","0.975")])

##
## make a plot of the hindcasts relative to proxy data
##


# pick which simulation set to display
bb <- "30"
dd <- "ct"

# don't want to plot proxy windows for time periods without data
idx_no_data <- which(windows$co2[,"low"]==0 | windows$co2[,"high"]==50000)
idx_data <- setdiff(1:n_time, idx_no_data)
ifirst <- idx_data[1]
idx_data <- c(ifirst-1, idx_data) # start time series 1 earlier for continuity in the figure


pdf('../figures/model_vs_proxy_supp.pdf',width=4,height=6,colormodel='cmyk', pointsize=11)

par(mfrow=c(2,1), mai=c(0.6,.9,.18,.15))
plot(-time, log10(model_quantiles_supp$co2[,"0.5"]), type='l', xlim=c(-425,0), ylim=c(-0.3,log10(40000)), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
grid()
polygon(-c(time[idx_data],rev(time[idx_data])), log10(c(windows$co2[idx_data,"high"],rev(windows$co2[idx_data,"low"]))), col=rgb(.5,.5,.5,.5), border=1, lty=1)
igood <- which(is.finite(model_quantiles[[bb]]$t$co2[,"0.025"])); lines(-time[igood], log10(model_quantiles[[bb]]$t$co2[igood,"0.025"]), lwd=1, lty=5)
igood <- which(is.finite(model_quantiles[[bb]]$t$co2[,"0.975"])); lines(-time[igood], log10(model_quantiles[[bb]]$t$co2[igood,"0.975"]), lwd=1, lty=5)
polygon(-c(time,rev(time)), log10(c(model_quantiles_supp$co2[,"0.05"],rev(model_quantiles_supp$co2[,"0.95"]))), col=rgb(0,.15,.7,.45), border=NA)
polygon(-c(time,rev(time)), log10(c(model_quantiles_supp$co2[,"0.25"],rev(model_quantiles_supp$co2[,"0.75"]))), col=rgb(0,.15,.7,.65), border=NA)
#polygon(-c(time,rev(time)), log10(c(windows$co2[,"high"],rev(windows$co2[,"low"]))), col=rgb(.5,.5,.5,.5), border=NA)
lines(-time, log10(model_quantiles_supp$co2[,"0.5"]), lwd=1, lty=1, col=rgb(0,.15,.7))
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

plot(-time, model_quantiles_supp$temp[,"0.5"], type='l', xlim=c(-425,0), ylim=c(0,51), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
grid()
igood <- which(is.finite(model_quantiles[[bb]]$c$temp[,"0.025"])); lines(-time[igood], model_quantiles[[bb]]$c$temp[igood,"0.025"], lwd=1, lty=5)
igood <- which(is.finite(model_quantiles[[bb]]$c$temp[,"0.975"])); lines(-time[igood], model_quantiles[[bb]]$c$temp[igood,"0.975"], lwd=1, lty=5)
polygon(-c(time,rev(time)), c(windows$temp_sol[,"high"],rev(windows$temp_sol[,"low"])), col=rgb(.5,.5,.5,.5), border=1, lty=1)
polygon(-c(time,rev(time)), c(model_quantiles_supp$temp[,"0.025"],rev(model_quantiles_supp$temp[,"0.975"])), col=rgb(.6,0,0,.25), border=NA)
polygon(-c(time,rev(time)), c(model_quantiles_supp$temp[,"0.05"],rev(model_quantiles_supp$temp[,"0.95"])), col=rgb(.6,0,0,.45), border=NA)
polygon(-c(time,rev(time)), c(model_quantiles_supp$temp[,"0.25"],rev(model_quantiles_supp$temp[,"0.75"])), col=rgb(.6,0,0,.65), border=NA)
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

dev.off()
