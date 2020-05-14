##==============================================================================
## figure1.R
##
## Plotting for Figure 1 in manuscript
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================



##==============================================================================
## MS FIGURE 1 -- past ESS estimates, and ours

#input <- read.csv("../input_data/Plot_ESS_Cenozoic_input.csv", header = TRUE)
#input <- read.csv("../input_data/ESS_SupplementaryTable.csv", header = TRUE)
input <- read.csv("../input_data/ESS_SupplementaryTable_TWnotes.csv", header = TRUE)

# from this study, to add to the previous ones in the input file
bb <- '30'
dd <- 'ct'
this_study_g <- quantile(par_calib[[bb]][[dd]][,"deltaT2Xglac"], c(.5,.16,.84))
this_study_ng <- quantile(par_calib[[bb]][[dd]][,"deltaT2X"], c(.5,.16,.84))

# Park and Royer 2011 glacial and non-glacial
pr2011_g <- c(7,6,8)
pr2011_ng <- c(3.778407, 3.778407-1.450958, 3.778407+2.181385)
# glacial period years
x_g1 <- c(260, 340)
x_g2 <- c(0, 40)

# separate pCO2 and method
age_Ma <- input[,1]
age_high <- input[,3]
age_low <- input[,2]
ESS <- input[,4]
ESS_high <- input[,5]
ESS_low <- input[,6]


fig1_names <- c("This work",
                "Park and Royer (2011)",
                "Krissansen-Totton and Catling (2017)",
                "Cramwinckel et al. (2018)",
                "Anagnostou et al. (2016)",
                "Rohling et al. (2012)",
                "Bijl et al. (2010)",
                "Lunt et al. (2010)",
                "Hargreaves and Annan (2016)",
                "Martínez-Botí et al. (2015)",
                "Cramwinckel et al. (2018)",
                "Haywood et al. (2013)",
                "Pagani et al. (2010)",
                "Shaffer et al. (2016) (PETM)",
                "Shaffer et al. (2016) (pre-PETM)",
                "Knobb and Schaller (2017)",
                "Farnsworth et al. (2019)",
                ""," (* denotes no confidence level given)")

fig1_lty <- c(NA, NA, NA, NA, NA, NA, NA, NA,
              NA, NA, NA, NA, NA, NA, NA, 1, NA,
              NA, NA )

fig1_pch <- c(15, 15, 15, 15,
              16, 17, 18, NA,
              NA, NA, NA, 15,
              15, 1, 1, NA, NA,
              NA, NA )

fig1_col <- c("firebrick2", rgb(.5,.5,.5,.4), rgb(.87,.87,.05,.4), rgb(.1,.6,.1,.4),
              "red", "darkorange2", "forestgreen", NA,
              NA, NA, NA, 'blueviolet',
              "dodgerblue2", "blue", "red", "magenta", NA,
              NA, "black" )
#
fig1_cex <- c(1.6, 1.6, 1.6, 1.6,
              1, 1, 1, 1,
              1, 1, 1, 1,
              1, 1, 1, 1, 1,
              0.25, 0.5 )


# figure
pdf('../figures/ESS_previous_work.pdf', width=7.5, height=3.2, pointsize=11, colormodel='cmyk')

layout(mat = matrix(c(1,2), nrow=1, ncol=2),
       heights = c(1), # Heights of the two rows
       widths = c(62, 38)) # Widths of the two columns

par(mai=c(.7,.7,.3,.1))

# start with plot for Park and Royer 2011
plot(-c(500, x_g1[2]), rep(pr2011_ng[1], 2), type='l', lwd=1.5, lty=2, xlim=c(-420,2), ylim=c(0,15.8), xlab='', ylab='', xaxt='n', yaxt='n', xaxs='i', yaxs='i', axes=FALSE)
grid()
lines(-c(x_g1[1], x_g2[2]), rep(pr2011_ng[1], 2), lwd=1.5, lty=2)
polygon(-c(x_g1[1], x_g2[2], x_g2[2], x_g1[1]), c(rep(pr2011_ng[2], 2), rep(pr2011_ng[3], 2)), col=rgb(.5,.5,.5,.4), border=FALSE)
polygon(-c(500, x_g1[2], x_g1[2], 500), c(rep(pr2011_ng[2], 2), rep(pr2011_ng[3], 2)), col=rgb(.5,.5,.5,.4), border=FALSE)
polygon(c(-x_g2, rev(-x_g2)), c(rep(pr2011_g[2], 2), rep(pr2011_g[3], 2)), col=rgb(.5,.5,.5,.4), border=FALSE)
polygon(c(-x_g1, rev(-x_g1)), c(rep(pr2011_g[2], 2), rep(pr2011_g[3], 2)), col=rgb(.5,.5,.5,.4), border=FALSE)
axis(1, at=seq(-500,0,100), labels=c('500','400','300','200','100','0'))
axis(2, at=seq(0,14,2), las=1)
minor.tick(nx=2, ny=5)
mtext("Time [Myr ago]", side=1, line=2.3)
mtext(expression("ESS ["*degree*"C]"), side=2, line=2.3)
# results from KTC 2017
idx <- 6
polygon(-c(age_high[idx], age_low[idx], age_low[idx], age_high[idx]), c(rep(ESS_high[idx], 2), rep(ESS_low[idx], 2)), col=rgb(.87,.87,.05,.4), border=FALSE)
idx <- 11 # results from Cramwinckel et al. (2018)
polygon(-c(age_low[idx], age_high[idx], age_high[idx], age_low[idx]), c(rep(ESS_low[idx],2), rep(ESS_high[idx],2)), col=rgb(.1,.6,.1,.4), border=FALSE)
# add results from this study
lines(-x_g1, rep(this_study_g[1],2), lwd=1.5, lty=1, col="firebrick2")
lines(-x_g2, rep(this_study_g[1],2), lwd=1.5, lty=1, col="firebrick2")
lines(-c(500, x_g1[2]), rep(this_study_ng[1], 2), lwd=1.5, lty=1, col="firebrick2")
lines(-c(x_g1[1], x_g2[2]), rep(this_study_ng[1], 2), lwd=1.5, lty=1, col="firebrick2")
polygon(c(-x_g2, rev(-x_g2)), c(rep(this_study_g[2], 2), rep(this_study_g[3], 2)), col=rgb(.9,.4,.4,.4), border=FALSE)
polygon(c(-x_g1, rev(-x_g1)), c(rep(this_study_g[2], 2), rep(this_study_g[3], 2)), col=rgb(.9,.4,.4,.4), border=FALSE)
polygon(-c(x_g1[1], x_g2[2], x_g2[2], x_g1[1]), c(rep(this_study_ng[2], 2), rep(this_study_ng[3], 2)), col=rgb(.9,.4,.4,.4), border=FALSE)
polygon(-c(500, x_g1[2], x_g1[2], 500), c(rep(this_study_ng[2], 2), rep(this_study_ng[3], 2)), col=rgb(.9,.4,.4,.4), border=FALSE)
# add designation for glacial periods
polygon(c(-500, 5, 5, -500), c(15,15,17,17), col=rgb(1,1,1,1), border=FALSE)
polygon(c(-x_g1, rev(-x_g1)), c(15.3,15.3,16,16), col=rgb(.1,.9,.96,.4), border=FALSE)
polygon(c(-x_g2, rev(-x_g2)), c(15.3,15.3,16,16), col=rgb(.1,.9,.96,.4), border=FALSE)
mtext("glacial", side=3, line=0.1, adj=0.992, cex=0.65)
mtext("glacial", side=3, line=0.1, adj=.265, cex=0.65)
# add the rest of the data points
idx <- c(18,19,20) # results from Pagani et al 2010
for (ii in idx) {
    arrows(x0=-age_Ma[ii], y0=ESS_low[ii], x1=-age_Ma[ii], y1=ESS_high[ii], length=0.03, angle=90, code=3, col="dodgerblue2", lwd=.9)
    points(-age_Ma[ii], ESS[ii], pch=15, col="dodgerblue2", cex=0.6)
}
idx <- c(14,15,16) # results from Rohling et al 2012
for (ii in idx) {
    arrows(x0=-age_Ma[ii], y0=ESS_low[ii], x1=-age_Ma[ii], y1=ESS_high[ii], length=0.03, angle=90, code=3, col="darkorange2", lwd=.9)
    points(-age_Ma[ii], ESS[ii], pch=17, col="darkorange3", cex=0.6)
}
idx <- c(12,13) # results from Anagnostou et al 2016
for (ii in idx) {
    arrows(x0=-age_Ma[ii], y0=ESS_low[ii], x1=-age_Ma[ii], y1=ESS_high[ii], length=0.03, angle=90, code=3, col="red", lwd=.9)
    points(-age_Ma[ii], ESS[ii], pch=16, col="red", cex=0.6)
}
idx <- 7 # results from Bijl et al. (2010)
arrows(x0=-age_Ma[idx], y0=ESS_low[idx], x1=-age_Ma[idx], y1=ESS_high[idx], length=0.03, angle=90, code=3, col="forestgreen", lwd=.9)
points(-age_Ma[idx], ESS[idx], pch=18, col="forestgreen", cex=0.85)
idx <- 21 # PETM from Shaffer et al 2016
points(-age_Ma[idx], ESS[idx], pch=1, cex=0.6, col='blue', lwd=0.8)
idx <- 22 # pre-PETM from Shaffer et al 2016
points(-age_Ma[idx], ESS[idx], pch=1, cex=0.6, col='red', lwd=0.8)
idx <- 23 # results from Knobbe and Schaller (2017)
arrows(x0=-age_low[idx], y0=ESS[idx], x1=-age_high[idx], y1=ESS[idx], length=0.03, angle=90, code=3, col="magenta", lwd=.9)
idx <- 17 # results from Haywood et al 2013
points(-age_Ma[idx], ESS[idx], pch=15, cex=0.6, col='blueviolet', lwd=0.8)

# HERE NOW

# add labels for periods and studies
# this study
text(-385, 6.3, "This study", srt=0, cex=0.65, col="firebrick3")
# park and royer 2011
text(-360, 1.6, "Park and Royer\n(2011)", srt=0, cex=0.65)
# Permo-carboniferous
text(-300, 9.5, "Permo-\n Carboniferous", srt=90, cex=0.65, adj=0)
# Triassic
text(-220, 5.25, "Triassic", srt=90, cex=0.65, adj=0)
# Eocene
text(-50, 0.2, "Eocene", srt=90, cex=0.65, adj=0)
# PETM
text(-56, 7, "PETM", srt=90, cex=0.65, adj=0)
# Pliocene
text(-5, 1.8, "Pliocene", srt=90, cex=0.65, adj=0)

par(mar=c(.3,.1,.3,.3))
plot(-1, -1, xlim=c(0,2), ylim=c(0,2), col="white", xlab='', ylab='', xaxt='n', yaxt='n', xaxs='i', yaxs='i', axes=FALSE)
legend(0, 2, fig1_names, lty=fig1_lty, pch=fig1_pch, col=fig1_col, pt.cex=fig1_cex, cex=0.75, bty='n')
dev.off()


##==============================================================================
## End
##==============================================================================
