##==============================================================================
## parameterSetup_tvq.R
##
## set up specific to sampling the time-varying parameters by CDF value
##
## Requires: filename.calibinput
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

l_Godderis <- TRUE  #set to "TRUE" to run time arrays of fA, fAw/fA, fD, and GEOG
                    #from Godderis et al, 2012; set to "FALSE" to run standard
                    #time arrays from GEOCARBSULF.
iteration_threshold <- 10 #maximum number of times the convergence equation for
                          #CO2 will iterate before signaling a failed run; in a
                          #test with reasonably-well-constrained input parameters
                          #(similar to simulations presented in Royer et al, 20XX),
                          #the number of iterations never exceeded 7; this variable
                          #is here mostly as a failsafe stop-gap.

# reading in the two input files

# input constant parameters and time-varying parameters (comes from calib driver)
input <- read.csv(filename.calibinput)

# time-varying parameters
#y=young; a=old; p=pyrite; s=sulfate; c=carbonate; si=silicates;
# g=organic matter; b=burial; m=degassing; w=weathering
#masses are in units of 10^18 mol
#fluxes ("F" prefix) are in units of 10^18 mol Myrs-1
#rates ("k" prefix) are in units of Myrs-1
#stable isotopic compositions ("d" prefix) are in per mil units
if(USE_LENTON_FSR) {
  time_arrays <- read.csv("../input_data/GEOCARB_input_arrays_Lenton.csv")
} else if (USE_DT2019_FSR) {
  time_arrays <- read.csv("../input_data/GEOCARB_input_arrays_NewDegassing.csv")
} else {
  time_arrays <- read.csv("../input_data/GEOCARB_input_arrays.csv")
}
ageN <- length(time_arrays$age) #number of time-steps
age <- matrix (time_arrays$age, nrow=ageN, ncol=1); colnames(age) <- "age (Myrs ago)"

# if using Godderis et al 2012 values for fA, fAw_FA, fD, and GEOG...
if (l_Godderis) {
  arrays_Godderis <- c('fA','fAw_fA', 'fD', 'GEOG')
  for (i in 1:length(arrays_Godderis)) {
    # reset the array to match Godderis values
	  time_arrays[,arrays_Godderis[i]] <- time_arrays[,paste(arrays_Godderis[i],'Godderis',sep='_')]
	  # also reset the error associated with the array
	  time_arrays[,paste('e',arrays_Godderis[i],sep='')] <- time_arrays[,paste(paste('e',arrays_Godderis[i],sep=''),'Godderis',sep='_')]
  }
}

# normalization of fAw_fA (as per Krause et al, 2018; doi: 10.1038/s41467-018-06383-y)
time_arrays[,"fAw_fA"] <- time_arrays[,"fAw_fA"]/time_arrays[nrow(time_arrays),"fAw_fA"]
time_arrays[,'fAw_fA_Godderis'] <- time_arrays[,'fAw_fA_Godderis']/time_arrays[58,'fAw_fA_Godderis']

ind_time <- which(input$type=='time array')
ind_const <- which(input$type=='constant')
parnames_time <- as.character(input$parameter[ind_time])
parnames_const <- as.character(input$parameter[ind_const])

ind_calib <- which(input$calib==1)
parnames_calib <- c(intersect(ind_calib, ind_time), intersect(ind_calib, ind_const))

# what are the columns within time_arrays of the calibration/fixed guys?
ind_time_calib <- intersect(ind_calib, ind_time)
ind_time_fixed <- setdiff(ind_time, ind_time_calib)

# what are the columns within input of the calibration/fixed guys?
ind_const_calib <- intersect(ind_calib, ind_const)
ind_const_fixed <- setdiff(ind_const, ind_const_calib)

# this needs to be a matrix initially, but for being passed around
# in the calibration, should be a vector (which matrix(..) can unwrap)
ind_time_mean <- NULL
ind_time_err <- NULL
for (i in 1:length(parnames_time)) {
  ind_time_mean <- c(ind_time_mean, match(parnames_time[i], colnames(time_arrays)))
  ind_time_err <- c(ind_time_err, match(paste('e',parnames_time[i],sep=''), colnames(time_arrays)))
}

par_time_center <- time_arrays[,ind_time_mean]
par_time_stdev <- time_arrays[,ind_time_err]
colnames(par_time_stdev) <- colnames(par_time_center)

par_const0 <- input$mean[ind_const]
par_const_fixed0 <- input$mean[ind_const_fixed]
par_const_calib0 <- input$mean[ind_const_calib]
par_time0 <- input$mean[ind_time]
par_time_fixed0 <- input$mean[ind_time_fixed]
par_time_calib0 <- input$mean[ind_time_calib]
par_calib0 <- c(par_const_calib0, par_time_calib0)
par_fixed0 <- c(par_const_fixed0, par_time_fixed0)

parnames_calib <- c( as.character(input$parameter[ind_const_calib]) , as.character(input$parameter[ind_time_calib]))
parnames_fixed <- c( as.character(input$parameter[ind_const_fixed]) , as.character(input$parameter[ind_time_fixed]))

# change the indices fed into the calibration to reflect which values within
# par_calib and par_fixed are constant parameters and which are time-varying
if (length(ind_const_calib)>0) {ind_const_calib <- 1:length(ind_const_calib)
} else                         {ind_const_calib <- NULL}
if (length(ind_const_fixed)>0) {ind_const_fixed <- 1:length(ind_const_fixed)
} else                         {ind_const_fixed <- NULL}

if (length(ind_time_calib)>0) {ind_time_calib <- (length(ind_const_calib)+1):length(par_calib0)
} else                        {ind_time_calib <- NULL}
if (length(ind_time_fixed)>0) {ind_time_fixed <- (length(ind_const_fixed)+1):length(par_fixed0)
} else                        {ind_time_fixed <- NULL}

# what are the orders that run_geocarb.f90 expects?
const_names_expected <- c('ACT',      #(1)
                          'ACTcarb',  #(2)
                          'VNV',      #(3)
                          'NV',       #(4)
                          'exp_NV',   #(5)
                          'LIFE',     #(6)
                          'GYM',      #(7)
                          'FERT',     #(8)
                          'exp_fnBb', #(9)
                          'deltaT2X', #(10)
                          'GLAC',     #(11)
                          'J',        #(12)
                          'n',        #(13)
                          'Ws',       #(14)
                          'exp_fD',   #(15)
                          'Fwpa_0',   #(16)
                          'Fwsa_0',   #(17)
                          'Fwga_0',   #(18)
                          'Fwca_0',   #(19)
                          'Fmg_0',    #(20)
                          'Fmc_0',    #(21)
                          'Fmp_0',    #(22)
                          'Fms_0',    #(23)
                          'Fwsi_0',   #(24)
                          'Xvolc_0',  #(25)
                          'CAPd13C_0',  #(26)
                          'CAPd34S_0',  #(27)
                          'oxygen_570', #(28)
                          'Gy_570',   #(29)
                          'Cy_570',   #(30)
                          'Ca_570',   #(31)
                          'Ssy_570',  #(32)
                          'Spy_570',  #(33)
                          'dlsy_570', #(34)
                          'dlcy_570', #(35)
                          'dlpy_570', #(36)
                          'dlpa_570', #(37)
                          'dlgy_570', #(38)
                          'dlga_570', #(39)
                          'Rcy_570',  #(40)
                          'Rca_570',  #(41)
                          'Rv_570',   #(42)
                          'Rg_570',   #(43)
                          'Fob',      #(44)
                          'COC',      #(45)
                          'Ga',       #(46)
                          'Ssa',      #(47)
                          'Spa',      #(48)
                          'ST',       #(49)
                          'dlst',     #(50)
                          'CT',       #(51)
                          'dlct',     #(52)
                          'kwpy',     #(53)
                          'kwsy',     #(54)
                          'kwgy',     #(55)
                          'kwcy')     #(56)

time_names_expected <- c('Sr',     #(:,1)
                         'd13C',   #(:,2)
                         'd34S',   #(:,3)
                         'fR',     #(:,4)
                         'fL',     #(:,5)
                         'fA',     #(:,6)
                         'fAw_fA', #(:,7)
                         'fD',     #(:,8)
                         'GEOG',   #(:,9)
                         'RT',     #(:,10)
                         'fSR',    #(:,11)
                         'fC')     #(:,12)
n_parameters_time <- length(time_names_expected)

# how do we rearrange Matrix_56 to feed into run_geocarb.f90 as it expects?
const_names_in <- c( parnames_calib[ind_const_calib],
                     parnames_fixed[ind_const_fixed] )
ind_expected_const <- rep(NA, length(const_names_in))
for (pp in 1:length(const_names_in)) {
  ind_expected_const[pp] <- match(const_names_expected[pp], const_names_in)
}

# how do we rearrange Matrix_12 to feed into run_geocarb.f90 as it expects?
time_names_in <- c(unique(parnames_calib[ind_time_calib]), unique(parnames_fixed[ind_time_fixed]))
ind_expected_time <- rep(NA, length(time_names_in))
for (pp in 1:length(ind_expected_time)) {
  ind_expected_time[pp] <- match(time_names_expected[pp], time_names_in)
}

# this is actually all we need, because setting the Matrix_12 time-varyign
# arrays equal to the distribution means initially sets the columns up in a way
# where `ind_time_rearr` below will rearrange them into `time_names_expected`
ind_time_rearr <- c(1,2,3,4,5,6,8,7,10,9,11,12)
ind_expected_time <- ind_time_rearr

##==============================================================================
## End
##==============================================================================
