##==============================================================================
## lhs_sampling_cret.R
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

idx_cret <- 49

# compute the degrees of freedom for the time series inverse Wishart prior on
# their covariance
source("time_series_df.R")

tbeg <- proc.time()
if (n_sample <= n_sample_per_chunk) {
  # business as usual

  set.seed(2021)
  parameters_lhs <- randomLHS(n_sample, n_parameters)
  par_calib <- parameters_lhs  # initialize

#TODO
print("here for some reason? Not used")
#do simulations... only return the results with at most 50% %outbound

} else {
  # divide the samples into chunks  <-- TODO: could try to distribute the load more evenly if parallelized
  n_chunk <- ceiling(n_sample/n_sample_per_chunk)
  n_sample_this_chunk <- rep(n_sample_per_chunk, n_chunk)
  n_sample_this_chunk[n_chunk] <- n_sample - (n_chunk-1)*n_sample_per_chunk

  # sample the constant parameters from one large LHS and divide
  set.seed(2021)
  parameters_lhs <- randomLHS(n_sample, n_parameters)

  # constant parameter sampling
  n_const_calib <- length(ind_const_calib)
  par_calib_tmp <- parameters_lhs  # initialize
  colnames(par_calib_tmp) <- parnames_const_calib
  for (ii in 1:n_const_calib) {
    row_num <- match(parnames_const_calib[ii],input$parameter)
    parname <- as.character(input[row_num,"parameter"])
    if(input[row_num, 'distribution_type']=='gaussian') {
      # scale the parameters so they won't be outside the bounds for the actual parameter values
      new_bounds <- c(0,1)
      if (!is.na(bounds[parname,"lower"])) {
        new_bounds[1] <- pnorm(q=bounds[parname, "lower"], mean=input[row_num,"mean"], sd=(0.5*input[row_num,"two_sigma"]))
      }
      if (!is.na(bounds[parname,"upper"])) {
        new_bounds[2] <- pnorm(q=bounds[parname, "upper"], mean=input[row_num,"mean"], sd=(0.5*input[row_num,"two_sigma"]))
      }
      # scale parameters
      parameters_lhs[,ind_const_calib[ii]] <- parameters_lhs[,ind_const_calib[ii]]*(new_bounds[2]-new_bounds[1]) + new_bounds[1]
      par_calib_tmp[,ii] <- qnorm(p=parameters_lhs[,ind_const_calib[ii]], mean=input[row_num,"mean"], sd=(0.5*input[row_num,"two_sigma"]))
    } else if(input[row_num, 'distribution_type']=='lognormal') {
      # scale the parameters so they won't be outside the bounds for the actual parameter values
      new_bounds <- c(0,1)
      if (!is.na(bounds[parname,"lower"])) {
        new_bounds[1] <- plnorm(q=bounds[parname, "lower"], meanlog=log(input[row_num,"mean"]), sdlog=log(0.5*input[row_num,"two_sigma"]))
      }
      if (!is.na(bounds[parname,"upper"])) {
        new_bounds[2] <- plnorm(q=bounds[parname, "upper"], meanlog=log(input[row_num,"mean"]), sdlog=log(0.5*input[row_num,"two_sigma"]))
      }
      # scale parameters
      parameters_lhs[,ind_const_calib[ii]] <- parameters_lhs[,ind_const_calib[ii]]*(new_bounds[2]-new_bounds[1]) + new_bounds[1]
      par_calib_tmp[,ii] <- qlnorm(p=parameters_lhs[,ind_const_calib[ii]], meanlog=log(input[row_num,"mean"]), sdlog=log(0.5*input[row_num,"two_sigma"]))
    } else {
      print('ERROR - unknown prior distribution type')
    }
  }

  # divide the scalar parameters across the chunks
  par_calib <- vector('list', n_chunk)
  for (cc in 1:(n_chunk-1)) {
    par_calib[[cc]] <- par_calib_tmp[((cc-1)*n_sample_per_chunk+1):(cc*n_sample_per_chunk),]
  }
  par_calib[[n_chunk]] <- par_calib_tmp[(cc*n_sample_per_chunk+1):n_sample,]
  rm(list=c("par_calib_tmp"))

  # set up arrays to hold the parameters and covariances with low enough %outbound
  par_calib_save <- NULL
  par_time_save <- NULL
  par_covar_save <- NULL

  for (cc in 1:n_chunk) {
    # do simulations... only return the results with at most 50% %outbound
    # within simulation loop, do the time series sampling
    covariance_samples <- array(dim=c(n_time,n_time,n_sample_this_chunk[cc],n_parameters_time))
    time_series_samples <- array(dim=c(n_time,n_sample_this_chunk[cc],n_parameters_time))
    for (pp in 1:n_parameters_time) {
      # these degrees of freedom give variances in line with the reduced variances
      # of Royer et al (2014, AJS), and set the sampled mean covariance on the
      # diagonal matrix of variances (df-(n_time+1))
      covariance_samples[,,,pp] <- rInvWishart(n_sample_this_chunk[cc], df[pp], (df[pp]-(n_time+1))*diag(par_time_stdev[,pp]^2))
      # now draw the actual time series
      for (ii in 1:n_sample_this_chunk[cc]) {
        time_series_samples[,ii,pp] <- mvrnorm(n=1, mu=par_time_center[,pp], Sigma=covariance_samples[,,ii,pp])
        # normalize the ones that must be normalized (first batch with 1 as present, second with 0 as present)
        if (pp %in% c("fR", "fL", "fA", "fAw_fA", "fC")) {
          time_series_samples[,ii,pp] <- time_series_samples[,ii,pp]/time_series_samples[n_time,ii,pp]
        } else if (pp %in% c("GEOG")) {
          time_series_samples[,ii,pp] <- time_series_samples[,ii,pp] - time_series_samples[n_time,ii,pp]
        }
      }
    }

    # run the simulations
    model_co2_this_chunk <- model_temp_this_chunk <- mat.or.vec(nr=n_time, nc=n_sample_this_chunk[cc])
    prcout_co2 <- prcout_temp <- rep(NA, n_sample_this_chunk[cc])
    for (ii in 1:n_sample_this_chunk[cc]) {
      model_out <- model_forMCMC(par_calib=par_calib[[cc]][ii,],
                                 #par_time=time_series_samples[,ii,],
                                 par_time=par_time_center,
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
      model_co2_this_chunk[,ii] <- model_out[,"co2"]
      model_temp_this_chunk[,ii] <- model_out[,"temp"] + 15 # Note: Berner (2004; Eq 2.8 and 2.28) assumes present CO2 is 280 ppmv and T is 15 deg C
      ##    # normalize relative to "present" (t=0 Mya)
      ##    for (ii in 1:ncol(model_temp)) {model_temp[,ii] <- model_temp[,ii] - model_temp[58,ii]}
      # use the percent-outbound approach of Mill et al 2019 (Gondwana Research)
      prcout_co2[ii] <- percout(model_co2_this_chunk[,ii], windows$co2)
      prcout_temp[ii] <- percout(model_temp_this_chunk[,ii], windows$temp)
      ## toss anything outside of the temperature window at 100 Mya
      if (any(model_temp_this_chunk[idx_cret,ii] < windows$temp[idx_cret,"low"]) |
          any(model_temp_this_chunk[idx_cret,ii] > windows$temp[idx_cret,"high"])) {
        prcout_temp[ii] <- 1
      }
    }

    # do the following for CO2-only, temperature-only, and CO2+temperature

    # combine the good ones with previous good estimates - and remove any runs with inf time steps
    idx_save_temp <- which(prcout_temp <= prcout_threshold)
    idx_save_co2 <- which(prcout_co2 <= prcout_threshold)
    idx_save_both <- intersect(idx_save_temp, idx_save_co2)

    # pick which to use based on which data sets
    if (use_temperature & use_co2) {
      idx_save <- idx_save_both
    } else if (use_temperature & !use_co2) {
      idx_save <- idx_save_temp
    } else if (!use_temperature & use_co2) {
      idx_save <- idx_save_co2
    } else {
      print("ERROR: need to use at least one of temperature and/or CO2 data")
    }

    if( (length(idx_save)+max(nrow(par_calib_save),0)) >= n_sample_min) {
      if (is.null(nrow(par_calib_save))) {
        idx_save <- idx_save[1:(n_sample_min)]
      } else {
        idx_save <- idx_save[1:(n_sample_min-nrow(par_calib_save))]
      }
      par_calib_save <- rbind(par_calib_save, par_calib[[cc]][idx_save,])
      par_time_save <- abind(par_time_save, time_series_samples[,idx_save,], along=2)
      par_covar_save <- abind(par_covar_save, covariance_samples[,,idx_save,], along=3)
      break
    } else {
      par_calib_save <- rbind(par_calib_save, par_calib[[cc]][idx_save,])
      par_time_save <- abind(par_time_save, time_series_samples[,idx_save,], along=2)
      par_covar_save <- abind(par_covar_save, covariance_samples[,,idx_save,], along=3)
    }
  }
}
tend <- proc.time()
print(paste('precalibration took ',(tend[3]-tbeg[3])/60,' minutes', sep=''))



##==============================================================================
## End
##==============================================================================
