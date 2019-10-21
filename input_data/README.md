Files in the `input_data` directory:

`CO2_Proxy_Foster2017_calib.csv`
* original CO2 proxy data set from Foster et al (2017, doi: 10.1038/ncomms14845)

`CO2_Proxy_Foster2017_calib_NM-co2_25Sep2018.csv`
* modified to include the estimated parameters if a Gaussian (normal) kernel is fitted to each data point, assuming the high/low uncertainties are a 68% range (16-84%)

`CO2_Proxy_Foster2017_calib_SN-co2_25Sep2018.csv`
* modified to include the estimated parameters if a skew-normal kernel is fitted to each data point, assuming the high/low uncertainties are a 68% range (16-84%)

`CO2_Proxy_PR2011_calib_SN-co2_28Aug2019.csv`
* CO2 proxy data set from Park and Royer (2011, doi: 10.2475/01.2011.01), for comparison (supplemental) experiments to evaluate the impacts of the calibration data

`GEOCARB_input_arrays.csv`
* Original data set of input time series, from Royer et al. (2014, doi: 10.2475/09.2014.01). Note that the fAw_fA time series needs to be normalized to have value of 1 in the last time step.

`GEOCARB_input_arrays_Lenton.csv`
* Same as above, but with the fSR (seafloor spreading rate) from Lenton et al. (2018, doi: 10.1016/j.earscirev.2017.12.004). Note that the fAw_fA time series needs to be normalized to have value of 1 in the last time step.

`GEOCARB_input_arrays_NewDegassing.csv`
* Same as above, but with the fSR (seafloor spreading rate) from Domeier and Torsvik (2019, doi: 10.1017/S0016756817001005). Note that the fAw_fA time series needs to be normalized to have value of 1 in the last time step.

`alternative fSR series.xlsx`
* Original data set used to create the above forcing data set.

`GEOCARB_input_summaries_calib_PR2011.csv`
* Only the 6 parameters of Park and Royer (2011): ACT, FERT, LIFE, GYM, deltaT2X and GLAC. Does not include the stdev parameter

`GEOCARB_input_summaries_calib_PR2011_stdev.csv`
* Only the 6 parameters of Park and Royer (2011): ACT, FERT, LIFE, GYM, deltaT2X and GLAC. Includes the stdev parameter

`GEOCARB_input_summaries_calib_all.csv`
* All 68 parameters are calibrated, not including the stdev parameter

`GEOCARB_input_summaries_calib_all_stdev.csv`
* All 69 parameters are calibrated, including the stdev parameter

`GEOCARB_input_summaries_calib_stdevOnly.csv`
* Only the stdev parameter is being calibrated, in order to run the preliminary warm-up chain for that parameter

`GEOCARB_output--both error envelopes_TW.xlsx`
* Model simulation data for CO2 hindcast of Royer et al. (2014). Used for comparison of the new model ensemble vs the old one (manuscript Figure 3).

`Park and Royer (2011) climate sensitivity (Figure 3).xlsx`
* The original full estimated deltaT2X (ESS) probability distribution data from Park and Royer (2011). Not used, except for grabbing the 85% variance reduction data for the next data set for comparison plotting.

`ParkRoyer2011_Fig3_85varred.csv`
* The deltaT2X distribution from Park and Royer (2011) that results from 85% reduction in variance

`co2proxy_2010.dat`
* CO2 proxy data set from Park and Royer (2011) for comparison plotting

`Plot_ESS_Cenozoic_input.csv`
* Data set of ESS estimates from previous work (used to generate manuscript Figure 1)
