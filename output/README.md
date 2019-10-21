Output files go in here.

The ones that should come with the model codes and results are detailed below. The naming convention for the calibration experiments is "d[data set identifier]p[parameter set identifier]s[seafloor spreading identifier]l[likelihood identifier]". The identifiers can take these values:
* data set is one of P (Park and Royer 2011) or F (Foster et al 2017)
* parameter set is one of P (Park and Royer 2011) or A (all parameters), followed by U if the uncertainty parameter (this work) is added
* seafloor spreading is one of O (original, from Royer et al 2014), L (Lenton et al 2018) or R (Domeier and Torsvik 2019)
* likelihood is one of U (unimodal, as in Park and Royer 2011) or M (mixture model, this work)
* the likelihood form is followed by either sn or nm, depending on whether the kernels in the mixture model are skew-normal (sn) or Gaussian (nm). If folks are interested, this form can also be modified to be log-normal, following Park and Royer (2011), or other forms too.

For the main calibration results:
* `geocarb_mcmcoutput_dFpAUsRlMsn_04Oct2019.RData`
  * The MCMC output for the main text experiments
* `chains_analysis_05Oct2019.rds`
  * The MCMC chains resulting from the thinning and burn-in processing. The main text chains are processed along with all of the supplemental experiments, which are all on this results file too.
* `analysis_06Oct2019.RData`
  * Analysis results like model hindcast quantiles and other numbers needed for generating plots and presentation in the manuscript.

For the main sensitivity experiment results:
* `precal_parameters_06Oct2019.rds`
  * Resulting parameter sets that pass the precalibration test
* `geocarb_sobol-1-tot_sensNS_dFpAUsRlMsn_08Oct2019.txt`
  * Sobol' results for the first-order and total sensitivity indices, and confidence intervals
* `geocarb_sobol-2_sensNS_dFpAUsRlMsn_08Oct2019.txt`
  * Sobol' results for the second-order sensitivity indices, and confidence intervals
* `sobol_sensNS_dFpAUsRlMsn_08Oct2019.rds`
  * The raw Sobol' analysis output object (comes out of `sobolTony.Rs`)

These files give the MCMC output for the supplemental experiments:
* `geocarb_mcmcoutput_dFpAUsRlMnm_03Oct2019.RData`
  * Same as main text, but uses Gaussian kernels instead of skew-normal
* `geocarb_mcmcoutput_dPpAUsRlMsn_03Oct2019.RData`
  * Same as main text, but uses the data set of Park and Royer (2011)
* `geocarb_mcmcoutput_dFpAUsRlUsn_03Oct2019.RData`
  * Same as main text, but uses unimodal likelihood for each time slice a la Park and Royer (2011), instead of the mixture model where each data point is represented by a kernel.
* `geocarb_mcmcoutput_dFpPUsRlMsn_03Oct2019.RData`
  * Same as main text, but only calibrates the 6 parameters of Park and Royer (2011) (ACT, FERT, GYM, LIFE, deltaT2X and GLAC), and the stdev uncertainty parameter.
* `geocarb_mcmcoutput_dFpPUsRlUsn_03Oct2019.RData`
  * Same as main text, but uses the unimodal likelihood AND only calibrates the 6 parameters of Park and Royer (2011) (ACT, FERT, GYM, LIFE, deltaT2X and GLAC), and the stdev uncertainty parameter. (so, a combination of the previous 2)
* `geocarb_mcmcoutput_dPpPUsRlUsn_04Oct2019.RData`
  * Using Park and Royer (2011) data, parameters and unimodal likelihood, but with the fSR time series of Domeier and Torsvik (2019)
* `geocarb_mcmcoutput_dPpPUsLlUsn_04Oct2019.RData`
  * Using Park and Royer (2011) data, parameters and unimodal likelihood, but with the fSR time series of Lenton et al (2018)
* `geocarb_mcmcoutput_dPpPUsOlUsn_04Oct2019.RData`
  * Using Park and Royer (2011) data, parameters, unimodal likelihood, and fSR time series. Main difference is the calibration approach (MCMC, this work)
* `geocarb_mcmcoutput_dPpPUsRlMsn_03Oct2019.RData`
  * Same as above, but using a mixture model likelihood (this work) instead of unimodal (Park and Royer 2011)

For the warm-up chain:
* First, running these chains:
  * `geocarb_mcmcoutput_mix_06Jul2019sn-mix.RData`
* Which yield these initial conditions files (the proposal covariance matrix and parameter estimates from the final Markov state):
  * `param_init_sn-mix_07Jul2019.rds`
  * `covar_init_sn-mix_07Jul2019.rds`
* Then run these chains, starting from the above initial conditions:
  * `geocarb_mcmcoutput_mix_07Jul2019sn-mix.RData`
* Which yield these initial conditions files:
  * `param_init_sn-mix_08Jul2019.rds`
  * `covar_init_sn-mix_08Jul2019.rds`
* Then run these chains, starting from the above initial conditions:
  * `geocarb_mcmcoutput_mix_08Jul2019sn-mix.RData`
* Which yield these initial conditions files:
  * `param_init_sn-mix_09Jul2019.rds`
  * `covar_init_sn-mix_09Jul2019.rds`
* Then run these chains, starting from the above initial conditions:
  * `geocarb_mcmcoutput_mix_12Jul2019sn-mix.RData`
* Which yield these initial conditions files:
  * `param_init_sn-mix_14Jul2019.rds`
  * `covar_init_sn-mix_14Jul2019.rds`
* Then run these chains, starting from the above initial conditions:
  * `geocarb_mcmcoutput_mix_15Jul2019sn-mix.RData`
* Which yield these initial conditions files:
  * `param_init_sn-mix_17Jul2019.rds`
  * `covar_init_sn-mix_17Jul2019.rds`
* Then run these chains, starting from the above initial conditions:
  * `geocarb_mcmcoutput_mix_20Jul2019sn-mix.RData`
* Which yield these initial conditions files:
  * `param_init_sn-mix_22Jul2019.rds`
  * `covar_init_sn-mix_22Jul2019.rds`
* Then run these chains, starting from the above initial conditions:
  * `geocarb_mcmcoutput_mix_10Aug2019sn-mix.RData`
* Which yield these initial conditions files:
  * `param_init_sn-mix_12Aug2019.rds`
  * `covar_init_sn-mix_12Aug2019.rds`
