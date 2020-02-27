rm(list=ls())

load("../output/lhs_param_c_out50.RData"); par_calib_c_50 <- par_calib_save
load("../output/lhs_param_t_out50.RData"); par_calib_t_50 <- par_calib_save
load("../output/lhs_param_ct_out50.RData"); par_calib_ct_50 <- par_calib_save

load("../output/lhs_param_c_out45.RData"); par_calib_c_45 <- par_calib_save
load("../output/lhs_param_t_out45.RData"); par_calib_t_45 <- par_calib_save
load("../output/lhs_param_ct_out45.RData"); par_calib_ct_45 <- par_calib_save

load("../output/lhs_param_c_out40.RData"); par_calib_c_40 <- par_calib_save
load("../output/lhs_param_t_out40.RData"); par_calib_t_40 <- par_calib_save
load("../output/lhs_param_ct_out40.RData"); par_calib_ct_40 <- par_calib_save

load("../output/lhs_param_c_out35.RData"); par_calib_c_35 <- par_calib_save
load("../output/lhs_param_t_out35.RData"); par_calib_t_35 <- par_calib_save
load("../output/lhs_param_ct_out35.RData"); par_calib_ct_35 <- par_calib_save

load("../output/lhs_param_c_out30.RData"); par_calib_c_30 <- par_calib_save
load("../output/lhs_param_t_out30.RData"); par_calib_t_30 <- par_calib_save
#load("../output/lhs_param_ct_out30.RData"); par_calib_ct_30 <- par_calib_save

par_calib <- par_calib_ct_50
idx <- sample(1:nrow(par_calib), size=10000, replace=FALSE)
quantile(par_calib[idx,10], c(.05,.5,.95))
par_calib <- par_calib_ct_45
idx <- sample(1:nrow(par_calib), size=10000, replace=FALSE)
quantile(par_calib[idx,10], c(.05,.5,.95))
par_calib <- par_calib_ct_40
idx <- sample(1:nrow(par_calib), size=10000, replace=FALSE)
quantile(par_calib[idx,10], c(.05,.5,.95))
par_calib <- par_calib_ct_35
idx <- sample(1:nrow(par_calib), size=10000, replace=FALSE)
quantile(par_calib[idx,10], c(.05,.5,.95))
#par_calib <- par_calib_ct_30
#idx <- sample(1:nrow(par_calib), size=10000, replace=FALSE)
#quantile(par_calib[idx,10], c(.05,.5,.95))

par_calib <- par_calib_c_50
idx <- sample(1:nrow(par_calib), size=10000, replace=FALSE)
quantile(par_calib[idx,10], c(.05,.5,.95))
par_calib <- par_calib_c_45
idx <- sample(1:nrow(par_calib), size=10000, replace=FALSE)
quantile(par_calib[idx,10], c(.05,.5,.95))
par_calib <- par_calib_c_40
idx <- sample(1:nrow(par_calib), size=10000, replace=FALSE)
quantile(par_calib[idx,10], c(.05,.5,.95))
par_calib <- par_calib_c_35
idx <- sample(1:nrow(par_calib), size=10000, replace=FALSE)
quantile(par_calib[idx,10], c(.05,.5,.95))
par_calib <- par_calib_c_30
idx <- sample(1:nrow(par_calib), size=10000, replace=FALSE)
quantile(par_calib[idx,10], c(.05,.5,.95))

par_calib <- par_calib_t_50
idx <- sample(1:nrow(par_calib), size=10000, replace=FALSE)
quantile(par_calib[idx,10], c(.05,.5,.95))
par_calib <- par_calib_t_45
idx <- sample(1:nrow(par_calib), size=10000, replace=FALSE)
quantile(par_calib[idx,10], c(.05,.5,.95))
par_calib <- par_calib_t_40
idx <- sample(1:nrow(par_calib), size=10000, replace=FALSE)
quantile(par_calib[idx,10], c(.05,.5,.95))
par_calib <- par_calib_t_35
idx <- sample(1:nrow(par_calib), size=10000, replace=FALSE)
quantile(par_calib[idx,10], c(.05,.5,.95))
par_calib <- par_calib_t_30
idx <- sample(1:nrow(par_calib), size=10000, replace=FALSE)
quantile(par_calib[idx,10], c(.05,.5,.95))
