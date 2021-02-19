#============================================================================#
# Simulation code for N=5000, 14 different simulation scenarios              #
# Author: Jiayi Ji                                                           #
# July 6, 2020                                                               #
#============================================================================#
library(survival)
library(BART)
library(randomForestSRC)
library(tidyverse)
library(AFTrees)
source("code/aft_bart_sp_HTE.R")
source("code/aft_bart_np_HTE.R")
source("code/aft_bart_np_SurvivalProb.R")
source("code/cox_HTE.R")
source("code/aft_weibull_HTE.R")
source("code/aft_lognormal_HTE.R")
source("code/rsf_HTE.R")
source("code/simulation_design")
# scenario 1: N = 5000, PH, 20% censor, heterogeneous setting (i), strong overlap-------------- 
for (i in 1:250){
  set.seed(1111+i)
  mydata <- data_gen_censor(n=5000, p=10,  PH=TRUE, censor = "20%", setting = 1, overlap = "strong")
  result_aft_bart_sp_HTE <- aft_bart_sp_HTE()
  result_cox_HTE <- cox_HTE()
  result_rsf_HTE <- rsf_HTE()
  result_aft_weibull_HTE <- aft_weibull_HTE()
  result_aft_lognormal_HTE <- aft_lognormal_HTE()
  result_aft_bart_np_HTE <- aft_bart_np_HTE()
  save(result_aft_bart_sp_HTE, result_cox_HTE, result_rsf_HTE, result_aft_weibull_HTE, result_aft_lognormal_HTE, result_aft_lognormal_HTE, result_aft_bart_np_HTE, file = paste0("result_scenario_1_",i,".Rdata"))
}

###### Testing ############################################################################
# i = 1
# set.seed(1111+i)
# mydata <- data_gen_censor(n=500, p=10,  PH=TRUE, censor = "20%", setting = 1, overlap = "strong")
# result_aft_bart_sp_HTE <- aft_bart_sp_HTE()
# result_cox_HTE <- cox_HTE()
# result_rsf_HTE <- rsf_HTE()
# result_aft_weibull_HTE <- aft_weibull_HTE()
# result_aft_lognormal_HTE <- aft_lognormal_HTE()
# result_aft_bart_np_HTE <- aft_bart_np_HTE()

# result_aft_bart_sp_HTE
# # # A tibble: 50 x 3
# # p1_category true_ate bias_ate
# # <int>    <dbl>    <dbl>
# #   1           1     14.6    0.348
# # 2           2     14.3    0.274
# # 3           3     14.3    0.241
# # 4           4     14.5    0.278
# # 5           5     13.9    0.296
# # 6           6     13.8    0.287
# # 7           7     13.3    0.301
# # 8           8     14.2    0.162
# # 9           9     13.4    0.269
# # 10          10     14.2    0.306
# # # … with 40 more rows
# result_cox_HTE
# # # A tibble: 50 x 3
# # p1_category true_ate bias_ate
# # <int>    <dbl>    <dbl>
# #   1           1     14.6     2.24
# # 2           2     14.0     2.33
# # 3           3     14.2     2.34
# # 4           4     14.6     2.17
# # 5           5     13.4     2.18
# # 6           6     14.5     2.21
# # 7           7     14.8     2.16
# # 8           8     14.0     2.24
# # 9           9     14.1     2.24
# # 10          10     14.0     2.27
# # # … with 40 more rows
# result_rsf_HTE
# # A tibble: 50 x 3
# # p1_category true_ate bias_ate
# # <int>    <dbl>    <dbl>
# #   1           1     13.8    0.365
# # 2           2     14.7    0.287
# # 3           3     14.5    0.361
# # 4           4     13.5    0.327
# # 5           5     13.9    0.395
# # 6           6     14.8    0.313
# # 7           7     14.3    0.213
# # 8           8     13.3    0.348
# # 9           9     14.1    0.417
# # 10          10     14.2    0.315
# # # … with 40 more rows
# result_aft_weibull_HTE
# # A tibble: 50 x 3
# # p1_category true_ate bias_ate
# # <int>    <dbl>    <dbl>
# #   1           1     13.2     2.51
# # 2           2     14.7     2.59
# # 3           3     14.2     2.59
# # 4           4     14.3     2.49
# # 5           5     15.4     2.60
# # 6           6     15.7     2.60
# # 7           7     15.4     2.62
# # 8           8     15.0     2.55
# # 9           9     14.1     2.50
# # 10          10     15.2     2.55
# # # … with 40 more rows
# result_aft_lognormal_HTE
# # # A tibble: 50 x 3
# # p1_category true_ate bias_ate
# # <int>    <dbl>    <dbl>
# #   1           1     13.7     3.50
# # 2           2     15.5     4.35
# # 3           3     14.9     2.30
# # 4           4     14.8     3.11
# # 5           5     14.3     2.72
# # 6           6     14.8     4.02
# # 7           7     14.5     4.29
# # 8           8     14.0     3.25
# # 9           9     15.2     3.63
# # 10          10     14.9     3.47
# # # … with 40 more rows
# result_aft_bart_sp_HTE
# # # A tibble: 50 x 3
# # p1_category true_ate bias_ate
# # <int>    <dbl>    <dbl>
# #   1           1     14.6    0.348
# # 2           2     14.3    0.274
# # 3           3     14.3    0.241
# # 4           4     14.5    0.278
# # 5           5     13.9    0.296
# # 6           6     13.8    0.287
# # 7           7     13.3    0.301
# # 8           8     14.2    0.162
# # 9           9     13.4    0.269
# # 10          10     14.2    0.306
# # # … with 40 more rows

# scenario 2: N = 5000, PH, 60% censor, heterogeneous setting (i), strong overlap-------------- 
for (i in 1:250){
  set.seed(1111+i)
  mydata <- data_gen_censor(n=5000, p=10,  PH=TRUE, censor = "60%", setting = 1, overlap = "strong")
  result_aft_bart_sp_HTE <- aft_bart_sp_HTE()
  result_cox_HTE <- cox_HTE()
  result_rsf_HTE <- rsf_HTE()
  result_aft_weibull_HTE <- aft_weibull_HTE()
  result_aft_lognormal_HTE <- aft_lognormal_HTE()
  result_aft_bart_np_HTE <- aft_bart_np_HTE()
  save(result_aft_bart_sp_HTE, result_cox_HTE, result_rsf_HTE, result_aft_weibull_HTE, result_aft_lognormal_HTE, result_aft_lognormal_HTE, result_aft_bart_np_HTE, file = paste0("result_scenario_2_",i,".Rdata"))
}

# scenario 3: N = 5000, nPH, 20% censor, heterogeneous setting (i), strong overlap-------------- 
for (i in 1:250){
  set.seed(1111+i)
  mydata <- data_gen_censor(n=5000, p=10,  PH=FALSE, censor = "20%", setting = 1, overlap = "strong")
  result_aft_bart_sp_HTE <- aft_bart_sp_HTE()
  result_cox_HTE <- cox_HTE()
  result_rsf_HTE <- rsf_HTE()
  result_aft_weibull_HTE <- aft_weibull_HTE()
  result_aft_lognormal_HTE <- aft_lognormal_HTE()
  result_aft_bart_np_HTE <- aft_bart_np_HTE()
  save(result_aft_bart_sp_HTE, result_cox_HTE, result_rsf_HTE, result_aft_weibull_HTE, result_aft_lognormal_HTE, result_aft_lognormal_HTE, result_aft_bart_np_HTE, file = paste0("result_scenario_3_",i,".Rdata"))
}

# scenario 4: N = 5000, nPH, 60% censor, heterogeneous setting (i), strong overlap-------------- 
for (i in 1:250){
  set.seed(1111+i)
  mydata <- data_gen_censor(n=5000, p=10,  PH=FALSE, censor = "60%", setting = 1, overlap = "strong")
  result_aft_bart_sp_HTE <- aft_bart_sp_HTE()
  result_cox_HTE <- cox_HTE()
  result_rsf_HTE <- rsf_HTE()
  result_aft_weibull_HTE <- aft_weibull_HTE()
  result_aft_lognormal_HTE <- aft_lognormal_HTE()
  result_aft_bart_np_HTE <- aft_bart_np_HTE()
  save(result_aft_bart_sp_HTE, result_cox_HTE, result_rsf_HTE, result_aft_weibull_HTE, result_aft_lognormal_HTE, result_aft_lognormal_HTE, result_aft_bart_np_HTE, file = paste0("result_scenario_4_",i,".Rdata"))
}

# scenario 5: N = 5000, PH, 20% censor, heterogeneous setting (ii), strong overlap-------------- 
for (i in 1:250){
  set.seed(1111+i)
  mydata <- data_gen_censor(n=5000, p=10,  PH=TRUE, censor = "20%", setting = 1, overlap = "strong")
  result_aft_bart_sp_HTE <- aft_bart_sp_HTE()
  result_cox_HTE <- cox_HTE()
  result_rsf_HTE <- rsf_HTE()
  result_aft_weibull_HTE <- aft_weibull_HTE()
  result_aft_lognormal_HTE <- aft_lognormal_HTE()
  result_aft_bart_np_HTE <- aft_bart_np_HTE()
  save(result_aft_bart_sp_HTE, result_cox_HTE, result_rsf_HTE, result_aft_weibull_HTE, result_aft_lognormal_HTE, result_aft_lognormal_HTE, result_aft_bart_np_HTE, file = paste0("result_scenario_5_",i,".Rdata"))
}

# scenario 6: N = 5000, PH, 60% censor, heterogeneous setting (ii), strong overlap-------------- 
for (i in 1:250){
  set.seed(1111+i)
  mydata <- data_gen_censor(n=5000, p=10,  PH=TRUE, censor = "60%", setting = 1, overlap = "strong")
  result_aft_bart_sp_HTE <- aft_bart_sp_HTE()
  result_cox_HTE <- cox_HTE()
  result_rsf_HTE <- rsf_HTE()
  result_aft_weibull_HTE <- aft_weibull_HTE()
  result_aft_lognormal_HTE <- aft_lognormal_HTE()
  result_aft_bart_np_HTE <- aft_bart_np_HTE()
  save(result_aft_bart_sp_HTE, result_cox_HTE, result_rsf_HTE, result_aft_weibull_HTE, result_aft_lognormal_HTE, result_aft_lognormal_HTE, result_aft_bart_np_HTE, file = paste0("result_scenario_6_",i,".Rdata"))
}

# scenario 7: N = 5000, nPH, 20% censor, heterogeneous setting (ii), strong overlap-------------- 
for (i in 1:250){
  set.seed(1111+i)
  mydata <- data_gen_censor(n=5000, p=10,  PH=FALSE, censor = "20%", setting = 1, overlap = "strong")
  result_aft_bart_sp_HTE <- aft_bart_sp_HTE()
  result_cox_HTE <- cox_HTE()
  result_rsf_HTE <- rsf_HTE()
  result_aft_weibull_HTE <- aft_weibull_HTE()
  result_aft_lognormal_HTE <- aft_lognormal_HTE()
  result_aft_bart_np_HTE <- aft_bart_np_HTE()
  save(result_aft_bart_sp_HTE, result_cox_HTE, result_rsf_HTE, result_aft_weibull_HTE, result_aft_lognormal_HTE, result_aft_lognormal_HTE, result_aft_bart_np_HTE, file = paste0("result_scenario_7_",i,".Rdata"))
}

# scenario 8: N = 5000, nPH, 60% censor, heterogeneous setting (ii), strong overlap-------------- 
for (i in 1:250){
  set.seed(1111+i)
  mydata <- data_gen_censor(n=5000, p=10,  PH=FALSE, censor = "60%", setting = 1, overlap = "strong")
  result_aft_bart_sp_HTE <- aft_bart_sp_HTE()
  result_cox_HTE <- cox_HTE()
  result_rsf_HTE <- rsf_HTE()
  result_aft_weibull_HTE <- aft_weibull_HTE()
  result_aft_lognormal_HTE <- aft_lognormal_HTE()
  result_aft_bart_np_HTE <- aft_bart_np_HTE()
  save(result_aft_bart_sp_HTE, result_cox_HTE, result_rsf_HTE, result_aft_weibull_HTE, result_aft_lognormal_HTE, result_aft_lognormal_HTE, result_aft_bart_np_HTE, file = paste0("result_scenario_8_",i,".Rdata"))
}

# scenario 9: N = 5000, PH, 20% censor, heterogeneous setting (iii), strong overlap-------------- 
for (i in 1:250){
  set.seed(1111+i)
  mydata <- data_gen_censor(n=5000, p=10,  PH=TRUE, censor = "20%", setting = 1, overlap = "strong")
  result_aft_bart_sp_HTE <- aft_bart_sp_HTE()
  result_cox_HTE <- cox_HTE()
  result_rsf_HTE <- rsf_HTE()
  result_aft_weibull_HTE <- aft_weibull_HTE()
  result_aft_lognormal_HTE <- aft_lognormal_HTE()
  result_aft_bart_np_HTE <- aft_bart_np_HTE()
  save(result_aft_bart_sp_HTE, result_cox_HTE, result_rsf_HTE, result_aft_weibull_HTE, result_aft_lognormal_HTE, result_aft_lognormal_HTE, result_aft_bart_np_HTE, file = paste0("result_scenario_9_",i,".Rdata"))
}

# scenario 10: N = 5000, PH, 60% censor, heterogeneous setting (iii), strong overlap-------------- 
for (i in 1:250){
  set.seed(1111+i)
  mydata <- data_gen_censor(n=5000, p=10,  PH=TRUE, censor = "60%", setting = 1, overlap = "strong")
  result_aft_bart_sp_HTE <- aft_bart_sp_HTE()
  result_cox_HTE <- cox_HTE()
  result_rsf_HTE <- rsf_HTE()
  result_aft_weibull_HTE <- aft_weibull_HTE()
  result_aft_lognormal_HTE <- aft_lognormal_HTE()
  result_aft_bart_np_HTE <- aft_bart_np_HTE()
  save(result_aft_bart_sp_HTE, result_cox_HTE, result_rsf_HTE, result_aft_weibull_HTE, result_aft_lognormal_HTE, result_aft_lognormal_HTE, result_aft_bart_np_HTE, file = paste0("result_scenario_10_",i,".Rdata"))
}

# scenario 11: N = 5000, nPH, 20% censor, heterogeneous setting (iii), strong overlap-------------- 
for (i in 1:250){
  set.seed(1111+i)
  mydata <- data_gen_censor(n=5000, p=10,  PH=FALSE, censor = "20%", setting = 1, overlap = "strong")
  result_aft_bart_sp_HTE <- aft_bart_sp_HTE()
  result_cox_HTE <- cox_HTE()
  result_rsf_HTE <- rsf_HTE()
  result_aft_weibull_HTE <- aft_weibull_HTE()
  result_aft_lognormal_HTE <- aft_lognormal_HTE()
  result_aft_bart_np_HTE <- aft_bart_np_HTE()
  save(result_aft_bart_sp_HTE, result_cox_HTE, result_rsf_HTE, result_aft_weibull_HTE, result_aft_lognormal_HTE, result_aft_lognormal_HTE, result_aft_bart_np_HTE, file = paste0("result_scenario_11_",i,".Rdata"))
}

# scenario 12: N = 5000, nPH, 60% censor, heterogeneous setting (iii), strong overlap-------------- 
for (i in 1:250){
  set.seed(1111+i)
  mydata <- data_gen_censor(n=5000, p=10,  PH=FALSE, censor = "60%", setting = 1, overlap = "strong")
  result_aft_bart_sp_HTE <- aft_bart_sp_HTE()
  result_cox_HTE <- cox_HTE()
  result_rsf_HTE <- rsf_HTE()
  result_aft_weibull_HTE <- aft_weibull_HTE()
  result_aft_lognormal_HTE <- aft_lognormal_HTE()
  result_aft_bart_np_HTE <- aft_bart_np_HTE()
  save(result_aft_bart_sp_HTE, result_cox_HTE, result_rsf_HTE, result_aft_weibull_HTE, result_aft_lognormal_HTE, result_aft_lognormal_HTE, result_aft_bart_np_HTE, file = paste0("result_scenario_12_",i,".Rdata"))
}

# scenario 13: N = 5000, nPH, 60% censor, heterogeneous setting (iii), medium overlap-------------- 
for (i in 1:250){
  set.seed(1111+i)
  mydata <- data_gen_censor(n=5000, p=10,  PH=FALSE, censor = "60%", setting = 1, overlap = "medium")
  result_aft_bart_sp_HTE <- aft_bart_sp_HTE()
  result_cox_HTE <- cox_HTE()
  result_rsf_HTE <- rsf_HTE()
  result_aft_weibull_HTE <- aft_weibull_HTE()
  result_aft_lognormal_HTE <- aft_lognormal_HTE()
  result_aft_bart_np_HTE <- aft_bart_np_HTE()
  save(result_aft_bart_sp_HTE, result_cox_HTE, result_rsf_HTE, result_aft_weibull_HTE, result_aft_lognormal_HTE, result_aft_lognormal_HTE, result_aft_bart_np_HTE, file = paste0("result_scenario_13_",i,".Rdata"))
}

# scenario 14: N = 5000, nPH, 60% censor, heterogeneous setting (iii), weak overlap-------------- 
for (i in 1:250){
  set.seed(1111+i)
  mydata <- data_gen_censor(n=5000, p=10,  PH=FALSE, censor = "60%", setting = 1, overlap = "weak")
  result_aft_bart_sp_HTE <- aft_bart_sp_HTE()
  result_cox_HTE <- cox_HTE()
  result_rsf_HTE <- rsf_HTE()
  result_aft_weibull_HTE <- aft_weibull_HTE()
  result_aft_lognormal_HTE <- aft_lognormal_HTE()
  result_aft_bart_np_HTE <- aft_bart_np_HTE()
  save(result_aft_bart_sp_HTE, result_cox_HTE, result_rsf_HTE, result_aft_weibull_HTE, result_aft_lognormal_HTE, result_aft_lognormal_HTE, result_aft_bart_np_HTE, file = paste0("result_scenario_14_",i,".Rdata"))
}

