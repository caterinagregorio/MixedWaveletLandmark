rm(list = ls())

Sys.setenv(JAGS_HOME = "C:\\Program Files\\JAGS\\JAGS-4.3.0\\")

# load data
load("data/data.RData")
load("data/long_res.RData")

library(plyr)
library(doParallel)
library(tidyverse)
library(magrittr)
library(rms)
library(rjags)
library(data.table)
library(splines)
library(JMbayes)

source("fun.R")

cl <- makeCluster(2)
registerDoParallel(cl)
opts <- list(preschedule=TRUE)
clusterSetRNGStream(cl, 123)

# Define landmark times
tLM <- seq(365, 365 * 5, by = 90)

# data for mixed and mixed+wavelets
LMdata_obj <- llply(
  tLM,
  w = 180,
  landmark_data_mixed_wavelet,
  longdataexp = potassioext,
  longdata = potassio,
  survdata = index_visit,
 fm=fp,
 kn=kn,
 formula.k=f.pot,
  .parallel=T,
  .paropts = list(.options.snow=opts))

LMdata <- as.data.frame(do.call(rbind,lapply(LMdata_obj,function(x)x$datasurv)))

LMdata %<>% dplyr::mutate_at(vars(contains("_d2")),as.factor) %>%
  dplyr::mutate_at(vars(contains("_d2")),relevel,ref="0") 

LMdata%<>% 
  mutate(over180_d=case_when(abs(over180)>0~1,
                             TRUE~0),
         b180_d=case_when(abs(b180)>0~1,
                          TRUE~0),
         b90_d=case_when(abs(b90)>0~1,
                         TRUE~0),
         b30_d=case_when(abs(b30)>0~1,
                         TRUE~0),
         b14_d=case_when(abs(b14)>0~1,
                         TRUE~0))

LMdata$Tstart <- LMdata$time

colnames(LMdata)

save(LMdata,LMdata_obj,file="data/landmark_data.RData")
