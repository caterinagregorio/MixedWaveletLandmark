rm(list = ls())

Sys.setenv(JAGS_HOME = "C:\\Program Files\\JAGS\\JAGS-4.3.0\\")

# load data
load("data/landmark_data.RData")
load("data/data.RData")

library(tidyverse)
library(magrittr)
library(survival)
library(JMbayes)
library(splines)
library(rms)
library(survival)
library(splines)
library(lme4)


source('fun.R')



####################
# MIXED & WAVELETS #
####################

LMdata$age_stand <- scale(LMdata$age)
LMdata$k.dtrend_stand <- scale(LMdata$k.dtrend*100)
LMdata$k.trend_stand <- scale(LMdata$k.trend)
LMdata$years <- cut(LMdata$Tstart,breaks = c(365,365*2,365*3,365*4,365*5),labels = c("second","third","fourth","fifth"),include.lowest = T)

summary(LMdata$years)

# Define landmark times
tLM <- seq(365, 365 * 5, by = 90)

LMdata %<>% filter(Tstart%in%tLM) 


## Simple (ipl) categorical covariates for oscillations
LMsupercox_mixed_wave1 <-
  coxph(
    Surv(Tstart, fup_days_w, dht) ~ age_stand + gender + nyha_34   + lvef_classi +  somma_comor_nc3 +
      ns(k.trend_stand, df = 4) +
      k.dtrend_stand +
      b14_d2 +
      b30_d2 +
      b90_d2 +
      b180_d2 +
      over180_d2 +
      strata(time) + cluster(id),
    data = LMdata,
    method = "breslow"
  )
summary(LMsupercox_mixed_wave1)
AIC(LMsupercox_mixed_wave1)

## Simple (ipl) continuous covariates for oscillations
LMsupercox_mixed_wave2 <-
  coxph(
    Surv(Tstart, fup_days_w, dht) ~ age_stand + gender + nyha_34   + lvef_classi +  somma_comor_nc3 +
      ns(k.trend_stand, df = 4) +
      k.dtrend_stand +
       ns(b14,df=4)+
       ns(b30,df=4)+
       ns(b90,df=4)+
       ns(b180,df=4)+
       ns(over180,df=4)+
      strata(time) + cluster(id),
    data = LMdata,
    method = "breslow"
  )

summary(LMsupercox_mixed_wave2)
AIC(LMsupercox_mixed_wave2)
AIC(LMsupercox_mixed_wave1)


# Create K=10 random folds
K <- 10
IDs <- unique(LMdata$id)
n <- length(IDs)
set.seed(123)
splits <- split(IDs, sample(rep(seq_len(K), length.out = n)))

cl <- makeCluster(2)
registerDoParallel(cl)
opts <- list(preschedule=TRUE)
clusterSetRNGStream(cl, 123)

perf <- foreach(k=1:K, .combine=rbind,
        .options.snow = opts) %dopar% {
          out_matrix <- cv_fun(k, splits, LMdata, potassio,index_visit,w=180)
          out_matrix}

# Stop cluster
stopCluster(cl)

colnames(perf)
summary(perf)

perf_data <- perf[,1]
perf_data <-rbind(perf_data$result.1,
                  perf_data$result.2,
                  perf_data$result.3,
                  perf_data$result.4,
                  perf_data$result.5,
                  perf_data$result.6,
                  perf_data$result.7,
                  perf_data$result.8,
                  perf_data$result.9,
                  perf_data$result.10)

colnames(perf_data)

perf_summ <- perf_data %>% 
  group_by(tLM,model) %>% 
  select(-k) %>% 
  summarise_all(mean)


perf_summ<- perf_data %>% 
  dplyr::group_by(tLM,model) %>% 
  dplyr::summarize(aucm=mean(auc),
            se.auc=sd(auc),
            bsm=mean(bs),
            se.bs=sd(bs))






save(LMdata,LMsupercox_mixed_wave1,LMsupercox_mixed_wave2,perf_summ,file="data/04_LM model_cv.RData")
