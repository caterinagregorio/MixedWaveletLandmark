Sys.setenv(JAGS_HOME = "C:\\Program Files\\JAGS\\JAGS-4.3.0\\")

library(tidyverse)
library(gtsummary)
library(magrittr)
library(purrr)
library(nlme)
library(splines)
library(tseries)
library(plyr)
library(rms)
library(data.table)
library(splines)
library(JMbayes)

source("fun.R")

# load data 
load("data/data.RData")


# set random seed
set.seed(12345)

# sample 200 id
sample_id <- sample(index_visit$id,200)
potassio_sample <- potassio %>% filter(id%in%sample_id)


#################
## MIXED MODEL ##
#################


# Formula for potassium mixed model
f.pot <- "k ~ gender+age_stand+irc+nyha_34"
kn <- as.numeric(quantile(potassio$time, probs=c(0.25,0.50,0.75)))
formula.k <- formula(paste0(f.pot,"+ns(time,knots=c(",paste(kn,collapse = ","),"))"))
formula.random.k <- formula(paste0("~ns(time,knots=c(",paste(kn,collapse = ","),"))| id"))


potassioext$z1 <- ns(potassioext$time,knots = kn)[,1]
potassioext$z2 <- ns(potassioext$time,knots = kn)[,2]
potassioext$z3 <- ns(potassioext$time,knots = kn)[,3]
potassioext$z4 <- ns(potassioext$time,knots = kn)[,4]

potassioext$dz1 <- dns(potassioext$time,knots = kn)[,1]
potassioext$dz2 <- dns(potassioext$time,knots = kn)[,2]
potassioext$dz3 <- dns(potassioext$time,knots = kn)[,3]
potassioext$dz4 <- dns(potassioext$time,knots = kn)[,4]


# Estimation of the linear mixed model
fp <- lme(formula.k,
          data = potassio,
          random = formula.random.k,
          control = lmeControl(opt = 'optim'))

f.pot <- formula(f.pot)

potassioext$k.trend <- predict(fp,newdata = potassioext)
potassioext$k.detrended <- potassioext$k-potassioext$k.trend

potassio_sample$k.trend <- predict(fp,newdata = potassio_sample)
potassio_sample$res <- potassio_sample$k.trend-potassio_sample$k

###################
# VARIOGRAM #
###################

var <- Variogram(fp,resType = "normalized")

#save
save(fp,kn,f.pot,potassio,potassio_sample,potassioext,sample_id,var,file="data//long_res.RData")
