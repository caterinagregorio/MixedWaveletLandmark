adf_test <- function(id,dataset){
  sel <- dataset[dataset$id==id,]
  res.int <- approx(sel$time,sel$res,xout=0:max(sel$time))$y
  skip <- F
  tryCatch(
    p <- adf.test(res.int,"stationary")$p.value
    ,
    error=function(e)skip <<-T )
  if(skip)p <- NA
  return(p)
}


fitwavelet <- function(x,data,lambda,deriv) {
  
  load("data//soglia significatività.RData")
  source("fun.R")
  require(tidyverse)
  require(magrittr)
  require(JMbayes2)
  require(WaveletComp)
  require(splines)
  require(plyr)
  require(doParallel)
  
  if(deriv){
    # I select data from subject x
    my.data <- data %>%
      group_by(id) %>%
      filter(id == x) %>%
      ungroup() %>%
      dplyr::select(time, k, k.trend,k.dtrend ,k.detrended, dist_m) %>%
      dplyr::rename(date = time)
    
  }else{
    # I select data from subject x
    my.data <- data %>%
      group_by(id) %>%
      filter(id == x) %>%
      ungroup() %>%
      dplyr::select(time, k, k.trend ,k.detrended, dist_m) %>%
      dplyr::rename(date = time)
    
  }

 
  
  
  skip <- FALSE
  
  tryCatch({
    
    # to scale data
    sd <- sd(my.data$k.detrended-mean(my.data$k.detrended))
    my.data$k.detrended.scaled <- scale(my.data$k.detrended,scale=T,center=T)
    
    time0 <- (max(my.data$date)+1):max(my.data$date)+100
    k.detrended0 <- rep(0,length(time0))
    my.data0 <- data.frame(date=time0,
                           k.detrended=k.detrended0,
                           dist_m=rep(1000,length(time0)))
    my.data1 <- bind_rows(my.data,my.data0)
    
    my.w <- analyze.wavelet(
      my.data1,
      "k.detrended",
      loess.span = 0,
      dt = 1,
      dj = 1 / 20,
      lowerPeriod = 2,
      upperPeriod = 365,
      make.pval = F,
      verbose = F
    )
  },
  error = function(e) {
    skip <<- T
  })
  
  if (skip) {
    my.w <- NULL
    my.data$k.est <- my.data$k.trend 
    my.data$id <-rep(x,nrow(my.data)) 
    my.data$b14 <- 0
    my.data$b30 <- 0
    my.data$b90 <- 0
    my.data$b180 <- 0
    my.data$over180 <- 0
    return(my.data)
  } else{
    max.period <- max(my.w$Period)
    
    
    my.w$Power.pval <-
      # Significant p-value only if power > corresponding thresh_quant 
      ifelse(my.w$Power > thresh_quant[1:nrow(my.w$Power), 1:ncol(my.w$Power)], 0.01, 0.99)
    # Significant p-value can exist only when period greater or equal than distance from last measurament
    dist_matrix <- outer(my.w$Period, my.data1$dist_m, ">=")
    my.w$Power.pval <- ifelse(dist_matrix == F, 0.99, my.w$Power.pval)
    
    invisible(capture.output({
      
      
      r14 <- reconstruct(
        my.w,
        sel.period = c(2:pmin(14, max.period)),
        only.sig = T,
        only.ridge = F,
        plot.rec = F,
        lvl=0.09,
        plot.waves = F,
        rescale = F
      )
      
      my.data$b14 <- r14$series[[3]][1:nrow(my.data)]
      my.data$b14 <- replace(my.data$b14, NaN, 0)
      
      
      if (max.period >= 15) {
        
        
        r30 <- reconstruct(
          my.w,
          sel.period = c(15:pmin(30, max.period)),
          only.sig = T,
          only.ridge = F,
          plot.rec = F,
          lvl=0.09,
          plot.waves = F,
          rescale = F
        )
        my.data$b30 <- r30$series[[3]][1:nrow(my.data)]
        
        my.data$b30 <- replace(my.data$b30, NaN, 0)
        
        
        if (max.period > 30) {
          
          
          r90 <- reconstruct(
            my.w,
            sel.period = c(31:pmin(90, max.period)),
            only.sig = T,
            only.ridge = F,
            lvl=0.09,
            plot.rec = F,
            plot.waves = F,
            rescale = F
          )
          
          my.data$b90 <- r90$series[[3]][1:nrow(my.data)]
          my.data$b90 <- replace(my.data$b90, NaN, 0)
          
          
          if (max.period > 90) {
            
            
            
            r180 <- reconstruct(
              my.w,
              sel.period = c(91:180),
              only.sig = T,
              only.ridge = F,
              lvl=0.09,
              plot.rec = F,
              plot.waves = F,
              rescale = F
            )
            
            my.data$b180 <- r180$series[[3]][1:nrow(my.data)]
            
            my.data$b180 <- replace(my.data$b180, NaN, 0)
            
            
            if (max.period > 180) {
              
              rover180 <- reconstruct(
                my.w,
                sel.period = c(181:365),
                only.sig = T,
                only.ridge = F,
                lvl=0.09,
                plot.rec = F,
                plot.waves = F,
                rescale = F
              )
              
              my.data$over180 <- rover180$series[[3]][1:nrow(my.data)]
              
              my.data$over180 <- replace(my.data$over180, NaN, 0)
            } else  {
              my.data$over180 <- 0
            }
            
            
          }else{
            my.data$b180 <- 0
          }}
        else{
          my.data$b90 <- 0
        }
      }else{
        my.data$b30 <- 0
      }
      
      
      
      
      r <- reconstruct(
        my.w,
        sel.period = c(2:365),
        only.sig = T,
        only.ridge = F,
        lvl = 0.09,
        plot.rec = F,
        plot.waves = F)
      
      
    }))
    
  
    my.data$k.est <-replace_na(r$series[[3]][1:nrow(my.data)], 0) + my.data$k.trend 
    my.data$id <-rep(x,nrow(my.data)) 
    
    
    my.data %<>% mutate_all(replace_na, 0)
    
    my.data %<>% mutate_all(replace_na, 0) %>% 
      mutate(
             k.mis=case_when(dist_m==0~k))
    
    
    
    return(my.data)
  }
}


f <- function(x)x[1]
fillData <- function(x){
  y <- ave(x, cumsum(!is.na(x)), FUN = f)
  
  return(y)
}


landmark_data_mixed_wavelet <- function(tLM,w,longdata,longdataexp,survdata,fm,kn,formula.k){
  source("fun.R")
  require(plyr)
  require(doParallel)
  require(tidyverse)
  require(magrittr)
  require(rms)
  require(JMbayes)
  require(data.table)
  require(splines)
  
  # One considers all data collected before the landmark time point for 
  # subjects still at risk at the landmark
  LMlong <- longdata %>% 
    arrange(id,time) %>% 
    filter(time<=tLM & fup_days>tLM)
  
  LMlongexp <- longdataexp %>% 
    arrange(id,time) %>% 
    filter(time<=tLM & fup_days>tLM)
  
  # ID at risk at landmark time
  Ri <- LMlong %>% 
    distinct(id) %>% 
    pull()
  
  LMsurv <- survdata %>% 
    arrange(id) %>% 
    filter(id%in%Ri) %>% 
    mutate(time=tLM)
  
  # BLUPs and parameters
  b <- ranef(fm)
  b <- b[which(survdata$id%in%Ri),]
  # sigma <- fm$sigma
  betas <- fixef(fm)
  formula.random.k <- formula(paste0("~ns(time,knots=c(",paste(kn,collapse = ","),"))| id"))
  
  Z.pot <- data.matrix(cbind(rep(1,nrow(LMlongexp)),LMlongexp$z1,LMlongexp$z2,LMlongexp$z3,LMlongexp$z4))
  Xlong.pot <- cbind(model.matrix(formula(paste0("~",gsub("k ","",as.character(formula.k))[3],collapse="")),data=LMlongexp),LMlongexp$z1,LMlongexp$z2,LMlongexp$z3,LMlongexp$z4)
  nreps <- LMlongexp %>% dplyr::group_by(id) %>% dplyr::summarize(n=n()) %>% dplyr::ungroup() %>%  pull(n) 
  b_rep <- data.matrix(b[rep(seq_len(nrow(b)), each = unique(nreps)), ])
  b_rep <- apply(b_rep, 2, as.numeric)
  LMlongexp$k.trend <- as.vector(c(Xlong.pot %*%betas) + rowSums(Z.pot*b_rep ))
  
  
  Z.pot <- data.matrix(cbind(LMlongexp$dz1,LMlongexp$dz2,LMlongexp$dz3,LMlongexp$dz4))
  Xlong.pot <- data.matrix(cbind(LMlongexp$dz1,LMlongexp$dz2,LMlongexp$dz3,LMlongexp$dz4))
  nreps <- LMlongexp %>% dplyr::group_by(id) %>% dplyr::summarize(n=n()) %>% dplyr::ungroup() %>%  pull(n) 
  b_rep <- data.matrix(b[rep(seq_len(nrow(b)), each = unique(nreps)),2:5 ])
  b_rep <- apply(b_rep, 2, as.numeric)
  LMlongexp$k.dtrend <- as.vector(c(Xlong.pot %*%betas[(length(betas)-3):length(betas)]) + rowSums(Z.pot*b_rep))
  
  LMlongexp$k.detrended <- LMlongexp$k-LMlongexp$k.trend
  LMlongexp$k.detrended <- fillData(LMlongexp$k.detrended)
  
  
  # Set up parallel computing
  cl <- makeCluster(2)
  registerDoParallel(cl)
  opts <- list(preschedule=TRUE)
  clusterSetRNGStream(cl, 123)
  LMwaveletfit <- llply(Ri,
                        fitwavelet,
                        data=LMlongexp,
                        deriv=T,
                        .parallel=T,
                        .paropts = list(.options.snow=opts))
  
  LMlongwave <- do.call(bind_rows,LMwaveletfit)
  LMlongwave %<>% dplyr::mutate(over180_d=case_when(abs(over180)>0~1,
                                                    TRUE~0),
                                b180_d=case_when(abs(b180)>0~1,
                                                 TRUE~0),
                                b90_d=case_when(abs(b90)>0~1,
                                                TRUE~0),
                                b30_d=case_when(abs(b30)>0~1,
                                                TRUE~0),
                                b14_d=case_when(abs(b14)>0~1,
                                                TRUE~0),
                                over180=case_when(abs(over180)>0~over180,
                                                  TRUE~0),
                                b180=case_when(abs(b180)>0~b180,
                                               TRUE~0),
                                b90=case_when(abs(b90)>0~b90,
                                              TRUE~0),
                                b30=case_when(abs(b30)>0~b30,
                                              TRUE~0),
                                b14=case_when(abs(b14)>0~b14,
                                              TRUE~0),
                                over180_d2=case_when(over180>0~1,
                                                     over180<0~2,
                                                     TRUE~0),
                                b180_d2=case_when(b180>0~1,
                                                  b180<0~2,
                                                  TRUE~0),
                                b90_d2=case_when(b90>0~1,
                                                 b90<0~2,
                                                 TRUE~0),
                                b30_d2=case_when(b30>0~1,
                                                 b30<0~2,
                                                 TRUE~0),
                                b14_d2=case_when(b14>0~1,
                                                 b14<0~2,
                                                 TRUE~0)) %>%
    dplyr::rename(time=date)
  LMlongwave$tLM <- tLM
  LMlongwave %<>%dplyr::select(id,time,k.trend,k.est,k.dtrend,contains("b14"),contains("b30"),contains("b90"),contains("b180"),contains("over180")) 
  
  
  # Last measurament
  LMsurv$k <- LMlongexp %>% 
    filter(time==tLM) %>% 
    dplyr::pull(k)
  
  
  # Long-term 
  LMsurv$k.trend <- LMlongwave %>% 
    filter(time==tLM) %>% 
    dplyr::pull(k.trend)
  
  # Derivative Long-term 
  LMsurv$k.dtrend <- LMlongwave %>% 
    filter(time==tLM) %>% 
    dplyr::pull(k.dtrend)
  
  
  # Long-term + Short term oscillation
  LMsurv$k.est <- LMlongwave %>% 
    filter(time==tLM) %>% 
    dplyr::pull(k.est)
  
  # Oscillations at different frequencies
  LMsurv$b14_d2 <- LMlongwave %>% 
    filter(time==tLM) %>% 
    dplyr::pull(b14_d2)
  
  LMsurv$b30_d2 <- LMlongwave %>% 
    filter(time==tLM) %>% 
    dplyr::pull(b30_d2)
  
  LMsurv$b90_d2 <- LMlongwave %>% 
    filter(time==tLM) %>% 
    dplyr::pull(b90_d2)
  
  LMsurv$b180_d2 <- LMlongwave %>% 
    filter(time==tLM) %>% 
    dplyr::pull(b180_d2)
  
  LMsurv$over180_d2 <- LMlongwave %>% 
    filter(time==tLM) %>% 
    dplyr::pull(over180_d2)
  
  LMsurv$over180 <- LMlongwave %>% 
    filter(time==tLM) %>% 
    dplyr::pull(over180)
  
  LMsurv$b180 <- LMlongwave %>% 
    filter(time==tLM) %>% 
    dplyr::pull(b180)
  
  LMsurv$b90 <- LMlongwave %>% 
    filter(time==tLM) %>% 
    dplyr::pull(b90)
  
  LMsurv$b30 <- LMlongwave %>% 
    filter(time==tLM) %>% 
    dplyr::pull(b30)
  
  LMsurv$b14 <- LMlongwave %>% 
    filter(time==tLM) %>% 
    dplyr::pull(b14)
  
  LMsurv %<>%mutate(dht=ifelse(fup_days>w+tLM,0,dht),
                    fup_days_w=pmin(w+tLM,fup_days)) 
  
  return(list(datasurv=LMsurv,
              knots=kn))
  
  
}


cv_fun <- function(k,index,datalm , datalong,datasurv,w) {
  source("fun.R")
  Sys.setenv(JAGS_HOME = "C:\\Program Files\\JAGS\\JAGS-4.3.0\\")
  require(data.table)
  require(plyr)
  require(survival)
  require(splines)
  require(doSNOW)
  require(foreach)
  require(JMbayes)
  library(tidyverse)
  library(magrittr)
  
  
  cl <- makeCluster(2)
  registerDoParallel(cl)
  opts <- list(preschedule=TRUE)
  clusterSetRNGStream(cl, 123)
  
  # Define landmark times
  tLM <- seq(365, 365 * 5, by = 90)
  
  datasurv_test <- datasurv[datasurv$id%in%index[[k]],]
  datalong_test <- datalong[datalong$id%in%index[[k]],]
  datalm_test <- datalm[datalm$id%in%index[[k]],]
  
  datasurv_train <- datasurv[!(datasurv$id%in%index[[k]]),]
  datalong_train <- datalong[!(datalong$id%in%index[[k]]),]
  datalm_train <- datalm[datalm$id%in%index[[k]],]
  
  datasurv_train$age_stand <- scale(datasurv_train$age)
  datasurv_test$age_stand <- scale(datasurv_test$age)
  
  
  datalm_train$age_stand <- scale(datalm_train$age)
  datalm_train$k.dtrend_stand <- scale(datalm_train$k.dtrend * 100)
  datalm_train$k.trend_stand <- scale(datalm_train$k.trend)
  datalm_train$k_stand <- scale(datalm_train$k)
  datalm_train$hypo <- ifelse(datalm_train$k<3.5,1,0)
  datalm_train$hyper <- ifelse(datalm_train$k>5,1,0)
  
  datalm_test$age_stand <- scale(datalm_test$age)
  datalm_test$k.dtrend_stand <- scale(datalm_test$k.dtrend * 100)
  datalm_test$k.trend_stand <- scale(datalm_test$k.trend)
  datalm_test$k_stand <- scale(datalm_test$k)
  datalm_test$hypo <- ifelse(datalm_test$k<3.5,1,0)
  datalm_test$hyper <- ifelse(datalm_test$k>5,1,0)
  
  datasurv_train %<>%arrange(id)
  datalong_train %<>%arrange(id,time)
  
  print("test")
  ########
  # LOCF #
  ########
  
  
  ## Simple (ipl) locf
  
  LMsupercox_locf <-
    coxph(
      Surv(Tstart, fup_days_w, dht) ~ age_stand + gender + nyha_34  + lvef_classi +  somma_comor_nc3 +
        ns(k_stand, df = 4) + strata(time) + cluster(id),
      data = datalm_train
    )
  
  
  ########
  # LOCF 2#
  ########
  
  
  ## Simple (ipl) locf
  
  LMsupercox_locf2 <-
    coxph(
      Surv(Tstart, fup_days_w, dht) ~ age_stand + gender + nyha_34  + lvef_classi +  somma_comor_nc3 +
        hypo+hyper + strata(time) + cluster(id),
      data = datalm_train
    )
  
  #########
  # MIXED #
  #########
  
  
  ## Simple (ipl) l
  LMsupercox_mixed <-
    coxph(
      Surv(Tstart, fup_days_w, dht) ~ age_stand + gender + nyha_34 + lvef_classi +  somma_comor_nc3 +
        ns(k.trend_stand,df=4)+ + k.dtrend_stand + strata(time) + cluster(id),
      data = datalm_train,
      method = "breslow"
    )
  
  
  ####################
  # MIXED & WAVELETS #
  ####################
  
  
  ## Simple (ipl)
  LMsupercox_mixed_wave <-
    coxph(
      Surv(Tstart, fup_days_w, dht) ~ age_stand + gender + nyha_34   + lvef_classi +  somma_comor_nc3 +
        ns(k.trend_stand,df=4)+
        k.dtrend_stand +
        b14_d2 +
        b30_d2 +
        b90_d2 +
        b180_d2 +
        over180_d2 +
        strata(time) + cluster(id),
      data = datalm_train,
      method = "breslow"
    )
  
  
  ###############
  # JOINT MODEL #
  ###############
  
  # Joint Model- Current value+ Slope
  # Cox model for the event death
  CoxFit <-
    coxph(
      Surv(fup_days, dht) ~ age_stand + gender + nyha_34  + lvef_classi +  somma_comor_nc3,
      data = datasurv_train,
      x = TRUE
    )
  
  # a linear mixed model for potassium
  fp <- lme(
    k ~ gender + age + irc + nyha_34 + ns(time, df = 4),
    data = datalong_train,
    random = ~ ns(time, df = 4) | id,
    control = lmeControl(opt = 'optim')
  )
  
  f1 <- function(x, data) {
    cbind(
      "1" = ns(x, df = 4)[, 1],
      "2" = ns(x, df = 4)[, 2],
      "3" = ns(x, df = 4)[, 3],
      "4" = ns(x, df = 4)[, 4]
    )
  }
  
  f2 <- function(x, data) {
    x
  }
  
  fns <- list("value" = f1,
              "extra" = f2)
  
  dform = list(
    fixed = ~ 0 + dns(time, df = 4),
    random = ~ 0 + dns(time, df = 4),
    indFixed = 6:9,
    indRandom = 2:5
  )
  
  # the joint model  with current value structure and slope
  fitJoint <- jointModelBayes(
    fp,
    CoxFit,
    timeVar = "time",
    param = "td-both",
    transFun = fns,
    extraForm = dform,
    control = list(n.iter = 10000,
                   n.burnin = 1000)
  )
  
  print("ok models")
  
  # Obtain performance measures
  
  # LOCF
  p_locf <- llply(
    tLM,
    perfLM,
    object = LMsupercox_locf,
    data = datalm_test,
    w = w,
    .parallel=T,
    .paropts = list(.options.snow=opts)
  )
  
  
  print("ok locf")
  
  
  # LOCF 2
  p_locf2 <- llply(
    tLM,
    perfLM,
    object = LMsupercox_locf2,
    data = datalm_test,
    w = w,
    .parallel=T,
    .paropts = list(.options.snow=opts)
  )
  
  
  print("ok locf2")
  
  #MIXED
  p_mixed <- llply(
    tLM,
    perfLM,
    object = LMsupercox_mixed,
    data = datalm_test,
    w = w,
    .parallel=T,
    .paropts = list(.options.snow=opts)
  )
  
  
  
  print("ok mixed")
  
  # MIXED-WAVELETS
  p_mixed_wave <- llply(
    tLM,
    perfLM,
    object = LMsupercox_mixed_wave,
    data = datalm_test,
    w =w,
    .parallel=T,
    .paropts = list(.options.snow=opts)
  )
  
  
  
  print("ok wave")
  
  # JOINT
  #debug(perfJM)
  p_joint <-
    llply(
      tLM,
      perfJM,
      object = fitJoint,
      data = datalm_test,
      w = w,
      .parallel=F
      #,
      #.paropts = list(.options.snow=opts)
    )
  
  print("ok joint")
  
  tranf_perf <- function(obj){
    extract<- function(landmark) lapply(obj, function(y) y[[landmark]])
    
    res<- sapply(1:9, function(my.year)
      do.call(rbind, extract(my.year)))   
    
    return(res)
    
  }
  
  p_locf_long <- tranf_perf(p_locf)
  p_locf2_long <- tranf_perf(p_locf2)
  p_mixed_long <- tranf_perf(p_mixed)
  p_mixed_wave_long <- tranf_perf(p_mixed_wave)
  p_joint_long <- tranf_perf(p_joint)
  
  
  data_perf <- data.frame(tLM=rep(tLM / 365, 5),
                          bs=c(p_locf_long[[1]],p_locf2_long[[1]],p_mixed_long[[1]],p_mixed_wave_long[[1]],p_joint_long[[1]]),
                          se.bs=c(p_locf_long[[2]],p_locf2_long[[2]],p_mixed_long[[2]],p_mixed_wave_long[[2]],p_joint_long[[2]]),
                          auc=c(p_locf_long[[4]],p_locf2_long[[4]],p_mixed_long[[4]],p_mixed_wave_long[[4]],p_joint_long[[4]]),
                          se.auc=c(p_locf_long[[5]],p_locf2_long[[5]],p_mixed_long[[5]],p_mixed_wave_long[[5]],p_joint_long[[5]]),
                          model=rep(1:5,each=17),
                          k=rep(k,17*5))
  
  
  c1 <- compare_perf(p_locf_long, p_mixed_wave_long,n=2981)
  datac1 <- as.data.frame(c1)
  colnames(datac1) <- names(c1)
  datac1$time <- tLM[1:17]
  datac1$k <- rep(k,17)
  
  c2 <- compare_perf(p_locf2_long, p_mixed_wave_long,n=2981)
  datac2 <- as.data.frame(c2)
  colnames(datac2) <- names(c2)
  datac2$time <- tLM[1:17]
  datac2$k <- rep(k,17)
  
  c3 <- compare_perf(p_mixed_long, p_mixed_wave_long,n=2981)
  datac3 <- as.data.frame(c3)
  colnames(datac3) <- names(c3)
  datac3$time <- tLM[1:17]
  datac3$k <- rep(k,17)
  
  
  c4 <- compare_perf(p_joint_long, p_mixed_wave_long,n=2981)
  datac4 <- as.data.frame(c4)
  colnames(datac3) <- names(c4)
  datac4$time <- tLM[1:17]
  datac4$k <- rep(k,17)
  
  
  return(list(data=data_perf,comp_locf=datac1,comp_locf2=datac2,comp_mixed=datac3,comp_joint=datac4))
  
}


perfJM <- function(Tstart,w,data,object){
  require(timeROC)
  require(JMbayes) 
  
  n <- length(unique(data$id))
  datasel <- data[data$time==Tstart,]
  datasel <- as.data.frame(datasel)
  Thoriz <- Tstart+w-0.5
  year <- as.character(unique(datasel$years))
  year2 <- setdiff(c("second", "third" , "fourth" ,"fifth" ),year)
  Tstart2 <- sample(unique(data$Tstart[data$years==year2[1]]),1)
  Tstart3 <- sample(unique(data$Tstart[data$years==year2[2]]),1)
  Tstart4 <- sample(unique(data$Tstart[data$years==year2[3]]),1)
  datasel2 <- data[data$Tstart%in%c(Tstart,Tstart2,Tstart3,Tstart4),]
  pred <- survfitJM(object,newdata = datasel,survTimes =Thoriz,last.time = Tstart,idVar="id",simulate=F)$summaries
  surv <- do.call(rbind,pred)[,2]
  #surv <- surv[which(datasel2$time==Tstart)]
  roc <- timeROC(T=datasel$fup_days_w,
                 delta=datasel$dht,
                 marker=1-surv,
                 cause=1,
                 weighting="marginal",
                 times=c(Thoriz),
                 iid=TRUE)
  
  bs <- BS(timepoints=c(Thoriz),
           times=datasel$fup_days_w,
           status=datasel$dht,
           pred=as.matrix(1-surv),
           cause=1)
  
  
  est.bs <- bs$BS # BS estimate
  se.bs <- bs$sd  # BS s.e. estimate
  Mat.iid.BS <-  c(rep(0,n-roc$n),as.vector(bs$iid))/(roc$n/n)  # BS iid decomposition
  est.auc <- roc$AUC[2] # AUC estimate
  sd.auc <- roc$inference$vect_sd_1[2]  # AUC s.e. estimate
  Mat.iid.AUC <- c( rep(0,n-roc$n),roc$inference$mat_iid_rep_1[,2]) # AUC iid decomposition   
  matrixstat <- roc$Stats[2,] # proportions of cases, controls, censored subjects within the prediction window and survivors at s+t
  CumInc <- roc$CumulativeIncidence[2]  # marginal cumulative incidence  = P(T<s+t,eta=1|T>s)
  Surv <- roc$survProb[2] # marginal survival probability = P(T>s+t|T>s)
  return(list(est.bs,se.bs,Mat.iid.BS,est.auc,sd.auc,Mat.iid.AUC,matrixstat,CumInc,Surv))
  
}



# Expected Brier score estimator iid representation 

# {{{ Inpus ;
# pred : prediction, e.g. P(T<t,cause=1|history) A MATRIX 
# timepoints : vector of time points for which we aim to compute the iid representations
# times : the vector of observed time-to-event
# status : 1=uncensored, 0=censored
# cause : cause for which we aim to compute the expected Brier score estimator
# }}}

# {{{ Outputs :
# matrix with iid representation for all time points
# }}}


# {{{ Functions 
# main function
BS <- function(timepoints,times,status,pred,cause=1,compute.iid=TRUE){ 
  
  require(pec)
  start_computation_time <- Sys.time()
  # define useful objects
  n <- length(times)
  n_times <- length(timepoints)
  timepoints <- timepoints[order(timepoints)]
  times_names <- paste("t=",timepoints,sep="")
  # output initialisation 
  BS <- rep(NA,n_times)
  CumInci <- rep(NA,n_times)
  surv <- rep(NA,n_times)
  Stats <- matrix(NA,nrow=n_times,ncol=4)
  hit1_all <- matrix(NA,nrow=n,ncol=n_times)
  hit2_all <- matrix(NA,nrow=n,ncol=n_times)
  epsilon_i <- matrix(NA,nrow=n,ncol=n_times)
  #adds name to outputs
  names(BS) <- times_names
  names(CumInci) <- times_names
  names(surv) <- times_names
  colnames(Stats) <- c("Cases","survivor at t","Other events at t","Censored at t")
  rownames(Stats) <- times_names
  colnames(epsilon_i) <- times_names
  colnames(hit1_all) <-  times_names
  colnames(hit2_all)  <- times_names 
  # we need to order to use the ipcw() function of the pec package
  #browser()
  order_T <- order(times)
  times <-  times[order_T]
  delta  <-  status[order_T]
  pred <-  pred[order_T,,drop=FALSE]
  #compute KM weights
  weights <- pec::ipcw(Surv(failure_time,status)~1,
                       data=data.frame(failure_time=times,status=as.numeric(delta!=0)),
                       method="marginal",times=timepoints,subjectTimes=times,subjectTimesLag=1)
  Mat_data <- cbind(times,delta,as.numeric(delta==0))
  colnames(Mat_data) <- c("T","delta","indic_Cens")
  # computate weights of cases
  Weights_cases_all <- 1/(weights$IPCW.subjectTimes*n)
  # compute KM censoring estimator iid representation
  if (compute.iid){ MatInt0TcidhatMCksurEff <- Compute.iid.KM(times,delta!=0)}
  # loop on all time points
  for(t in 1:n_times){
    Cases <- (Mat_data[,"T"]<= timepoints[t] &  Mat_data[,"delta"]==cause)
    Controls_1 <- (Mat_data[,"T"]> timepoints[t] )
    Controls_2 <- (Mat_data[,"T"]<= timepoints[t] &  Mat_data[,"delta"]!=cause & Mat_data[,"delta"]!=0)  
    # compute weights
    Weights_controls_1 <- rep(1/(weights$IPCW.times[t]*n),times=n)
    Weights_cases <- Weights_cases_all
    Weights_controls_2 <- Weights_cases_all
    Weights_cases[!Cases] <- 0
    Weights_controls_1[!Controls_1] <- 0
    Weights_controls_2[!Controls_2] <- 0   
    #compute outputs
    CumInci[t] <- c(sum(Weights_cases))
    surv[t] <- c(sum(Weights_controls_1))
    Stats[t,] <- c(sum(Cases),sum(Controls_1),sum(Controls_2),n-sum(Cases)-sum(Controls_1)-sum(Controls_2)) 
    hit1_all[,t] <- (Weights_controls_1*((pred[,t])^2))*n
    hit2_all[,t] <- (Weights_cases*((1-pred[,t])^2) + Weights_controls_2*((pred[,t])^2))*n
    BS[t] <- (sum(hit1_all[,t]) +sum(hit2_all[,t]))/n
    if (compute.iid){
      # compute 
      Int0tdMCsurEffARisk <- MatInt0TcidhatMCksurEff[max(which(Mat_data[,"T"]<=timepoints[t])),]   
      #browser()
      epsilon_i[,t] <- hit1_all[,t]+hit2_all[,t]-BS[t] + mean(hit1_all[,t])*Int0tdMCsurEffARisk +  colMeans(MatInt0TcidhatMCksurEff*hit2_all[,t])
    }
  } 
  #compute mean and sd of iid representation
  sd_all <- rep(NA,n_times)
  mean_all <- rep(NA,n_times)
  if (compute.iid){sd_all <- apply(epsilon_i,2,sd)/sqrt(n)
  mean_all <- apply(epsilon_i,2,mean)}
  #compute a table to print 
  print.tab <- cbind(Stats,BS,sd_all,mean_all)
  colnames(print.tab) <- c(colnames(Stats),"BS","sd","mean_iid")
  #compute the computation time
  stop_computation_time <- Sys.time()
  computation_time=difftime(stop_computation_time,start_computation_time,units="secs")
  
  out <- list(BS=BS,iid=epsilon_i,sd=sd_all,res=(hit1_all+hit2_all),
              CumulativeIncidence=CumInci,survProb=surv,n=n,Stats=Stats,print.tab=print.tab,timepoints=timepoints,
              computation_time=difftime(stop_computation_time,start_computation_time,units="secs"))
  class(out) <- "ipcwEBS"
  out 
}
# print function
print.ipcwEBS <- function(x,No.lines=5,...){
  tab_ou_print <- round(cbind(x$Stats,x$BS*100,x$sd*100),2)
  colnames(tab_ou_print) <- c("Cases","survivors","Other events","Censored","EBS(*100)","se(*100)")
  l <- length(x$timepoints)
  if(l<=No.lines){  print(tab_ou_print) }
  else{print(tab_ou_print[unique(round(quantile(1:length(x$timepoints),probs=seq(0,1,length.out=No.lines)),0)),])}
  cat("\n")
  cat("Total computation time :",round(x$computation_time,2)," secs. \n")
}
# }}}

# function to compute iid decomposition of KM estimator of the censoring survival distribution

# {{{ Input
#times : observed time
#status : 1 if non censored, 0 if censored
# }}}

# {{{ Output
#iid.mat : matrix of the iid representation of KM for all time in vector times
# }}}

# {{{ Function
Compute.iid.KM <- function(times,status){
  #browser()
  times <- times[order(times)]
  status <- status[order(times)] 
  n <- length(times)
  mat.data<-cbind(times,as.numeric(status==0))
  colnames(mat.data)<-c("T","indic.Cens")
  # compute the empirical survival function corresponding to the counting process 1(\tilde{eta}=0, \tilde{T}<=t)
  hatSdeltaCensTc<-1-cumsum(mat.data[,c("indic.Cens")])/n  
  # Build the matrix required for computing  dM_C(u) for all time u (all observed times \tilde{T}_i)
  temp1 <- cbind(mat.data[,c("T","indic.Cens")],1-(1:n)/n,hatSdeltaCensTc)
  temp1 <- rbind(c(0,0,1,1),temp1) # Add the first row corresponding to time t=0
  colnames(temp1)<-c("T","indic.Cens","hatSTc","hatSdeltaCensTc")
  # compute hazard function of the censoring
  lambdaC<-(temp1[-1,"indic.Cens"])/(n:1)  
  # Add the column of the hazard function of the censoring (equal to 0 at time t=0)
  temp1<-cbind(temp1,c(0,lambdaC))
  colnames(temp1)[ncol(temp1)]<-"lambdaC"
  # Cumulative hazard of censoring
  LambdaC<-cumsum(lambdaC)         
  # Add the column of the cumulative hazard function of the censoring (equal to 0 at time t=0)
  temp1 <- cbind(temp1,c(0,LambdaC))
  colnames(temp1)[ncol(temp1)]<-"LambdaC"
  temp2<-temp1[-1,]
  # compute  martingale of censoring \hat{M}_{C_i}(u) for all time u (all observed times \tilde{T}_i) using previous matrix
  # We obtain a matrix. Each column contains the vector of M_{C_i}(\tilde{T}_j) for  all j.
  hatMC<-matrix(NA,n,n)
  for (i in 1:n){
    hatMC[,i] <-temp2[i,2]*as.numeric(temp2[i,1]<=temp2[,"T"])- c(temp2[0:i,"LambdaC"], rep(temp2[i,6],(n-i)))
  }  
  # In order to draw martingale paths
  #matplot(mat.data[,"T"],hatMC,type="l")
  #lines(mat.data[,"T"],rowMeans(hatMC),lwd=5)  
  # Compute d \hat{M}_{C_i} (u) for all time u (all observed times \tilde{T}_i)
  dhatMC<-rbind(hatMC[1,],hatMC[-1,]-hatMC[-nrow(hatMC),])
  # Compute d \hat{M}_{C_i} (u)/(S_{\tilde{T}}(u)) for all time u (all observed times \tilde{T}_i)
  # We need this for integrals in the martingale representation of the Kaplan-Meier estimator of the censoring survival function
  # function to divide d \hat{M}_{C_i} (u) by (S_{\tilde{T}}(u))
  MulhatSTc<-function(v){
    n <- length(v)
    v/c(1,1-(1:(n-1))/n)      # c(1,1-(1:(n-1))/n) is the at risk probability (S_{\tilde{T}}(u))
  }
  # apply the function for each column (corresponding to the
  # vector M_{C_i}(u)  for all time u (all observed times \tilde{T}_i), 
  # time \tilde{T}_i corresponds to the i-th row of the matrix)
  dhatMCdivST<-apply(dhatMC,2,MulhatSTc)
  # Compute \int_0^{\tilde{T}_j} d{ \hat{M}_{C_l} (u) } / (S_{\tilde{T}}(u)) for each subject l, we compute for all time \tilde{T}_j.
  # l=column, j=row
  MatInt0TcidhatMCksurEff<-apply(dhatMCdivST,2,cumsum)  # (Remark : on of the row corresponds to the previous step...) 
  colnames(MatInt0TcidhatMCksurEff)<-paste("M_{C_",1:length(times),"}",sep="")
  rownames(MatInt0TcidhatMCksurEff)<-times  
  return(MatInt0TcidhatMCksurEff)  
}




perfLM <- function(Tstart,w,data,object){
  require(timeROC)
  require(survival)
  require(splines)
  
  n <- length(unique(data$id))
  datasel <- data[data$time==Tstart,]
  Thoriz <- Tstart+w-0.5
  year <- as.character(unique(datasel$years))
  year2 <- setdiff(c("second", "third" , "fourth" ,"fifth" ),year)
  Tstart2 <- sample(unique(data$Tstart[data$years==year2[1]]),1)
  Tstart3 <- sample(unique(data$Tstart[data$years==year2[2]]),1)
  Tstart4 <- sample(unique(data$Tstart[data$years==year2[3]]),1)
  datasel2 <- data[data$Tstart%in%c(Tstart,Tstart2,Tstart3,Tstart4),]
  sfit <- survfit(object, newdata = datasel2)
  pred <- summary(sfit, times = Thoriz,extend=T)
  surv <- pred$surv[which(datasel2$time==Tstart)]
  roc <- timeROC(T=datasel$fup_days_w,
                 delta=datasel$dht,
                 marker=1-surv,
                 cause=1,
                 weighting="marginal",
                 times=c(Thoriz),
                 iid=TRUE)
  
  bs <- BS(timepoints=c(Thoriz),
           times=datasel$fup_days_w,
           status=datasel$dht,
           pred=as.matrix(1-surv),
           cause=1)
  
  
  est.bs <- bs$BS # BS estimate
  se.bs <- bs$sd  # BS s.e. estimate
  Mat.iid.BS <-  c(rep(0,n-roc$n),as.vector(bs$iid))/(roc$n/n)  # BS iid decomposition
  est.auc <- roc$AUC[2] # AUC estimate
  sd.auc <- roc$inference$vect_sd_1[2]  # AUC s.e. estimate
  Mat.iid.AUC <- c( rep(0,n-roc$n),roc$inference$mat_iid_rep_1[,2]) # AUC iid decomposition 
  matrixstat <- roc$Stats[2,] # proportions of cases, controls, censored subjects within the prediction window and survivors at s+t
  CumInc <- roc$CumulativeIncidence[2]  # marginal cumulative incidence  = P(T<s+t,eta=1|T>s)
  Surv <- roc$survProb[2] # marginal survival probability = P(T>s+t|T>s)
  return(list(est.bs,se.bs,Mat.iid.BS,est.auc,sd.auc,Mat.iid.AUC,matrixstat,CumInc,Surv))
}


# {{{ function to compute confidence band
PerturbIIDQuantConfBand <- function(Estimate,iid,level=0.95,N.Monte.Carlo=2000){
  #browser()
  n.vis <- length(Estimate)
  n.subj <- nrow(iid)
  vect.Delta <- rep(NA,n.vis)
  Conf.band <- matrix(NA,nrow=2,ncol=n.vis)
  s.e <- apply(iid,2,sd)/sqrt(n.subj)
  # genrate the quantile
  for (i in 1:N.Monte.Carlo){
    temp1 <- iid*rnorm(n.subj)
    temp2 <- t(t(temp1)/s.e)
    vect.Delta[i] <- max(abs(colMeans(temp2)))   
  }
  C.alpha <- quantile(vect.Delta,level)
  return(C.alpha=C.alpha)
}


fix_res <- function(index){
  perf_data <- perf[,index]
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
  return(perf_data)
}


