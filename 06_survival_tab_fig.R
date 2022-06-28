rm(list=ls())

library(tidyverse)
library(gtsummary)
library(ggplot2)
library(ggpubr)
library(magrittr)
library(purrr)
library(xtable)
library(data.table)
library(survival)
library(splines)
library(flextable)

load("data/04_LM model_cv.RData")
load("data/long_res.RData")
source('utils.R')


############
# TABLE 3  #
############

table3 <- rownames_to_column(tab_surv(LMsupercox_mixed_wave1,digits = 2),"Variable")
table3 %>% flextable()
print(xtable(table3,digits=4),include.rownames=F)



#############
# FIGURE 6  #
#############

perf_summ %<>%mutate(Model=case_when(model==1~"CM1",
                                     model==2~"CM2",
                                     model==3~"CM3",
                                     model==5~"CM4",
                                     model==4~"M")) 

gg_disc <- ggplot(perf_summ)+
  geom_line(aes(tLM,aucm,col=as.factor(Model)),size=1)+
  geom_point(aes(tLM,aucm,col=as.factor(Model)),size=2)+
  scale_color_manual("",values=c("Forestgreen","grey60","black","Darkblue","Indianred"))+
  scale_fill_manual("",values=c("Forestgreen","grey60","black","Darkblue","Indianred"))+
  scale_x_continuous("Time point (Years)",breaks = seq(1,5,by=0.5))+
  scale_y_continuous("AUC(t) - 6 months horizon")+
  theme_light()+
  theme(text = element_text(face = "bold",
                            size=16))
gg_disc
 

gg_cal <-ggplot(perf_summ)+
  geom_line(aes(tLM,bsm*100,col=as.factor(Model)),size=1)+
  geom_point(aes(tLM,bsm*100,col=as.factor(Model)),size=2)+
  scale_color_manual("",values=c("Forestgreen","grey60","black","Darkblue","Indianred"))+
   scale_x_continuous("Time point (Years)",breaks = seq(1,5,by=0.5))+
   scale_y_continuous("Brier Score (%) - 6 months horizon",limits = c(3,6.5))+
  theme_light()+
  theme(text = element_text(face = "bold",
                             size=16))


fig6 <- ggarrange(gg_disc,gg_cal,ncol=2,nrow=1,common.legend = T)
fig6

ggsave(plot = fig6,filename = "fig/Fig6.jpg",device = "jpeg",dpi=300,width = 7,height = 4.5)



#################
#  FIGURE 7 #####
#################

LMsupercox_mixed_wavefin <-
  coxph(
    Surv(Tstart, fup_days_w, dht) ~ age_stand + gender + nyha_34   + lvef_classi +  somma_comor_nc3 +
      b14_d2 +
      b30_d2 +
      b90_d2 +
      b180_d2 +
      over180_d2+
      strata(time) + cluster(id),
    data = LMdata,
    method = "breslow"
  )

summary(LMsupercox_mixed_wavefin)
LMdata_HFKS <- calc_dynhfdks_red(LMdata,LMsupercox_mixed_wavefin,c(6:10))
table(LMdata_HFKS$HFDKS)

potassiosel <- potassio %>% select(id,time,k)
gg_predsurv <- gg_dynsurv(data=LMdata,
                          w=180,
                          datalong = potassiosel,
                          fit=LMsupercox_mixed_wavefin,
                          id=1164807   ,
                          times = c(365,455,635, 815,1715))
gg_predscore <- gg_dynhfdks_cat(data=LMdata_HFKS,
                                id=1164807   ,
                                times = c(365,455,635, 815,1715))


fig7 <- ggarrange(gg_predscore,gg_predsurv$plot,ncol=1,nrow=2,heights = c(0.5,2),
                  common.legend=T,legend = "bottom",align = "v")
fig7
ggsave(plot = fig7,filename = "fig/Fig7.jpg",device = "jpeg",dpi=300,width = 10,height = 4.5)


