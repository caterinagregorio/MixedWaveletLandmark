
rm(list=ls())

library(tidyverse)
library(gtsummary)
library(ggplot2)
library(ggpubr)
library(magrittr)
library(purrr)
library(xtable)
library(data.table)
library(nlme)


load("data\\long_res.RData")
source('utils.R')
source('fun.R')


#########
# FIG 1 #
#########


dht_lab <- c("Censored","Death")
names(dht_lab) <- c(0,1)

potassio_sample %<>% mutate(time_rev=(time-fup_days)/365) 

fig1 <- ggplot(potassio_sample)+
  geom_point(aes(time_rev,k),alpha=0.2,col="grey30")+
  geom_line(aes(time_rev,k,group=as.factor(id)),alpha=0.3,col="grey30",size=1)+
  geom_smooth(aes(time_rev,k,group=as.factor(dht)),col="red",size=1,method = "loess",se=F)+
  facet_grid(~dht,labeller = labeller(dht=dht_lab))+
  scale_y_continuous("K",limits = c(2,8))+
  scale_x_continuous("Time until end of observation period (Years)")+
  theme_light()+
  theme(text = element_text(face = "bold",
                            size=16))

fig1
ggsave(plot = fig1,filename = "fig/Fig1.jpg",device = "jpeg",dpi=300,
       width = 7,height = 4.5)



#############
# TABLE 2 ##
############

table2 <- tab_long(fp,2)
table2 %>% flextable()
print(xtable(table2,digits=4), include.rownames=FALSE)

#########
# FIG 2 #
#########

dht_lab <- c("Censored","Death")
names(dht_lab) <- c(0,1)

fig2 <- ggplot(potassio_sample)+
  geom_point(aes(time_rev,k),alpha=0.2,col="grey30")+
  geom_line(aes(time_rev,k.trend,group=as.factor(id)),alpha=0.5,col="darkblue",size=1)+
  facet_grid(~dht,labeller = labeller(dht=dht_lab))+
  scale_y_continuous(expression(hat("K")),limits = c(2,8))+
  scale_x_continuous("Time until end of observation period (Years)")+
  theme_light()+
  theme(text = element_text(face = "bold",
                            size=16))
fig2
ggsave(plot = fig2,filename = "fig/Fig2.jpg",device = "jpeg",dpi=300,
       width = 7,height = 4.5)

#########
# FIG 3 #
#########


fig3 <- ggplot(potassio_sample)+
  geom_point(aes(time_rev,res),alpha=0.2,col="grey30")+
  geom_line(aes(time_rev,res,group=as.factor(id)),alpha=0.5,col="darkblue",size=1)+
  facet_grid(~dht,labeller = labeller(dht=dht_lab))+
  scale_y_continuous("Residuals from LME")+
  scale_x_continuous("Time until end of observation period (Years)")+
  theme_light()+
  theme(text = element_text(face = "bold",
                            size=16))
fig3
ggsave(plot = fig3,filename = "fig/Fig3.jpg",device = "jpeg",dpi=300,width = 7,height = 4.5)

#########
# FIG 4 #
#########


fig4 <- ggplot(var)+
  geom_smooth(aes(dist,variog),col="grey30",span=0.6,se=F,fullrange=F)+
  geom_point(aes(dist,variog),col="grey30")+
  scale_x_continuous("Time distance (Days)",breaks = c(1,7,14,30,90,180),limits = c(0,180))+
  scale_y_continuous("Semi-variogram",limits = c(0,1.5))+
  theme_light()+
  theme(text = element_text(face = "bold",
                            size=16))


fig4

ggsave(plot = fig4,filename = "fig/Fig4.jpg",device = "jpeg",dpi=300,width = 7,height = 4.5)

#########
# FIG 5 #
#########

gg51 <- gg_periodogram( 1164807,data=potassioext)
gg52 <- plot_osc(gg51$data)
fig5 <- egg::ggarrange(gg52,gg51$plot,ncol=1,nrow=2)
fig5
ggsave(plot = fig5,filename = "fig/Fig5.jpg",device = "jpeg",dpi=300,width = 7,height = 4.5)




      