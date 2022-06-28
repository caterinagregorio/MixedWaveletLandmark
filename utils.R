

tab_long <- function(fit,digits){
  s <- summary(fit)
  ci <- intervals(fit,which = 'fixed')$fixed
  table <- rownames_to_column( as.data.frame(s$tTable[,c(1,2,5)]), "Variable")
  table$lower <- ci[,1]
  table$upper <- ci[,3]
  table$ci <- paste(round(table$lower,digits),round(table$upper,2),sep=";")
  table$Value <- round(table$Value,digits)
  table$Std.Error <- round(table$Std.Error,digits)
  table$`p-value` <- round(table$`p-value`,4)
  table <- table[,-c(5,6)]
  return(table)
}


tab_surv <- function(fit,digits){
  s <- summary(fit)
  tab <- s$coefficients[,c(1,4,6)]
  tab <- data.frame(tab)
  tab$lower <- round(tab[,1]-1.96*tab[,2],digits)
  tab$upper <- round(tab[,1]+1.96*tab[,2],digits)
  tab$coef <- round(tab$coef,digits)
  tab$robust.se <- round(tab$robust.se,digits)
  tab$p.value <- round(tab$Pr...z..,4)
  tab$ci <- paste(tab$lower,tab$upper,sep = ";")
  tab <- tab[,-c(3,4,5)]
  return(tab)
 
}


gg_periodogram <- function(id,data){
  require(WaveletComp)
  
  datasel <- data[data$id==id,]
  datasel$date <- datasel$time
  my.w <- analyze.wavelet(datasel, 
                          "k.detrended",
                          loess.span = 0,
                          dt = 1, 
                          dj=1/20,
                          lowerPeriod = 2,
                          upperPeriod = 365,
                          make.pval = F,
                          verbose=F)
  
  im <- wt.image(my.w)
  
  
  datapower <- data.frame(time=rep(my.w$axis.1,each=151),
                          period=rep(my.w$axis.2,times=length(my.w$axis.1)),
                          power=c(my.w$Power))
  
  
  
  newcol <- colorRampPalette(c("white","#225560","#7EA8BE","#FCC05F","#BB4430"))
  ncols <- 1000
  col <- newcol(ncols)#apply the function to get 100 colours
  
  gg <- ggplot(datapower)+
    geom_tile(aes(time/365,period,fill=power))+
    scale_fill_gradientn("Power",colours = col)+
    scale_x_continuous("Time (Years)",breaks = seq(0,max(data$time/365),by=1))+
    scale_y_continuous("Period",breaks = log2(c(2,4,7,14,30,90,180,365)),labels =c(2,4,7,14,30,90,180,365) )+
    #scale_height_continuous(guide=F)+
    ggtitle("Periodogram")+
    theme_classic()+
    theme(text = element_text(face="bold",
                              size=16),
          legend.position = "bottom")
  
  return(list(plot=gg,data=datasel))
  
}


plot_osc <- function(data){
  
  id <- unique(data$id)
  data <- fitwavelet(id,data,deriv = F)
  

    gg<- ggplot(data)+
      geom_line(aes(date/365,k.est,group=as.factor(id)),size=1,col="grey60")+
      geom_line(aes(date/365,k.trend,group=as.factor(id)),size=1,col="black",linetype="dashed")+
      geom_point(aes(date/365,k.mis),shape=8)+
      scale_y_continuous(expression(italic("f")~hat("K")))+
      scale_x_continuous("Time (Years)",breaks = seq(0,max(data$date/365),by=1))+
      theme_light()+
      theme(text = element_text(face = "bold",
                                size=16))
  return(gg)
    
  
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

gg_dynsurv <- function(data,fit,id,times,datalong,w){
  
  datasellong <- datalong[datalong$id==id,]
  datapredall <- data.frame(id=NA,time=NA,k=NA,surv=NA,lower=NA,upper=NA,times=NA)
  for (i in 1:length(times)){
    datasel <- data[data$time==times[i] & data$id==id,]
    Thoriz <- seq(times[i],times[i]+w-0.5,by=7)
    datapred <- datasellong[datasellong$time<=times[i],]
    
    pred <- summary(survfit(fit, newdata = datasel), times = Thoriz,extend=T)
    datapred <- bind_rows(datapred,data.frame(id=rep(id,length(Thoriz)),
                                              time=Thoriz,
                                              surv=pred$surv[pred$strata==paste0("time=",times[i])],
                                              lower=pred$lower[pred$strata==paste0("time=",times[i])],
                                              upper=pred$upper[pred$strata==paste0("time=",times[i])]))
    datapred$times <- times[i]
    datapredall <- bind_rows(datapredall,datapred)
  }
  
  
  datapredall <- as.data.frame(datapredall[!is.na(datapredall$times),])
  
  times_lab <- paste0("t=",round(times/365,1))
  names(times_lab) <- times
  gg <- ggplot(datapredall) +
    geom_line(aes(time/365,k),col="grey60")+
    geom_point(aes(time/365,k),shape=8)+
    geom_ribbon(aes(time/365,ymin=lower/0.19,ymax=upper/0.19),fill="ivory4",alpha=0.5)+
    geom_line(aes(time/365,surv/0.19),size=1)+
    geom_vline(aes(xintercept=times/365),linetype="dashed")+
    scale_x_continuous("Follow-up (Years)")+
    facet_grid(~times,labeller = labeller(times=times_lab),scales = "free_x")+
    scale_y_continuous("K (mmol/L)",
                       sec.axis = sec_axis(~.*0.19, name = "Survival Probability "))+
    theme_bw()+
    theme(legend.position = "null",
          plot.margin = unit(c(0.05,0.05,0.05,0.05), "cm"),
          text = element_text(face = "bold",
                              size=16)
    )
  
  return(list(plot=gg,data=datapredall))
}


gg_dynhfdks <- function(data,fit,times,id){
  
  times_lab <- paste0("t=",times)
  names(times_lab) <- times
  
  pred <- predict(fit,type="terms",newdata=data,terms=c(6:11))
  score <-  rowSums(pred)
  datasel <- data[data$Tstart%in%times & data$id==id,]
  score <- scale(score)
  HFDKS <- pmax(pmin(10,round(10*(score-median(score)+IQR(score)*1.5)/(median(score)+IQR(score)*1.5-median(score)+IQR(score)*1.5),0)),0)
  HFDKS2 <- 10*((score-min(score))/(max(score)-min(score)))
  datasel$HFDKS <-HFDKS[data$Tstart%in%times & data$id==id]
  
  gg <- ggplot(datasel)+
    geom_line(aes(time/365,HFDKS),size=1,col="grey90")+
    geom_point(aes(time/365,HFDKS,col=HFDKS),size=3)+
    scale_x_continuous("Time (Years)")+
    scale_y_continuous("HFDKS")+
    ggtitle("Heart Functional Dynamic potassium(K) Score")+
    scale_color_gradientn("",colours = c("#001219","#005F73","#0A9396","#94D2BD",
                                         "#E9D8A6","#EE9B00","#CA6702","#BB3E03",
                                         "#AE2012","#9B2226"),limits=c(0,10))+
    theme_classic2()+
    theme(text = element_text(face = "bold",
                              size=16),
          legend.position = "right")
  
  prepGradient <- function(x,y,spacing=max(y)/100){
    stopifnot(length(x)==length(y))
    df <- data.frame(x=x,y=y)
    newDf = data.frame(x=NULL,y=NULL,z=NULL)
    for (r in 1:nrow(df)){
      n <- floor(df[r,"y"]/spacing)
      for (s in c(1:n)){
        tmp <- data.frame(x=df[r,"x"],y=spacing,z=s*spacing)
        newDf <- rbind(newDf,tmp)
      }
      tmp <- data.frame(x=df[r,"x"],y=df[r,"y"]%%spacing,z=df[r,"y"])
      newDf <- rbind(newDf,tmp)
    }
    return(newDf)
  }
  
  df2 <- prepGradient(datasel$time,datasel$HFDKS)
  df2$w <-1 
  
  gg2 <- ggplot(df2) + 
    geom_bar(aes(x=w,y=y,fill=z),stat="identity") +
    scale_x_continuous("",breaks = c())+
    scale_y_continuous("HFDKS",breaks = c(0,5,10))+
    ggtitle("Heart Functional Dynamic potassium(K) Score")+
    scale_fill_gradientn("",colours = c("#001219","#005F73","#0A9396","#94D2BD",
                                        "#E9D8A6","#EE9B00","#CA6702","#BB3E03",
                                        "#AE2012","#9B2226"),
                         limits=c(0,10))+
    coord_flip()+
    facet_grid(~x,labeller = labeller(x=times_lab))+
    theme_classic2()+
    theme(text = element_text(face = "bold",
                              size=16),
          legend.position = "right",
          plot.margin = unit(c(0.05,0.05,0.05,0.05), "cm"))
  
  
  return(list(gg,gg2,df2))
  
}




gg_dynhfdks_red <- function(data,times,id){
  
  times_lab <- paste0("t=",times)
  names(times_lab) <- times
  
  datasel <- data[data$Tstart%in%times & data$id==id,]
  
  
  gg <- ggplot(datasel)+
    geom_line(aes(time/365,HFDKS),size=1,col="grey90")+
    geom_point(aes(time/365,HFDKS,col=HFDKS),size=3)+
    scale_x_continuous("Time (Years)")+
    scale_y_continuous("HFDKS")+
    ggtitle("Heart Functional Dynamic potassium(K) Score")+
    scale_color_gradientn("",colours = c("#001219","#005F73","#0A9396","#94D2BD",
                                         "#E9D8A6","#EE9B00","#CA6702","#BB3E03",
                                         "#AE2012","#9B2226"),limits=c(0,10))+
    theme_classic2()+
    theme(text = element_text(face = "bold",
                              size=16),
          legend.position = "right")
  
  prepGradient <- function(x,y,spacing=max(y)/100){
    stopifnot(length(x)==length(y))
    df <- data.frame(x=x,y=y)
    newDf = data.frame(x=NULL,y=NULL,z=NULL)
    for (r in 1:nrow(df)){
      n <- floor(df[r,"y"]/spacing)
      for (s in c(1:n)){
        tmp <- data.frame(x=df[r,"x"],y=spacing,z=s*spacing)
        newDf <- rbind(newDf,tmp)
      }
      tmp <- data.frame(x=df[r,"x"],y=df[r,"y"]%%spacing,z=df[r,"y"])
      newDf <- rbind(newDf,tmp)
    }
    return(newDf)
  }
  
  df2 <- prepGradient(datasel$time,datasel$HFDKS)
  df2$w <-1 
  
  gg2 <- ggplot(df2) + 
    geom_bar(aes(x=w,y=y,fill=z),stat="identity") +
    scale_x_continuous("",breaks = c())+
    scale_y_continuous("HFDKS",breaks = c(0,5,10))+
    ggtitle("Heart Functional Dynamic potassium(K) Score")+
    scale_fill_gradientn("",colours = c("#001219","#005F73","#0A9396","#94D2BD",
                                        "#E9D8A6","#EE9B00","#CA6702","#BB3E03",
                                        "#AE2012","#9B2226"),
                         limits=c(0,10))+
    coord_flip()+
    facet_grid(~x,labeller = labeller(x=times_lab))+
    theme_classic2()+
    theme(text = element_text(face = "bold",
                              size=16),
          legend.position = "right",
          plot.margin = unit(c(0.05,0.05,0.05,0.05), "cm"))
  
  
  return(list(gg,gg2,df2,data))
  
}





calc_dynhfdks_red <- function(data,fit,indexes){
  
  
  pred <- predict(fit,type="terms",newdata=data,terms=c(indexes),reference="zero")
  score <-  rowSums(pred)
  score <- exp(score)
  quart <- quantile(score[score>1],probs=c(0.25,0.5,0.75))
  print(quart)
  HFDKS <- case_when(score==1~0,
                     score>1 & score<=quart[1]~ 1,
                     score>quart[1] & score<=quart[2]~2,
                     score>quart[2] & score<=quart[3]~3,
                     score>quart[3]~4)
  
  
  
  data$HFDKS <-HFDKS
  data$HFDKS_cat <- factor( data$HFDKS,labels=c("Vey Low","Low","Intermediate","High","Very High"))
  
  
  data %<>% select(id,time, contains("k"), contains("d2"),dht) 
  return(data)
  
}



gg_dynhfdks_cat <- function(data,times,id){
  
  times_lab <- paste0("t=",times)
  names(times_lab) <- times
  
  datasel <- data[data$time%in%times & data$id==id,]
  
  datasel$HFDKS <- factor( datasel$HFDKS,levels = c("0","1","2","3","4"))
  gg <- ggplot(datasel)+
    geom_point(aes(x=0.5,y=1,col=HFDKS),size=10)+
    facet_grid(~time)+
    scale_x_continuous("",limits = c(0,1))+
    scale_y_continuous("HFDKS",breaks = NULL,limits = c(0.9,1.2))+
    ggtitle("Heart Functional Dynamic potassium(K) Score")+
    scale_color_manual("",values =  c("#3a8ebb","#86A861","#ffba08","#c95203","#a10500"),
                       labels=c("Very Low","Low","Intermediate","High","Very High"),
                       drop=F)+
    theme_void()+
    theme(text = element_text(face = "bold",
                              size=16),
          legend.position = "bottom",
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          plot.margin = unit(c(0.05,0.05,0.05,0.05), "cm"))
  
  
  return(gg)
  
}





