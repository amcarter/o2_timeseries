# read in all powell center data and create an RDS of the good stuff
library(tidyverse)
library(plotfunctions)

setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/powell_center_analyses")
pc_metab <- readRDS("../data/gapPhilled_data.rds")
pc_metab$regionID <- substr(pc_metab$sitecode, 1,2)
pc_metab$siteID <- substr(pc_metab$sitecode, 4,nchar(pc_metab$sitecode))
pc_metab$DO_pctsat <- pc_metab$DO_obs/pc_metab$DO_sat
sites<- unique(pc_metab$sitecode)

pc_metab$ER<- -pc_metab$ER/quantile(pc_metab$ER,.5, na.rm=T)
pc_metab$discharge<- pc_metab$discharge/quantile(pc_metab$discharge,.5, na.rm=T)

metadat <- read_csv("../data/all_basic_site_data.csv")

# Calculate exceedence curve for a site

plot_exceedence <- function(DO, col="black", lty=1, plot=FALSE){
  DO <- DO[which(!is.na(DO))]
  DO<- sort(DO, decreasing=F)
  index <- seq(1:length(DO))
  freq <- 100*(index/(1+length(DO)))
  dat <- data.frame(freq=freq,
                    return_int_days = 100/freq,
                    DO=DO)
  if(plot){
    lines(freq,DO, lwd = 2, col = col, lty=lty)
  }else return(dat)
}


par(mfrow = c(1,2), mar = c(4,4,3,0))
DOrange <- range(pc_metab$DO_obs, na.rm=T)
Qrange <- range(pc_metab$discharge, na.rm=T)


# Find which sites experience hypoxia
freq_hypoxia  <-  pc_metab %>%
      filter(!is.na(DO_pctsat)) %>%
      select(DO_pctsat, sitecode, Year) %>%      
      group_by(sitecode, Year) %>%
      summarize(p_hx50= sum(DO_pctsat <= .50 ) / length(DO_pctsat),
                p_hx25= sum(DO_pctsat <= .25) / length(DO_pctsat),
                p_anx= sum(DO_pctsat <= 0) / length(DO_pctsat)) 


hpx_sites <- freq_hypoxia[freq_hypoxia$p_hx50!=0,]
sites <- unique(hpx_sites$sitecode)


return_ints <- data.frame()
for(site in sites){
  sdat <- pc_metab[pc_metab$sitecode==site,]
  yrs<-unique(sdat$Year) 
 
 # plot(c(0,100),c(0,130), type="n", xlab = "exceedence", ylab="DO (mg/L)", bty="n")
  tmp <- data.frame(site=site,
                    year=yrs,
                    T_hx50=NA,
                    T_hx25=NA,
                    T_anx=NA,
                    Q_hx50=NA,
                    Q_hx25=NA,
                    Q_anx=NA,
                    Qmin=NA,
                    Qmax=NA)
  for(y in 1:length(yrs)){
    sydat <- sdat[sdat$Year==yrs[y],]
    tmp[y,"Qmin"]<- min(sydat$discharge, na.rm=T)
    tmp[y,"Qmax"]<- max(sydat$discharge, na.rm=T)
    exc <- plot_exceedence(sydat$DO_pctsat*100, lty=y)
    minDO <- exc$DO[1]
    if(minDO<=50){
      w<- which(exc$DO<=50)
      tmp[y,"T_hx50"]<- exc$return_int_days[w[length(w)]]
    }
    if(minDO<=25){
      w<- which(exc$DO<=25)
      tmp[y,"T_hx25"]<- exc$return_int_days[w[length(w)]]
    }
    if(minDO<=0){
      w<- which(exc$DO<=0)
      tmp[y,"T_anx"]<- exc$return_int_days[w[length(w)]]
    }
    
  }
#  plot(c(0,100),c(min(sdat$discharge, na.rm=T), max(sdat$discharge, na.rm=T)),log="y", type="n", xlab = "exceedence", ylab="discharge (m3s)", bty="n")
  for(y in 1:length(yrs)){
    sydat <- sdat[sdat$Year==yrs[y],]
    exc<-  plot_exceedence(sydat$discharge, lty=y)
    if(!is.na(tmp[y,"T_hx50"])){
      tmp[y,"Q_hx50"]<-exc$DO[exc$return_int_days==tmp[y,"T_hx50"]]
    }
    if(!is.na(tmp[y,"T_hx25"])){
      tmp[y,"Q_hx25"]<-exc$DO[exc$return_int_days==tmp[y,"T_hx25"]]
    }
    if(!is.na(tmp[y,"T_anx"])){
      tmp[y,"Q_anx"]<-exc$DO[exc$return_int_days==tmp[y,"T_anx"]]
    }
  }
  # legend("topright",legend=yrs,lty=seq(1,length(yrs)) , bty="n")
  # mtext(site, 3, xpd=NA, outer=T, line=-2.5)
  return_ints <- bind_rows(return_ints, tmp)
}

summary_ints <- return_ints %>% group_by(site) %>%
  summarize(n_years=length(year),
            cv_Q=sd(Q_hx50, na.rm=T)/mean(Q_hx50, na.rm=T),
            cv_hx_ri=sd(T_hx50, na.rm=T)/mean(T_hx50, na.rm=T))
plot(summary_ints$cv_hx_ri, summary_ints$cv_Q, cex=summary_ints$n_years / 3)

for(site in sites){
  a<- c(return_ints$T_anx[return_ints$site==site],return_ints$T_hx25[return_ints$site==site],return_ints$T_hx50[return_ints$site==site])
  if(length(which(!is.na(a)))==1) next
  a<- max(a,na.rm=T)
  plot(return_ints$Q_hx50[return_ints$site==site],return_ints$T_hx50[return_ints$site==site],
       cex=2, log="x", xlab="Q (m3s)", ylab="Hypoxia Return Interval (days)",
       main=paste0(site,", n = ",nrow(return_ints[return_ints$site==site,])," years"),
       xlim=c(min(return_ints$Qmin[return_ints$site==site],na.rm=T),max(return_ints$Qmax[return_ints$site==site],na.rm=T)),
       ylim=c(0,a+1))
  points(return_ints$Q_hx25[return_ints$site==site],return_ints$T_hx25[return_ints$site==site],
         pch=19, col="grey50", cex=2)
  points(return_ints$Q_anx[return_ints$site==site],return_ints$T_anx[return_ints$site==site],
         pch=19, cex=2)
  legend("topright", c("<50%","<25%","0%"), col=c("black","grey50","black"), pch=c(1,19,19),  cex=1.5)
}


calc_Qeff<- function(dat, ox){
  dat <- dat[!is.na(dat$discharge),]
  hh<-hist(log(dat$discharge),plot=F)
  bins <- hh$breaks
  dat$bin <- cut(log(dat$discharge),bins, labels=1:(length(bins)-1))
  tmp <- dat %>% group_by(bin)%>%
    summarise(pr_hx=length(which(DO_pctsat<=ox))/length(DO_pctsat))
  tmp$bin <- as.integer(tmp$bin)
  tmp<- left_join(data.frame(bin=seq(1:(length(bins)-1))),tmp)
  tmp$Q <- exp(hh$mids)
  tmp$p_Q <- hh$density
  tmp$p_Qhx <- tmp$p_Q*tmp$pr_hx
 # plot(tmp$Q, tmp$p_Qhx, log="x")
  if(sum(tmp$pr_hx, na.rm=T)==0){Qeff=NA} else{
  Qeff <- tmp$Q[which.max(tmp$p_Qhx)]}
  return(Qeff)
}

calc_QeffvsDO <- function(sdat){  
  Qeff <- data.frame(ox=seq(0,1.2,by=.05),
                     Qeff=NA)
  for(i in 1:nrow(Qeff)){
    Qeff[i,"Qeff"]<-calc_Qeff(sdat, Qeff$ox[i])  
  }
  return(Qeff)
}


#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c("white","black"))
pc_metab$K600[pc_metab$K600 >= 25]<- 25
pc_metab$color <- rbPal(10)[as.numeric(cut(pc_metab$K600,breaks = 10))]

#This adds a column of color values
# based on the y collection times
par(mar = c(4,4,2,1), mfrow=c(2,2))
for(site in sites){
  sdat<- pc_metab[pc_metab$sitecode==site,]
 
  Qeff<- calc_QeffvsDO(sdat)
  plot( Qeff$Qeff,Qeff$ox, main=site, xlim = range(sdat$discharge, na.rm=T),ylim=c(0,1.2), type="l", log="x",
        xlab = "effective discharge (m3s)", ylab = "dissolved oxygen threshold (% sat)", )
  points(sdat$discharge, sdat$DO_pctsat, pch=20, col = alpha(1,.2))
  lines(Qeff$Qeff, Qeff$ox, lwd=2)
  # gradientLegend(valRange=c(0,120),
  #                color = rbPal(10), 
  #                pos = c(max(sdat$discharge, na.rm=T)*5/4-min(sdat$discharge, na.rm=T)/4,.3,
  #                        max(sdat$discharge, na.rm=T)*26/20-min(sdat$discharge, na.rm=T)*6/20,.7),side=4,length = 0.5, 
  #                depth = 0.1, inside = FALSE, coords=T,fit.margin=F,
  #                pos.num = NULL, n.seg =4, border.col = "black", dec = 0)
  plot(sdat$discharge, sdat$ER, pch=20, log="x", main="ER")
  plot(sdat$discharge, sdat$GPP, pch=20, log="x", main="GPP")
  plot(sdat$discharge, sdat$K600, pch=20, log="x", main="K")
}


for(site in sites){
  sdat<- pc_metab[pc_metab$sitecode==site,]
  print(ggplot(sdat, aes(GPP, ER, col = Date))+
    geom_point()+
    ggtitle(site)+
    geom_abline(slope = -1, intercept = 0)
    # scale_color_gradientn(colors = rainbow(5))
 )
}

