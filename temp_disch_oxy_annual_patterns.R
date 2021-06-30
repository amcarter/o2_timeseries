# Patterns in Temp, DO, Hydrology across Powell Center sites
library(lubridate)
library(tidyverse)

dat <- readRDS("../data/gapPhilled_data.rds")
sites <- unique(dat$sitecode)

meds <- dat %>% 
  group_by(sitecode) %>%
  filter(length(!is.na(DO_obs)) > 3000,
         median(discharge, na.rm = T) >= 50,
         median(discharge, na.rm = T) <= 500) 


for (site %in% unique(sites)) {
  sdat <- dat[dat$sitecode==site,]
  yrs<-unique(sdat$Year) 
  if yrs < 3 next
  for(y %in% yrs){
    sydat <- sdat[sdat$Year==y,]
    if(length(!is.na(sydat$DO_obs))/365<0.6) next
    dat2 <- bind_rows(dat2, sydat)
  }
}

#annual temp patterns
sites <- c( "NHC","PM","CBP", "WB", "WBP", "UNHC" )

pdf(width=7, height=2, file="tempswings_byMonth_NHC.pdf", onefile=T)
par(mfrow = c(3,2),mar=c(4,4,1,4))
for(s in sites) {
  dat <- readRDS(paste0("modeled/",s,".rds"))

  avg <- dat$data %>% group_by(date)%>%
    summarise(max.temp = max(temp.water, na.rm=T),
              min.temp = min(temp.water, na.rm=T),
              mean.temp=mean(temp.water, na.rm=T))
  avg$amp <- avg$max.temp-avg$min.temp
  avg$month <- month(avg$date)
  avg[is.infinite(avg$amp),2:5]<- NA
  mons <- avg %>% group_by(month)%>%
    summarise(temp.amp=mean(amp, na.rm=T),
              temp.max = mean(max.temp, na.rm=T),
              temp.min=mean(min.temp, na.rm=T),
              temp.mean=mean(mean.temp, na.rm=T))
  
  plot(mons$month, mons$temp.mean, type="n", axes=F, ylab="", xlab="", ylim = c(5,30))
  polygon(c(mons$month, 12,1), c(mons$temp.mean,0,0), col=alpha("darkblue",.4),border=NA)
  axis(4, col="darkblue", col.axis="darkblue")
  mtext("mean temp C", side=4, line=2.2, col="darkblue", cex=.8)
  mtext(s,side=3, adj=.9, line=-1.5)
  par(new=T)
  plot(mons$month, mons$temp.amp, type="l", lwd=2, axes=F, xlab="", ylab = "delta temp C", xpd=T, ylim = c(1,4))
  axis(2)
  axis(1,labels=F )
  t <- seq(2,12, by=2)
  text(x = t, y = .8 , labels = as.character(month(t, label=T, abbr=T)), 
       pos = 1, xpd=NA,  cex=.9)
}
  
dev.off()
