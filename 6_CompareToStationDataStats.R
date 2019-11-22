#########
# Author: Alex Crawford
# Date Modified: 20 May 2019 (added figure creation)
# Purpose: Compare MERRA-2 precipitation data to data from various stations.
#########
library(gridExtra)
library(plotrix)
library(ggplot2)
library(reshape2) # for melt

theme1 <- theme(
  panel.background = element_rect(fill="white",color="grey20"),
  plot.background = element_rect(fill="white"),
  plot.title = element_text(color="black",size=12, hjust=0.5),
  axis.text = element_text(size=7),
  axis.title =element_text(size=9)
)

### LOAD DATA ###
path = "/Volumes/Miranda/RainOnSnow/StationComparisons/"
figpath = "/Volumes/Miranda/RainOnSnow/Figures/Histogram_DailyPrecip"
vs = c(11,12,13,15,16,17,19,20,21)
pt = 2.54 #(in mm)
ra = "MERRA2"
names = c("Aniak","Pargon Creek","Kelly Station","Fort Yukon Other","Teuchet Creek",
          "Indian Pass","Monument Creek","Cooper Lake","Kenai Moose Pens","Mt Ryan",
          "Bethel","Nome","Anchorage","Juneau","Fairbanks","Utqiagvik",
          "Prudhoe Bay","Coleville Village","Kotzebue","Fort Yukon","Galena")

### CALCULATE COMPARISONS ########
ndf <- data.frame(Name=names,P_r0=NA,P_p0=NA,P_sM20=NA,P_sST0=NA,P_RMSE0=NA,P_r1=NA,P_p1=NA,P_sM21=NA,P_sST1=NA,P_RMSE1=NA,
                  T_r0=NA,T_p0=NA,T_sM20=NA,T_sST0=NA,T_RMSE0=NA,T_r1=NA,T_p1=NA,T_sM21=NA,T_sST1=NA,T_RMSE1=NA,
                  perRA0=NA,perST0=NA,perRA1=NA,perST1=NA,nPrcpRA=NA,nPrcpST=NA,PrcpRA=NA,PrcpST=NA)
for (v in vs){
  # Read in Data
  df <- read.csv(paste0(path,names[v],"_",ra,".csv"))
  df <- df[is.finite(df$PRCP),]
  
  # Compare Counts
  ndf[ndf$Name == names[v],"nPrcpST"] <- sum(df$PRCP > 0)
  ndf[ndf$Name == names[v],"nPrcpRA"] <- sum(df$MERRA2Precip > 2.54)
  ndf[ndf$Name == names[v],"PrcpST"] <- mean(df$PRCP)*365.25
  ndf[ndf$Name == names[v],"PrcpRA"] <- sum(df[df$MERRA2Precip > 2.54,"MERRA2Precip"])/39
  
  # Compare correlation of total precip
  ct0 <- cor.test(df$PRCP,df$MERRA2Precip,method = "spearman")
  
  ndf[ndf$Name == names[v],"P_r0"] <- ct0$estimate
  ndf[ndf$Name == names[v],"P_p0"] <- ct0$p.value
  ndf[ndf$Name == names[v],"P_sM20"] <- sd(df$MERRA2Precip,na.rm = T) 
  ndf[ndf$Name == names[v],"P_sST0"] <- sd(df$PRCP,na.rm = T)
  ndf[ndf$Name == names[v],"P_RMSE0"] <- sqrt(mean((df$PRCP-df$MERRA2Precip)^2))
  
  # Compare correlation of total precip when at least 0.1 in. are recorded in Reanalysis
  ct1 <- cor.test(df[df$MERRA2Precip >= pt,]$PRCP,df[df$MERRA2Precip >= pt,]$MERRA2Precip,method = "spearman")
  
  ndf[ndf$Name == names[v],"P_r1"] <- ct1$estimate
  ndf[ndf$Name == names[v],"P_p1"] <- ct1$p.value
  ndf[ndf$Name == names[v],"P_sM21"] <- sd(df[df$MERRA2Precip >= pt,]$MERRA2Precip,na.rm = T) 
  ndf[ndf$Name == names[v],"P_sST1"] <- sd(df[df$MERRA2Precip >= pt,]$PRCP,na.rm = T)
  ndf[ndf$Name == names[v],"P_RMSE1"] <- sqrt(mean((df[df$MERRA2Precip >= pt,]$PRCP-df[df$MERRA2Precip >= pt,]$MERRA2Precip)^2))
  
  # Compare correlation of TMAX
  ct3 <- cor.test(df$TMAX,df$MERRA2TMax,method = "spearman")
  
  ndf[ndf$Name == names[v],"T_r0"] <- ct3$estimate
  ndf[ndf$Name == names[v],"T_p0"] <- ct3$p.value
  ndf[ndf$Name == names[v],"T_sM20"] <- sd(df$MERRA2TMax,na.rm = T) 
  ndf[ndf$Name == names[v],"T_sST0"] <- sd(df$TMAX,na.rm = T)
  ndf[ndf$Name == names[v],"T_RMSE0"] <- sqrt(mean((df$TMAX-df$MERRA2TMax)^2))
  
  # Compare correlation of TMAX when at least 0.1 in. are recorded in Reanalysis
  ct4 <- cor.test(df[df$MERRA2Precip >= pt,]$TMAX,df[df$MERRA2Precip >= pt,]$MERRA2TMax,method = "spearman")
  
  ndf[ndf$Name == names[v],"T_r1"] <- ct4$estimate
  ndf[ndf$Name == names[v],"T_p1"] <- ct4$p.value  
  ndf[ndf$Name == names[v],"T_sM21"] <- sd(df[df$MERRA2Precip >= pt,]$MERRA2TMax,na.rm = T) 
  ndf[ndf$Name == names[v],"T_sST1"] <- sd(df[df$MERRA2Precip >= pt,]$TMAX,na.rm = T)
  ndf[ndf$Name == names[v],"T_RMSE1"] <- sqrt(mean((df[df$MERRA2Precip >= pt,]$TMAX-df[df$MERRA2Precip >= pt,]$MERRA2TMax)^2))
  
  # Compare Comission/Omission
  rasp0 <- sum( (df$PRCP > 0) & (df$MERRA2Precip > 0)) # Number of precip events shared by both sources
  rasp1 <- sum( (df$PRCP >= pt) & (df$MERRA2Precip >= pt)) # Number of notable precip events shared
  
  ndf[ndf$Name == names[v],"perRA0"] <- rasp0 / sum(df$MERRA2Precip > 0) # % of MERRA 2 events shared
  ndf[ndf$Name == names[v],"perST0"] <- rasp0 / sum(df$PRCP > 0) # % of station events shared
  
  ndf[ndf$Name == names[v],"perRA1"] <- rasp1 / sum(df$MERRA2Precip >= pt) # % of MERRA 2 events shared
  ndf[ndf$Name == names[v],"perST1"] <- rasp1 / sum(df$PRCP >= pt) # % of station events shared
}

##### PLOTTING TAYLOR DIAGRAMS #####
colorlist = c("red","blue","green3","orange","black","magenta","brown","salmon","green4","grey50","gold")
shapelist = c(15,15,16,16,18,17,17,17,15,18,18)

# Precipitation - Total
i = 0
for (v in vs){
  i = i+1
  df <- read.csv(paste0(path,names[v],"_",ra,".csv"))
  df <- df[is.finite(df$PRCP),]
  
  if (v == vmin){
    taylor.diagram(ref=df$PRCP, model=df$MERRA2Precip, normalize=T, ref.sd=F, sd.arcs = T, show.gamma = T,
              main = "MERRA2 v. GHCN-D Precipitation", ylab= "Normalized Standard Deviation",cex=0.8,cex.lab=0.8,
              ngamma = 5, pch=shapelist[i], col=colorlist[i], grad.corr.lines = c(0.2,0.4,0.6,0.8,0.9,0.99))
  }
  if (v != vmin){
    taylor.diagram(ref=df$PRCP, model=df$MERRA2Precip, normalize=T, add = T, pch=shapelist[i], col=colorlist[i])
  } }

# Precipitation - > 0.1 in.
i = 0
for (v in vs){
  i = i+1
  df <- read.csv(paste0(path,names[v],"_",ra,".csv"))
  df <- df[is.finite(df$PRCP),]
  
  if (v == vs[1]){
    taylor.diagram(ref=df[df$MERRA2Precip >= pt,]$PRCP, model=df[df$MERRA2Precip >= pt,]$MERRA2Precip, normalize=T, ref.sd=F, sd.arcs = T, show.gamma = T,
    main = "MERRA2 v. GHCN-D > 0.1 in. Precip.", ylab= "Normalized Standard Deviation",cex=0.8,cex.lab=0.8,
    ngamma = 5, pch=shapelist[i], col=colorlist[i], grad.corr.lines = c(0.2,0.4,0.6,0.8,0.9,0.99))
  }
  if (v != vs[1]){
    taylor.diagram(ref=df[df$MERRA2Precip >= pt,]$PRCP, model=df[df$MERRA2Precip >= pt,]$MERRA2Precip, normalize=T, add = T, pch=shapelist[i], col=colorlist[i])
    
  } }

# Temperature
i = 0
for (v in vs){
  i = i+1
  df <- read.csv(paste0(path,names[v],"_",ra,".csv"))
  df <- df[is.finite(df$PRCP),]
  
  if (v == vs[1]){
    taylor.diagram(ref=df$TMAX, model=df$MERRA2TMax, normalize=T, ref.sd=F, sd.arcs = T, show.gamma = T,
                   main = "MERRA2 v. GHCN-D TMax", ylab= "Normalized Standard Deviation",cex=0.8,cex.lab=0.8,
                   ngamma = 5, pch=shapelist[i], col=colorlist[i], grad.corr.lines = c(0.2,0.4,0.6,0.8,0.9,0.99))
  }
  if (v != vs[1]){
    taylor.diagram(ref=df$TMAX, model=df$MERRA2TMax, normalize=T, add = T, pch=shapelist[i], col=colorlist[i])
    
  } }

######## Taylor Diagram Logged ########
# Precipitation
i = 0
for (v in vs){
  i = i+1
  df <- read.csv(paste0(path,names[v],"_",ra,".csv"))
  df <- df[is.finite(df$PRCP),]
  
  df$PRCPLOG <- ifelse(df$PRCP > 0, log(df$PRCP, base = 10), NA)
  df$M2PRCPLOG <- ifelse(df$MERRA2Precip > 0, log(df$MERRA2Precip, base = 10), NA)
  
  
  if (v == vmin){
    taylor.diagram(ref=df$PRCPLOG, model=df$M2PRCPLOG, normalize=T, ref.sd=F, sd.arcs = T, show.gamma = T,
                   main = "MERRA2 v. GHCN-D Precip.", ylab= "Normalized Standard Deviation",cex=0.8,cex.lab=0.8,
                   ngamma = 5, pch=shapelist[i], col=colorlist[i], grad.corr.lines = c(0.2,0.4,0.6,0.8,0.9,0.99))
  }
  if (v != vmin){
    taylor.diagram(ref=df$PRCPLOG, model=df$M2PRCPLOG, normalize=T, add = T, pch=shapelist[i], col=colorlist[i])
    
  } }

# Precipitation - > 0.1 in.
i = 0
for (v in vs){
  i = i+1
  df <- read.csv(paste0(path,names[v],"_",ra,".csv"))
  df <- df[is.finite(df$PRCP),]
  
  df$PRCPLOG <- ifelse(df$PRCP > pt, log(df$PRCP, base = 10), NA)
  df$M2PRCPLOG <- ifelse(df$MERRA2Precip > pt, log(df$MERRA2Precip, base = 10), NA)
  
  
  if (v == vmin){
    taylor.diagram(ref=df[df$MERRA2Precip >= pt,]$PRCPLOG, model=df[df$MERRA2Precip >= pt,]$M2PRCPLOG, normalize=T, ref.sd=F, sd.arcs = T, show.gamma = T,
                   main = "MERRA2 v. GHCN-D > 0.1 in. Precip.", ylab= "Normalized Standard Deviation",cex=0.8,cex.lab=0.8,
                   ngamma = 5, pch=shapelist[i], col=colorlist[i], grad.corr.lines = c(0.2,0.4,0.6,0.8,0.9,0.99))
  }
  if (v != vmin){
    taylor.diagram(ref=df[df$MERRA2Precip >= pt,]$PRCPLOG, model=df[df$MERRA2Precip >= pt,]$M2PRCPLOG, normalize=T, add = T, pch=shapelist[i], col=colorlist[i])
    
  } }

##### PLOT HISTOGRAMS ########
# Kotzebue
df <- read.csv(paste0(path,names[19],"_",ra,".csv"))
df <- df[is.finite(df$PRCP),]
df$PRCPLOG <- ifelse(df$PRCP > 0, log(df$PRCP, base = 10), NA)
df$M2PRCPLOG <- ifelse(df$MERRA2Precip > 0, log(df$MERRA2Precip, base = 10), NA)
p19 = ggplot(df, aes(x=M2PRCPLOG)) + geom_histogram(binwidth = 0.2, fill="blue",alpha=0.4) + theme1 +
  annotate('text', x=-3.5, y = 1650, label = names[19], hjust='left',size=3) + xlab("") +
  scale_x_continuous(limits=c(-3.5,2),breaks=seq(-3,2,1),labels = c(0.001,0.01,0.1,1,10,100)) +
  scale_y_continuous(limits=c(0,1700),breaks=seq(0,2000,400)) +
  geom_histogram(data=df, aes(x=PRCPLOG), binwidth= 0.2,fill="darkorange",alpha=0.8)

# Utqiagvik
df <- read.csv(paste0(path,names[16],"_",ra,".csv"))
df <- df[is.finite(df$PRCP),]
df$PRCPLOG <- ifelse(df$PRCP > 0, log(df$PRCP, base = 10), NA)
df$M2PRCPLOG <- ifelse(df$MERRA2Precip > 0, log(df$MERRA2Precip, base = 10), NA)
p16 = ggplot(df, aes(x=M2PRCPLOG)) + geom_histogram(binwidth = 0.2, fill="blue",alpha=0.4) + 
  theme1 + ylab("") + xlab("") +
  annotate('text', x=-3.5, y = 2500, label = names[16], hjust='left',size=3) +
  scale_x_continuous(limits=c(-3.5,2),breaks=seq(-3,2,1),labels = c(0.001,0.01,0.1,1,10,100)) +
  scale_y_continuous(limits=c(0,2600),breaks=seq(0,2600,600)) +
  geom_histogram(data=df, aes(x=PRCPLOG), binwidth= 0.2,fill="darkorange",alpha=0.8)

# Prudhoe Bay
df <- read.csv(paste0(path,names[17],"_",ra,".csv"))
df <- df[is.finite(df$PRCP),]
df$PRCPLOG <- ifelse(df$PRCP > 0, log(df$PRCP, base = 10), NA)
df$M2PRCPLOG <- ifelse(df$MERRA2Precip > 0, log(df$MERRA2Precip, base = 10), NA)
p17 = ggplot(df, aes(x=M2PRCPLOG)) + geom_histogram(binwidth = 0.2, fill="blue",alpha=0.4) + 
  theme1 + ylab("") + xlab("") +
  annotate('text', x=-3.5, y = 1670, label = names[17], hjust='left',size=3) +
  scale_x_continuous(limits=c(-3.5,2),breaks=seq(-3,2,1),labels = c(0.001,0.01,0.1,1,10,100)) +
  scale_y_continuous(limits=c(0,1730),breaks=seq(0,2000,400)) +
  geom_histogram(data=df, aes(x=PRCPLOG), binwidth= 0.2,fill="darkorange",alpha=0.8)

# Nome
df <- read.csv(paste0(path,names[12],"_",ra,".csv"))
df <- df[is.finite(df$PRCP),]
df$PRCPLOG <- ifelse(df$PRCP > 0, log(df$PRCP, base = 10), NA)
df$M2PRCPLOG <- ifelse(df$MERRA2Precip > 0, log(df$MERRA2Precip, base = 10), NA)
p12 = ggplot(df, aes(x=M2PRCPLOG)) + geom_histogram(binwidth = 0.2, fill="blue",alpha=0.4) + theme1 +
  annotate('text', x=-3.5, y = 1650, label = names[12], hjust='left',size=3) + xlab("") +
  scale_x_continuous(limits=c(-3.5,2),breaks=seq(-3,2,1),labels = c(0.001,0.01,0.1,1,10,100)) +
  scale_y_continuous(limits=c(0,1700),breaks=seq(0,2000,400)) +
  geom_histogram(data=df, aes(x=PRCPLOG), binwidth= 0.2,fill="darkorange",alpha=0.8)

# Fairbanks
df <- read.csv(paste0(path,names[15],"_",ra,".csv"))
df <- df[is.finite(df$PRCP),]
df$PRCPLOG <- ifelse(df$PRCP > 0, log(df$PRCP, base = 10), NA)
df$M2PRCPLOG <- ifelse(df$MERRA2Precip > 0, log(df$MERRA2Precip, base = 10), NA)
p15 = ggplot(df, aes(x=M2PRCPLOG)) + geom_histogram(binwidth = 0.2, fill="blue",alpha=0.4) + 
  theme1 + ylab("") + xlab("") +
  annotate('text', x=-3.5, y = 1650, label = names[15], hjust='left',size=3) +
  scale_x_continuous(limits=c(-3.5,2),breaks=seq(-3,2,1),labels = c(0.001,0.01,0.1,1,10,100)) +
  scale_y_continuous(limits=c(0,1700),breaks=seq(0,2000,400)) +
  geom_histogram(data=df, aes(x=PRCPLOG), binwidth= 0.2,fill="darkorange",alpha=0.8)

# Fort Yukon
df <- read.csv(paste0(path,names[20],"_",ra,".csv"))
df <- df[is.finite(df$PRCP),]
df$PRCPLOG <- ifelse(df$PRCP > 0, log(df$PRCP, base = 10), NA)
df$M2PRCPLOG <- ifelse(df$MERRA2Precip > 0, log(df$MERRA2Precip, base = 10), NA)
p20 = ggplot(df, aes(x=M2PRCPLOG)) + geom_histogram(binwidth = 0.2, fill="blue",alpha=0.4) + 
  theme1 + ylab("") + xlab("") +
  annotate('text', x=-3.5, y = 1650, label = names[20], hjust='left',size=3) +
  scale_x_continuous(limits=c(-3.5,2),breaks=seq(-3,2,1),labels = c(0.001,0.01,0.1,1,10,100)) +
  scale_y_continuous(limits=c(0,1700),breaks=seq(0,2000,400)) +
  geom_histogram(data=df, aes(x=PRCPLOG), binwidth= 0.2,fill="darkorange",alpha=0.8)

# Bethel
df <- read.csv(paste0(path,names[11],"_",ra,".csv"))
df <- df[is.finite(df$PRCP),]
df$PRCPLOG <- ifelse(df$PRCP > 0, log(df$PRCP, base = 10), NA)
df$M2PRCPLOG <- ifelse(df$MERRA2Precip > 0, log(df$MERRA2Precip, base = 10), NA)
p11 = ggplot(df, aes(x=M2PRCPLOG)) + geom_histogram(binwidth = 0.2, fill="blue",alpha=0.4) + 
  theme1 + xlab("Daily Precipitation (mm)") + 
  annotate('text', x=-3.5, y = 1650, label = names[11], hjust='left',size=3) +
  scale_x_continuous(limits=c(-3.5,2),breaks=seq(-3,2,1),labels = c(0.001,0.01,0.1,1,10,100)) +
  scale_y_continuous(limits=c(0,1700),breaks=seq(0,2000,400)) +
  geom_histogram(data=df, aes(x=PRCPLOG), binwidth= 0.2,fill="darkorange",alpha=0.8)

# Galena
df <- read.csv(paste0(path,names[21],"_",ra,".csv"))
df <- df[is.finite(df$PRCP),]
df$PRCPLOG <- ifelse(df$PRCP > 0, log(df$PRCP, base = 10), NA)
df$M2PRCPLOG <- ifelse(df$MERRA2Precip > 0, log(df$MERRA2Precip, base = 10), NA)
p21 = ggplot(df, aes(x=M2PRCPLOG)) + geom_histogram(binwidth = 0.2, fill="blue",alpha=0.4) + 
  theme1 + xlab("Daily Precipitation (mm)") + ylab("") +
  annotate('text', x=-3.5, y = 1650, label = names[21], hjust='left',size=3) +
  scale_x_continuous(limits=c(-3.5,2),breaks=seq(-3,2,1),labels = c(0.001,0.01,0.1,1,10,100)) +
  scale_y_continuous(limits=c(0,1700),breaks=seq(0,2000,400)) +
  geom_histogram(data=df, aes(x=PRCPLOG), binwidth= 0.2,fill="darkorange",alpha=0.8)

# Anchorage
df <- read.csv(paste0(path,names[13],"_",ra,".csv"))
df <- df[is.finite(df$PRCP),]
df$PRCPLOG <- ifelse(df$PRCP > 0, log(df$PRCP, base = 10), NA)
df$M2PRCPLOG <- ifelse(df$MERRA2Precip > 0, log(df$MERRA2Precip, base = 10), NA)
p13 = ggplot(df, aes(x=M2PRCPLOG)) + geom_histogram(binwidth = 0.2, fill="blue",alpha=0.4) + 
  theme1 + xlab("Daily Precipitation (mm)") + ylab("") +
  annotate('text', x=-3.5, y = 1650, label = names[13], hjust='left',size=3) +
  scale_y_continuous(limits=c(0,1700),breaks=seq(0,2000,400)) +
  scale_x_continuous(limits=c(-3.5,2),breaks=seq(-3,2,1),labels = c(0.001,0.01,0.1,1,10,100)) +
  geom_histogram(data=df, aes(x=PRCPLOG), binwidth= 0.2,fill="darkorange",alpha=0.8)

ptotal <- grid.arrange(p19,p16,p17,p12,p15,p20,p11,p21,p13,nrow=3)
ggsave("Histogram_DailyPrecip_Panel9.png",ptotal,'png',figpath,dpi=300,width=6.5,height=5.5,units="in")

##### PLOT RAINFALL V. TEMP ########
for (v in vmin:vmax){
  df <- read.csv(paste0(path,names[v],"_",ra,".csv"))
  df <- df[is.finite(df$PRCP),]
  
  # ggplot(df, aes(x=TMAX,y=PRCP-(SNOW/10))) + geom_point()
  # ggplot(df, aes(x=MERRA2TMax,y=MERRA2Precip-MERRA2Snowfall)) + geom_point(size=0.3) + 
  #   scale_x_continuous(limits=c(-30,0)) + scale_y_continuous(limits = c(0,15))

  ddf <- data.frame(MERRA2TMax = seq(-40,40,2), MERRA2Precip=0, MERRA2Snowfall=0, MERRA2Rainfall=0,MERRA2TMaxCount=0,MERRA2PrecipCount=0,MERRA2SnowCount=0,MERRA2RainCount=0)
  for (i in 1:nrow(ddf)){
    ddf[i,"MERRA2Precip"] <- mean(df[df$MERRA2TMax > ddf[i,"MERRA2TMax"]-1 & df$MERRA2TMax <= ddf[i,"MERRA2TMax"]+1,"MERRA2Precip"])
    ddf[i,"MERRA2Snowfall"] <- mean(df[df$MERRA2TMax > ddf[i,"MERRA2TMax"]-1 & df$MERRA2TMax <= ddf[i,"MERRA2TMax"]+1,"MERRA2Snowfall"])
    ddf[i,"MERRA2Rainfall"] <- ddf[i,"MERRA2Precip"] - ddf[i,"MERRA2Snowfall"] 
    ddf[i,"MERRA2TMaxCount"] <- nrow(df[df$MERRA2TMax > ddf[i,"MERRA2TMax"]-1 & df$MERRA2TMax <= ddf[i,"MERRA2TMax"]+1,])
    ddf[i,"MERRA2PrecipCount"] <- nrow(df[df$MERRA2TMax > ddf[i,"MERRA2TMax"]-1 & df$MERRA2TMax <= ddf[i,"MERRA2TMax"]+1 & df$MERRA2Precip > 0,])
    ddf[i,"MERRA2SnowCount"] <- nrow(df[df$MERRA2TMax > ddf[i,"MERRA2TMax"]-1 & df$MERRA2TMax <= ddf[i,"MERRA2TMax"]+1 & df$MERRA2Snowfall > 0,])
    ddf[i,"MERRA2RainCount"] <- nrow(df[df$MERRA2TMax > ddf[i,"MERRA2TMax"]-1 & df$MERRA2TMax <= ddf[i,"MERRA2TMax"]+1 & (df$MERRA2Precip-df$MERRA2Snowfall) > 0,])
  }
  ggplot(ddf, aes(x=MERRA2TMax, y=MERRA2Precip)) + geom_col(fill="blue") +
    geom_col(data = ddf, aes(x=MERRA2TMax, y=MERRA2Rainfall),fill="darkgreen") + 
    scale_x_continuous(limits=c(-40,28),breaks=seq(-40,40,10)) + theme1 + ggtitle(paste0("",names[v])) + 
    xlab("Maximum Daily Temperature (°C)") + ylab("Average Daily Precipitation (mm)")
  ggsave(paste0(figpath,"/Histogram_DailyPrecip_v_TMax_",names[v],".png"),dpi = 300,units = "in",width=3,height=3)

  ggplot(ddf, aes(x=MERRA2TMax, y=MERRA2PrecipCount)) + geom_col(fill="blue") +
    geom_col(data = ddf, aes(x=MERRA2TMax, y=MERRA2RainCount),fill="darkgreen") + 
    scale_x_continuous(limits=c(-40,28),breaks=seq(-40,40,10)) + theme1 + ggtitle(paste0("",names[v])) + 
    xlab("Maximum Daily Temperature (°C)") + ylab("# of Precip Days")
  ggsave(paste0(figpath,"/Histogram_PrecipDays_v_TMax_",names[v],".png"),dpi = 300,units = "in",width=3,height=3)
  
  # ggplot(ddf, aes(x=MERRA2TMax, y=MERRA2PrecipCount/MERRA2TMaxCount)) + geom_col(fill="blue") +
  #   geom_col(data = ddf, aes(x=MERRA2TMax, y=MERRA2RainCount/MERRA2TMaxCount),fill="darkgreen") + 
  #   scale_x_continuous(limits=c(-40,28),breaks=seq(-40,40,10)) + theme1 + ggtitle(paste0("",names[v])) + 
  #   xlab("Maximum Daily Temperature (°C)") + ylab("# of Precip Days")
  # ggsave(paste0(figpath,"/Histogram_PrecipDays_v_TMaxDays_",names[v],".png"),dpi = 300,units = "in",width=3,height=3)
  
  # ggplot(ddf, aes(x=MERRA2TMax, y=MERRA2Rainfall/MERRA2Precip)) + geom_col(fill="darkgreen") +
    # scale_x_continuous(limits=c(-28,28)) + theme1 + ggtitle(paste0("",names[v])) + 
    # xlab("Maximum Daily Temperature (°C)") + ylab("% of Precip as Rainfall (mm)")
}

# MERRA 2 #####
# Kotzebue
df <- read.csv(paste0(path,names[19],"_",ra,".csv"))
df$MERRA2Rainfall <- df$MERRA2Precip-df$MERRA2Snowfall
ddf <- data.frame(MERRA2TMax = seq(-40,40,2),MERRA2SnowCount=0,MERRA2MixedCount=0,MERRA2RainCount=0)
for (i in 1:nrow(ddf)){
  ddf[i,"MERRA2SnowCount"] <- nrow(df[df$MERRA2TMax > ddf[i,"MERRA2TMax"]-1 & df$MERRA2TMax <= ddf[i,"MERRA2TMax"]+1 & df$MERRA2Snowfall >= 0.254 & df$MERRA2Rainfall < 0.254,])
  ddf[i,"MERRA2MixedCount"] <- nrow(df[df$MERRA2TMax > ddf[i,"MERRA2TMax"]-1 & df$MERRA2TMax <= ddf[i,"MERRA2TMax"]+1 & df$MERRA2Snowfall >= 0.254 & df$MERRA2Rainfall >= 0.254,])
  ddf[i,"MERRA2RainCount"] <- nrow(df[df$MERRA2TMax > ddf[i,"MERRA2TMax"]-1 & df$MERRA2TMax <= ddf[i,"MERRA2TMax"]+1 & df$MERRA2Snowfall < 0.254 & df$MERRA2Rainfall >= 0.254,])
}
edf <- melt(ddf,"MERRA2TMax")
edf$value2 <- edf$value/nrow(df)*365.25 # convert to a #/year measure
p19 <- ggplot() + geom_bar(data=edf,aes(x=MERRA2TMax,y=value2,fill=variable),stat='identity',position='stack',show.legend=FALSE) + 
  scale_fill_manual(values=c("#0000FF", "#00DCDC","#027C02")) + 
  scale_x_continuous(limits=c(-40,28),breaks=seq(-40,40,10)) + scale_y_continuous(limits=c(0,32),breaks=seq(0,40,8)) +
  annotate('text', x=-38, y = 30.5, label = names[19], hjust='left',size=3) +
  ylab("Precip Days/Year") + theme1 + xlab("")

# Utqiagvik
df <- read.csv(paste0(path,names[16],"_",ra,".csv"))
df$MERRA2Rainfall <- df$MERRA2Precip-df$MERRA2Snowfall
ddf <- data.frame(MERRA2TMax = seq(-40,40,2),MERRA2SnowCount=0,MERRA2MixedCount=0,MERRA2RainCount=0)
for (i in 1:nrow(ddf)){
  ddf[i,"MERRA2SnowCount"] <- nrow(df[df$MERRA2TMax > ddf[i,"MERRA2TMax"]-1 & df$MERRA2TMax <= ddf[i,"MERRA2TMax"]+1 & df$MERRA2Snowfall >= 0.254 & df$MERRA2Rainfall < 0.254,])
  ddf[i,"MERRA2MixedCount"] <- nrow(df[df$MERRA2TMax > ddf[i,"MERRA2TMax"]-1 & df$MERRA2TMax <= ddf[i,"MERRA2TMax"]+1 & df$MERRA2Snowfall >= 0.254 & df$MERRA2Rainfall >= 0.254,])
  ddf[i,"MERRA2RainCount"] <- nrow(df[df$MERRA2TMax > ddf[i,"MERRA2TMax"]-1 & df$MERRA2TMax <= ddf[i,"MERRA2TMax"]+1 & df$MERRA2Snowfall < 0.254 & df$MERRA2Rainfall >= 0.254,])
}
edf <- melt(ddf,"MERRA2TMax")
edf$value2 <- edf$value/nrow(df)*365.25 # convert to a #/year measure
p16 <- ggplot() + geom_bar(data=edf,aes(x=MERRA2TMax,y=value2,fill=variable),stat='identity',position='stack',show.legend=FALSE) + 
  scale_fill_manual(values=c("#0000FF", "#00DCDC","#027C02")) + 
  scale_x_continuous(limits=c(-40,28),breaks=seq(-40,40,10)) + scale_y_continuous(limits=c(0,32),breaks=seq(0,40,8)) +
  annotate('text', x=-38, y = 30.5, label = names[16], hjust='left',size=3) +
  ylab("") + theme1 + xlab("")

# Prudhoe Bay
df <- read.csv(paste0(path,names[17],"_",ra,".csv"))
df$MERRA2Rainfall <- df$MERRA2Precip-df$MERRA2Snowfall
ddf <- data.frame(MERRA2TMax = seq(-40,40,2),MERRA2SnowCount=0,MERRA2MixedCount=0,MERRA2RainCount=0)
for (i in 1:nrow(ddf)){
  ddf[i,"MERRA2SnowCount"] <- nrow(df[df$MERRA2TMax > ddf[i,"MERRA2TMax"]-1 & df$MERRA2TMax <= ddf[i,"MERRA2TMax"]+1 & df$MERRA2Snowfall >= 0.254 & df$MERRA2Rainfall < 0.254,])
  ddf[i,"MERRA2MixedCount"] <- nrow(df[df$MERRA2TMax > ddf[i,"MERRA2TMax"]-1 & df$MERRA2TMax <= ddf[i,"MERRA2TMax"]+1 & df$MERRA2Snowfall >= 0.254 & df$MERRA2Rainfall >= 0.254,])
  ddf[i,"MERRA2RainCount"] <- nrow(df[df$MERRA2TMax > ddf[i,"MERRA2TMax"]-1 & df$MERRA2TMax <= ddf[i,"MERRA2TMax"]+1 & df$MERRA2Snowfall < 0.254 & df$MERRA2Rainfall >= 0.254,])
}
edf <- melt(ddf,"MERRA2TMax")
edf$value2 <- edf$value/nrow(df)*365.25 # convert to a #/year measure
p17 <- ggplot() + geom_bar(data=edf,aes(x=MERRA2TMax,y=value2,fill=variable),stat='identity',position='stack',show.legend=FALSE) + 
  scale_fill_manual(values=c("#0000FF", "#00DCDC","#027C02")) + 
  scale_x_continuous(limits=c(-40,28),breaks=seq(-40,40,10)) + scale_y_continuous(limits=c(0,32),breaks=seq(0,40,8)) +
  annotate('text', x=-38, y = 30.5, label = names[17], hjust='left',size=3) +
  ylab("") + theme1 + xlab("")

# Nome
df <- read.csv(paste0(path,names[12],"_",ra,".csv"))
df$MERRA2Rainfall <- df$MERRA2Precip-df$MERRA2Snowfall
ddf <- data.frame(MERRA2TMax = seq(-40,40,2),MERRA2SnowCount=0,MERRA2MixedCount=0,MERRA2RainCount=0)
for (i in 1:nrow(ddf)){
  ddf[i,"MERRA2SnowCount"] <- nrow(df[df$MERRA2TMax > ddf[i,"MERRA2TMax"]-1 & df$MERRA2TMax <= ddf[i,"MERRA2TMax"]+1 & df$MERRA2Snowfall >= 0.254 & df$MERRA2Rainfall < 0.254,])
  ddf[i,"MERRA2MixedCount"] <- nrow(df[df$MERRA2TMax > ddf[i,"MERRA2TMax"]-1 & df$MERRA2TMax <= ddf[i,"MERRA2TMax"]+1 & df$MERRA2Snowfall >= 0.254 & df$MERRA2Rainfall >= 0.254,])
  ddf[i,"MERRA2RainCount"] <- nrow(df[df$MERRA2TMax > ddf[i,"MERRA2TMax"]-1 & df$MERRA2TMax <= ddf[i,"MERRA2TMax"]+1 & df$MERRA2Snowfall < 0.254 & df$MERRA2Rainfall >= 0.254,])
}
edf <- melt(ddf,"MERRA2TMax")
edf$value2 <- edf$value/nrow(df)*365.25 # convert to a #/year measure
p12 <- ggplot() + geom_bar(data=edf,aes(x=MERRA2TMax,y=value2,fill=variable),stat='identity',position='stack',show.legend=FALSE) + 
  scale_fill_manual(values=c("#0000FF", "#00DCDC","#027C02")) + 
  scale_x_continuous(limits=c(-40,28),breaks=seq(-40,40,10)) + scale_y_continuous(limits=c(0,32),breaks=seq(0,40,8)) +
  annotate('text', x=-38, y = 30.5, label = names[12], hjust='left',size=3) +
  ylab("Precip Days/Year") + theme1 + xlab("")

# Fairbanks
df <- read.csv(paste0(path,names[15],"_",ra,".csv"))
df$MERRA2Rainfall <- df$MERRA2Precip-df$MERRA2Snowfall
ddf <- data.frame(MERRA2TMax = seq(-40,40,2),MERRA2SnowCount=0,MERRA2MixedCount=0,MERRA2RainCount=0)
for (i in 1:nrow(ddf)){
  ddf[i,"MERRA2SnowCount"] <- nrow(df[df$MERRA2TMax > ddf[i,"MERRA2TMax"]-1 & df$MERRA2TMax <= ddf[i,"MERRA2TMax"]+1 & df$MERRA2Snowfall >= 0.254 & df$MERRA2Rainfall < 0.254,])
  ddf[i,"MERRA2MixedCount"] <- nrow(df[df$MERRA2TMax > ddf[i,"MERRA2TMax"]-1 & df$MERRA2TMax <= ddf[i,"MERRA2TMax"]+1 & df$MERRA2Snowfall >= 0.254 & df$MERRA2Rainfall >= 0.254,])
  ddf[i,"MERRA2RainCount"] <- nrow(df[df$MERRA2TMax > ddf[i,"MERRA2TMax"]-1 & df$MERRA2TMax <= ddf[i,"MERRA2TMax"]+1 & df$MERRA2Snowfall < 0.254 & df$MERRA2Rainfall >= 0.254,])
}
edf <- melt(ddf,"MERRA2TMax")
edf$value2 <- edf$value/nrow(df)*365.25 # convert to a #/year measure
p15 <- ggplot() + geom_bar(data=edf,aes(x=MERRA2TMax,y=value2,fill=variable),stat='identity',position='stack',show.legend=FALSE) + 
  scale_fill_manual(values=c("#0000FF", "#00DCDC","#027C02")) + 
  scale_x_continuous(limits=c(-40,28),breaks=seq(-40,40,10)) + scale_y_continuous(limits=c(0,32),breaks=seq(0,40,8)) +
  annotate('text', x=-38, y = 30.5, label = names[15], hjust='left',size=3) +
  ylab("") + theme1 + xlab("")

# Fort Yukon
df <- read.csv(paste0(path,names[20],"_",ra,".csv"))
df$MERRA2Rainfall <- df$MERRA2Precip-df$MERRA2Snowfall
ddf <- data.frame(MERRA2TMax = seq(-40,40,2),MERRA2SnowCount=0,MERRA2MixedCount=0,MERRA2RainCount=0)
for (i in 1:nrow(ddf)){
  ddf[i,"MERRA2SnowCount"] <- nrow(df[df$MERRA2TMax > ddf[i,"MERRA2TMax"]-1 & df$MERRA2TMax <= ddf[i,"MERRA2TMax"]+1 & df$MERRA2Snowfall >= 0.254 & df$MERRA2Rainfall < 0.254,])
  ddf[i,"MERRA2MixedCount"] <- nrow(df[df$MERRA2TMax > ddf[i,"MERRA2TMax"]-1 & df$MERRA2TMax <= ddf[i,"MERRA2TMax"]+1 & df$MERRA2Snowfall >= 0.254 & df$MERRA2Rainfall >= 0.254,])
  ddf[i,"MERRA2RainCount"] <- nrow(df[df$MERRA2TMax > ddf[i,"MERRA2TMax"]-1 & df$MERRA2TMax <= ddf[i,"MERRA2TMax"]+1 & df$MERRA2Snowfall < 0.254 & df$MERRA2Rainfall >= 0.254,])
}
edf <- melt(ddf,"MERRA2TMax")
edf$value2 <- edf$value/nrow(df)*365.25 # convert to a #/year measure
p20 <- ggplot() + geom_bar(data=edf,aes(x=MERRA2TMax,y=value2,fill=variable),stat='identity',position='stack',show.legend=FALSE) + 
  scale_fill_manual(values=c("#0000FF", "#00DCDC","#027C02")) + 
  scale_x_continuous(limits=c(-40,28),breaks=seq(-40,40,10)) + scale_y_continuous(limits=c(0,32),breaks=seq(0,40,8)) +
  annotate('text', x=-38, y = 30.5, label = names[20], hjust='left',size=3) +
  ylab("") + theme1 + xlab("")

# Bethel
df <- read.csv(paste0(path,names[11],"_",ra,".csv"))
df$MERRA2Rainfall <- df$MERRA2Precip-df$MERRA2Snowfall
ddf <- data.frame(MERRA2TMax = seq(-40,40,2),MERRA2SnowCount=0,MERRA2MixedCount=0,MERRA2RainCount=0)
for (i in 1:nrow(ddf)){
  ddf[i,"MERRA2SnowCount"] <- nrow(df[df$MERRA2TMax > ddf[i,"MERRA2TMax"]-1 & df$MERRA2TMax <= ddf[i,"MERRA2TMax"]+1 & df$MERRA2Snowfall >= 0.254 & df$MERRA2Rainfall < 0.254,])
  ddf[i,"MERRA2MixedCount"] <- nrow(df[df$MERRA2TMax > ddf[i,"MERRA2TMax"]-1 & df$MERRA2TMax <= ddf[i,"MERRA2TMax"]+1 & df$MERRA2Snowfall >= 0.254 & df$MERRA2Rainfall >= 0.254,])
  ddf[i,"MERRA2RainCount"] <- nrow(df[df$MERRA2TMax > ddf[i,"MERRA2TMax"]-1 & df$MERRA2TMax <= ddf[i,"MERRA2TMax"]+1 & df$MERRA2Snowfall < 0.254 & df$MERRA2Rainfall >= 0.254,])
}
edf <- melt(ddf,"MERRA2TMax")
edf$value2 <- edf$value/nrow(df)*365.25 # convert to a #/year measure
p11 <- ggplot() + geom_bar(data=edf,aes(x=MERRA2TMax,y=value2,fill=variable),stat='identity',position='stack',show.legend=FALSE) + 
  scale_fill_manual(values=c("#0000FF", "#00DCDC","#027C02")) + 
  scale_x_continuous(limits=c(-40,28),breaks=seq(-40,40,10)) + scale_y_continuous(limits=c(0,32),breaks=seq(0,40,8)) +
  annotate('text', x=-38, y = 30.5, label = names[11], hjust='left',size=3) +
  ylab("Precip Days/Year") + theme1 + xlab("Max Daily Temperature (°C)")

# Galena
df <- read.csv(paste0(path,names[21],"_",ra,".csv"))
df$MERRA2Rainfall <- df$MERRA2Precip-df$MERRA2Snowfall
ddf <- data.frame(MERRA2TMax = seq(-40,40,2),MERRA2SnowCount=0,MERRA2MixedCount=0,MERRA2RainCount=0)
for (i in 1:nrow(ddf)){
  ddf[i,"MERRA2SnowCount"] <- nrow(df[df$MERRA2TMax > ddf[i,"MERRA2TMax"]-1 & df$MERRA2TMax <= ddf[i,"MERRA2TMax"]+1 & df$MERRA2Snowfall >= 0.254 & df$MERRA2Rainfall < 0.254,])
  ddf[i,"MERRA2MixedCount"] <- nrow(df[df$MERRA2TMax > ddf[i,"MERRA2TMax"]-1 & df$MERRA2TMax <= ddf[i,"MERRA2TMax"]+1 & df$MERRA2Snowfall >= 0.254 & df$MERRA2Rainfall >= 0.254,])
  ddf[i,"MERRA2RainCount"] <- nrow(df[df$MERRA2TMax > ddf[i,"MERRA2TMax"]-1 & df$MERRA2TMax <= ddf[i,"MERRA2TMax"]+1 & df$MERRA2Snowfall < 0.254 & df$MERRA2Rainfall >= 0.254,])
}
edf <- melt(ddf,"MERRA2TMax")
edf$value2 <- edf$value/nrow(df)*365.25 # convert to a #/year measure
p21 <- ggplot() + geom_bar(data=edf,aes(x=MERRA2TMax,y=value2,fill=variable),stat='identity',position='stack',show.legend=FALSE) + 
  scale_fill_manual(values=c("#0000FF", "#00DCDC","#027C02")) + 
  scale_x_continuous(limits=c(-40,28),breaks=seq(-40,40,10)) + scale_y_continuous(limits=c(0,32),breaks=seq(0,40,8)) +
  annotate('text', x=-38, y = 30.5, label = names[21], hjust='left',size=3) +
  ylab("") + theme1 + xlab("Max Daily Temperature (°C)")

# Anchorage
df <- read.csv(paste0(path,names[13],"_",ra,".csv"))
df$MERRA2Rainfall <- df$MERRA2Precip-df$MERRA2Snowfall
ddf <- data.frame(MERRA2TMax = seq(-40,40,2),MERRA2SnowCount=0,MERRA2MixedCount=0,MERRA2RainCount=0)
for (i in 1:nrow(ddf)){
  ddf[i,"MERRA2SnowCount"] <- nrow(df[df$MERRA2TMax > ddf[i,"MERRA2TMax"]-1 & df$MERRA2TMax <= ddf[i,"MERRA2TMax"]+1 & df$MERRA2Snowfall >= 0.254 & df$MERRA2Rainfall < 0.254,])
  ddf[i,"MERRA2MixedCount"] <- nrow(df[df$MERRA2TMax > ddf[i,"MERRA2TMax"]-1 & df$MERRA2TMax <= ddf[i,"MERRA2TMax"]+1 & df$MERRA2Snowfall >= 0.254 & df$MERRA2Rainfall >= 0.254,])
  ddf[i,"MERRA2RainCount"] <- nrow(df[df$MERRA2TMax > ddf[i,"MERRA2TMax"]-1 & df$MERRA2TMax <= ddf[i,"MERRA2TMax"]+1 & df$MERRA2Snowfall < 0.254 & df$MERRA2Rainfall >= 0.254,])
}
edf <- melt(ddf,"MERRA2TMax")
edf$value2 <- edf$value/nrow(df)*365.25 # convert to a #/year measure
p13 <- ggplot() + geom_bar(data=edf,aes(x=MERRA2TMax,y=value2,fill=variable),stat='identity',position='stack',show.legend=FALSE) + 
  scale_fill_manual(values=c("#0000FF", "#00DCDC","#027C02")) + 
  scale_x_continuous(limits=c(-40,28),breaks=seq(-40,40,10)) + scale_y_continuous(limits=c(0,32),breaks=seq(0,40,8)) +
  annotate('text', x=-38, y = 30.5, label = names[13], hjust='left',size=3) +
  ylab("") + theme1 + xlab("Max Daily Temperature (°C)")

ptotal <- grid.arrange(p19,p16,p17,p12,p15,p20,p11,p21,p13,nrow=3)
ggsave("Histogram_PrecipDaysvTMax_Panel9_M2.png",ptotal,'png',figpath,dpi=300,width=6.5,height=5.5,units="in")

##### SNOW COVER ######
# ANNUAL NUMBER OF SNOW COVER DAYS #
sc <- data.frame(Station=names,SnowNM2=0,SnowM2N=0,SnowN=0,SnowM2=0)
for (v in vmin:vmax){
  # Read in Data
  df <- read.csv(paste0(path,names[v],"_",ra,".csv"))
  df$SNOWN <- ifelse((df$SNOW > 0) | (df$SNWD > 0), 1, 0)
  df$SNOWM2 <- ifelse((df$MERRA2SnowDepth > 25.4), 1, 0)
  
  df <- df[is.finite(df$SNOWN) == 1,]
  
  sc[sc$Station == names[v],]$SnowN <- (sum(df$SNOWN) - sum(df[df$SNOWN == 1,]$SNOWM2))/nrow(df)*365.25
  sc[sc$Station == names[v],]$SnowM2 <- (sum(df$SNOWM2) - sum(df[df$SNOWM2 == 1,]$SNOWN))/nrow(df)*365.25
  sc[sc$Station == names[v],]$SnowNM2 <- (sum(df[df$SNOWN == 1,]$SNOWM2))/nrow(df)*365.25
  sc[sc$Station == names[v],]$SnowM2N <- (sum(df[df$SNOWM2 == 1,]$SNOWN))/nrow(df)*365.25
}

library(ggplot2)
sc2 <- melt(sc,"Station")
sc2$cat1 <- ifelse(sc2$variable == "SnowN" | sc2$variable == "SnowNM2", 0, 1)
sc2$cat2 <- ifelse(sc2$variable == "SnowM2N" | sc2$variable == "SnowNM2", 0, 1)

ggplot() +
  geom_bar(data=sc2[sc2$value > 0, ], aes(y = value, x = cat1, fill = cat2), stat="identity", position='stack') +
  theme_bw() + facet_grid( ~ Station) + ylab("Days Per Year") + xlab('') + scale_x_continuous(breaks=c(0,1),labels=c("GHCN","M2")) + 
  theme(legend.position="none")

# PERCENT SNOW COVER #
sc <- data.frame(Station=names,SnowNM2=0,SnowM2N=0,SnowN=0,SnowM2=0,Order=0)
for (v in c(11:21)){
  # Read in Data
  df <- read.csv(paste0(path,names[v],"_",ra,".csv"))
  df$SNOWN <- ifelse((df$SNOW > 0) | (df$SNWD > 0), 1, 0)
  df$SNOWM2 <- ifelse((df$MERRA2SnowDepth > 25.4), 1, 0)
  
  df <- df[is.finite(df$SNOWN) == 1,]
  
  SNOWN = sum(df$SNOWN)
  SNOWM2 = sum(df$SNOWM2)
  
  sc[sc$Station == names[v],]$SnowN <- (SNOWN - sum(df[df$SNOWN == 1,]$SNOWM2))/SNOWN
  sc[sc$Station == names[v],]$SnowM2 <- (SNOWM2 - sum(df[df$SNOWM2 == 1,]$SNOWN))/SNOWM2
  sc[sc$Station == names[v],]$SnowNM2 <- (sum(df[df$SNOWN == 1,]$SNOWM2))/SNOWN
  sc[sc$Station == names[v],]$SnowM2N <- (sum(df[df$SNOWM2 == 1,]$SNOWN))/SNOWM2
}

sc2 <- melt(sc,"Station")
sc2$cat1 <- ifelse(sc2$variable == "SnowN" | sc2$variable == "SnowNM2", 0, 1)
sc2$cat2 <- ifelse(sc2$variable == "SnowM2N" | sc2$variable == "SnowNM2", 0, 1)

# VERSION 1
newnames <- c("1"=names[19],"2"=names[16],"3"=names[17],"4"=names[12],"5"=names[15])
sc2$Order <- ifelse(sc2$Station == names[19],"1",ifelse(sc2$Station == names[16],"2",ifelse(sc2$Station == names[17],"3",ifelse(sc2$Station == names[12],"4",ifelse(sc2$Station == names[15],"5","0")))))

ggplot() +
  geom_bar(data=sc2[sc2$value > 0,], aes(y = value, x = cat1, fill = cat2), stat="identity", position='stack') +
  theme_bw() + facet_grid( ~ Order,labeller=as_labeller(newnames)) + ylab("Fraction") + xlab('') + scale_x_continuous(breaks=c(0,1),labels=c("GHCN","M2")) + 
  theme(legend.position="none")

# VERSION 2
newnames <- c("1"=names[20],"2"=names[11],"3"=names[21],"4"=names[13],"5"=names[14])
sc2$Order <- ifelse(sc2$Station == names[20],"1",ifelse(sc2$Station == names[11],"2",ifelse(sc2$Station == names[21],"3",ifelse(sc2$Station == names[13],"4",ifelse(sc2$Station == names[14],"5","0")))))

ggplot() +
  geom_bar(data=sc2[sc2$Order != "0",], aes(y = value, x = cat1, fill = cat2), stat="identity", position='stack') +
  theme_bw() + facet_grid( ~ Order,labeller=as_labeller(newnames)) + ylab("Fraction") + xlab('') + scale_x_continuous(breaks=c(0,1),labels=c("GHCN","M2")) + 
  theme(legend.position="none")