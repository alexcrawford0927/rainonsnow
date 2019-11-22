#########
# Author: Alex Crawford
# Date Created: 20 May 2019
# Date Modified: 17 Jun 2019 --> Added SOS Events and removed some of the ROSF Counts
# Purpose: Create a table of counts for precip & ROS events at various locations
#########

### DECLARE VARIABLES #####
path = "/Volumes/Miranda/RainOnSnow/PrecipDetection_ByPoint_withStorms/V2/"
csvpath = "/Volumes/Miranda/RainOnSnow/PrecipDetection_Aggregation/V5/"
figpath = "/Volumes/Miranda/RainOnSnow/Figures/"

time = "19800101_20190101"
hthresh = 2 # in hr
rthresh = 6.096/24. # in mm/hr
pthresh = 0.254 # in mm
ppthresh = 2.54 # 0.254 # in mm 
sthresh = 0.0254 # in m
nyrs = 39

vmin = 11
vmax = 21
ra = "MERRA2"
names = c("Aniak","Pargon Creek","Kelly Station","Fort Yukon","Teuchet Creek",
          "Indian Pass","Monument Creek","Cooper Lake","Kenai Moose Pens","Mt Ryan",
          "Bethel","Nome","Anchorage","Juneau","Fairbanks","Utqiagvik",
          "Prudhoe Bay","Coleville Village","Kotzebue","Fort Yukon","Galena")
ys = c(122,129,135,132,129,121,129,120,120,NA,
       121,128,121,116,129,142,140,140,133,132,128) # Y of station
xs = c(33,27,28,56,55,49,55,48,47,NA,
       29,23,48,73,51,37,50,47,28,56,37) # X of station

### CREATE COUNTS ########
count <- data.frame(Location=names[vmin:vmax],PrecipCT=0, RainCT=0, ROSCT=0, ROSf01CT=0, ROSf10CT=0, ROSf30CT=0, ROSf90CT=0, SOSCT=0, PrecipTotal=0, RainTotal=0, AvgLength=0)
storm <- data.frame(Location=names[vmin:vmax],PrecipCT=0, RainCT=0, ROSCT=0, ROSf01CT=0, ROSf10CT=0, ROSf30CT=0, ROSf90CT=0, SOSCT=0, PrecipTotal=0, RainTotal=0, AvgLength=0)

for (v in vmin:vmax){
  # Read in Data
  df <- read.csv(paste0(path,names[v],"_",time,"_Gap",hthresh,"_Rate",rthresh,"_Total",pthresh,"_X",xs[v],"_Y",ys[v],".csv"))
  df <- df[df$Precip > ppthresh,]
  ros <- df[df$Ice == 3 & df$Rain > ppthresh,]
  
  
  # Make Basic Counts
  count[count$Location == names[v],"PrecipCT"] <- nrow(df)/nyrs
  count[count$Location == names[v],"RainCT"] <- sum(df$Rain > ppthresh)/nyrs
  count[count$Location == names[v],"ROSCT"] <- nrow(ros)/nyrs
  count[count$Location == names[v],"ROSf01CT"] <- sum(ros$Frz01 >= 0.25)/nyrs
  count[count$Location == names[v],"ROSf10CT"] <- sum(ros$Frz10 >= 0.95)/nyrs
  count[count$Location == names[v],"ROSf30CT"] <- sum(ros$Frz30 >= 0.95)/nyrs
  count[count$Location == names[v],"ROSf90CT"] <- sum(ros$Frz90 >= 0.95)/nyrs
  count[count$Location == names[v],"SOSCT"] <- nrow(df[df$EdSnoDep > sthresh & df$StSnoDep > sthresh & df$Rain < pthresh & df$Precip > ppthresh,])/nyrs

  count[count$Location == names[v],"PrecipTotal"] <- sum(df$Precip)/nyrs
  count[count$Location == names[v],"RainTotal"] <- sum(df$Rain)/nyrs
  count[count$Location == names[v],"AvgLength"] <- mean(df$Length)
  
  # Make Cyc-Associated Counts
  stdf <- df[is.finite(df$SID) == 1,]
  stros <- ros[is.finite(ros$SID) == 1,]
  storm[storm$Location == names[v],"PrecipCT"] <- nrow(stdf)/nyrs
  storm[storm$Location == names[v],"RainCT"] <- sum(stdf$Rain > ppthresh)/nyrs
  storm[storm$Location == names[v],"ROSCT"] <- nrow(stros)/nyrs
  storm[storm$Location == names[v],"ROSf01CT"] <- sum(stros$Frz01 >= 0.25)/nyrs
  storm[storm$Location == names[v],"ROSf10CT"] <- sum(stros$Frz10 >= 0.95)/nyrs
  storm[storm$Location == names[v],"ROSf30CT"] <- sum(stros$Frz30 >= 0.95)/nyrs
  storm[storm$Location == names[v],"ROSf90CT"] <- sum(stros$Frz90 >= 0.95)/nyrs
  storm[storm$Location == names[v],"SOSCT"] <- nrow(stdf[stdf$EdSnoDep > sthresh & stdf$StSnoDep > sthresh & stdf$Rain < pthresh & stdf$Precip > ppthresh,])/nyrs
  
  storm[storm$Location == names[v],"PrecipTotal"] <- sum(stdf$Precip)/nyrs
  storm[storm$Location == names[v],"RainTotal"] <- sum(stdf$Rain)/nyrs
  storm[storm$Location == names[v],"AvgLength"] <- mean(stdf$Length)
}

write.csv(count,paste0(csvpath,"EventCounts_",time,"_Gap",hthresh,"_Rate",rthresh,"_Total",ppthresh,".csv"),row.names = F)
write.csv(storm,paste0(csvpath,"CycEventCounts_",time,"_Gap",hthresh,"_Rate",rthresh,"_Total",ppthresh,".csv"),row.names = F)

### CREATE MONTHLY COUNTS ########
count <- data.frame(Location=NA,Month=rep(1:12,length(vmin:vmax)),PrecipCT=0, RainCT=0, ROSCT=0, ROSf01CT=0, ROSf10CT=0, ROSf30CT=0, ROSf90CT=0, SOSCT=0, PrecipTotal=0, RainTotal=0, AvgLength=0)
storm <- data.frame(Location=NA,Month=rep(1:12,length(vmin:vmax)),PrecipCT=0, RainCT=0, ROSCT=0, ROSf01CT=0, ROSf10CT=0, ROSf30CT=0, ROSf90CT=0, SOSCT=0, PrecipTotal=0, RainTotal=0, AvgLength=0)

row = 0
for (v in vmin:vmax){
  # Read in Data
  df1 <- read.csv(paste0(path,names[v],"_",time,"_Gap",hthresh,"_Rate",rthresh,"_Total",pthresh,"_X",xs[v],"_Y",ys[v],".csv"))
  df1 <- df1[df1$Precip > ppthresh,]

  # Subset for ROS events & Cyc events
  ros1 <- df1[df1$Ice == 3 & df1$Rain > ppthresh,]
  stdf1 <- df1[is.finite(df1$SID) == 1,]
  stros1 <- ros1[is.finite(ros1$SID) == 1,]
  
  # For Each Month...
  for (m in 1:12){
    row = row + 1
    count[row,"Location"] <- names[v]
    storm[row,"Location"] <- names[v]
    
    # Subset by month
    df <- df1[df1$Month == m,]
    ros <- ros1[ros1$Month == m,]
    stdf <- stdf1[stdf1$Month == m,]
    stros <- stros1[stros1$Month == m,]
    
    # Make counts
    count[row,"PrecipCT"] <- nrow(df)/nyrs
    count[row,"RainCT"] <- sum(df$Rain > ppthresh)/nyrs
    count[row,"ROSCT"] <- nrow(ros)/nyrs
    count[row,"ROSf01CT"] <- sum(ros$Frz01 >= 0.25)/nyrs
    count[row,"ROSf10CT"] <- sum(ros$Frz10 >= 0.95)/nyrs
    count[row,"ROSf30CT"] <- sum(ros$Frz30 >= 0.95)/nyrs
    count[row,"ROSf90CT"] <- sum(ros$Frz90 >= 0.95)/nyrs
    count[row,"SOSCT"] <- nrow(df[df$EdSnoDep > sthresh & df$StSnoDep > sthresh & df$Rain < pthresh & df$Precip > ppthresh,])/nyrs
    
    count[row,"PrecipTotal"] <- sum(df$Precip)/nyrs
    count[row,"RainTotal"] <- sum(df$Rain)/nyrs
    count[row,"AvgLength"] <- mean(df$Length)
    
    # Make counts for Cyc events
    storm[row,"PrecipCT"] <- nrow(stdf)/nyrs
    storm[row,"RainCT"] <- sum(stdf$Rain > ppthresh)/nyrs
    storm[row,"ROSCT"] <- nrow(stros)/nyrs
    storm[row,"ROSf01CT"] <- sum(stros$Frz01 >= 0.25)/nyrs
    storm[row,"ROSf10CT"] <- sum(stros$Frz10 >= 0.95)/nyrs
    storm[row,"ROSf30CT"] <- sum(stros$Frz30 >= 0.95)/nyrs
    storm[row,"ROSf90CT"] <- sum(stros$Frz90 >= 0.95)/nyrs
    storm[row,"SOSCT"] <- nrow(stdf[stdf$EdSnoDep > sthresh & stdf$StSnoDep > sthresh & stdf$Rain < pthresh & stdf$Precip > ppthresh,])/nyrs
    
    storm[row,"PrecipTotal"] <- sum(stdf$Precip)/nyrs
    storm[row,"RainTotal"] <- sum(stdf$Rain)/nyrs
    storm[row,"AvgLength"] <- mean(stdf$Length)
  }}

write.csv(count,paste0(csvpath,"MonthlyEventCounts_",time,"_Gap",hthresh,"_Rate",rthresh,"_Total",ppthresh,".csv"),row.names = F)
write.csv(storm,paste0(csvpath,"MonthlyCycEventCounts_",time,"_Gap",hthresh,"_Rate",rthresh,"_Total",ppthresh,".csv"),row.names = F)

############################################
# OLD VERSION #############
# ### CREATE COUNTS ########
# count <- data.frame(Location=names[vmin:vmax],PrecipCT=0, RainCT=0, ROSCT=0, ROSh120CT=0, ROSh240CT=0, ROSf01CT=0, ROSf10CT=0, ROSf30CT=0, ROSf90CT=0, ROSf01f30CT=0, ROSh240f30CT=0, ROSh240f90CT=0, PrecipTotal=0, RainTotal=0, AvgLength=0)
# storm <- data.frame(Location=names[vmin:vmax],PrecipCT=0, RainCT=0, ROSCT=0, ROSh120CT=0, ROSh240CT=0, ROSf01CT=0, ROSf10CT=0, ROSf30CT=0, ROSf90CT=0, ROSf01f30CT=0, ROSh240f30CT=0, ROSh240f90CT=0, PrecipTotal=0, RainTotal=0, AvgLength=0)
# 
# for (v in vmin:vmax){
#   # Read in Data
#   df <- read.csv(paste0(path,names[v],"_",time,"_Gap",hthresh,"_Rate",rthresh,"_Total",pthresh,"_X",xs[v],"_Y",ys[v],".csv"))
#   df <- df[df$Precip > ppthresh,]
#   ros <- df[df$Ice == 3 & df$Rain > ppthresh,]
#   
#   # Make Basic Counts
#   count[count$Location == names[v],"PrecipCT"] <- nrow(df)
#   count[count$Location == names[v],"RainCT"] <- sum(df$Rain > ppthresh)
#   count[count$Location == names[v],"ROSCT"] <- nrow(ros)
#   count[count$Location == names[v],"ROSh120CT"] <- sum(ros$FrzHrs >= 120)
#   count[count$Location == names[v],"ROSh240CT"] <- sum(ros$FrzHrs >= 240)
#   count[count$Location == names[v],"ROSf01CT"] <- sum(ros$Frz01 >= 0.25)
#   count[count$Location == names[v],"ROSf10CT"] <- sum(ros$Frz10 >= 0.95)
#   count[count$Location == names[v],"ROSf30CT"] <- sum(ros$Frz30 >= 0.95)
#   count[count$Location == names[v],"ROSf90CT"] <- sum(ros$Frz90 >= 0.95)
#   count[count$Location == names[v],"ROSf01f30CT"] <- sum((ros$Frz30 >= 0.95) & (ros$Frz01 >= 0.25))
#   count[count$Location == names[v],"ROSh240f30CT"] <- sum((ros$Frz30 >= 0.95) & (ros$FrzHrs >= 240))
#   count[count$Location == names[v],"ROSh240f90CT"] <- sum((ros$Frz90 >= 0.95) & (ros$FrzHrs >= 240))
#   
#   
#   count[count$Location == names[v],"PrecipTotal"] <- sum(df$Precip)
#   count[count$Location == names[v],"RainTotal"] <- sum(df$Rain)
#   count[count$Location == names[v],"AvgLength"] <- mean(df$Length)
#   
#   # Make Cyc-Associated Counts
#   stdf <- df[is.finite(df$SID) == 1,]
#   stros <- ros[is.finite(ros$SID) == 1,]
#   storm[storm$Location == names[v],"PrecipCT"] <- nrow(stdf)
#   storm[storm$Location == names[v],"RainCT"] <- sum(stdf$Rain > ppthresh)
#   storm[storm$Location == names[v],"ROSCT"] <- nrow(stros)
#   storm[storm$Location == names[v],"ROSh120CT"] <- sum(stros$FrzHrs >= 120)
#   storm[storm$Location == names[v],"ROSh240CT"] <- sum(stros$FrzHrs >= 240)
#   storm[storm$Location == names[v],"ROSf01CT"] <- sum(stros$Frz01 >= 0.25)
#   storm[storm$Location == names[v],"ROSf10CT"] <- sum(stros$Frz10 >= 0.95)
#   storm[storm$Location == names[v],"ROSf30CT"] <- sum(stros$Frz30 >= 0.95)
#   storm[storm$Location == names[v],"ROSf90CT"] <- sum(stros$Frz90 >= 0.95)
#   storm[storm$Location == names[v],"ROSf01f30CT"] <- sum((stros$Frz30 >= 0.95) & (stros$Frz01 >= 0.25))
#   storm[storm$Location == names[v],"ROSh240f30CT"] <- sum((stros$Frz30 >= 0.95) & (stros$FrzHrs >= 240))
#   storm[storm$Location == names[v],"ROSh240f90CT"] <- sum((stros$Frz90 >= 0.95) & (stros$FrzHrs >= 240))
#   
#   storm[storm$Location == names[v],"PrecipTotal"] <- sum(stdf$Precip)
#   storm[storm$Location == names[v],"RainTotal"] <- sum(stdf$Rain)
#   storm[storm$Location == names[v],"AvgLength"] <- mean(stdf$Length)
# }
# 
# write.csv(count,paste0(csvpath,"EventCounts_",time,"_Gap",hthresh,"_Rate",rthresh,"_Total",ppthresh,".csv"),row.names = F)
# write.csv(storm,paste0(csvpath,"CycEventCounts_",time,"_Gap",hthresh,"_Rate",rthresh,"_Total",ppthresh,".csv"),row.names = F)
# 
# ### CREATE MONTHLY COUNTS ########
# count <- data.frame(Location=NA,Month=rep(1:12,length(vmin:vmax)),PrecipCT=0, RainCT=0, ROSCT=0, ROSh120CT=0, ROSh240CT=0, ROSf01CT=0, ROSf10CT=0, ROSf30CT=0, ROSf90CT=0, ROSf01f30CT=0, ROSh240f30CT=0, ROSh240f90CT=0, PrecipTotal=0, RainTotal=0, AvgLength=0)
# storm <- data.frame(Location=NA,Month=rep(1:12,length(vmin:vmax)),PrecipCT=0, RainCT=0, ROSCT=0, ROSh120CT=0, ROSh240CT=0, ROSf01CT=0, ROSf10CT=0, ROSf30CT=0, ROSf90CT=0, ROSf01f30CT=0, ROSh240f30CT=0, ROSh240f90CT=0, PrecipTotal=0, RainTotal=0, AvgLength=0)
# 
# row = 0
# for (v in vmin:vmax){
#   # Read in Data
#   df1 <- read.csv(paste0(path,names[v],"_",time,"_Gap",hthresh,"_Rate",rthresh,"_Total",pthresh,"_X",xs[v],"_Y",ys[v],".csv"))
#   df1 <- df1[df1$Precip > ppthresh,]
#   
#   # Subset for ROS events & Cyc events
#   ros1 <- df1[df1$Ice == 3 & df1$Rain > ppthresh,]
#   stdf1 <- df1[is.finite(df1$SID) == 1,]
#   stros1 <- ros1[is.finite(ros1$SID) == 1,]
#   
#   # For Each Month...
#   for (m in 1:12){
#     row = row + 1
#     count[row,"Location"] <- names[v]
#     storm[row,"Location"] <- names[v]
#     
#     # Subset by month
#     df <- df1[df1$Month == m,]
#     ros <- ros1[ros1$Month == m,]
#     stdf <- stdf1[stdf1$Month == m,]
#     stros <- stros1[stros1$Month == m,]
#   
#     # Make counts
#     count[row,"PrecipCT"] <- nrow(df)
#     count[row,"RainCT"] <- sum(df$Rain > ppthresh)
#     count[row,"ROSCT"] <- nrow(ros)
#     count[row,"ROSh120CT"] <- sum(ros$FrzHrs >= 120)
#     count[row,"ROSh240CT"] <- sum(ros$FrzHrs >= 240)
#     count[row,"ROSf01CT"] <- sum(ros$Frz01 >= 0.25)
#     count[row,"ROSf10CT"] <- sum(ros$Frz10 >= 0.95)
#     count[row,"ROSf30CT"] <- sum(ros$Frz30 >= 0.95)
#     count[row,"ROSf90CT"] <- sum(ros$Frz90 >= 0.95)
#     count[row,"ROSf01f30CT"] <- sum((ros$Frz30 >= 0.95) & (ros$Frz01 >= 0.25))
#     count[row,"ROSh240f30CT"] <- sum((ros$Frz30 >= 0.95) & (ros$FrzHrs >= 240))
#     count[row,"ROSh240f90CT"] <- sum((ros$Frz90 >= 0.95) & (ros$FrzHrs >= 240))
#     
#     count[row,"PrecipTotal"] <- sum(df$Precip)
#     count[row,"RainTotal"] <- sum(df$Rain)
#     count[row,"AvgLength"] <- mean(df$Length)
#     
#     # Make counts for Cyc events
#     storm[row,"PrecipCT"] <- nrow(stdf)
#     storm[row,"RainCT"] <- sum(stdf$Rain > ppthresh)
#     storm[row,"ROSCT"] <- nrow(stros)
#     storm[row,"ROSh120CT"] <- sum(stros$FrzHrs >= 120)
#     storm[row,"ROSh240CT"] <- sum(stros$FrzHrs >= 240)
#     storm[row,"ROSf01CT"] <- sum(stros$Frz01 >= 0.25)
#     storm[row,"ROSf10CT"] <- sum(stros$Frz10 >= 0.95)
#     storm[row,"ROSf30CT"] <- sum(stros$Frz30 >= 0.95)
#     storm[row,"ROSf90CT"] <- sum(stros$Frz90 >= 0.95)
#     storm[row,"ROShf01f30CT"] <- sum((stros$Frz30 >= 0.95) & (stros$Frz01 >= 0.25))
#     storm[row,"ROSh240f30CT"] <- sum((stros$Frz30 >= 0.95) & (stros$FrzHrs >= 240))
#     storm[row,"ROSh240f90CT"] <- sum((stros$Frz90 >= 0.95) & (stros$FrzHrs >= 240))
#     
#     storm[row,"PrecipTotal"] <- sum(stdf$Precip)
#     storm[row,"RainTotal"] <- sum(stdf$Rain)
#     storm[row,"AvgLength"] <- mean(stdf$Length)
# }}
# 
# write.csv(count,paste0(csvpath,"MonthlyEventCounts_",time,"_Gap",hthresh,"_Rate",rthresh,"_Total",ppthresh,".csv"),row.names = F)
# write.csv(storm,paste0(csvpath,"MonthlyCycEventCounts_",time,"_Gap",hthresh,"_Rate",rthresh,"_Total",ppthresh,".csv"),row.names = F)
