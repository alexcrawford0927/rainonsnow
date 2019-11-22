#########
# Author: Alex Crawford
# Date Created: 2 July 2019
# Date Modified: 
#########

path = "/Volumes/Miranda/RainOnSnow/PrecipDetection_Aggregation/"
V = "V5"
  
time = "19800101_20190101"
minls = 1
mintl = 100

cs = c(1,2,3,4,5,6,7,8,10,14)

vs = c(4,11,12,13,15,16,17,19,20)
names = c("Aniak","Pargon Creek","Kelly Station","Fort Yukon","Teuchet Creek",
          "Indian Pass","Monument Creek","Cooper Lake","Kenai Moose Pens","Mt Ryan",
          "Bethel","Nome","Anchorage","Juneau","Fairbanks","Utqiagvik",
          "Prudhoe Bay","Coleville Village","Kotzebue","Galena")
ys = c(122,129,135,132,129,121,129,120,120,NA,
       121,128,121,116,129,142,140,140,133,128) # Y of station
xs = c(33,27,28,56,55,49,55,48,47,NA,
       29,23,48,73,51,37,50,47,28,37) # X of station

tdf = data.frame("Location"=NA,"Variable"=NA,"MedianROS"=NA,"MedianSOS"=NA,"Diff"=NA,"pvalue"=NA)

i = 0
for (v in vs){
  # Read in Data
  ros = read.csv( paste0(path,names[v],"_",V,"/ROS_AggregatedStats",time,".csv") )
  sos = read.csv( paste0(path,names[v],"_",V,"/SOS_AggregatedStats",time,".csv") )
  
  ros <- ros[ros$lifespan >= minls & ros$trlen >= mintl,]
  sos <- sos[sos$lifespan >= minls & sos$trlen >= mintl,]
  
  ros$genLon <- ifelse(ros$genLon < 0, ros$genLon+360, ros$genLon)
  sos$genLon <- ifelse(sos$genLon < 0, sos$genLon+360, sos$genLon)
  
  cnames <- colnames(ros)
  
  for (c in cs){
    i = i+1
    tdf[i,"Location"] <- names[v]
    tdf[i,"Variable"] <- cnames[c]
    tdf[i,"MedianROS"] <- median(ros[,c])
    tdf[i,"MedianSOS"] <- median(sos[,c])
    tdf[i,"Diff"] <- tdf[i,"MedianROS"] - tdf[i,"MedianSOS"]
    tdf[i,"pvalue"] <- wilcox.test(ros[,c],sos[,c])$p.value
  }
}

write.csv(tdf,paste0(path,V,"/CycloneStats_WilcoxTests_",time,"_minls",minls,"_mintl",mintl,".csv"),row.names = FALSE)
