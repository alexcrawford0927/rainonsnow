'''
Author: Alex Crawford
Date Created: 19 Mar 2019
Date Modified: 29 May 2019 --> Update for Python 3
Purpose: Aggregate storm characteristics (e.g., genesis and track location) for
all storms that relate to precipitation at a given location.
'''

'''********************
Import Modules
********************'''
import os
import numpy as np
import pandas as pd
from osgeo import gdal, gdalnumeric, gdalconst
import CycloneModule_11_1 as md

'''*******************************************
Define Variables
*******************************************'''
# Path Variables
csvpath = "/Volumes/Miranda/RainOnSnow/PrecipDetection_ByPoint_withStorms/V2B"
outpath = "/Volumes/Miranda/RainOnSnow/PrecipDetection_Aggregation"
ctpath = "/Volumes/Ferdinand/ArcticCyclone/detection11_3AM2/SystemTracks"
suppath = "/Volumes/Ferdinand/Projections"
latN = "EASE2_N0_100km_Lats.tif"
lonN = "EASE2_N0_100km_Lons.tif"

### Location Variables ###
v = 14 # (0-4, 5-9, 10-15, 16-19)
names = ["Aniak","Pargon Creek","Kelly Station","Fort Yukon","Teuchet Creek",\
         "Indian Pass","Monument Creek","Cooper Lake","Kenai Moose Pens","Mt Ryan",\
         "Bethel","Nome","Anchorage","Juneau","Fairbanks","Utqiagvik",\
         "Prudhoe Bay","Coleville Village","Kotzebue","Galena"]
rows = [122, 129, 135, 132, 129, 121, 129, 120, 120, 129,\
        121, 128, 121, 116, 129, 142, 140, 140, 133, 128]
cols = [33, 27, 28, 56, 55, 49, 55, 48, 47, 54,\
        29, 23, 48, 73, 51, 37, 50, 47, 28, 37]

### Physical Variables ###
rthresh = 6.096/24 # in mm/hr (equiva. to a rate of 0.01 in/day)
pthresh = 0.254 # total mm (equiv. to 0.1 in/event)
minrain = 2.54 # minimum total mm of rain needed to count as ROS
vmons = [10,11,12,1,2] # valid months to compare (should be based on which
## months actually contain events of interest)
vice = [3] # valid icing conditions
V = "V7"

### Time Variables ###
starttime = [1980,1,1,0,0,0]
endtime = [2019,1,1,0,0,0]
reftime = [1900,1,1,0,0,0]
hrs = 1 # Temporal Resolution in hours of the ROS
t = 3 # Temporal Resolution in hours of the cyclone tracking
hthresh = 2 # Number of hours w/out measurable precip needed to mark separate
# precip events

months = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
mons = ["01","02","03","04","05","06","07","08","09","10","11","12"]
days = ["01","02","03","04","05","06","07","08","09","10","11","12","13",\
    "14","15","16","17","18","19","20","21","22","23","24","25","26","27",\
    "28","29","30","31"]
hours = ["0000","0100","0200","0300","0400","0500","0600","0700","0800",\
    "0900","1000","1100","1200","1300","1400","1500","1600","1700","1800",\
    "1900","2000","2100","2200","2300"]

# Output Variables
dtype = gdal.GDT_Float64
kSize = 2

vNames = ["gen","lys","spl","mrg",\
"maxuv","maxdpdt","maxdep","minp","maxdsqp",\
"trkden","countC","countP","countU","countR","mcc",\
"uvAb","uvDi","radius","area",\
"depth","dpdr","dpdt","pcent","dsqp"]
varsi = range(0,len(vNames)) # range(0,1) #

'''*******************************************
Main Analysis
*******************************************'''

########## READ IN ROS DATA ############
# Load data frame
name = names[v]
row = rows[v]
col = cols[v]

SDAY = str(starttime[0]) + mons[starttime[1]-1] + days[starttime[2]-1]
EDAY = str(endtime[0]) + mons[endtime[1]-1] + days[endtime[2]-1]
pdf = pd.read_csv(csvpath+"/"+name+"_"+SDAY+"_"+EDAY+"_Gap"+str(hthresh)+"_Rate"+str(np.round(rthresh,4))+"_Total"+str(pthresh)+"_X"+str(col)+"_Y"+str(row)+".csv")

# Limit the events based on ROS criteria
qdf = pdf[(pdf['Rain'] >= minrain) & (np.isfinite(pdf['SID'])) & \
          (np.apply_along_axis(np.sum,0,np.array([pdf['Ice'] == m for m in vice]))) & \
          (np.apply_along_axis(np.sum,0,np.array([pdf['Month'] == m for m in vmons])))]

# Create Output Directory
try:
    os.chdir(outpath+"/"+name+"_"+V)
except:
    os.mkdir(outpath+"/"+name+"_"+V)
    os.chdir(outpath+"/"+name+"_"+V)
    
######### READ IN ALL STORM TRACKS #######
# Identify the unique months to parse through
YM = np.unique([str(int(qdf.iloc[i]['CYear']))+mons[int(qdf.iloc[i]['CMonth'])-1] for i in range(len(qdf))])

# Prep empty pdf & list
sdf = pd.DataFrame(columns=["tid","year","month","genLat","genLon",\
"maxuv","maxdpdt","maxdepth","minp","lifespan","trlen","avgarea","mcc"])

st = [] # Empty list to store all storm tracks

for i in YM: # For each valid month...
    print(i)
    
    # Load the storm tracks & SIDs
    ct = pd.read_pickle(ctpath+"/"+i[0:4]+"/systemtracks"+i+".pkl")
    sids = np.array([c.sid for c in ct])
    
    # Subset events that occurred
    rdf = qdf[(qdf['CYear'] == int(i[0:4])) & (qdf['CMonth'] == int(i[4:6]))]
    
    # Isolate each storm in turn, then record stats
    for j in range(len(rdf)):
        tr = ct[np.where(sids == int(rdf.iloc[j].SID))[0][0]]
        
        st.append(tr)
        
        # Store a row for each storm in a pandas dataframe
        row = pd.DataFrame([{"sid":tr.sid, "year":int(i[0:4]), \
        "month":int(i[4:6]), "maxuv":tr.maxUV()[0], "maxdpdt":tr.maxDpDt()[0], \
        "maxdepth":tr.maxDepth()[0], "maxdsqp":tr.maxDsqP()[0], "minp":tr.minP()[0], \
        "lifespan":tr.lifespan(), "trlen":tr.trackLength(), "avgarea":tr.avgArea(), \
        "mcc":tr.mcc(), "genLat":tr.data.loc[tr.data['type'] != 0]['lat'].iloc[0], \
        "genLon":tr.data.loc[tr.data['type'] != 0]['long'].iloc[0]},])
        sdf = sdf.append(row,ignore_index=True,sort=True)

# Write results to file
pd.to_pickle(st,outpath+"/"+name+"_"+V+"/ROS_systemtracks_"+SDAY+"_"+EDAY+".pkl")
sdf.to_csv(outpath+"/"+name+"_"+V+"/ROS_AggregatedStats"+SDAY+"_"+EDAY+".csv",index=0)

########### Calculate aggregate stats for these storms ############
# Read in attributes of reference files
ref = gdal.Open(suppath+"/"+latN,gdalconst.GA_ReadOnly)
refA = gdalnumeric.LoadFile(suppath+"/"+latN)
cellsize = ref.GetGeoTransform()[1] # cell size
n = md.daysBetweenDates(starttime,endtime)*24/hrs # Number of valid observations

print("events")
fields1 = md.aggregateEvents([st,[],[]],'system',0,refA.shape)
print("tracks")
fields2 = md.aggregateTrackWiseStats(st,[0,0,0],refA.shape)
print("points")
try:
    fields3 = md.aggregatePointWiseStats(st,n,refA.shape)
except:
    fields3 = [np.zeros_like(refA), np.zeros_like(refA), np.zeros_like(refA), \
    np.zeros_like(refA), np.zeros_like(refA), np.zeros_like(refA), \
    np.zeros_like(refA)*np.nan, np.zeros_like(refA)*np.nan, np.zeros_like(refA)*np.nan, \
    np.zeros_like(refA)*np.nan, np.zeros_like(refA)*np.nan, np.zeros_like(refA)*np.nan, \
    np.zeros_like(refA)*np.nan, np.zeros_like(refA)*np.nan, np.zeros_like(refA)*np.nan]

fields = fields1 + fields2[1:] + fields3

# Write results to file
print("smoothing and writing")
for v in varsi:
    varFieldsm = md.smoothField(fields[v],kSize) # Smooth
    md.writeNumpy_gdalObj(varFieldsm,outpath+"/"+name+"_"+V+"/ROS_"+vNames[v]+"Field"+SDAY+"_"+EDAY+".tif",ref,dtype) # Write File