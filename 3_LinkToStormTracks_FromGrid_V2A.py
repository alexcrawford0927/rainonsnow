'''
Author: Alex Crawford
Date Created: 15 Mar 2019
Date Modified: 3 Jun 2019 --> modified to work with a "grid-based" detection of precip events
Purpose: Identify which storms from a Lagrangian cyclone algorithm are 
associated with a particular class of precipitation event at a given location.
'''

'''********************
Import Modules
********************'''
import netCDF4 as nc
import os
import numpy as np
import pandas as pd
import scipy
from copy import deepcopy
from osgeo import gdal, gdalnumeric
import MERRA_Module as md

'''********************
Define Variables
********************'''
ra = "MERRA2"
region = "Alaska"
V = "V2"

ncpath = "/Volumes/Miranda/"+ra+"_nc/Hourly/MERRA-LND"
csvpath = "/Volumes/Miranda/RainOnSnow/PrecipDetection_ByGrid/"+region+"/ROSIdentified_"+V
csvpath2 = "/Volumes/Miranda/RainOnSnow/PrecipDetection_ByPoint_withStorms/"+V+"A"

cfpath = "/Volumes/Ferdinand/ArcticCyclone/detection11_3AM2/CycloneFields"
ctpath = "/Volumes/Ferdinand/ArcticCyclone/detection11_3AM2/SystemTracks"
suppath = "/Volumes/Ferdinand/Projections"
latN = "EASE2_N0_100km_Lats.tif"
lonN = "EASE2_N0_100km_Lons.tif"

### Location Variables ###
v = 16 # (0-4, 5-9, 10-15, 16-19)
names = ["Aniak","Pargon Creek","Kelly Station","Fort Yukon","Teuchet Creek",\
         "Indian Pass","Monument Creek","Cooper Lake","Kenai Moose Pens","Mt Ryan",\
         "Bethel","Nome","Anchorage","Juneau","Fairbanks","Utqiagvik",\
         "Prudhoe Bay","Coleville Village","Kotzebue","Galena"]
xs = [-159.58, -163.1, -162.28, -145.25, -145.52,\
      -149.48, -145.87, -149.69, -150.48, -146.15,\
      -161.8293, -165.44, -149.7833, -134.564, -147.8761, -156.7815,\
      -148.57, -150.4094, -162.6333,-156.7855] # longitude of station
ys = [61.58, 64.99, 67.93, 66.57, 64.95,\
      61.07, 65.08, 60.39, 60.73, 65.25,\
      60.785, 64.511, 61.2, 58.357, 64.8039, 71.2834,\
      70.27, 70.4322, 66.8667,64.7408] # latitude of station

### Physical Variables ###
rthresh = 6.096/24 # in mm/hr (equiva. to a rate of 0.01 in/day)
pthresh = 0.254 # total mm (equiv. to 0.01 in/event)
dthresh = 1200000 # max distance between a point and it's assoc. sys center

### Time Variables ###
starttime = [1980,1,1,0,0,0]
endtime = [2019,1,1,0,0,0]
reftime = [1900,1,1,0,0,0]
monthstep = [0,1,0,0,0,0]
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

'''*******************************************
Main Analysis
*******************************************'''
print("Load Data")
########## READ IN ROS DATA ############
files1 = os.listdir(ncpath)
files1 = [f for f in files1 if f.startswith(ra)]

# Load latitude and longitude arrays from a sample file
ncf = nc.Dataset(ncpath+"/"+files1[0])
lats = ncf.variables['lat'][:]
lons = ncf.variables['lon'][:]
ncf.close()

# Find nearest grid cell to the desired location
name = names[v]
x = xs[v] 
y = ys[v] 

row = md.findNearest(lats,y)[1]
col = md.findNearest(lons,x)[1]

# Load precip events for the location of choice
pdf = pd.DataFrame()

time = deepcopy(starttime)
while time != endtime:
    # Open monthly PDF
    tdf = pd.read_csv(csvpath+"/PrecipEvents_"+region+"_Gap"+str(hthresh)+"_Rate"+str(np.round(rthresh,4))+"_Total"+str(pthresh)+"_"+str(time[0])+mons[time[1]-1]+".csv")

    # Subset to the location of interest
    ldf = tdf[(tdf["X"] == col) & (tdf["Y"] == row)]
    
    # Append to main pdf
    pdf = pdf.append(ldf, sort=True, ignore_index=True)
    
    # Advance to next time
    time = md.timeAdd(time,monthstep)

########## READ IN BASIC CYCLONE INFO ############
# Load lats and lons
elats = gdalnumeric.LoadFile(suppath+"/"+latN)
elons = gdalnumeric.LoadFile(suppath+"/"+lonN)

# Find the appropriate location in the cyclone tracking grid
eloc = md.findNearest_latlong(elats,elons,y,x)

pdf['CMonth'] = 0
pdf['CYear'] = 0
pdf['SID'] = 0
SDAY = str(starttime[0]) + mons[starttime[1]-1] + days[starttime[2]-1]
EDAY = str(endtime[0]) + mons[endtime[1]-1] + days[endtime[2]-1]

for i in range(0,len(pdf)):
    link = 0

    # Find the tracking hour that most closely follows to the median time of 
    ## the precip event
    MedHrs = pdf['EndTime'][i]-(pdf['Length'][i]/2.)
    resid = MedHrs%t
    if resid < t/2.:
        MedHrs = MedHrs - (resid)
    else:
        MedHrs = MedHrs + (t-resid)
    MedTime = md.timeAdd(reftime,(0,0,0,MedHrs,0,0))
    MedTime1 = md.timeAdd(MedTime,[0,1,0,0,0,0])
        
    yr = MedTime[0]
    m = MedTime[1]
    d = MedTime[2]
    h = int(MedTime[3])
    DATE = str(yr)+mons[m-1]+days[d-1]+"_"+hours[h]
    
    # Load Cyclone Field & Tracks for that time
    ct0 = pd.read_pickle(ctpath+"/"+str(yr)+"/systemtracks"+str(yr)+mons[m-1]+".pkl")
    ct1 = pd.read_pickle(ctpath+"/"+str(MedTime1[0])+"/systemtracks"+str(MedTime1[0])+mons[MedTime1[1]-1]+".pkl")
    cf = pd.read_pickle(cfpath+"/"+str(yr)+"/"+months[m-1]+"/CF"+DATE+".pkl")
    cAreas, nC = scipy.ndimage.measurements.label(cf.fieldAreas)
    
    # Identify which storms existed at this time
    cta = [ct for ct in ct0 if (MedHrs/24. in list(ct.data.time))] + [ct for ct in ct1 if (MedHrs/24. in list(ct.data.time))]
    
    # Limit to type 1 centers (not ones that have just died or split)
    cttyp = np.array([int(ct.data['type'].loc[ct.data['time'] == MedHrs/24.]) for ct in cta])
    cta = [cta[j] for j in range(len(cta)) if cttyp[j] != 0]
    
    ctaxs = [int(ct.data['x'].loc[ct.data['time'] == MedHrs/24.]) for ct in cta]
    ctays = [int(ct.data['y'].loc[ct.data['time'] == MedHrs/24.]) for ct in cta]
    ctp = np.array([float(ct.data['p_cent'].loc[ct.data['time'] == MedHrs/24.]) for ct in cta])

    # Option 1: Location lies within a system's area:
    if cAreas[eloc[0][0],eloc[1][0]] > 0:
        ctareas = [cAreas[ctays[j],ctaxs[j]] for j in range(len(cta))]
        
        try:
            ctb = cta[int(np.where(cAreas[eloc[0][0],eloc[1][0]] == ctareas)[0])]
            ctbtime = md.timeAdd(reftime,[0,0,list(ctb.data.time)[-1],0,0,0])
            pdf.loc[i,'SID'] = ctb.sid
            pdf.loc[i,'CMonth'] = ctbtime[1]
            pdf.loc[i,'CYear'] = ctbtime[0]
            
            link = 1
        
        except:
            link = 0            
        
    # Option 2: Location lies within some distance of a cyclone system
    ## Only comes into play if there isn't a more direct area overlap.
    if link == 0:
        # Calculate distance to each system center
        dista = np.array([md.haversine(y,elats[ctays[j],ctaxs[j]],x,elons[ctays[j],ctaxs[j]]) for j in range(len(cta))])
        
        # If there's 1 or more possible storm(s)
        if sum(dista < dthresh) > 0:
            # Identify the closest storm
            ctb = cta[np.where(dista == np.min(dista))[0][0]]
            # Record its info in the main data frame
            ctbtime = md.timeAdd(reftime,[0,0,list(ctb.data.time)[-1],0,0,0])
            pdf.loc[i,'SID'] = ctb.sid
            pdf.loc[i,'CMonth'] = ctbtime[1]
            pdf.loc[i,'CYear'] = ctbtime[0]
        
        # If there are no possible storms:
        else:
            # Record NANs so that it's clear that no good link exists
            pdf.loc[i,'SID'] = np.nan
            pdf.loc[i,'CMonth'] = np.nan
            pdf.loc[i,'CYear'] = np.nan 
    
    if i%100 == 0:
        print(str(i) + " of " + str(len(pdf)))
        pdf.to_csv(csvpath2+"/"+name+"_"+SDAY+"_"+EDAY+"_Gap"+str(hthresh)+"_Rate"+str(np.round(rthresh,4))+"_Total"+str(pthresh)+"_X"+str(col)+"_Y"+str(row)+".csv",index=0)
