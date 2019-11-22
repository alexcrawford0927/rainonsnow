'''
Author: Alex Crawford
Date Created: 21 Mar 2019
Date Modified: ???
Purpose: Compares snow presence, precip occurrence, and temperature between a 
reanalysis and station data. Assumes the original station data is in inches and
F and the reanlysis is in C.
'''

'''********************
Import Modules
********************'''
# Import clock:
from time import perf_counter
# Start script stopwatch. The clock starts running when time is imported
start = perf_counter()

import netCDF4 as nc
import os
import numpy as np
import pandas as pd
from osgeo import gdal, gdalnumeric, gdalconst
import CycloneModule_11_1 as md

'''*******************************************
Define Variables
*******************************************'''
# Path Variables
ra = "MERRA2"
ncpath = "/Volumes/Miranda/"+ra+"_nc"
ncpath1 = ncpath+"/Hourly/MERRA-LND"
ncpath2 = ncpath+"/Hourly/T2M"
stationpath = "/Volumes/Miranda/RainOnSnow/NOAA Data"#SNOTEL Data"
outpath = "/Volumes/Miranda/RainOnSnow/StationComparisons"

### Location Variables ###
v = 19 # (0-4, 5-9, 10-15, 16-19)
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
      70.27, 70.4322, 66.8667,64,7408] # latitude of station

### Physical Variables ###
pconversion = 3600 # to convert precip to a value of mm. day^-1
# In MERRA2 --> kg m^-2 s^-1 --> mm hr^-1 is *3600/1000*1000... then sum
sconversion =  1000 # convert from m to mm --> *1000

### Time Variables ###
starttime1 = [1980,1,1,0,0,0]
endtime1 = [2019,1,1,0,0,0]
reftime = [1900,1,1,0,0,0]
daystep = [0,0,1,0,0,0]
hrs = 1 # Temporal Resolution in hours of the ROS
dt = -9 # Time difference in hours between UTC and Local Time; -9 for Alaska

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

name = names[v]
x = xs[v]
y = ys[v]

# Load Station Data
sfiles = os.listdir(stationpath)
sf = [s for s in sfiles if names[v] in s]
sdf = pd.read_csv(stationpath+"/"+sf[0])
#sdf = pd.read_csv(stationpath+"/"+names[0]+".csv")
sdf = sdf[np.isfinite(sdf['Year'])]
starttime2 = [int(sdf['Year'][0]),int(sdf['Month'][0]),int(sdf['Day'][0]),0,0,0]
endtime2 = [int(sdf['Year'][len(sdf)-1]),int(sdf['Month'][len(sdf)-1]),int(sdf['Day'][len(sdf)-1]),0,0,0]

# Set Start and End Times
if md.daysBetweenDates(reftime,starttime1) > md.daysBetweenDates(reftime,starttime2):
    starttime = md.timeAdd(starttime1,daystep)
else:
    starttime = starttime2
    
if md.daysBetweenDates(reftime,endtime1) < md.daysBetweenDates(reftime,endtime2):
    endtime = endtime1
else:
    endtime = endtime2
    
# Subset the data frame to just these times:
sdf['Days'] = [md.daysBetweenDates(reftime,[int(sdf['Year'][i]),int(sdf['Month'][i]),int(sdf['Day'][i]),0,0,0]) for i in range(len(sdf))]
sdf2 = sdf[( sdf['Days'] >= md.daysBetweenDates(reftime,starttime) ) & ( sdf['Days'] < md.daysBetweenDates(reftime,endtime) )]

# Convert form F to C and in. to mm
sdf2['TMAX'] = ((sdf2['TMAX']-32)*5/9)
sdf2['TMIN'] = ((sdf2['TMIN']-32)*5/9)
sdf2['SNWD'] = sdf2['SNWD']*25.4
sdf2['SNOW'] = sdf2['SNOW']*25.4
sdf2['PRCP'] = sdf2['PRCP']*25.4

# Load latitude and longitude arrays from a sample file
files1 = os.listdir(ncpath1)
files2 = os.listdir(ncpath2)

ncf = nc.Dataset(ncpath1+"/"+files1[0])
lats = ncf.variables['lat'][:]
lons = ncf.variables['lon'][:]
ncf.close()

# Find nearest grid cell to the desired location
row = md.findNearest(lats,y)[1]
col = md.findNearest(lons,x)[1]

# Run loop
p, sf, sd, tmax, tmin = [], [], [], [], []
for i in np.arange(len(sdf2)):
    ii = sdf2.index[i]
    
    # Identify time
    time = [int(sdf2['Year'].iloc[i]),int(sdf2['Month'].iloc[i]),int(sdf2['Day'].iloc[i]),0,0,0]
    if time[2] == 1:
        print(time)
    
    Y = str(time[0])
    M = mons[time[1]-1]
    D = days[int(time[2]-1)]
    
    # Identify Time for Day Before
    timeA = md.timeAdd(time,[0,0,-1,0,0,0])
    YA = str(time[0])
    MA = mons[time[1]-1]
    DA = days[int(time[2]-1)]
    
    # Load NC files
    nc1 = nc.Dataset(ncpath1+"/"+[f for f in files1 if Y+M+D in f][0])
    nc2 = nc.Dataset(ncpath2+"/"+[f for f in files2 if Y+M+D in f][0])
    nc1a = nc.Dataset(ncpath1+"/"+[f for f in files1 if YA+MA+DA in f][0])
    nc2a = nc.Dataset(ncpath2+"/"+[f for f in files2 if YA+MA+DA in f][0])
    
    # Extract from the given point
    p = np.sum(np.hstack((nc1a.variables['PRECTOTLAND'][dt:,row,col],nc1.variables['PRECTOTLAND'][:dt,row,col])))
    sf = np.sum(np.hstack((nc1a.variables['PRECSNOLAND'][dt:,row,col],nc1.variables['PRECSNOLAND'][:dt,row,col]))) 
    sd = np.mean(np.hstack((nc1a.variables['SNODP'][dt:,row,col],nc1.variables['SNODP'][:dt,row,col]))) 
    tmax = np.max(np.hstack((nc2a.variables['T2M'][dt:,row,col],nc2.variables['T2M'][:dt,row,col])))
    tmin = np.min(np.hstack((nc2a.variables['T2M'][dt:,row,col],nc2.variables['T2M'][:dt,row,col])))
    
    sdf2.loc[ii,ra+'Precip'] = p
    sdf2.loc[ii,ra+'Snowfall'] = sf
    sdf2.loc[ii,ra+'SnowDepth'] = sd
    sdf2.loc[ii,ra+'TMax'] = tmax
    sdf2.loc[ii,ra+'TMin'] = tmin
    
    # Transition to Next Day
    nc1a.close(), nc2a.close()
    nc1.close(), nc2.close()

# Write to File
sdf2[ra+'Precip'] = sdf2[ra+'Precip'] * pconversion
sdf2[ra+'Snowfall'] = sdf2[ra+'Snowfall'] * pconversion
sdf2[ra+'SnowDepth'] = sdf2[ra+'SnowDepth'] * sconversion
sdf2[ra+'TMax'] = sdf2[ra+'TMax'] - 273.15
sdf2[ra+'TMin'] = sdf2[ra+'TMin'] - 273.15
    
#sdf2 = sdf2[[sdf2.columns]]
sdf2.to_csv(outpath+"/"+names[v]+"_"+ra+".csv",index=False)

print(perf_counter()-start)