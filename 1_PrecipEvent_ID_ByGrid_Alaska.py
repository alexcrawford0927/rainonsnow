'''
Author: Alex Crawford
Date Created: 29 May 2019
Date Modified: 29 May 2019
Purpose: Identify precipitation events over a grid of locations using reanalysis 
inputs. Save a separate file for each month of data.

Default units: convert to mm for precip/snow, hours for time, and deg C for T
'''

'''********************
Import Modules
********************'''
import netCDF4 as nc
import os
import numpy as np
import pandas as pd
from copy import deepcopy
import MERRA_Module as md

'''********************
Define Variables
********************'''

### File Path Variables ###
ra = "MERRA2"
ncpath = "/Volumes/Miranda/"+ra+"_nc"
ncpath1 = ncpath+"/Hourly/MERRA-LND"
ncpath2 = ncpath+"/Hourly/T2M"
outpath = "/Volumes/Miranda/RainOnSnow/PrecipDetection_ByGrid/Alaska/PrecipIdentified"

### Physical Variables ###
pconversion = 3600 # to convert precip to a value of mm
# In MERRA2 --> kg m^-2 s^-1 --> mm hr^-1 is *3600/1000*1000
rthresh = 6.096/24 # in mm/hr (equiva. to a rate of 0.254/24 = 0.01 in/day)
pthresh = 0.254 # total mm (equiv. to 0.01 in/event)

### Time Variables ###
starttime = [1980,5,1,0,0,0]
endtime = [2019,1,1,0,0,0]
daystep = [0,0,1,0,0,0]
hrs = 1 # Temporal Resolution in hours
hthresh = 2 # Number of hours w/out measurable precip needed to mark separate
# precip events

init = 1 # 0 = no initialization file; 1 = initialization file present from prior month

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
files1 = os.listdir(ncpath1)
files2 = os.listdir(ncpath2)

# Load latitude and longitude arrays from a sample file
ncf = nc.Dataset(ncpath1+"/"+files1[0])
lats = ncf.variables['lat'][:]
lons = ncf.variables['lon'][:]
examp = ncf.variables['PRECTOTLAND'][0,:,:]
ncf.close()

# Identify Alaska Grid Cells
longrid, latgrid = np.meshgrid(lons, lats)
reg1 = np.where( ( ((longrid <= -141) & (longrid >= -170) & (latgrid >= 50) & (latgrid <= 75)) |
        ((longrid <= -130) & (longrid >= -135) & (latgrid >= 54.5) & (latgrid < 56.5)) |
        ((longrid <= -131.5) & (longrid >= -137) & (latgrid >= 56.5) & (latgrid < 57.5)) |
        ((longrid <= -132.5) & (longrid >= -139) & (latgrid >= 57.5) & (latgrid < 59)) |
        ((longrid <= -134.5) & (longrid >= -141) & (latgrid >= 59) & (latgrid < 61)) )
         & (np.isfinite(examp) == 1), 1, np.nan)
rows, cols = np.where(reg1 == 1)

del examp

# Set up initial conditions
if init == 1:
    evprec, evsnof, evhrs, evtsrf, evt2m, st, stnod = pd.read_pickle(outpath+"/ActiveEvents/Active_"+str(starttime[0])+mons[starttime[1]-1]+".pkl")
else: 
    evprec = np.zeros(rows.shape)
    evsnof = np.zeros(rows.shape)
    evhrs = np.zeros(rows.shape)
    evtsrf = [[] for i in range(len(rows))]
    evt2m = [[] for i in range(len(rows))]
    st = [[] for i in range(len(rows))]
    stsnod = np.zeros(rows.shape)

time = deepcopy(starttime)
pdf = pd.DataFrame(columns = ["Year","Month","Day","Hour","X","Y","Length",\
                "Precip","Snowfall","StSnoDep","EdSnoDep",\
                "MaxTSurf","MinTSurf","MaxT2m","MinT2m"])

while time != endtime:
    # Load data for first day of analysis
    Y, M, D = str(time[0]), mons[time[1]-1], days[time[2]-1]
    
    # Load data for current day
    fN1 = [fN for fN in files1 if Y+M+D in fN][0]
    fN2 = [fN for fN in files2 if Y+M+D in fN][0]
    
    ncf1 = nc.Dataset(ncpath1+"/"+fN1)
    ncf2 = nc.Dataset(ncpath2+"/"+fN2)
    
    for h in range(24): # For each hour of the day...
        prec = ncf1['PRECTOTLAND'][h,:,:][(rows,cols)].flatten()*pconversion
        snof = ncf1['PRECSNOLAND'][h,:,:][(rows,cols)].flatten()*pconversion
        snod = ncf1['SNODP'][h,:,:][(rows,cols)].flatten()
        tsrf = ncf1['TSURF'][h,:,:][(rows,cols)].flatten()
        t2m = ncf2['T2M'][h,:,:][(rows,cols)].flatten()
        
        # For any grid cell with measurable precip...
        pis = np.where(prec >= rthresh)[0]
    
        for i in pis:
            # If there is measurable precip, and a NEW event...
            if evprec[i] == 0:
                ### Initiate a new event: ###
                st[i] = deepcopy(time)
                st[i][3] = h
                # Go back 1 hr b/c precip is accumlative over the past hour
                st[i] = md.timeAdd(starttime,[0,0,md.daysBetweenDates(starttime,st[i])-(1/24.),0,0,0])
                stsnod[i] = snod[i]
                evhrs[i] = 0
            
            ### Add the new precip, snowfall, and temp ###
            evprec[i] = evprec[i] + prec[i]
            evsnof[i] = evsnof[i] + snof[i]
            evtsrf[i].append(tsrf[i])
            evt2m[i].append(t2m[i])
        
        # If there is no measurable precipitation...
        nis = np.where((prec < rthresh) & (evprec > 0))[0]
        
        for i in nis:
            # Add one to the tally
            evhrs[i] = evhrs[i] + hrs
            
            # If the threshold for no-precip observations has past...
            if evhrs[i] >= hthresh:
                ### Terminate the event: ###
                # If there was sufficient precip...
                if evprec[i] >= pthresh:
                    # Record the event
                    ed = deepcopy(time)
                    ed[3] = h                    

                    rdf = pd.DataFrame([{"Year":st[i][0], "Month":st[i][1], \
                            "Day":st[i][2],"Hour":st[i][3],"X":cols[i],"Y":rows[i],\
                            "Length":md.daysBetweenDates(st[i],ed)*24-hthresh,\
                            "Precip":evprec[i],"Snowfall":evsnof[i],\
                            "StSnoDep":stsnod[i],"EdSnoDep":snod[i],\
                            "MaxTSurf":np.max(evtsrf[i]),"MinTSurf":np.min(evtsrf[i]),\
                            "MaxT2m":np.max(evt2m[i]),"MinT2m":np.min(evt2m[i])},])
                    
                    pdf = pdf.append(rdf, ignore_index=1, sort=1)
                
                # Reset values
                evprec[i] = 0
                evsnof[i] = 0
                evhrs[i] = 0
                evtsrf[i] = []
                evt2m[i] = []
                st[i] = []
                stsnod[i] = 0
    
    # Advance day step
    ncf1.close(), ncf2.close()
    time = md.timeAdd(time,daystep)
    
    if time[2] == 1:
        # Unit Conversions
        pdf['MaxTSurf'] = pdf['MaxTSurf'] - 273.15
        pdf['MinTSurf'] = pdf['MinTSurf'] - 273.15
        pdf['MaxT2m'] = pdf['MaxT2m'] - 273.15
        pdf['MinT2m'] = pdf['MinT2m'] - 273.15
        
        # Write to File
        pdf.to_csv(outpath+"/PrecipEvents_Alaska_Gap"+str(hthresh)+"_Rate"+str(np.round(rthresh,4))+"_Total"+str(pthresh)+"_"+Y+M+".csv",index=0)
        pd.to_pickle([evprec,evsnof,evhrs,evtsrf,evt2m,st,stsnod],outpath+"/ActiveEvents/Active_"+str(time[0])+mons[time[1]-1]+".pkl")
        
        print("Completed " + Y + M)
        
        # Restart PDF
        pdf = pd.DataFrame(columns = ["Year","Month","Day","Hour","X","Y","Length",\
                "Precip","Snowfall","StSnoDep","EdSnoDep",\
                "MaxTSurf","MinTSurf","MaxT2m","MinT2m"])
