'''
Author: Alex Crawford
Date Created: 13 Mar 2019
Date Modified: 30 May 2019 --> Added 24-hr slot
Purpose: Given a precipitation event, identify whether it counts as an icing
event and for how long the ice is likely to persist.
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
ra = "MERRA2"

ncpath = "/Volumes/Miranda/"+ra+"_nc/Hourly/MERRA-LND"
csvpath = "/Volumes/Miranda/RainOnSnow/PrecipDetection_ByGrid/Alaska/PrecipIdentified"
outpath = "/Volumes/Miranda/RainOnSnow/PrecipDetection_ByGrid/Alaska/ROSIdentified_V2"

### Physical Variables ###
pconversion = 3600 # to convert precip to a value of mm
# In MERRA2 --> kg m^-2 s^-1 --> mm hr^-1 is *3600/1000*1000
rthresh = 6.096/24 # in mm/hr (equiva. to a rate of 0.01 in/day)
pthresh = 0.254 # total mm (equiv. to 0.01 in/event)
tthresh = -10 # minimum temperature (deg C) for which rain is allowed to be detected
sthresh = 0.0254 # minimum snow depth in meters

### Time Variables ###
starttime = [2000,1,1,0,0,0]
endtime = [2019,1,1,0,0,0]
reftime = [1900,1,1,0,0,0]
daystep = [0,0,1,0,0,0]
hrs = 1 # Temporal Resolution in hours
hthresh = 2 # Number of hours w/out measurable precip needed to mark separate
# precip events (also used as gap for re-freezing check)

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
########## READ IN INITAL DATA ############
files1 = os.listdir(ncpath)
files1 = [f for f in files1 if f.startswith(ra)]

df = pd.DataFrame()

# Start Time Loop
time1 = deepcopy(starttime)
t = 24*md.daysBetweenDates(reftime,time1)
tO = []
while time1 != endtime:
    # Identify date
    Y, M, D = str(time1[0]), mons[time1[1]-1], days[time1[2]-1]
    
    # Load data for current day
    fN = [fN for fN in files1 if Y+M+D in fN][0]
    ncf = nc.Dataset(ncpath+"/"+fN)
    
    # If the 1st day of the month, load new CSV file
    if time1[2] == 1:
        print("Starting " + Y + " " + M)
        pdf = pd.read_csv(csvpath+"/PrecipEvents_Alaska_Gap"+str(hthresh)+"_Rate"+str(np.round(rthresh,4))+"_Total"+str(pthresh)+"_"+Y+M+".csv")

        ########## Identify Icing Events ############
        pdf['Rain'] = pdf['Precip'] - pdf['Snowfall']
        pdf['Ice'] = 0
        pdf.loc[(pdf['MaxT2m'] > tthresh) & (pdf['Rain'] > pthresh) & (pdf['StSnoDep'] > sthresh) & (pdf['EdSnoDep'] <= sthresh),'Ice'] = 1
        pdf.loc[(pdf['MaxT2m'] > tthresh) & (pdf['Rain'] > pthresh) & (pdf['StSnoDep'] <= sthresh) & (pdf['MaxTSurf'] <= 0),'Ice'] = 2
        pdf.loc[(pdf['MaxT2m'] > tthresh) & (pdf['Rain'] > pthresh) & (pdf['StSnoDep'] > sthresh) & (pdf['EdSnoDep'] > sthresh),'Ice'] = 3
        
        ########## Identify Number of Frozen Days that Follow Icing ############
        pdf['EndTime'] = [24*md.daysBetweenDates(reftime,(pdf['Year'][i],pdf['Month'][i],pdf['Day'][i],pdf['Hour'][i],0,0)) for i in range(len(pdf))] + pdf['Length']
        pdf['Frz01'] = 0
        pdf['Frz10'] = 0
        pdf['Frz30'] = 0
        pdf['Frz90'] = 0
        
        ########## Append Input PDF to Main DF #########
        df = df.append(pdf, ignore_index=1)
        tO.append( 24*md.daysBetweenDates(reftime,md.timeAdd(time1,[0,1,0,0,0,0])) )
    
    # Measure Total Number of Frozen Hours
    for h in range(24): # For each hour of the day...
        # Extract temperature
        tsurf = ncf['TSURF'][h,:,:] <= 273.15 # units are Kelvin

        # Identify Relevant Events Based on Time
        indices01 = df.loc[((t - df['EndTime']) <= 24) & ((t - df['EndTime']) > 0)].index
        indices10 = df.loc[((t - df['EndTime']) <= 240) & ((t - df['EndTime']) > 0)].index
        indices30 = df.loc[((t - df['EndTime']) <= 720) & ((t - df['EndTime']) > 0)].index
        indices90 = df.loc[((t - df['EndTime']) <= 2160) & ((t - df['EndTime']) > 0)].index

        # Add Info For Fixed-Time Examinations        
        df.loc[indices01,'Frz01'] += tsurf[df.loc[indices01,'Y'],df.loc[indices01,'X']]
        df.loc[indices10,'Frz10'] += tsurf[df.loc[indices10,'Y'],df.loc[indices10,'X']]
        df.loc[indices30,'Frz30'] += tsurf[df.loc[indices30,'Y'],df.loc[indices30,'X']]
        df.loc[indices90,'Frz90'] += tsurf[df.loc[indices90,'Y'],df.loc[indices90,'X']]
        
        t =  t + 1
    
    # Advance day step
    ncf.close()
    time1 = md.timeAdd(time1,daystep)

    # If it has been long enough, write to file by month...
    if t-tO[0] >= 2160:
        timeO = md.timeAdd(reftime,[0,0,0,tO[0]-1,0,0])
        
        YO = str(timeO[0])
        MO = mons[timeO[1]-1]

        # Subset Data Frames
        odf = df[df['EndTime'] < tO[0]]
        df = df[df['EndTime'] >= tO[0]]
        
        tO = tO[1:] # removes month from consideration
        
        # Unit Conversions (to a %)
        odf['Frz01'] = odf['Frz01']/24.
        odf['Frz10'] = odf['Frz10']/240.
        odf['Frz30'] = odf['Frz30']/720.
        odf['Frz90'] = odf['Frz90']/2160.
        
        # Write to file
        odf.to_csv(csvpath+"/PrecipEvents_Alaska_Gap"+str(hthresh)+"_Rate"+str(np.round(rthresh,4))+"_Total"+str(pthresh)+"_"+YO+MO+"m.csv",index=0)
        

# Write Preliminary Files for all remaining Months
for i in range(len(tO)):
    timeO = md.timeAdd(reftime,[0,0,0,tO[i]-1,0,0])
    
    YO = str(timeO[0])
    MO = mons[timeO[1]-1]

    # Subset Data Frames
    odf = df[df['EndTime'] < tO[i]]
    df = df[df['EndTime'] >= tO[i]]
    
    # Unit Conversions (to a %)
    odf.index = range(len(odf))
    odf['Frz01'] = [odf['Frz01'][j]/np.min([t-odf['EndTime'][j],24.]) for j in range(len(odf))]
    odf['Frz10'] = [odf['Frz10'][j]/np.min([t-odf['EndTime'][j],240.]) for j in range(len(odf))]
    odf['Frz30'] = [odf['Frz30'][j]/np.min([t-odf['EndTime'][j],720.]) for j in range(len(odf))]
    odf['Frz90'] = [odf['Frz90'][j]/np.min([t-odf['EndTime'][j],2160.]) for j in range(len(odf))]
    
    odf.to_csv(outpath+"/PrecipEvents_Alaska_Gap"+str(hthresh)+"_Rate"+str(np.round(rthresh,4))+"_Total"+str(pthresh)+"_"+YO+MO+".csv",index=0)