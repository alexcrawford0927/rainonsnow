'''*********************************************
Author: Alex Crawford
Date Created: 31 May 2019
Date Modified: 7 Jun 2019
Purpose: Generates composites of atmospheric fields for Rain events, Non-Rain Events, 
ROS Events, Non-ROS Events, ROSF Events, and ROSNonF Events; includes tests of difference
to assess statistical significance. Must be done by point location.
*********************************************'''

'''**********
Import Modules
**********'''
import numpy as np
import pandas as pd
import os
from osgeo import gdalnumeric, gdal
import netCDF4 as nc
from scipy import stats
from copy import deepcopy
import CycloneModule_11_1 as md

'''*******************************************
Define Variables
*******************************************'''
### Location Variables ###
names = ["Aniak","Pargon Creek","Kelly Station","Fort Yukon","Teuchet Creek",\
         "Indian Pass","Monument Creek","Cooper Lake","Kenai Moose Pens","Mt Ryan",\
         "Bethel","Nome","Anchorage","Juneau","Fairbanks","Utqiagvik",\
         "Prudhoe Bay","Coleville Village","Kotzebue","Galena"]
rows = [122, 129, 135, 132, 129, 121, 129, 120, 120, 129,\
        121, 128, 121, 116, 129, 142, 140, 140, 133, 128]
cols = [33, 27, 28, 56, 55, 49, 55, 48, 47, 54,\
        29, 23, 48, 73, 51, 37, 50, 47, 28, 37]
mos = [[],[],[],[9,10,11,12,2,4],[],\
       [],[],[],[],[],\
       [10,11,12,1,2,3],[10,11,12,1,2],[10,11,12,1,2],[10,11,12,1,2,3,4],[10,11,12,1,2],[9,10],\
       [9,10],[9,10],[10,11,12,2],[10,11,12,1,2]]

### Physical Variables ###
rthresh = 6.096/24 # in mm/hr (equiva. to a rate of 0.01 in/day)
pthresh = 0.254 # total mm (equiv. to 0.1 in/event)
minrain = 2.54 # minimum total mm of rain needed to count as ROS
f30thresh = 0.95
sthresh = 0.0254 # minimum snow depth in m
vice = [3] # valid icing conditions
V = "V6"
hshift = 72 # Number of hours before (give a positive number for "before")

# Path Variables
path = "/Volumes/Miranda"
csvpath = path+"/RainOnSnow/PrecipDetection_ByPoint_withStorms/V2"
ncpaths = [path+"/MERRA2_nc/Hourly/",path+"/MERRA2_nc/Hourly/",path+"/MERRA2_nc/Hourly/","/Volumes/Ferdinand/MERRA2_nc/Hourly/"] 
suppath = "/Volumes/Miranda/Projections"
latN = "MERRA2_Lats.tif"
lonN = "MERRA2_Lons.tif"
dtype = gdal.GDT_Float64

vi = 0
v1s = ["T2M","TQV","GPH","SLP"]
v2s = ["","","500",""]
v3s = ["MERRA2.inst1_2d_asm_Nx.","MERRA2.inst1_2d_asm_Nx.","MERRA2.tavg1_2d_slv_Nx.","MERRA2.inst1_2d_asm_Nx."]
ncvars = ["T2M","TQV","H500","SLP"]

### Time Variables ###
YYYY = "19800101_20190101"
YY = "1980_2018"
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
# Load Files
ref = gdal.Open(suppath+"/"+latN)

v1, v2, v3, ncvar, ncpath = v1s[vi], v2s[vi], v3s[vi], ncvars[vi], ncpaths[vi]

#####################
##### ROS v. SOS ####
#####################
for v in [3,18]:#[3,10,11,12,14,15,16,17,18,19]:
    outpath = path+"/RainOnSnow/Composites/Pre"+str(hshift)+"hr/"+names[v]+"_"+V
    
    try:
        os.chdir(outpath)
            
    except:
        os.mkdir(outpath)
        os.chdir(outpath)
  
    pdf = pd.read_csv(csvpath+"/"+names[v]+"_"+YYYY+"_Gap"+str(hthresh)+"_Rate"+str(rthresh)+"_Total"+str(pthresh)+"_X"+str(cols[v])+"_Y"+str(rows[v])+".csv")

    print("Load all instances for Case 3 " + names[v])
    df1 = pdf[np.array(pdf["Rain"] >= minrain) & np.array(pdf["Ice"] == 3) & np.isin(pdf["Month"], mos[v])]
    
    # Loop through each hour of each precip event, extract the variable
    listD = []    
    for i in range(len(df1)): 
        t0 = [int(df1.iloc[i].Year),int(df1.iloc[i].Month),int(df1.iloc[i].Day),int(df1.iloc[i].Hour),0,0]
        t = md.timeAdd(t0,[0,0,0,-1*hshift,0,0])
        
        YDM = str(t[0])+mons[t[1]-1]+days[t[2]-1]
        
        ncfD = nc.Dataset(ncpath+"/"+v1+"/"+v2+"/"+v3+YDM+".SUB.nc")
        
        listD.append(ncfD.variables[ncvar][t[3],:,:])
        
        # Close NetCDF
        ncfD.close()
    
    # Convert to Arrays
    var1 = np.array(listD)
    
    print("Load all instances for Case 4")
    df2 = pdf[np.array(pdf["Precip"] >= minrain) & np.array(pdf["Rain"] < minrain) & \
              np.array(pdf["StSnoDep"] > sthresh) & np.array(pdf["EdSnoDep"] > sthresh) & np.isin(pdf["Month"], mos[v])]
    
    # Loop through each hour of each precip event, extract the variable
    listD = []    
    for i in range(len(df2)): 
        t0 = [int(df2.iloc[i].Year),int(df2.iloc[i].Month),int(df2.iloc[i].Day),int(df2.iloc[i].Hour),0,0]
        t = md.timeAdd(t0,[0,0,0,-1*hshift,0,0])
        
        YDM = str(t[0])+mons[t[1]-1]+days[t[2]-1]
        
        ncfD = nc.Dataset(ncpath+"/"+v1+"/"+v2+"/"+v3+YDM+".SUB.nc")
        
        listD.append(ncfD.variables[ncvar][t[3],:,:])
        
        # Close NetCDF
        ncfD.close()
    
    # Convert to Arrays
    var2 = np.array(listD)
    
    print("Calculate Averages for 3 & 4")
    var1avg = np.apply_along_axis(np.mean,0,var1)
    
    md.writeNumpy_gdalObj(np.flipud(var1avg),outpath+"/ROSEvents_"+v1+v2+"_Avg"+YYYY+".tif",ref,dtype)
    
    var2avg = np.apply_along_axis(np.mean,0,var2)
    
    md.writeNumpy_gdalObj(np.flipud(var2avg),outpath+"/SOSEvents_"+v1+v2+"_Avg"+YYYY+".tif",ref,dtype)
    
    print("Calculate Difference Tests 3 & 4")
    t12D, p12D = np.zeros_like(var1[0,:,:]), np.zeros_like(var1[0,:,:])
    
    for ri in range(t12D.shape[0]):
        if ri%10 == 0:
            print(ri)
        for ci in range(t12D.shape[1]):
            t12D[ri,ci], p12D[ri,ci] = stats.mannwhitneyu(var1[:,ri,ci],var2[:,ri,ci],True,"two-sided")
    
    md.writeNumpy_gdalObj(np.flipud(var1avg-var2avg),outpath+"/ROS_v_SOS_"+v1+"_Diff"+YYYY+".tif",ref,dtype)
    md.writeNumpy_gdalObj(np.flipud(p12D),outpath+"/ROS_v_SOS_"+v1+"_MWUp"+YYYY+".tif",ref,dtype)
    
    del var1, var2, t12D, p12D, df1, df2
    del var1avg, var2avg
    
    ##########################
    ##### ROSF v. ROSNonF ####
    ##########################
    #print("Load all instances for Case 5")
    #df1 = pdf[np.array(pdf["Rain"] >= minrain) & np.array(pdf["Ice"] == 3) & np.array(pdf["Frz30"] >= f30thresh) & np.isin(pdf["Month"], mos[v])]
    ## Loop through each hour of each precip event, record the time
    #ts = []
    #for i in range(len(df1)): 
    #    t0 = [int(df1.iloc[i].Year),int(df1.iloc[i].Month),int(df1.iloc[i].Day),int(df1.iloc[i].Hour),0,0]
    #    ts.append( md.timeAdd(t0,[0,0,0,-1*hshift,0,0]) )
    #
    ## For each of those month/day/hour combos, load every relevant year
    ## Extract values per hour and place in an array stack
    ## Do this for EVERY variable of interest
    #ts = [t for t in ts if ((t[1] != 2) | (t[2] != 29))]
    #MDHs = [mons[t[1]-1]+days[t[2]-1]+"_"+hours[t[3]] for t in ts]
    #
    #listD = []     
    #for MDH in np.unique(MDHs):
    #    ncfD = nc.Dataset(ncpath+"/"+v1+"_"+YY+"/"+v2+"/"+v1+v2+"_"+YY+"_"+MDH+".nc")
    #
    #    tis = np.where(np.array(MDHs) == MDH)[0]
    #    
    #    for t in tis:
    #        y = ts[t][0] - int(YYYY[:4])
    #        listD.append(ncfD.variables[v1+'det'][y,:,:])
    #    
    #    # Close NetCDF
    #    ncfD.close()
    #
    ## Convert to Arrays
    #var1 = np.array(listD)
    #
    #print("Load all instances for Case 6")
    #df2 = pdf[np.array(pdf["Rain"] >= minrain) & np.array(pdf["Ice"] == 3) & np.array(pdf["Frz30"] < f30thresh) & np.isin(pdf["Month"], mos[v])]
    ## Loop through each hour of each precip event, record the time
    #ts = []
    #for i in range(len(df2)): 
    #    t0 = [int(df2.iloc[i].Year),int(df2.iloc[i].Month),int(df2.iloc[i].Day),int(df2.iloc[i].Hour),0,0]
    #    ts.append( md.timeAdd(t0,[0,0,0,-1*hshift,0,0]) )
    
    #
    ## For each of those month/day/hour combos, load every relevant year
    ## Extract values per hour and place in an array stack
    ## Do this for EVERY variable of interest
    #ts = [t for t in ts if ((t[1] != 2) | (t[2] != 29))]
    #MDHs = [mons[t[1]-1]+days[t[2]-1]+"_"+hours[t[3]] for t in ts]
    #
    #listD = []   
    #for MDH in np.unique(MDHs):
    #    ncfD = nc.Dataset(ncpath+"/"+v1+"_"+YY+"/"+v2+"/"+v1+v2+"_"+YY+"_"+MDH+".nc")
    #
    #    tis = np.where(np.array(MDHs) == MDH)[0]
    #    
    #    for t in tis:
    #        y = ts[t][0] - int(YYYY[:4])
    #        listD.append(ncfD.variables[v1+'det'][y,:,:])
    #    
    #    # Close NetCDF
    #    ncfD.close()
    #
    ## Convert to Arrays
    #var2 = np.array(listD)
    #
    #print("Calculate Averages for 5 & 6")
    #var1avg = np.apply_along_axis(np.mean,0,var1)
    #
    #md.writeNumpy_gdalObj(np.flipud(var1avg),outpath+"/ROSFEvents_"+v1+v2+"_Avg"+YYYY+".tif",ref,dtype)
    #
    #var2avg = np.apply_along_axis(np.mean,0,var2)
    #
    #md.writeNumpy_gdalObj(np.flipud(var2avg),outpath+"/ROSNonFEvents_"+v1+v2+"_Avg"+YYYY+".tif",ref,dtype)
    #
    #print("Calculate Difference Tests 5 & 6")
    #t12D, p12D = np.zeros_like(var1[0,:,:]), np.zeros_like(var1[0,:,:])
    #
    #for ri in range(t12D.shape[0]):
    #    if ri%10 == 0:
    #        print(ri)
    #    for ci in range(t12D.shape[1]):
    #        t12D[ri,ci], p12D[ri,ci] = stats.mannwhitneyu(var1[:,ri,ci],var2[:,ri,ci],True,"two-sided")
    #
    #md.writeNumpy_gdalObj(np.flipud(var1avg-var2avg),outpath+"/ROSF_v_ROSNonF_"+v1+"_Diff"+YYYY+".tif",ref,dtype)
    #md.writeNumpy_gdalObj(np.flipud(p12D),outpath+"/ROSF_v_ROSNonF_"+v1+"_MWUp"+YYYY+".tif",ref,dtype)
    #
    #del var1, var2, df1, df2, t12D, p12D
    #del var1avg, var2avg
    #
    #########################
    ##### Rain v. NoRain ####
    #########################
    #print("Load all instances for Case 1")
    #df1 = pdf[np.array(pdf["Rain"] >= minrain) & np.isin(pdf["Month"],mos[v])]
    #
    ## Loop through each hour of each precip event, record the time
    #ts = []
    #for i in range(len(df1)): 
    #    t0 = [int(df1.iloc[i].Year),int(df1.iloc[i].Month),int(df1.iloc[i].Day),int(df1.iloc[i].Hour),0,0]
    #    ts.append( md.timeAdd(t0,[0,0,0,-1*hshift,0,0]) )
    #
    ## For each of those month/day/hour combos, load every relevant year
    ## Extract values per hour and place in an array stack
    ## Do this for EVERY variable of interest
    #ts = [t for t in ts if ((t[1] != 2) | (t[2] != 29))]
    #MDHs = [mons[t[1]-1]+days[t[2]-1]+"_"+hours[t[3]] for t in ts]
    #
    #listD = []    
    #for MDH in np.unique(MDHs):
    #    ncfD = nc.Dataset(ncpath+"/"+v1+"_"+YY+"/"+v2+"/"+v1+v2+"_"+YY+"_"+MDH+".nc")
    #
    #    tis = np.where(np.array(MDHs) == MDH)[0]
    #    
    #    for t in tis:
    #        y = ts[t][0] - int(YYYY[:4])
    #        listD.append(ncfD.variables[v1+'det'][y,:,:])
    #    
    #    # Close NetCDF
    #    ncfD.close()
    #
    ## Convert to Arrays
    #var1 = np.array(listD)
    #
    #print("Load all instances for Case 2")
    #df2 = pdf[np.array(pdf["Rain"] < pthresh) & np.array(pdf["Precip"] >= minrain) & np.isin(pdf["Month"],mos[v])]
    #
    ## Loop through each hour of each precip event, record the time
    #ts = []
    #for i in range(len(df2)): 
    #    t0 = [int(df2.iloc[i].Year),int(df2.iloc[i].Month),int(df2.iloc[i].Day),int(df2.iloc[i].Hour),0,0]
    #    ts.append( md.timeAdd(t0,[0,0,0,-1*hshift,0,0]) )
    #
    ## For each of those month/day/hour combos, load every relevant year
    ## Extract values per hour and place in an array stack
    ## Do this for EVERY variable of interest
    #ts = [t for t in ts if ((t[1] != 2) | (t[2] != 29))]
    #MDHs = [mons[t[1]-1]+days[t[2]-1]+"_"+hours[t[3]] for t in ts]
    #
    #listD = []     
    #for MDH in np.unique(MDHs):
    #    ncfD = nc.Dataset(ncpath+"/"+v1+"_"+YY+"/"+v2+"/"+v1+v2+"_"+YY+"_"+MDH+".nc")
    #
    #    tis = np.where(np.array(MDHs) == MDH)[0]
    #    
    #    for t in tis:
    #        y = ts[t][0] - int(YYYY[:4])
    #        listD.append(ncfD.variables[v1+'det'][y,:,:])
    #    
    #    # Close NetCDF
    #    ncfD.close()
    #
    ## Convert to Arrays
    #var2 = np.array(listD)
    #
    #print("Calculate Averages for 1 & 2")
    #var1avg = np.apply_along_axis(np.mean,0,var1)
    #
    #md.writeNumpy_gdalObj(np.flipud(var1avg),outpath+"/RainEvents_"+v1+v2+"_Avg"+YYYY+".tif",ref,dtype)
    #
    #var2avg = np.apply_along_axis(np.mean,0,var2)
    #
    #md.writeNumpy_gdalObj(np.flipud(var2avg),outpath+"/NoRainEvents_"+v1+v2+"_Avg"+YYYY+".tif",ref,dtype)
    #
    #print("Calculate Difference Tests 1 & 2")
    #t12D, p12D = np.zeros_like(var1[0,:,:]), np.zeros_like(var1[0,:,:])
    #
    #for ri in range(t12D.shape[0]):
    #    if ri%10 == 0:
    #        print(ri)
    #    for ci in range(t12D.shape[1]):
    #        t12D[ri,ci], p12D[ri,ci] = stats.mannwhitneyu(var1[:,ri,ci],var2[:,ri,ci],True,"two-sided")
    #
    #md.writeNumpy_gdalObj(np.flipud(var1avg-var2avg),outpath+"/Rain_v_NoRain_"+v1+"_Diff"+YYYY+".tif",ref,dtype)
    #md.writeNumpy_gdalObj(np.flipud(p12D),outpath+"/Rain_v_NoRain_"+v1+"_MWUp"+YYYY+".tif",ref,dtype)
    #
    #del var1, var2, t12D, p12D, df1, df2
    #del var1avg, var2avg