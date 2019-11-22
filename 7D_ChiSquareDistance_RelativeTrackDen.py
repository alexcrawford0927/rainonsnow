'''
Author: Alex Crawford
Date Created: 24 Apr 2019
Date Modified: 31 May 2019 --> Added more parameterizaton
                11 Jun 2019 --> Switch from absolute to relative measure of storm density for determining area of interest
                4 Jul 2019 --> Make it relative
Purpose: Runs a chi-square distance test for grid cells around Alaska for storm
track presence/absence.
'''

'''********************
Import Modules
********************'''
import numpy as np
import pandas as pd
from osgeo import gdal, gdalnumeric
import CycloneModule_11_1 as md

def chidist(inArr):
    '''Calculates the chi-square distance for an array of counts with rows 
    representing profiles for each observation and columns representing 
    the measured parameters.  The input array should be discrete count data.
    Warning: This is not the best distance metric for binary data or for
    continuous data, but the function will run on any numerical data.
    
    Returns a square 2-D numpy array with dimensions equal to the number of 
    rows (i.e., observations or sites) of the input. NaNs are filled in for 
    the diagonal and upper-right section of the matrix because removing the 
    redundancy makes the process run faster. 
    '''
    # Step 0: Identify number of rows and columns
    Nr, Nc = inArr.shape
    
    # Step 1: Calculate relative proportions for each row
    tots = [float(np.sum(inArr[row,:])) for row in range(Nr)]
    tots = np.repeat(tots,Nc).reshape(Nr,Nc)
    arr = inArr/tots

    # Step 2: Establish mean proporation for each parameter (each column)
    means = [np.mean(arr[:,col]) for col in range(Nc)]
    
    # Step 3: Calculate chi distance for each pair of of rows
    chi = np.zeros((Nr,Nr))*np.nan
    for ii in range(Nr):
        for jj in range(ii):
            chi[ii,jj] = np.nansum([(((arr[ii,col]-arr[jj,col]))**2/means[col]) for col in range(len(means))])**0.5

    return chi

'''*******************************************
Define Variables
*******************************************'''

### Location Variables ###
v = 19 # (0-4, 5-9, 10-15, 16-19)
names = ["Aniak","Pargon Creek","Kelly Station","Fort Yukon","Teuchet Creek",\
         "Indian Pass","Monument Creek","Cooper Lake","Kenai Moose Pens","Mt Ryan",\
         "Bethel","Nome","Anchorage","Juneau","Fairbanks","Utqiagvik",\
         "Prudhoe Bay","Coleville Village","Kotzebue","Galena"]
rows = [122, 129, 135, 132, 129, 121, 129, 120, 120, 129,\
        121, 128, 121, 116, 129, 142, 140, 140, 133, 127]
cols = [33, 27, 28, 56, 55, 49, 55, 48, 47, 54,\
        29, 23, 48, 73, 51, 37, 50, 47, 28, 37]

kSize = 2 # smoothing radius of cyclone tracks
sigma = 1000. # std dev for Gaussian; larger means slower fall-off
td_thresh = 5 # Relative Threshold for area of interest (between 0 and 100) -- exclusive
ct_thresh = 2 # Absolute Threshold for area of interest (between 0 and positive infinity) -- inclusive
nMC = 1000 # Number of Monte Carlo Simulations
V = "V6"
T1, T2 = "ROS", "SOS"

### Path Variables ###
path = "/Volumes/Miranda/RainOnSnow"
inpath = path+"/PrecipDetection_Aggregation/"
cpath = "/Volumes/Ferdinand/ArcticCyclone/detection11_3AM2/SystemTracks"

dtype = gdal.GDT_Float64
suppath = "/Volumes/Ferdinand/Projections"
latN = "EASE2_N0_100km_Lats.tif"

### Time Variables ###
starttime = [1980,1,1,0,0,0]
endtime = [2019,1,1,0,0,0]
reftime = [1900,1,1,0,0,0]

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
########## READ IN ROS DATA ############
print("Load Data")
SDAY = str(starttime[0]) + mons[starttime[1]-1] + days[starttime[2]-1]
EDAY = str(endtime[0]) + mons[endtime[1]-1] + days[endtime[2]-1]

name = names[v]
pdf1 = pd.read_csv(inpath+"/"+name+"_"+V+"/"+T1+"_AggregatedStats"+SDAY+"_"+EDAY+".csv")
pdf2 = pd.read_csv(inpath+"/"+name+"_"+V+"/"+T2+"_AggregatedStats"+SDAY+"_"+EDAY+".csv")

# Identify the number of storms in each category
n1 = pdf1.shape[0]
n2 = pdf2.shape[0]

# Create Gaussian Kernel
x = np.arange(-1*kSize,kSize+1,1)
y = np.arange(-1*kSize,kSize+1,1)
xx, yy = np.meshgrid(x,y)
k = np.exp(-1/(2*sigma**2) * (xx**2 + yy**2))

# Identify Locations of Interest
td1 = gdalnumeric.LoadFile(inpath+"/"+name+"_"+V+"/"+T1+"_trkdenField"+SDAY+"_"+EDAY+".tif")
td2 = gdalnumeric.LoadFile(inpath+"/"+name+"_"+V+"/"+T2+"_trkdenField"+SDAY+"_"+EDAY+".tif")
ref = gdal.Open(suppath+"/"+latN)

########### AGGREGATE CYCLONE LOCATIONS ############
print("Calculate Observed Dissimilarity")
# create an list to store the counts
counts1 = [] # For First Category
counts2 = [] # For Second Category

# Find each storm for PDF #1
print (T1 + " Count: " + str(n1))

for i in range(n1):
    if i%10 == 0:
        print(i)
    
    sid = pdf1.iloc[i]['sid']
    y = int(pdf1.iloc[i]['year'])
    m = int(pdf1.iloc[i]['month'])
     
    # Read in Cyclone data, extract the data frame for locations    
    cycs = pd.read_pickle(cpath+"/"+str(y)+"/systemtracks"+str(y)+mons[m-1]+".pkl")
    cyc = [c for c in cycs if c.sid == sid][0]
    cdata = cyc.data.loc[cyc.data.type != 0]
    
    counts1a = np.zeros_like(td1) # Create empty count field for storm

    # For each observation...
    for j in range(cdata.shape[0]):
        col = int(cdata.iloc[j].x)
        row = int(cdata.iloc[j].y)
        
        # Add to the count
        counts1a[(row-kSize):(row+kSize+1),(col-kSize):(col+kSize+1)] = k + counts1a[(row-kSize):(row+kSize+1),(col-kSize):(col+kSize+1)]

    # Append to list
    counts1.append(counts1a > 0)
    
# Find each storm for PDF #2
print (T2 + " Count: " + str(n2))
for i in range(n2):
    if i%10 == 0:
        print(i)
        
    sid = pdf2.iloc[i]['sid']
    y = int(pdf2.iloc[i]['year'])
    m = int(pdf2.iloc[i]['month'])
    
    # Read in Cyclone data, extract the data frame for locations
    cycs = pd.read_pickle(cpath+"/"+str(y)+"/systemtracks"+str(y)+mons[m-1]+".pkl")
    cyc = [c for c in cycs if c.sid == sid][0]
    cdata = cyc.data.loc[cyc.data.type != 0]

    counts2a = np.zeros_like(td2) # Create empty count field for storm

    # For each observation...
    for j in range(cdata.shape[0]):
        col = int(cdata.iloc[j].x)
        row = int(cdata.iloc[j].y)
        
        # Add to the count
        counts2a[(row-kSize):(row+kSize+1),(col-kSize):(col+kSize+1)] = k + counts2a[(row-kSize):(row+kSize+1),(col-kSize):(col+kSize+1)]

    # Append to list
    counts2.append(counts2a > 0)
    
# Find total counts for both categories of storm
sum1 = np.apply_along_axis(np.sum,0,np.array(counts1))
sum2 = np.apply_along_axis(np.sum,0,np.array(counts2))

# Write to File
md.writeNumpy_gdalObj(sum1,inpath+"/"+name+"_"+V+"/"+T1+"_trkdenRAW_"+SDAY+"_"+EDAY+".tif",ref,dtype)
md.writeNumpy_gdalObj(sum2,inpath+"/"+name+"_"+V+"/"+T2+"_trkdenRAW_"+SDAY+"_"+EDAY+".tif",ref,dtype)

aoi = np.where( ( (sum1+sum2)/(n1+n2)*100 > td_thresh) & (sum1+sum2 >= ct_thresh) )

sums = np.vstack((sum1[aoi],sum2[aoi]))

# Calculate Dissimilarity of Real Data
chi = chidist(sums)[1,0]

# Remove unnecessary files
#del sum1, sum2, sums, col, row, i, j, cdata, cycs, sid, counts1a, counts2a

############ MONTE CARLO ############
print("Start Monte Carlo Simulation")
# Combine all storms into one list
counts = np.array(counts1 + counts2)

# Create array to fill with chi distance results from Monte Carlo
chiMC = np.zeros((nMC))

# Perform Monte Carlo Simulation
for mc in range(nMC): 
    if mc%10 == 0:
        print(mc)
    
    # Generate a set of random integers to subsample the full cyclone population
    i1 = np.random.choice(n1+n2,size=n1,replace=0)
    i2 = np.delete(np.arange(n1+n2),i1)
    
    # Subset the total population
    mcounts1 = counts[i1]
    mcounts2 = counts[i2]
    
    # Recount storms from MC population 1
    sum1 = np.apply_along_axis(np.sum,0,mcounts1[:,aoi[0],aoi[1]])
    sum2 = np.apply_along_axis(np.sum,0,mcounts2[:,aoi[0],aoi[1]])
    
    sums = np.vstack((sum1,sum2))
    
    # Calculate Dissimilarity of MC Data
    chiMC[mc] = chidist(sums)[1,0]

# Write output to file, making the observations the first row of the file
output = pd.DataFrame(data=np.hstack((np.array(chi),chiMC)),columns=["chi"])
output.to_csv(inpath+"/"+name+"_"+V+"/MonteCarloRelativeTrackDen"+str(td_thresh)+"_Tracks_"+T1+"v"+T2+"_"+SDAY+"_"+EDAY+".csv",index=False)
print(names[v]+": " + T1 + " v. " + T2)
print("Observed Chi: " + str(round(chi,2)))
print("Percentile: " + str(np.sum(chiMC <= chi)/10))
print("Max Monte Carlo Chi: " + str(np.round(np.max(chiMC),2)))
print("TrkDen (" + str(td_thresh) + "), AOI: " + str(aoi[0].shape[0]))