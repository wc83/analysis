#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 16:10:15 2019

@author: william
"""
# THIS SEEMS TO ONLY APPEAR TO WORK DUE TO ANOMALIES CAUSED WHEN LOCATION IS NEXT TO A STATION. 
# REMOVAL OF LB04 SHOWS NO INDICATION OF THE LOW MISFIT SEEN WHEN IT IS USED IN ANALYSIS


import numpy as np
from numpy import genfromtxt
from obspy import UTCDateTime
import numpy as np
import matplotlib.pyplot as plt
import obspy
import io
from obspy import Stream
from obspy.clients.earthworm import Client

#%%


#Network distance to vent

# station = (d_N, d_E, d_Z, horizontal_distance) positive is north, east and up. all in metres
LB01 = (-4432, -112.344332, -1362, 4430)
LB02 = (-3198, -232.666434, -1100, 3210)
LB03 = (1287, -887.031446, -601, 2280)
LB04 = (486,	 212.22201, -1196, 530)
LB05 = (-46, 809.99919, -229, 810)
LB06 = (-7383, -1390.120832, -1665, 7510)
LS01 = (5327,	 2187.031146, -49, 5760)
LS02 = (2348,	 5122.2171, -300, 5630)
LS03 = (3262,	 -2229.664437, -24, 3950)
LS04 = (-1744, -5254.316968, -1256, 5540)
LS05 = (-1636, 4004.929328, -842, 4330)
LS06 = (1728, 1770.364896, 967, 2470)
Caliente = (0, 0, 0, 0)

#%% all stations and vent
station = np.zeros(shape=(13,4))
station[:][0]=LB01[:]
station[:][1]=LB02[:]
station[:][2]=LB03[:]
station[:][3]=LB04[:]
station[:][4]=LB05[:]
station[:][5]=LB06[:]
station[:][6]=LS01[:]
station[:][7]=LS02[:]
station[:][8]=LS03[:]
station[:][9]=LS04[:]
station[:][10]=LS05[:]
station[:][11]=LS06[:]
station[:][10]=Caliente[:]

#%% stations active for test event
a=11
stations = np.zeros(shape=(a,4))
stations[:][0]=LB01[:]
stations[:][1]=LB02[:]
stations[:][2]=LB03[:]
stations[:][3]=LB04[:]
stations[:][4]=LB05[:]
stations[:][5]=LS01[:]
stations[:][6]=LS02[:]
stations[:][7]=LS03[:]
stations[:][8]=LS04[:]
stations[:][9]=LS05[:]
stations[:][10]=LS06[:]


#%%
#constants   
# what units do these need to be in?!

f=2
Q=150
V=500
p=3.14159
B=(p*f)/(Q*V)        
        
n=0.5   #body waves, surface =0.5      
z=-500
#%% read in event, calculate true amplitude ratios
t1 = UTCDateTime(2014, 12, 2, 11, 4, 40)

magnitudes=event_mag(t1,a)

ratio_obs_grid = np.zeros(shape=(a,a))
for i in range(0,a):
    for j in range(i+1,a):
        ratio_obs=magnitudes[i]/magnitudes[j]
        ratio_obs_grid[i][j]=ratio_obs



#%%

# grid of 2D points (ignoring vertical differeces)


grid_pairs = np.zeros(shape=(641601,3))
count = 0
for x in range(-4000,4010,10):
    for y in range(-4000,4010,10):
        grid_pairs[count][0]=x
        grid_pairs[count][1]=y
        
        
        #calculate theroretical magnitudes ratio for all station pairs at all grid points
       
        #source-receiver distances for new source location
        srd = np.zeros(shape=(a,1))
        for d in range(0,a):
            srd[d]=np.sqrt(((x-stations[d][1])**2)+((y-stations[d][0])**2)+((z-stations[d][2])**2))
        # ratio of distances
        dist_rat_grid = np.zeros(shape=(a,a))
        for i in range(0,a):
            for j in range(0,a):
                dist_rat = (srd[j]/srd[i])
                dist_rat_grid[i][j] = dist_rat
        #difference of distances
        dist_diff_grid = np.zeros(shape=(a,a))
        for i in range(0,a):
            for j in range(0,a):
                dist_diff = (srd[i]-srd[j])
                dist_diff_grid[i][j] = dist_diff                
        #theoretical ratio of magnitudes
        Ratio_grid = np.zeros(shape=(a,a))
        for i in range(0,a):
            for j in range(i+1,a):
                Ratio = ((dist_rat_grid[i][j])**n)*np.exp(-B*(dist_diff_grid[i][j]))
                Ratio_grid[i][j]=Ratio
                
        #Calculate misfit for location ij
        sum_misf=0
        misf_ij_grid=np.zeros(shape=(a,a))
        for i in range(0,a):
            for j in range(i+1,a):
                misf_ij=(Ratio_grid[i][j]-ratio_obs_grid[i][j])**2
                sum_misf += misf_ij
        
        misfit=np.sqrt(sum_misf)
        
        if misfit > 50:
            misfit=50
            
        grid_pairs[count][2]=misfit        
        

#        
        count+=1



#%%      
min_misf = np.argmin(grid_pairs[:,2])              
location = grid_pairs[min_misf]
print("location of source = ", location[0],'metres East, and', location[1],'metres North of Caliente Vent') 
print("misfit =", location[2])  


#%%


sc = plt.scatter(x=grid_pairs[:,0], y=grid_pairs[:,1], c=grid_pairs[:,2], cmap="rainbow")
plt.colorbar(sc)

plt.plot(LB01[1],LB01[0],'kx')
plt.plot(LB02[1],LB02[0],'kx')
plt.plot(LB03[1],LB03[0],'kx')
plt.plot(LB04[1],LB04[0],'kx')
plt.plot(LB05[1],LB05[0],'kx')
##plt.plot(LB06[1],LB06[0],'kx')
plt.plot(LS01[1],LS01[0],'kx')
plt.plot(LS02[1],LS02[0],'kx')
plt.plot(LS03[1],LS03[0],'kx')
plt.plot(LS04[1],LS04[0],'kx')
plt.plot(LS05[1],LS05[0],'kx')
plt.plot(LS06[1],LS06[0],'kx')
plt.plot(Caliente[1],Caliente[0],'r^')       
plt.plot(location[0],location[1],'rx')        
        
plt.title("Explosion origin - Seismic Amplitude Ratio Method")
plt.xlabel("Distance East (m)")
plt.ylabel("Distance North (m)")        
    
        
       
        
        
        
        
        