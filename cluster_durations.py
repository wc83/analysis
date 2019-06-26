#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 11:38:50 2019

@author: william
"""


import numpy as np
import obspy
import matplotlib.pyplot as plt 
from obspy.clients.earthworm import Client
from obspy import UTCDateTime
from obspy.signal.trigger import classic_sta_lta, recursive_sta_lta
from obspy.signal.trigger import plot_trigger, trigger_onset
from obspy import Stream
from numpy import genfromtxt

cat = genfromtxt("/Users/william/Documents/scanner/all_stations/Final_Catalogue_2014_2018.csv", delimiter=',',skip_header=1)
eventt = cat[:,0]

c1 = genfromtxt("/Users/william/Documents/scanner/analysis/LB01_cluster1.csv", delimiter=',',skip_header=1)
c2 = genfromtxt("/Users/william/Documents/scanner/analysis/LB01_cluster2.csv", delimiter=',',skip_header=1)
c3 = genfromtxt("/Users/william/Documents/scanner/analysis/LB01_cluster3.csv", delimiter=',',skip_header=1)
c4 = genfromtxt("/Users/william/Documents/scanner/analysis/LB01_cluster4.csv", delimiter=',',skip_header=1)
c5 = genfromtxt("/Users/william/Documents/scanner/analysis/LB01_cluster5.csv", delimiter=',',skip_header=1)
c6 = genfromtxt("/Users/william/Documents/scanner/analysis/LB01_cluster6.csv", delimiter=',',skip_header=1)
c7 = genfromtxt("/Users/william/Documents/scanner/analysis/LB01_cluster7.csv", delimiter=',',skip_header=1)
c8 = genfromtxt("/Users/william/Documents/scanner/analysis/LB01_cluster8.csv", delimiter=',',skip_header=1)
c9 = genfromtxt("/Users/william/Documents/scanner/analysis/LB01_cluster9.csv", delimiter=',',skip_header=1)
c10 = genfromtxt("/Users/william/Documents/scanner/analysis/LB01_cluster10.csv", delimiter=',',skip_header=1)

#%% c1


c1_d = np.zeros(shape=(0,1))
num1 = 0

for x in range(0,len(c1)):
    wt = UTCDateTime(int(c1[x,0]),int(c1[x,1]),int(c1[x,2]),int(c1[x,3]),int(c1[x,4]),int(c1[x,5])).timestamp

    evt,idx = find_nearest(eventt,wt)
    
    if abs(evt - wt) < 20:

        c1_d = np.lib.pad(c1_d, ((0,1),(0,0)), 'constant', constant_values=(0))
        c1_d[num1][0] = cat[idx,24]
        num1+=1


c1_maxd = max(c1_d)
c1_mind = min(c1_d)
c1_medd = np.median(c1_d)
c1_meand = np.mean(c1_d)

print('min =', c1_mind, 'max =', c1_maxd, 'mean =', c1_meand, 'median =', c1_medd)


#%% c2


c2_d = np.zeros(shape=(0,1))
num2 = 0

for x in range(0,len(c2)):
    wt = UTCDateTime(int(c2[x,0]),int(c2[x,1]),int(c2[x,2]),int(c2[x,3]),int(c2[x,4]),int(c2[x,5])).timestamp
    
    evt,idx = find_nearest(eventt,wt)
    
    if abs(evt - wt) < 20:

        c2_d = np.lib.pad(c2_d, ((0,1),(0,0)), 'constant', constant_values=(0))
        c2_d[num2][0] = cat[idx,24]
        num2+=1


c2_maxd = max(c2_d)
c2_mind = min(c2_d)
c2_medd = np.median(c2_d)
c2_meand = np.mean(c2_d)

print('min =', c2_mind, 'max =', c2_maxd, 'mean =', c2_meand, 'median =', c2_medd)


#%% c3


c3_d = np.zeros(shape=(0,1))
num3 = 0

for x in range(0,len(c3)):
    wt = UTCDateTime(int(c3[x,0]),int(c3[x,1]),int(c3[x,2]),int(c3[x,3]),int(c3[x,4]),int(c3[x,5])).timestamp

    evt,idx = find_nearest(eventt,wt)
    
    if abs(evt - wt) < 20:

        c3_d = np.lib.pad(c3_d, ((0,1),(0,0)), 'constant', constant_values=(0))
        c3_d[num3][0] = cat[idx,24]
        num3+=1


c3_maxd = max(c3_d)
c3_mind = min(c3_d)
c3_medd = np.median(c3_d)
c3_meand = np.mean(c3_d)

print('min =', c3_mind, 'max =', c3_maxd, 'mean =', c3_meand, 'median =', c3_medd)

#%% c4


c4_d = np.zeros(shape=(0,1))
num4 = 0

for x in range(0,len(c4)):
    wt = UTCDateTime(int(c4[x,0]),int(c4[x,1]),int(c4[x,2]),int(c4[x,3]),int(c4[x,4]),int(c4[x,5])).timestamp

    evt,idx = find_nearest(eventt,wt)
    
    if abs(evt - wt) < 20:

        c4_d = np.lib.pad(c4_d, ((0,1),(0,0)), 'constant', constant_values=(0))
        c4_d[num4][0] = cat[idx,24]
        num4+=1


c4_maxd = max(c4_d)
c4_mind = min(c4_d)
c4_medd = np.median(c4_d)
c4_meand = np.mean(c4_d)

print('min =', c4_mind, 'max =', c4_maxd, 'mean =', c4_meand, 'median =', c4_medd)

#%% c5


c5_d = np.zeros(shape=(0,1))
num5 = 0

for x in range(0,len(c5)):
    wt = UTCDateTime(int(c5[x,0]),int(c5[x,1]),int(c5[x,2]),int(c5[x,3]),int(c5[x,4]),int(c5[x,5])).timestamp

    evt,idx = find_nearest(eventt,wt)
    
    if abs(evt - wt) < 20:

        c5_d = np.lib.pad(c5_d, ((0,1),(0,0)), 'constant', constant_values=(0))
        c5_d[num5][0] = cat[idx,24]
        num5+=1


c5_maxd = max(c5_d)
c5_mind = min(c5_d)
c5_medd = np.median(c5_d)
c5_meand = np.mean(c5_d)

print('min =', c5_mind, 'max =', c5_maxd, 'mean =', c5_meand, 'median =', c5_medd)

#%% c6


c6_d = np.zeros(shape=(0,1))
num6 = 0

for x in range(0,len(c6)):
    wt = UTCDateTime(int(c6[x,0]),int(c6[x,1]),int(c6[x,2]),int(c6[x,3]),int(c6[x,4]),int(c6[x,5])).timestamp

    evt,idx = find_nearest(eventt,wt)
    
    if abs(evt - wt) < 20:

        c6_d = np.lib.pad(c6_d, ((0,1),(0,0)), 'constant', constant_values=(0))
        c6_d[num6][0] = cat[idx,24]
        num6+=1


c6_maxd = max(c6_d)
c6_mind = min(c6_d)
c6_medd = np.median(c6_d)
c6_meand = np.mean(c6_d)

print('min =', c6_mind, 'max =', c6_maxd, 'mean =', c6_meand, 'median =', c6_medd)

#%% c7


c7_d = np.zeros(shape=(0,1))
num7 = 0

for x in range(0,len(c7)):
    wt = UTCDateTime(int(c7[x,0]),int(c7[x,1]),int(c7[x,2]),int(c7[x,3]),int(c7[x,4]),int(c7[x,5])).timestamp

    evt,idx = find_nearest(eventt,wt)
    
    if abs(evt - wt) < 20:

        c7_d = np.lib.pad(c7_d, ((0,1),(0,0)), 'constant', constant_values=(0))
        c7_d[num7][0] = cat[idx,24]
        num7+=1


c7_maxd = max(c7_d)
c7_mind = min(c7_d)
c7_medd = np.median(c7_d)
c7_meand = np.mean(c7_d)

print('min =', c7_mind, 'max =', c7_maxd, 'mean =', c7_meand, 'median =', c7_medd)

#%% c8


c8_d = np.zeros(shape=(0,1))
num8 = 0

for x in range(0,len(c8)):
    wt = UTCDateTime(int(c8[x,0]),int(c8[x,1]),int(c8[x,2]),int(c8[x,3]),int(c8[x,4]),int(c8[x,5])).timestamp

    evt,idx = find_nearest(eventt,wt)
    
    if abs(evt - wt) < 20:

        c8_d = np.lib.pad(c8_d, ((0,1),(0,0)), 'constant', constant_values=(0))
        c8_d[num8][0] = cat[idx,24]
        num8+=1


c8_maxd = max(c8_d)
c8_mind = min(c8_d)
c8_medd = np.median(c8_d)
c8_meand = np.mean(c8_d)



print('min =', c8_mind, 'max =', c8_maxd, 'mean =', c8_meand, 'median =', c8_medd)


#%% c9


c9_d = np.zeros(shape=(0,1))
num9 = 0

for x in range(0,len(c9)):
    wt = UTCDateTime(int(c9[x,0]),int(c9[x,1]),int(c9[x,2]),int(c9[x,3]),int(c9[x,4]),int(c9[x,5])).timestamp

    evt,idx = find_nearest(eventt,wt)
    
    if abs(evt - wt) < 20:

        c9_d = np.lib.pad(c9_d, ((0,1),(0,0)), 'constant', constant_values=(0))
        c9_d[num9][0] = cat[idx,24]
        num9+=1


c9_maxd = max(c9_d)
c9_mind = min(c9_d)
c9_medd = np.median(c9_d)
c9_meand = np.mean(c9_d)


print('min =', c9_mind, 'max =', c9_maxd, 'mean =', c9_meand, 'median =', c9_medd)



#%% c10


c10_d = np.zeros(shape=(0,1))
num10 = 0

for x in range(0,len(c10)):
    wt = UTCDateTime(int(c10[x,0]),int(c10[x,1]),int(c10[x,2]),int(c10[x,3]),int(c10[x,4]),int(c10[x,5])).timestamp

    evt,idx = find_nearest(eventt,wt)
    
    if abs(evt - wt) < 20:

        c10_d = np.lib.pad(c10_d, ((0,1),(0,0)), 'constant', constant_values=(0))
        c10_d[num10][0] = cat[idx,24]
        num10+=1


c10_maxd = max(c10_d)
c10_mind = min(c10_d)
c10_medd = np.median(c10_d)
c10_meand = np.mean(c10_d)

print('min =', c10_mind, 'max =', c10_maxd, 'mean =', c10_meand, 'median =', c10_medd)





