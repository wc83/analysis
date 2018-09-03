#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  3 15:23:28 2018

@author: william
"""
import obspy
from obspy import read
import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
from obspy import Stream
from obspy.signal.trigger import classic_sta_lta, recursive_sta_lta
from obspy.signal.trigger import plot_trigger, trigger_onset
from obspy import UTCDateTime
from matplotlib import pyplot
import scipy
from obspy import UTCDateTime

cat= genfromtxt("/Users/william/Documents/scanner/all_stations/Explosion_catalogue_v2.csv", delimiter=',',skip_header=1)

pands=np.zeros(shape=(1,3))
pands[0,0]=cat[0,0]
pands[0,1]=1
pands[0,2]=0
pands = np.lib.pad(pands, ((0,1),(0,0)), 'constant', constant_values=(0))

last_p = cat[0,0]

p_tot=1
s_tot=0

for x in range(1,len(cat)):
    pands[x,0]=cat[x,0]
    
    if cat[x,0] - last_p < 10*60:
        pands[x,1]=2
        s_tot += 1
    else:
        pands[x,1]=1
        last_p=cat[x,0]
        p_tot +=1
    pands[x,2]=(cat[x,0]-cat[x-1,0] )/60 
    if x < len(cat)-1:   
        pands = np.lib.pad(pands, ((0,1),(0,0)), 'constant', constant_values=(0))
        
        
#np.savetxt("/Users/william/Documents/scanner/all_stations/primary_secondary_v1.csv", pands,delimiter=",",header="time_stanp,p_or_s,repose time")


print('number of Events =',p_tot)
print("Number of secondary EXP's =",s_tot)