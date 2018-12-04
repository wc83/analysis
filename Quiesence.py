#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 09:10:37 2018

@author: william
"""

import obspy
from obspy import read
import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
from obspy import Stream
from obspy import UTCDateTime
from matplotlib import pyplot
from obspy import UTCDateTime
from obspy.clients.earthworm import Client

cat= genfromtxt("/Users/william/Documents/scanner/all_stations/Explosion_Catalogue_V4c.csv", delimiter=',',skip_header=1)

TS = cat[:,0]
day_no=cat[:,2]
year=cat[:,3]
month=cat[:,4]
m_day=cat[:,5]

sta = 'LB01' # STATION 
cha = 'HHZ' # CHANNEL
net = 'Z4'  # 
loc = ''    # location, it depends mostly of which network you are in. 

NL=np.zeros(shape=(0,2))
num=0
for x in range(100,160):
#    print('EVENT', x)
    if TS[x]-TS[x-1] > 3600:
        gap=TS[x]-TS[x-1]
        start = TS[x] + 600 - gap 
        
        while TS[x] - start > 300: 
        
            client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
            t1 = UTCDateTime(start) #the format is year:day_of_the_year:month
            st = Stream()
            st = client.get_waveforms(net, sta, '', cha, t1 , t1+300)
            st.detrend(type='linear')
            st.detrend(type='demean')
            
            noise = np.mean(abs(st[0].data))
                                
            NL = np.lib.pad(NL, ((0,1),(0,0)), 'constant', constant_values=(0))
            NL[num][0] = -1*(TS[x] - start)
            NL[num][1] = noise
            
            
#            print(np.mean(abs(st[0].data)))
            
            start += 60
            num +=1
    
    if len(NL) > 2:      
        plt.figure()
        plt.plot(NL[:,0],NL[:,1],'-')
        plt.ylim([0, max(NL[:,1])+50])
        plt.xlabel("Time Before Explosion [s]")
        plt.ylabel("Amplitude (5 min average) [counts]")
        plt.title("Noise Level Pre-Explosion")
    NL=np.zeros(shape=(0,2))
    num=0






