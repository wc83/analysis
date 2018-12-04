#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 16 10:49:37 2018

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


cat= genfromtxt("/Users/william/Documents/scanner/all_stations/Explosion_Catalogue_V4c.csv", delimiter=',',skip_header=1)
#cat= genfromtxt("/Users/william/Documents/scanner/all_stations/Earthquake_Catalogue_v1.csv", delimiter=',',skip_header=1)

day_no=cat[:,1]
year=cat[:,3]
#month=cat[:,4]
#m_day=cat[:,5]


epd=np.zeros(shape=(0,3))

exp_no = 0


for x in range(0,904):
    exp_pd = 0
    epd = np.lib.pad(epd, ((0,1),(0,0)), 'constant', constant_values=(0))
    epd[x][0] = x+1
    epd[x][2] = 1416787200 + x*24*60*60
    for e in range(exp_no,exp_no+100):
        if e < 2420:
            if day_no[e]==x+1:
                exp_pd += 1
                exp_no += 1
            
    
    epd[x][1]=exp_pd
    
    
plt.figure(10001)
plt.plot(epd[:,0],epd[:,1])
#plt.ylim([0,max(epd[:,1])+10])
plt.xlabel('Day Number')
plt.ylabel('Number of Explosions')
plt.title('Explosions Detected per Day')
#
##np.savetxt("/Users/william/Documents/scanner/output_data/earthquakes_per_day_v1.csv", epd ,delimiter=",",header="day,detections,time")
#
#epw=np.zeros(shape=(0,3))
#week=1
#event_c =0
#day_c =0
#
#for x in range(0,len(epd)):
#    day_c += 1
#    event_c += epd[x][1]
#    
#    if day_c == 7:
#        epw = np.lib.pad(epw, ((0,1),(0,0)), 'constant', constant_values=(0))
#        epw[week-1][0]=week
#        epw[week-1][1]=event_c
#        epw[week-1][2] = 1416787200 + (week-1)*7*24*60*60
#        week += 1
#        event_c = 0
#        day_c =0
#
#
#
#plt.figure(10002)
#plt.plot(epw[:,0],epw[:,1])
##plt.ylim([0,max(epw[:,1])+10])
#plt.xlabel('Week Number')
#plt.ylabel('Number of Explosions')
#plt.title('Explosions Detected per Week')

#np.savetxt("/Users/william/Documents/scanner/output_data/earthquakes_per_week_v1.csv", epw ,delimiter=",",header="week,detections,time")