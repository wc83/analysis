#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 09:56:29 2018

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


#%% import data

cat= genfromtxt("/Users/william/Documents/scanner/all_stations/Explosion_Catalogue_V3.csv", delimiter=',',skip_header=1)
#%% Events 
start_day = 1450785600
d2n = 1450785600+(14*60*60)
end_night = 1450785600+(24*60*60)
day_count=0
night_count=0
day_n=1
for x in range(100,14000):
    tim = UTCDateTime(cat[x,0])
#    print(tim)
    if 1450785600 < tim < 1450785600+(25*24*60*60) :
        if start_day < tim < end_night:
            if start_day < tim < d2n:
                day_count += 1
            if d2n < tim < end_night:
                night_count += 1
        if tim > end_night:
            start_day += 24*60*60
            d2n += 24*60*60
            end_night += 24*60*60
            day_n += 1
            if start_day < tim < d2n:
                day_count += 1
            if d2n < tim < end_night:
                night_count += 1
            if tim > end_night:
                print('more than one day gap')
            

print('Number of days sampled = ',day_n)
#print(day_count)
print('Events per hour in day =',day_count/(day_n*14))
#print(night_count)
print('Events per hour in night =',night_count/(day_n*10))
#    print(tim)