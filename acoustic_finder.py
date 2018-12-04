#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  3 16:05:20 2018

@author: william
"""


import os
from collections import OrderedDict
import numpy as np
import obspy
import scipy.signal as sgn
import matplotlib.pyplot as plt 
import matplotlib.mlab as mlab

from obspy.clients.earthworm import Client
from obspy import UTCDateTime
from obspy.signal.trigger import trigger_onset

from obspy.clients.earthworm import Client
from obspy import UTCDateTime
#from scipy.signal import welch
from obspy import Stream


cat= genfromtxt("/Users/william/Documents/scanner/all_stations/Explosion_Catalogue_V3.csv", delimiter=',',skip_header=1)


###########  ADD IN CALIBRATIONS ##############

#LB01sc1=0.000001/750            # before 2015-12-05T00:00:01.000000Z  
#
#LB01ac = 0.000001/0.0250
#LB03ac = 0.000001/(0.0250*256)


#for x in range(0,len(cat)) : 
for x in range(0,10) :
    try:
        # STATION, CHANNEL (DDF --> 400 Hz), NETWWORK AND LOCATION CODES 
        sta = 'LB01' # STATION 
        cha = 'HHZ' # CHANNEL
        net = 'Z4'  # 
        loc = ''    # location, it depends mostly of which network you are in. 
        
        # Corner frequency for high-pass filter
        #hp_corner = 0.05
        
        # t1. and t2 are in hours:minutes:seconds
        # Get data from (Liverpool Winston default) wave server between times t1 and t2 for all stations in stalist      
        client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
        t1 = UTCDateTime(cat[x,0]-30) #the format is year:day_of_the_year:month
        t2 = t1 + 1.5*60
        st = Stream()
        st = client.get_waveforms(net, sta, '', cha, t1 , t2)
        
        
        st.filter(type='bandpass',freqmin=0.3, freqmax=15)
        st.detrend(type='linear')
        st.detrend(type='demean')
#        st[0].data=st[0].data*LB01sc1
        st.plot(color='r',starttime=t1+20, endtime=t2)
        
        #peak,cf,bwid50,bwid25=freq_info(st[0].data,t1,t2)
        #print(peak,cf,bwid50)
        
    except:
        print("LB01 Seismic not found")
    
    
    
    try:
        # STATION, CHANNEL (DDF --> 400 Hz), NETWWORK AND LOCATION CODES 
        sta = 'LB01' # STATION 
        cha = 'HDF' # CHANNEL
        net = 'Z4'  # 
        loc = ''    # location, it depends mostly of which network you are in. 
        
        # Corner frequency for high-pass filter
        #hp_corner = 0.05
        
        # t1. and t2 are in hours:minutes:seconds
        # Get data from (Liverpool Winston default) wave server between times t1 and t2 for all stations in stalist      
        client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
        t1 = UTCDateTime(cat[x,0]-30) #the format is year:day_of_the_year:month
        t2 = t1 + 1.5*60
        st = Stream()
        st = client.get_waveforms(net, sta, '', cha, t1 , t2)
        
        
        st.filter(type='bandpass',freqmin=0.3, freqmax=15)
        st.detrend(type='linear')
        st.detrend(type='demean')
#        st[0].data=st[0].data*LB01ac
        st.plot(color='b',starttime=t1+20, endtime=t2)
        
        #peak,cf,bwid50,bwid25=freq_info(st[0].data,t1,t2)
        #print(peak,cf,bwid50)
        
    except:
        print("LB01 Acoustic not found")

    try:
        # STATION, CHANNEL (DDF --> 400 Hz), NETWWORK AND LOCATION CODES 
        sta = 'LB03' # STATION 
        cha = 'HDF' # CHANNEL
        net = 'Z4'  # 
        loc = ''    # location, it depends mostly of which network you are in. 
        
        # Corner frequency for high-pass filter
        #hp_corner = 0.05
        
        # t1. and t2 are in hours:minutes:seconds
        # Get data from (Liverpool Winston default) wave server between times t1 and t2 for all stations in stalist      
        client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
        t1 = UTCDateTime(cat[x,0]-30) #the format is year:day_of_the_year:month
        t2 = t1 + 1.5*60
        st = Stream()
        st = client.get_waveforms(net, sta, '', cha, t1 , t2)
        
        
        st.filter(type='bandpass',freqmin=0.3, freqmax=15)
        st.detrend(type='linear')
        st.detrend(type='demean')
#        st[0].data=st[0].data*LB03ac
        st.plot(color='g',starttime=t1+20, endtime=t2)
        
        #peak,cf,bwid50,bwid25=freq_info(st[0].data,t1,t2)
        #print(peak,cf,bwid50)
        
    except:
        print("LB03 Acoustic not found")

