#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 10:35:36 2018

@author: william
"""


import obspy
from obspy import read
import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
from obspy import Stream
from obspy import UTCDateTime
from obspy.clients.earthworm import Client


cha = 'HHZ' # CHANNEL
net = 'Z4'  # 
loc = ''    # location, it depends mostly of which network you are in. 
client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23

def second(sta,last_p,cat1):
    
    
    arrival1=last_p
    arrival2=cat1
    print("arrival1",arrival1)
    
    t1 = UTCDateTime(arrival1)-30 #the format is year:day_of_the_year:month
    t2 = t1 + 120
    st1= Stream()
    st1 = client.get_waveforms(net, sta, '', cha, t1 , t2)
    st1.detrend(type='linear')
    st1.detrend(type='demean')
    peak1=max(st1[0].data)
    

    print("arrival2")
    t1 = UTCDateTime(arrival2)-30 #the format is year:day_of_the_year:month
    t2 = t1 + 120
    st2 = Stream()
    st2 = client.get_waveforms(net, sta, '', cha, t1 , t2)
    st2.detrend(type='linear')
    st2.detrend(type='demean')
    peak2=max(st2[0].data)
    
    
    print(st1,st2)
    

    

    if peak1 > peak2:

        sec=1

    
#        arrival=last_p
#        t1 = UTCDateTime(arrival)-60 #the format is year:day_of_the_year:month
#        t2 = t1 + (cat1 - last_p) + 120
#        st = Stream()
#        st = client.get_waveforms(net, sta, '', cha, t1 , t2)
#        st.detrend(type='linear')
#        st.detrend(type='demean')
#        st.filter(type='bandpass',freqmin=0.5, freqmax=6)
#        st.plot(color='b',starttime=t1, endtime=t2)
        print(cat1 - last_p)
        
    else: 
        sec=0
    

    
    
    return(sec,st1,st2)