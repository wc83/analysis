#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 21 17:30:18 2018

@author: william
"""

from obspy.clients.earthworm import Client
from obspy import UTCDateTime
from obspy.signal.trigger import trigger_onset

from obspy.clients.earthworm import Client
from obspy import UTCDateTime
#from scipy.signal import welch
from obspy import Stream

# STATION, CHANNEL (DDF --> 400 Hz), NETWWORK AND LOCATION CODES 
sta = 'LB03' # STATION 
cha = 'HHZ' # CHANNEL
net = 'Z4'  # 
loc = ''    # location, it depends mostly of which network you are in. 

# Corner frequency for high-pass filter
#hp_corner = 0.05

# t1. and t2 are in hours:minutes:seconds
# Get data from (Liverpool Winston default) wave server between times t1 and t2 for all stations in stalist      
client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
t1 = UTCDateTime(2014, 11, 24 ,5 ,40, 0) #the format is year:day_of_the_year:month
t2 = t1 + 3*60
st = Stream()
st = client.get_waveforms(net, sta, '', cha, t1 , t2)


st.filter(type='bandpass',freqmin=0.5, freqmax=15)
st.detrend(type='linear')
st.detrend(type='demean')
st.plot(color='b',starttime=t1+20, endtime=t2)

#peak,cf,bwid50,bwid25=freq_info(st[0].data,t1,t2)
#print(peak,cf,bwid50)

tr=st[0]
st_c = calibrate(tr)

st_c.plot(color='r',starttime=t1+20, endtime=t2)