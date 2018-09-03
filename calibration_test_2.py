#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 14:41:50 2018

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

#%%    Seismic calibrations


LB01c1=0.000001/750            # before 2015-12-05T00:00:01.000000Z
LB01c2=(10*0.000001)/(750*256) # 2015-12-05T00:00:01 - 2016-06-15T00:00:01 + maybe till end
LB01c3= 0.000001/750  

LB02c=0.000001/750
LB03c=0.000001/(750*256)
LB04c=0.000001/(750*256)
LB05c=0.000001/(750*256)
LB06c=(10*0.000001)/750
#
#LSsc=0.000000122/800

# Alternative factors
LB01c1=0.000001/750            # before 2015-12-05T00:00:01.000000Z
LB01c2=(10*0.0000041)/(750*256) # 2015-12-05T00:00:01 - 2016-06-15T00:00:01 + maybe till end
LB01c3 = 0.000001/750 

LB02c=0.0000000268/(750*256)
LB03c=0.0000000268/(750*256)
LB04c=0.0000000268/(750*256)
LB05c=0.0000000268/(750*256)
LB06c=(10*0.0000041)/750
#%% first station
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
t1 = UTCDateTime(2014, 11 ,29, 7, 56 , 0) #the format is year:day_of_the_year:month
t2 = t1 + 3*60
st = Stream()
st = client.get_waveforms(net, sta, '', cha, t1 , t2)

st.detrend(type='linear')
st.detrend(type='demean')
#st.plot(color='r',starttime=t1, endtime=t2)

st[0].data=st[0].data*LB01c1

st.filter("bandpass",freqmin=0.8, freqmax=20)
st.detrend(type='linear')
st.detrend(type='demean')
st.plot(color='b',starttime=t1, endtime=t2)
#
#
#
#%% second station
# STATION, CHANNEL (DDF --> 400 Hz), NETWWORK AND LOCATION CODES 
sta = 'LB02' # STATION 
cha = 'HHZ' # CHANNEL
net = 'Z4'  # 
loc = ''    # location, it depends mostly of which network you are in. 


# Get data from (Liverpool Winston default) wave server between times t1 and t2 for all stations in stalist      
client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
st2 = Stream()
st2 = client.get_waveforms(net, sta, '', cha, t1 , t2)

st2.detrend(type='linear')
st2.detrend(type='demean')
#st2.plot(color='r',starttime=t1, endtime=t2)

st2[0].data=st2[0].data*LB02c

st2.filter("bandpass",freqmin=0.8, freqmax=20)
st2.detrend(type='linear')
st2.detrend(type='demean')
st2.plot(color='g',starttime=t1, endtime=t2)

#%% third station
# STATION, CHANNEL (DDF --> 400 Hz), NETWWORK AND LOCATION CODES 
sta = 'LB03' # STATION 
cha = 'HHZ' # CHANNEL
net = 'Z4'  # 
loc = ''    # location, it depends mostly of which network you are in. 


# Get data from (Liverpool Winston default) wave server between times t1 and t2 for all stations in stalist      
client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
st3 = Stream()
st3 = client.get_waveforms(net, sta, '', cha, t1 , t2)

st3.detrend(type='linear')
st3.detrend(type='demean')
#st3.plot(color='r',starttime=t1, endtime=t2)

st3[0].data=st3[0].data*LB03c

#st.filter("bandpass",freqmin=0.01, freqmax=50)
st3.detrend(type='linear')
st3.detrend(type='demean')
st3.plot(color='r',starttime=t1, endtime=t2)

##%% fourth station
## STATION, CHANNEL (DDF --> 400 Hz), NETWWORK AND LOCATION CODES 
#sta = 'LB04' # STATION 
#cha = 'HHZ' # CHANNEL
#net = 'Z4'  # 
#loc = ''    # location, it depends mostly of which network you are in. 
#
## Corner frequency for high-pass filter
##hp_corner = 0.05
#
## t1. and t2 are in hours:minutes:seconds
## Get data from (Liverpool Winston default) wave server between times t1 and t2 for all stations in stalist      
#client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
#t1 = UTCDateTime(2014, 11 ,30, 11, 16 , 0) #the format is year:day_of_the_year:month
#t2 = t1 + 2*60
#st3 = Stream()
#st3 = client.get_waveforms(net, sta, '', cha, t1 , t2)
#
#st3.detrend(type='linear')
#st3.detrend(type='demean')
##st3.plot(color='r',starttime=t1, endtime=t2)
#
#st3[0].data=st3[0].data*LB04c
#
##st.filter("bandpass",freqmin=0.01, freqmax=50)
#st3.detrend(type='linear')
#st3.detrend(type='demean')
#st3.plot(color='b',starttime=t1, endtime=t2)
#
##%% fifth station
## STATION, CHANNEL (DDF --> 400 Hz), NETWWORK AND LOCATION CODES 
#sta = 'LB05' # STATION 
#cha = 'HHZ' # CHANNEL
#net = 'Z4'  # 
#loc = ''    # location, it depends mostly of which network you are in. 
#
## Corner frequency for high-pass filter
##hp_corner = 0.05
#
## t1. and t2 are in hours:minutes:seconds
## Get data from (Liverpool Winston default) wave server between times t1 and t2 for all stations in stalist      
#client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
#t1 = UTCDateTime(2014, 11 ,30, 11, 16 , 0) #the format is year:day_of_the_year:month
#t2 = t1 + 2*60
#st3 = Stream()
#st3 = client.get_waveforms(net, sta, '', cha, t1 , t2)
#
#st3.detrend(type='linear')
#st3.detrend(type='demean')
##st3.plot(color='r',starttime=t1, endtime=t2)
#
#st3[0].data=st3[0].data*LB05c
#
##st.filter("bandpass",freqmin=0.01, freqmax=50)
#st3.detrend(type='linear')
#st3.detrend(type='demean')
#st3.plot(color='b',starttime=t1, endtime=t2)
#


##%% sixth station
## STATION, CHANNEL (DDF --> 400 Hz), NETWWORK AND LOCATION CODES 
#sta = 'LB06' # STATION 
#cha = 'HHZ' # CHANNEL
#net = 'Z4'  # 
#loc = ''    # location, it depends mostly of which network you are in. 
#
## Corner frequency for high-pass filter
##hp_corner = 0.05
#
## t1. and t2 are in hours:minutes:seconds
## Get data from (Liverpool Winston default) wave server between times t1 and t2 for all stations in stalist      
#client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
#t1 = UTCDateTime(2016, 10 ,27, 4, 24 , 0) #the format is year:day_of_the_year:month
#t2 = t1 + 2*60
#st3 = Stream()
#st3 = client.get_waveforms(net, sta, '', cha, t1 , t2)
#
#st3.detrend(type='linear')
#st3.detrend(type='demean')
##st3.plot(color='r',starttime=t1, endtime=t2)
#
#st3[0].data=st3[0].data*LB06c
#
##st.filter("bandpass",freqmin=0.01, freqmax=50)
#st3.detrend(type='linear')
#st3.detrend(type='demean')
#st3.plot(color='b',starttime=t1, endtime=t2)
#









    