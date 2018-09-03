#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  3 14:26:04 2018

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



#%% import data
stre = read("/Users/william/Documents/scanner/output_data/EXP_all_data_stream_month_1.mseed")
per_day = genfromtxt("/Users/william/Documents/scanner/all_stations/Explosions_per_day_month_1.csv", delimiter=',',skip_header=1,skip_footer=1)
data = genfromtxt("/Users/william/Documents/scanner/all_stations/EXP_all_coincidence_month_1.csv", delimiter=',',skip_header=1)


#%%
day_one = 1416787200.0

LB01c1=7.50E+08
LB01c2=3.07E+09
LB01c3=7.50E-04
LB02c=2.01E+11
LB03c=2.01E+11
LB04c=2.01E+11
LB05c=2.01E+11
LB06c=3.07E+09


# STATION, CHANNEL (DDF --> 400 Hz), NETWWORK AND LOCATION CODES 
sta = 'LB01' # STATION 
cha = 'HHZ' # CHANNEL
net = 'Z4'  # 
loc = ''    # location, it depends mostly of which network you are in. 

# t1. and t2 are in hours:minutes:seconds
# Get data from (Liverpool Winston default) wave server between times t1 and t2 for all stations in stalist      
client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
t1 = UTCDateTime(2014, 12 ,25, 6, 45, 22) #the format is year:day_of_the_year:month
t2 = t1 + 20*60
st = Stream()
st = client.get_waveforms(net, sta, '', cha, t1 , t2)


st.filter(type='bandpass',freqmin=1, freqmax=10)
st.detrend(type='linear')
st.detrend(type='demean')
st[0].data=st[0].data/LB01c1
st[0].plot(type='relative',color='b',starttime=t1+20, endtime=t2)

##%%
## STATION, CHANNEL (DDF --> 400 Hz), NETWWORK AND LOCATION CODES 
#sta = 'LB02' # STATION 
#cha = 'HHZ' # CHANNEL
#net = 'Z4'  # 
#loc = ''    # location, it depends mostly of which network you are in. 
#
## t1. and t2 are in hours:minutes:seconds
## Get data from (Liverpool Winston default) wave server between times t1 and t2 for all stations in stalist      
#client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
#t1 = UTCDateTime(2015, 2 ,17, 14, 35, 0) #the format is year:day_of_the_year:month
#t2 = t1 + 60*60
#st = Stream()
#st = client.get_waveforms(net, sta, '', cha, t1 , t2)
#
#
#st.filter(type='bandpass',freqmin=1, freqmax=10)
#st.detrend(type='linear')
#st.detrend(type='demean')
#st[0].data=st[0].data/LB02c
#st[0].plot(color='b',starttime=t1+20, endtime=t2)














