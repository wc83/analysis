#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 13 11:19:55 2018

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
t1 = UTCDateTime(2016, 6, 19, 10, 28, 0) #the format is year:day_of_the_year:month
t2 = t1 + 120
st = Stream()
st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
#print(st)

#st.filter(type='bandpass',freqmin=0.1, freqmax=10)
st.detrend(type='linear')
st.detrend(type='demean')
st.plot(color='b',starttime=t1, endtime=t2)

#tr=st[0].slice(starttime=t1, endtime=t2)
#
#st_c = calibrate(tr)
#
#st_c.plot(color='b',starttime=t1+20, endtime=t2)


#plt.figure(10)
#plt.plot(st[0].data)
#plt.xlim((1000, 35000))
#plt.ylim((-12000,12000))

#st.plot(type="dayplot", interval=60,
#        vertical_scaling_range=5e3, one_tick_per_line=True,
#        color=['k', 'r', 'b', 'g'],show_y_UTC_label=False, 
#        title='48 hour plot of LS02 from 08/01/2015',)



#peak,cf, bwid50, bwid25 = freq_info(tr,t1 ,t2)

#print('peak=',peak,'cf=',cf, 'bwid50=',bwid50)
#if 0.2<peak<2.75 and 0.75<cf<3 and 0.3<bwid50<2.5:
#    print('peak=',peak,'cf=',cf,'bwid=', bwid50)
#    print('freq ok')
#else:
#    print('wrong freq')
#    print('peak=',peak,'cf=',cf,'bwid=', bwid50)














    