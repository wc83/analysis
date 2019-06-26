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
from scipy import integrate
from obspy.clients.earthworm import Client
from obspy import UTCDateTime
from obspy.signal.trigger import trigger_onset
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
t1 = UTCDateTime(1420044444)-10 #the format is year:day_of_the_year:month
t2 = t1 + 60
st = Stream()
st = client.get_waveforms(net, sta, '', cha, t1 , t2)
#print(st)

#2014-11-25 20:20:16 (UTC)

st.detrend(type='linear')
st.detrend(type='demean')
st.filter(type='bandpass',freqmin=0.2, freqmax=10)
#st.plot(color='b',starttime=t1, endtime=t2)

tr=st[0].slice(starttime=t1, endtime=t2)
#
st_c = calibrate1(tr)

st_c.plot(color='b',starttime=t1, endtime=t2)



#%%
A=1
rhoE=2500
cE=2000
pi=3.14159
r1=4630

r=r1
B=2*pi*rhoE*cE*(1/A)



dl=len(st_c[0].data)
p = np.linspace(0,dl/100, num=dl)
y= st_c[0].data
y2= np.square(st_c[0].data)

y_int = integrate.cumtrapz(y, p, initial=0)
y_int2 = integrate.cumtrapz(y2, p, initial=0)
y_int3 = B*(r*r)*integrate.cumtrapz(y2, p, initial=0)

EI = y_int2[-1]
EE= B*(r*r)*EI

E25=EE*0.025
E975=EE*0.975


#plt.figure()
#plt.plot(p, y, 'k-')
#plt.figure()
#plt.plot(p, y_int, 'r-')
#plt.figure()
#plt.plot(p, y_int2, 'b-')
plt.figure()
plt.plot(p, y_int3, 'r-')
plt.plot([0,60],[E25,E25],'k--')
plt.plot([0,60],[E975,E975],'k--')
plt.xlabel('Time [s]')
plt.ylabel('Cumulative Energy [s]')


# find closest point to 5% and 95% - get timems from them, then get duration 
n25, ind25 = find_nearest(y_int3, E25)
n975, ind975 = find_nearest(y_int3, E975)
duration = (ind975/100) - (ind25/100)
print('duration of event =', duration)

plt.plot([ind25/100,ind25/100],[0,EE],'b--')
plt.plot([ind975/100,ind975/100],[0,EE],'b--')

plt.figure()
plt.plot(p, y, 'k-')
plt.plot([ind25/100,ind25/100],[-0.000008,0.000008],'b--')
plt.plot([ind975/100,ind975/100],[-0.000008,0.000008],'b--')

#startt=UTCDateTime(st_c[0].stats.starttime.timestamp)
#st_c.plot(color='r',starttime=startt+(ind5/100)-5, endtime=t2)

#%%

#from obspy.signal.trigger import classic_sta_lta, recursive_sta_lta
#from obspy.signal.trigger import plot_trigger, trigger_onset
#
#sr = 100
#nsta=int(1*sr)                                      
#nlta=int(20*sr)                                   
#trig_on=5
#trig_off=0.2 
#
##st_c.plot(color='b',starttime=t1, endtime=t2)
#
#cft=recursive_sta_lta(tr, nsta, nlta)
#on_off = trigger_onset(cft,trig_on,trig_off)
#plot_trigger(tr, cft, trig_on, trig_off) 
#
#print(t1-20+on_off[0,0]/100)


#%%




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














    