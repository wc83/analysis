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
sta1 = 'LB01' # STATION 
sta2 = 'LB02' # STATION 
sta3 = 'LB03' # STATION 
sta4 = 'LB04' # STATION 
sta5 = 'LB05' # STATION 
sta6 = 'LB06'
cha = 'HHZ' # CHANNEL
net = 'Z4'  # 
loc = ''    # location, it depends mostly of which network you are in. 

client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
t1 = UTCDateTime(2016, 6, 25, 17, 7, 30) #the format is year:day_of_the_year:month
t2 = t1 + 2*60 #+0.01
st1 = Stream()
st1 = client.get_waveforms(net, sta1, '', cha, t1-20 , t2)

#st1.filter(type='bandpass',freqmin=0.1, freqmax=10)
st1.detrend(type='linear')
st1.detrend(type='demean')
#st.plot(color='b',starttime=t1, endtime=t2)

tr1=st1[0].slice(starttime=t1, endtime=t2)
st_c1 = calibrate(tr1)
trace1=st_c1[0].slice(starttime=t1+20, endtime=t2)
trace1.detrend(type='linear')
trace1.detrend(type='demean')
#st_c1.plot(color='b',starttime=t1+20, endtime=t2)



#st2 = Stream()
#st2 = client.get_waveforms(net, sta2, '', cha, t1-20 , t2)
#
##st2.filter(type='bandpass',freqmin=0.1, freqmax=10)
#st2.detrend(type='linear')
#st2.detrend(type='demean')
##st.plot(color='b',starttime=t1, endtime=t2)
#
#tr2=st2[0].slice(starttime=t1, endtime=t2)
#st_c2 = calibrate(tr2)
#trace2=st_c2[0].slice(starttime=t1+20, endtime=t2)
#trace2.detrend(type='linear')
#trace2.detrend(type='demean')
##st_c2.plot(color='b',starttime=t1+20, endtime=t2)



st3 = Stream()
st3 = client.get_waveforms(net, sta3, '', cha, t1-20 , t2)

#st3.filter(type='bandpass',freqmin=0.1, freqmax=10)
st3.detrend(type='linear')
st3.detrend(type='demean')
#st.plot(color='b',starttime=t1, endtime=t2)

tr3=st3[0].slice(starttime=t1, endtime=t2)
st_c3 = calibrate(tr3)
trace3=st_c3[0].slice(starttime=t1+20, endtime=t2)
trace3.detrend(type='linear')
trace3.detrend(type='demean')
#st_c3.plot(color='b',starttime=t1+20, endtime=t2)


#
#st4 = Stream()
#st4 = client.get_waveforms(net, sta4, '', cha, t1-20 , t2)
#
##st4.filter(type='bandpass',freqmin=0.1, freqmax=10)
#st4.detrend(type='linear')
#st4.detrend(type='demean')
##st.plot(color='b',starttime=t1, endtime=t2)
#
#tr4=st4[0].slice(starttime=t1, endtime=t2)
#st_c4 = calibrate(tr4)
#trace4=st_c4[0].slice(starttime=t1+20, endtime=t2)
#trace4.detrend(type='linear')
#trace4.detrend(type='demean')
##st_c4.plot(color='b',starttime=t1+20, endtime=t2)
#
#st5 = Stream()
#st5 = client.get_waveforms(net, sta5, '', cha, t1-20 , t2)
#
##st5.filter(type='bandpass',freqmin=0.1, freqmax=10)
#st5.detrend(type='linear')
#st5.detrend(type='demean')
##st.plot(color='b',starttime=t1, endtime=t2)
#
#tr5=st5[0].slice(starttime=t1, endtime=t2)
#st_c5 = calibrate(tr5)
#trace5=st_c5[0].slice(starttime=t1+20, endtime=t2)
#trace5.detrend(type='linear')
#trace5.detrend(type='demean')

st6 = Stream()
st6 = client.get_waveforms(net, sta6, '', cha, t1-20 , t2)

#st6.filter(type='bandpass',freqmin=0.1, freqmax=10)
st6.detrend(type='linear')
st6.detrend(type='demean')
#st.plot(color='b',starttime=t1, endtime=t2)

tr6=st6[0].slice(starttime=t1, endtime=t2)
st_c6 = calibrate(tr6)
trace6=st_c6[0].slice(starttime=t1+20, endtime=t2)
trace6.detrend(type='linear')
trace6.detrend(type='demean')

#%%

x = np.linspace(0,t2-t1-20,(t2-t1-20)*100)

fig, axes = plt.subplots(3,1)
fig.set_figheight(8)
plt.tight_layout()
ti = fig.suptitle("         Explosion Seismograms at "  + str(t1) , fontsize=13)
ti.set_y(0.97)
axes[0].plot(x,trace1.data*1000)
axes[1].plot(x,trace3.data*1000)
axes[2].plot(x,trace6.data*1000)
#axes[3].plot(x,trace4.data*1000)
#axes[4].plot(x,trace5.data*1000)
axes[0].set_xlabel('Time [s]')
axes[0].set_ylabel('LB01 ground velocity [mm/s]')
axes[1].set_xlabel('Time [s]')
axes[1].set_ylabel('LB03 ground velocity [mm/s]')
axes[2].set_xlabel('Time [s]')
axes[2].set_ylabel('LB06 ground velocity [mm/s]')
#axes[3].set_xlabel('Time [s]')
#axes[3].set_ylabel('LB04 ground velocity [mm/s]')
#axes[4].set_xlabel('Time [s]')
#axes[4].set_ylabel('LB05 ground velocity [m/s]')

fig.subplots_adjust(hspace = 0.3)
fig.subplots_adjust(top=0.92)
fig.subplots_adjust(left=0.175)
fig.subplots_adjust(bottom=0.075)
fig


fig.savefig('/Users/william/Documents/Figures_and_Tables/EXP_V3/EXP_stations_5.png')







