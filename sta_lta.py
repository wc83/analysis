#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 15:02:54 2018

@author: william
"""

from obspy.core import read
from obspy.signal.cross_correlation import correlate
import glob
import time
import obspy
import numpy as np
from numpy import argmax
from obspy.signal.trigger import classic_sta_lta, recursive_sta_lta
from obspy.signal.trigger import plot_trigger, trigger_onset

   #%% constants    
shift=100
step=2
#s_window=7*60*60 
#window=40
fmin=0.1
fmax=10
#%% Read in waveforms      
sample = read('/Users/william/Documents/lb01/14_330z.mseed')
trace=sample[0]
trace.detrend(type='demean')
trace.filter("bandpass", freqmin=3,freqmax=10)

#window endpoints
start= sample[0].stats.starttime +18*60*60 +20*60  #time window start 
end= start + 3*60
#end=sample[0].stats.endtime 
trs = trace.slice(starttime = start  , endtime= end) #cut out sample waveform with same window length as chosen event
trs.filter("bandpass", freqmin=fmin,freqmax=fmax)
#trs_e = obspy.signal.filter.envelope(trs.data)
#print('reference waveform')


trs.plot(type='relative',color='b')#, starttime=start , endtime=end)

sr = trace.stats.sampling_rate
nsta=int(1*sr)                                      #2
nlta=int(20*sr)                                     #20
stream=trs.data
cft=recursive_sta_lta(stream, nsta, nlta)
trig_on=5                                            #8
trig_off=1                                        #0.2
plot_trigger(trs, cft, trig_on, trig_off) 

on_off = trigger_onset(cft,trig_on,trig_off)

#for x in range(0,len(on_off)):
#    tr = trace.slice(starttime=start+(on_off[x,0]/sr)-10 , endtime=start+(on_off[x,1]/sr)+10)
#    tr.filter("bandpass", freqmin=fmin,freqmax=fmax)
#    tr_e = obspy.signal.filter.envelope(tr.data)
#    
#    window = tr.stats.endtime - tr.stats.starttime
#
#    peak,cf,bwid,bwid25 = freq_info(tr.data,tr.stats.starttime,tr.stats.endtime)        
#    
##%% plot if not EXP  
#    crite=0.6
#    
#    trc2_e1,trc2_e2,trc2_e3,trc2_e4,trc2_e5,trc2_e6=stacked_LBz_exp(fmin,fmax)
#    top_v_eq,top_eq,corell_eq = corel(trc2_e2,tr_e,shift)
##    if abs(top_v) < crite: #only allow waveforms that are not EXPs
#    if (0.75 < cf < 2.75) and (0.2 < peak < 2.5) and (0.3 < bwid < 2.5):
#        if abs(top_v_eq) > crite:
##        print(cf)
#            trs.plot(type='relative',color='b', starttime=start+(on_off[x,0]/sr)-10 , endtime=start+(on_off[x,1]/sr)+10)









# DO NOT DELETE
