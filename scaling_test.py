#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 10:29:26 2018

@author: william
"""


from numpy import argmax
import numpy as np
from obspy.clients.earthworm import Client
from obspy import UTCDateTime
from obspy.signal.trigger import trigger_onset
from obspy import Stream
from obspy import Trace
import matplotlib.pyplot as plt
from obspy import read, read_inventory
import io
import obspy 

reclen = 512
chunksize = 100000 * reclen # Around 50 MB

scale=np.zeros(shape=(0,9))

event_num = 0

A=1
rhoE=2500
cE=2000
pi=3.14159
r1=4630
r2=3370
r3=2310
r4=1300
r5=810
r6=7660
rs1=5710
rs2=5490
rs3=3900
rs4=5520
rs5=4290
rs6=2610

first=0

with io.open("/Users/william/Documents/scanner/output_data/m30.mseed", "rb") as fh:
        # just month 2
#with io.open("/Users/william/Documents/scanner/output_data/EXP_all_data_stream_2_month_1.mseed", "rb") as fh:
    while True:
        with io.BytesIO() as buf:
            c = fh.read(chunksize)
            if not c:
                break
            buf.write(c)
            buf.seek(0, 0)
            st = obspy.read(buf)
            
            for x in range(0,len(st)):
                if st[x].stats.station == "LB01":
                    if 1418580590 > st[x].stats.starttime.timestamp > 1418580580:
                        
                        tr2 = st[x]
                        t1 = st[x].stats.starttime
                        t2 = st[x].stats.endtime
                        window=t2-t1
                        tr2_data=st[x].data
                        m=np.mean(tr2_data)
                        tr2_data = tr2_data-m
                        
                        st_c2 = calibrate(tr2)
                        st_c2.plot(color='g')
                    
                        famp2 = abs(np.fft.rfft(tr2_data))
                        
                    if  1438145500 > st[x].stats.starttime.timestamp > 1438145480:
                        tr = st[x]
                        t1 = st[x].stats.starttime
                        t2 = st[x].stats.endtime
                        window=t2-t1
                        tr_data=st[x].data
                        m=np.mean(tr_data)
                        tr_data = tr_data-m
                        
                        st_c1 = calibrate(tr)
                        st_c1.plot(color='r')
                    
                        famp = abs(np.fft.rfft(tr_data))
                        
                        DC = famp2[2:-2]/famp[2:-2]
                        st_dc = np.fft.ifft(DC)
                        st_r = st_dc[10:-10].real
                        plt.plot(st_r)
                        
                        #%% get into proper obspy trace in proper time
                        
                        trace_obs = Trace()
                        trace_obs.data = st_r
                        st_obs = Stream(trace_obs)
                        st_obs.filter(type='lowpass',freq=0.5/100)
                        st_obs.plot(color='k')
                        
                        print('greatest amplitude = ', abs(min(st_obs[0].data)))
                        
                        
                    
#                        st_c1.remove_response(st_c2,water_level=60, plot=True)
                       
                
                        
                        
                        
                        # Dominant, Central, Bamdwidth50
    #                    dom1,cent1, bw, bwid25 = freq_info(tr,t1 ,t2)
    #                    
    #                    r=r1
            
#                    st_c1 = calibrate(tr)
##                    B=2*pi*rhoE*cE*(1/A)
##                    EI = sum(np.square(st_c1[0].data))
##                    Energy1 = B*(r*r)*EI
##            
#                    Amp1 = max(st_c1[0].data)
#                    
#                    
#            
#                    if 1418580590 > st[x].stats.starttime.timestamp > 1418580580:
#                        st_c1.plot(color='g')
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            