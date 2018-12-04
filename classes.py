#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 12:02:23 2018

@author: william
"""


import obspy
import io
from obspy import read
from obspy.clients.earthworm import Client
from obspy import UTCDateTime
from obspy import Stream
import numpy as np
import matplotlib.pyplot as plt


reclen = 512
chunksize = 100000 * reclen # Around 50 MB

A=1
rhoE=2500
cE=2000
pi=3.14159
r1=4630

with io.open("/Users/william/Documents/scanner/output_data/m30.mseed", "rb") as fh:
        # just month 2
#with io.open("/Users/william/Documents/scanner/output_data/EXP_all_data_stream_2_month_2.mseed", "rb") as fh:
    while True:
        with io.BytesIO() as buf:
            c = fh.read(chunksize)
            if not c:
                break
            buf.write(c)
            buf.seek(0, 0)
            st = obspy.read(buf)
        
            
        # For each chunck of time, do analysis
#        print(st)
        for x in range (0,int(len(st))):
            if st[x].stats.station == "LB01":
#                print(st[x].stats.starttime.timestamp)
                if 1418580590 > st[x].stats.starttime.timestamp > 1418580580:
                        tr = st[x]
                        tr.detrend(type='linear')
                        tr.detrend(type='demean')
                        t1 = st[x].stats.starttime
                        t2 = st[x].stats.endtime
                        window=t2-t1
                        tr_data=st[x].data
                        m=np.mean(tr_data)
                        tr_data = tr_data-m
                        
                        st_c1 = calibrate(tr)
                        st_c1.plot(color='r')
                        
                        B=2*pi*rhoE*cE*(1/A)                    
                        EI = sum(np.square(tr.data))
                        EE1= B*(r1*r1)*EI
                        print("Energy =", EE1)
                        
                        peak,cf, bwid50, bwid25 = freq_info(tr,t1,t2)
                        print("central f = ",cf)
                
                if  1438145500 > st[x].stats.starttime.timestamp > 1438145480:
                        tr = st[x]
                        tr.detrend(type='linear')
                        tr.detrend(type='demean')
                        t1 = st[x].stats.starttime
                        t2 = st[x].stats.endtime
                        window=t2-t1
                        tr_data=st[x].data
                        m=np.mean(tr_data)
                        tr_data = tr_data-m
                        
                        st_c1 = calibrate(tr)
                        st_c1.plot(color='b')
                        
                        B=2*pi*rhoE*cE*(1/A)                    
                        EI = sum(np.square(tr.data))
                        EE2= B*(r1*r1)*EI
                        print("Energy =", EE2)
                        
                        peak,cf, bwid50, bwid25 = freq_info(tr,t1,t2)
                        print("central f = ",cf)
                        
                if  1466176600 < st[x].stats.starttime.timestamp < 1466176620:
                        tr = st[x]
                        tr.detrend(type='linear')
                        tr.detrend(type='demean')
                        t1 = st[x].stats.starttime
                        t2 = st[x].stats.endtime
                        window=t2-t1
                        tr_data=st[x].data
                        m=np.mean(tr_data)
                        tr_data = tr_data-m
                        
                        st_c1 = calibrate(tr)
                        st_c1.plot(color='g')
                        
                        B=2*pi*rhoE*cE*(1/A)                    
                        EI = sum(np.square(tr.data))
                        EE3= B*(r1*r1)*EI
                        print("Energy =", EE3)
                            
                        peak,cf, bwid50, bwid25 = freq_info(tr,t1,t2)
                        print("central f = ",cf)






















































