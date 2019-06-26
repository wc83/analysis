#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 11:16:46 2018

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


sum_events = 0
day_start = 1416787200
day_end = 1416787200 + 24*60*60 
day =0
daily_tot_e=np.zeros(shape=(0,5))
daily_av_e=np.zeros(shape=(0,5))



Energy=np.zeros(shape=(0,2))

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
       
with io.open("/Users/william/Documents/scanner/output_data/m32.mseed", "rb") as fh:
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
            
            for x in range (0,int(len(st))):
                if st[x].stats.station == "LB01":
#                    if 1451607000 < st[x].stats.starttime.timestamp: #just 2014 and 2015 explosions
                        tr = st[x]
                        st_c = calibrate1(tr)
    #                    print(st_c[0])
                        B=2*pi*rhoE*cE*(1/A)
                        
                        EI = sum(np.square(st_c[0].data))
                        EE= B*(r1*r1)*EI
                        if EE < 5e8:
                            Energy = np.lib.pad(Energy, ((0,1),(0,0)), 'constant', constant_values=(0))
                            Energy[sum_events][0]=x
                            Energy[sum_events][1]=EE
                            
                            sum_events += 1

#plt.figure(10) 
#plt.semilogy(Energy[:,0],Energy[:,1],'bx')


Emax = max(Energy[:,1])
Emin = min(Energy[:,1])

print('Energy from', Emin, 'to', Emax)

    
# 10000000000

#%%

plt.figure(4001)
plt.hist(Energy[:,1],bins=30, histtype='step')
plt.yscale('log')
plt.xscale('log')
plt.xlim([1e5,1e9])
plt.xlabel('Energy')
plt.ylabel('Occurance [#]')
#plt.title('Mag-Freq: Power-law Relationship')

            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            