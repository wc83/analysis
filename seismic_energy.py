#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  3 13:39:16 2018

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

##%% single events
#
#reclen = 512
#chunksize = 100000 * reclen # Around 50 MB
#
#
#sum_events = 0
#day_start = 1416787200
#day_end = 1416787200 + 24*60*60 
#day =0
#daily_tot_e=np.zeros(shape=(0,5))
#daily_av_e=np.zeros(shape=(0,5))
#
#
#
#Energy=np.zeros(shape=(0,2))
#
#A=1
#rhoE=2500
#cE=2000
#pi=3.14159
#r1=4630
#r2=3370
#r3=2310
#r4=1300
#r5=810
#r6=7660
#       
##with io.open("/Users/william/Documents/scanner/output_data/m32.mseed", "rb") as fh:
#        # just month 2
#with io.open("/Users/william/Documents/scanner/output_data/EXP_all_data_stream_2_month_2.mseed", "rb") as fh:
#    while True:
#        with io.BytesIO() as buf:
#            c = fh.read(chunksize)
#            if not c:
#                break
#            buf.write(c)
#            buf.seek(0, 0)
#            st = obspy.read(buf)
#            
#            for x in range (0,int(len(st))):
#                if st[x].stats.station == "LB01":
#                    tr = st[x]
#                    st_c = calibrate(tr)
##                    print(st_c[0])
#                    B=2*pi*rhoE*cE*(1/A)
#                    
#                    EI = sum(np.square(st_c[0].data))
#                    EE= B*(r1*r1)*EI
#                    if EE < 1e10:
#                        Energy = np.lib.pad(Energy, ((0,1),(0,0)), 'constant', constant_values=(0))
#                        Energy[sum_events][0]=x
#                        Energy[sum_events][1]=EE
#                        
#                        sum_events += 1
#
#plt.figure(10) 
#plt.semilogy(Energy[:,0],Energy[:,1],'bx')
#
#

#%% daily events
#
#reclen = 512
#chunksize = 100000 * reclen # Around 50 MB
#
#day_events=0
#sum_events = 0
#sum_EE =0
#ad=0
#
#day_start = 1416787200
#day_end = 1416787200 + 24*60*60 
#day =0
#Energy_day=np.zeros(shape=(0,5))
#Energy=np.zeros(shape=(0,2))
#
#A=1
#rhoE=2500
#cE=2000
#pi=3.14159
#r1=4630
#r2=3370
#r3=2310
#r4=1300
#r5=810
#r6=7660
#
#E_max=1e14
#
#with io.open("/Users/william/Documents/scanner/output_data/m32.mseed", "rb") as fh:
#        # just month 2
##with io.open("/Users/william/Documents/scanner/output_data/EXP_all_data_stream_2_month_2.mseed", "rb") as fh:
#    while True:
#        with io.BytesIO() as buf:
#            c = fh.read(chunksize)
#            if not c:
#                break
#            buf.write(c)
#            buf.seek(0, 0)
#            st = obspy.read(buf)
#        
#            
#        # For each chunck of time, do analysis
##        print(st)
#        for x in range (0,int(len(st))):
#            if st[x].stats.station == "LB06":# or "LB02" or "LB03" or "LB04" or "LB05" or "LB06" :
#                if st[x].stats.station == "LB01":
#                    r=r1
#                if st[x].stats.station == "LB02":
#                    r=r2
#                if st[x].stats.station == "LB03":
#                    r=r3
#                if st[x].stats.station == "LB04":
#                    r=r4
#                if st[x].stats.station == "LB05":
#                    r=r5
#                if st[x].stats.station == "LB06":
#                    r=r6
#                if day_start < st[x].stats.starttime.timestamp < day_end:
#                    
#                    tr = st[x]
#                    st_c = calibrate(tr)
##                    print(st_c[0])
#                    B=2*pi*rhoE*cE*(1/A)
#                    
#                    EI = sum(np.square(st_c[0].data))
#                    EE= B*(r*r)*EI
#                    if EE < E_max:
#                        Energy = np.lib.pad(Energy, ((0,1),(0,0)), 'constant', constant_values=(0))
#                        Energy[sum_events][0]=x
#                        Energy[sum_events][1]=EE
#                        
#                        day_events += 1
#                        sum_events += 1
#                        sum_EE += EE
#                    
#                    
#                if st[x].stats.starttime.timestamp > day_end:
#                    if day_events > 0:
#                        
#                        av_EE = sum_EE/sum_events 
#                        
#                        Energy_day = np.lib.pad(Energy_day, ((0,1),(0,0)), 'constant', constant_values=(0))
#                        
#                        Energy_day[ad][0]=day + 1
#                        Energy_day[ad][1]=sum_EE
#                        Energy_day[ad][2]=sum_EE/day_events
#                        Energy_week[ad][3]=day_start
#                        ad += 1
#   
#                    else:
#                        print("LB inactive in day", day+1)
#                    
#                    sum_EE =0
#                    day_events =0
#                    
#                    for p in range(0,500):
#
#                        day_start=day_end
#                        day_end += 24*60*60
#                        day += 1
#                        if day_start < st[x].stats.starttime.timestamp < day_end:
#                            tr = st[x]
#                            st_c = calibrate(tr)
#        #                    print(st_c[0])
#                            B=2*pi*rhoE*cE*(1/A)
#                            
#                            EI = sum(np.square(st_c[0].data))
#                            EE= B*(r*r)*EI
#                            if EE < E_max:
#                                Energy = np.lib.pad(Energy, ((0,1),(0,0)), 'constant', constant_values=(0))
#                                Energy[sum_events][0]=x
#                                Energy[sum_events][1]=EE
#                                
#                                day_events += 1
#                                sum_events += 1
#                                sum_EE += EE
#                                break  
#                
#plt.figure(20)                                   
#plt.plot(Energy_day[:,0],Energy_day[:,1],'bx')
#plt.title('Total Energy per day')
#plt.xlabel('Day')
#plt.ylabel('Energy')
#plt.xlim([0,920])
#plt.ylim([0,2e12])
#
#plt.figure(30) 
#plt.plot(Energy_day[:,0],Energy_day[:,2],'bx')
#plt.title('Average Energy per day')
#plt.xlabel('Day')
#plt.xlabel('Day')
#plt.ylabel('Energy')
#plt.xlim([0,920])
#plt.ylim([0,1e12])



#%% weekly events

reclen = 512
chunksize = 100000 * reclen # Around 50 MB

week_events=0
sum_events = 0
sum_EE =0
ad=0

week_start = 1416787200
week_end = 1416787200 + 7*24*60*60 
day =0
Energy_week=np.zeros(shape=(0,4))
Energy=np.zeros(shape=(0,3))

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

E_max=1e16

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
        
            
        # For each chunck of time, do analysis
#        print(st)
        for x in range (0,int(len(st))):
            if st[x].stats.station == "LB06":# or "LB02" or "LB03" or "LB04" or "LB05" or "LB06" :
                if st[x].stats.station == "LB01":
                    r=r1
                if st[x].stats.station == "LB02":
                    r=r2
                if st[x].stats.station == "LB03":
                    r=r3
                if st[x].stats.station == "LB04":
                    r=r4
                if st[x].stats.station == "LB05":
                    r=r5
                if st[x].stats.station == "LB06":
                    r=r6
                if week_start < st[x].stats.starttime.timestamp < week_end:
                    
                    tr = st[x]
                    st_c = calibrate(tr)
#                    print(st_c[0])
                    B=2*pi*rhoE*cE*(1/A)
                    
                    EI = sum(np.square(st_c[0].data))
                    EE= B*(r*r)*EI
                    if EE < E_max:
                        Energy = np.lib.pad(Energy, ((0,1),(0,0)), 'constant', constant_values=(0))
                        Energy[sum_events][0]=sum_events + 1
                        Energy[sum_events][1]=EE
                        Energy[sum_events][2]=st[x].stats.starttime
                        
                        week_events += 1
                        sum_events += 1
                        sum_EE += EE
                    
                    
                if st[x].stats.starttime.timestamp > week_end:
                    if week_events > 1:
                        
                        av_EE = sum_EE/sum_events 
                        
                        Energy_week = np.lib.pad(Energy_week, ((0,1),(0,0)), 'constant', constant_values=(0))
                        
                        Energy_week[ad][0]=day + 1
                        Energy_week[ad][1]=sum_EE
                        Energy_week[ad][2]=sum_EE/week_events
                        Energy_week[ad][3]=week_start
                        ad += 1
   
                    else:
                        print("LB inactive in day", day+1)
                    
                    sum_EE =0
                    day_events =0
                    
                    for p in range(0,500):

                        week_start=week_end
                        week_end += 7*24*60*60
                        day += 1
                        if week_start < st[x].stats.starttime.timestamp < week_end:
                            tr = st[x]
                            st_c = calibrate(tr)
        #                    print(st_c[0])
                            B=2*pi*rhoE*cE*(1/A)
                            
                            EI = sum(np.square(st_c[0].data))
                            EE= B*(r*r)*EI
                            if EE < E_max:
                                Energy = np.lib.pad(Energy, ((0,1),(0,0)), 'constant', constant_values=(0))
                                Energy[sum_events][0]=x
                                Energy[sum_events][1]=EE
                                
                                week_events += 1
                                sum_events += 1
                                sum_EE += EE
                                break  
                
plt.figure(20)                                   
plt.semilogy(Energy_week[:,0],Energy_week[:,1],'bx')
plt.title('Total Energy per week')
plt.xlabel('week')
plt.ylabel('Energy')
plt.xlim([0,140])
plt.ylim([1e1,1e14])

plt.figure(30) 
plt.semilogy(Energy_week[:,0],Energy_week[:,2],'bx')
plt.title('Average Energy per week')
plt.xlabel('week')
plt.ylabel('Energy')
plt.xlim([0,140])
plt.ylim([1e1,1e12])


plt.figure(40)
plt.semilogy(Energy[:,0],Energy[:,1],'bx')
plt.title('Energy per Explosion')
plt.xlabel('Explosion')
plt.ylabel('Energy')
plt.ylim([1e1,1e14])



