#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 12:51:04 2019

@author: william
"""

import os
from collections import OrderedDict
import numpy as np
import obspy
import scipy.signal as sgn
import matplotlib.pyplot as plt 
import matplotlib.mlab as mlab
from obspy.core import read
from obspy.clients.earthworm import Client
from obspy import UTCDateTime
from obspy.signal.trigger import trigger_onset
from numpy import genfromtxt
#from scipy.signal import welch
from obspy import Stream
from scipy import integrate

#%%
f=open("/Users/william/Documents/scanner/Olivers_Final_exp_2014_2017.txt", "r")

cat =f.read()

events=np.zeros(shape=(6104,1))
num=0

for x in range(0,len(cat),20):  
    
    year = int(float(cat[x])*1000 + float(cat[x+1])*100 + float(cat[x+2])*10 + float(cat[x+3]))
    month = int(float(cat[x+5])*10+float(cat[x+6]))
    day = int(float(cat[x+8])*10 + float(cat[x+9])) 

    
    hour = int(float(cat[x+11])*10 + float(cat[x+12])) 
    minute = int(float(cat[x+14])*10 + float(cat[x+15]))
    second = int(float(cat[x+17])*10 + float(cat[x+18]))


    events[num] = UTCDateTime(year, month, day, hour, minute, second)
    num+=1
    
#%%  
    
#Repose = np.zeros(shape=(0,1))
#num=0
#for x in range(1,len(events)):
#    rep = events[x]-events[x-1]
#    
#    if rep < 6*60*60:
#        Repose = np.lib.pad(Repose, ((0,1),(0,0)), 'constant', constant_values=(0))             
#        Repose[num][0]=rep/60
#        num+=1
#        
#        
#        
#plt.figure(2)
#plt.hist(Repose,bins=180)
#plt.xlabel('Repose Time [mins]')
#plt.ylabel('Occurance [#]')
#plt.title('Repose Times - Olivers Catalogue')

#%% constants

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
  
        
color1='k'
net = 'Z4'  # 
loc = ''    # location, it depends mostly of which network you are in. 
client = Client('138.253.112.23', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23       

num2=0

Energy=np.zeros(shape=(len(events),2))


for x in range(0,len(events)):
            print(x)
            try:
                sta = 'LB01' # STATION 
                cha = 'HHZ' # CHANNEL
                t1 = UTCDateTime(events[x]) #the format is year:day_of_the_year:month
                t2 = t1 + 80
                st = Stream()
                st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
            
                st.detrend(type='linear')
                st.detrend(type='demean')
                st.filter(type='bandpass',freqmin=0.2, freqmax=10)
                
                r=r1
                tr=st[0]
                st_c = calibrate1(tr)
        
                B=2*pi*rhoE*cE*(1/A)
                dl=len(st_c[0].data)
                p = np.linspace(0,dl/100, num=dl)
                
                y= st_c[0].data
                y2= np.square(st_c[0].data)
                
                
                y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
                EI = y_int2[-1]
                EE= B*(r*r)*EI
                
                Energy[x][0]=t1.timestamp
                Energy[x][1]=EE
                
                
                
            except:
                try:
                    sta = 'LB02' # STATION 
                    cha = 'HHZ' # CHANNEL
                    t1 = UTCDateTime(events[x]) #the format is year:day_of_the_year:month
                    t2 = t1 + 80
                    st = Stream()
                    st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
                
                    st.detrend(type='linear')
                    st.detrend(type='demean')
                    st.filter(type='bandpass',freqmin=0.2, freqmax=10)
                    
                    r=r2
                    tr=st[0]
                    st_c = calibrate1(tr)
            
                    B=2*pi*rhoE*cE*(1/A)
                    dl=len(st_c[0].data)
                    p = np.linspace(0,dl/100, num=dl)
                    
                    y= st_c[0].data
                    y2= np.square(st_c[0].data)
                    
                    
                    y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
                    EI = y_int2[-1]
                    EE= B*(r*r)*EI
                    
                    Energy[x][0]=t1.timestamp
                    Energy[x][1]=EE
                
                except:
                    try:
                        sta = 'LB03' # STATION 
                        cha = 'HHZ' # CHANNEL
                        t1 = UTCDateTime(events[x]) #the format is year:day_of_the_year:month
                        t2 = t1 + 80
                        st = Stream()
                        st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
                    
                        st.detrend(type='linear')
                        st.detrend(type='demean')
                        st.filter(type='bandpass',freqmin=0.2, freqmax=10)
                        
                        r=r3
                        tr=st[0]
                        st_c = calibrate1(tr)
                
                        B=2*pi*rhoE*cE*(1/A)
                        dl=len(st_c[0].data)
                        p = np.linspace(0,dl/100, num=dl)
                        
                        y= st_c[0].data
                        y2= np.square(st_c[0].data)
                        
                        
                        y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
                        EI = y_int2[-1]
                        EE= B*(r*r)*EI
                        
                        Energy[x][0]=t1.timestamp
                        Energy[x][1]=EE
                
                    except:
                        try:
                            sta = 'LB04' # STATION 
                            cha = 'HHZ' # CHANNEL
                            t1 = UTCDateTime(events[x]) #the format is year:day_of_the_year:month
                            t2 = t1 + 80
                            st = Stream()
                            st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
                        
                            st.detrend(type='linear')
                            st.detrend(type='demean')
                            st.filter(type='bandpass',freqmin=0.2, freqmax=10)
                            
                            r=r4
                            tr=st[0]
                            st_c = calibrate1(tr)
                    
                            B=2*pi*rhoE*cE*(1/A)
                            dl=len(st_c[0].data)
                            p = np.linspace(0,dl/100, num=dl)
                            
                            y= st_c[0].data
                            y2= np.square(st_c[0].data)
                            
                            
                            y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
                            EI = y_int2[-1]
                            EE= B*(r*r)*EI
                            
                            Energy[x][0]=t1.timestamp
                            Energy[x][1]=EE
                
                        except:
                            try:
                                sta = 'LB05' # STATION 
                                cha = 'HHZ' # CHANNEL
                                t1 = UTCDateTime(events[x]) #the format is year:day_of_the_year:month
                                t2 = t1 + 80
                                st = Stream()
                                st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
                            
                                st.detrend(type='linear')
                                st.detrend(type='demean')
                                st.filter(type='bandpass',freqmin=0.2, freqmax=10)
                                
                                r=r5
                                tr=st[0]
                                st_c = calibrate1(tr)
                        
                                B=2*pi*rhoE*cE*(1/A)
                                dl=len(st_c[0].data)
                                p = np.linspace(0,dl/100, num=dl)
                                
                                y= st_c[0].data
                                y2= np.square(st_c[0].data)
                                
                                
                                y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
                                EI = y_int2[-1]
                                EE= B*(r*r)*EI
                                
                                Energy[x][0]=t1.timestamp
                                Energy[x][1]=EE
                
                            except:
                                try:
                                    sta = 'LB06' # STATION 
                                    cha = 'HHZ' # CHANNEL
                                    t1 = UTCDateTime(events[x]) #the format is year:day_of_the_year:month
                                    t2 = t1 + 80
                                    st = Stream()
                                    st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
                                
                                    st.detrend(type='linear')
                                    st.detrend(type='demean')
                                    st.filter(type='bandpass',freqmin=0.2, freqmax=10)
                                    
                                    r=r6
                                    tr=st[0]
                                    st_c = calibrate1(tr)
                            
                                    B=2*pi*rhoE*cE*(1/A)
                                    dl=len(st_c[0].data)
                                    p = np.linspace(0,dl/100, num=dl)
                                    
                                    y= st_c[0].data
                                    y2= np.square(st_c[0].data)
                                    
                                    
                                    y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
                                    EI = y_int2[-1]
                                    EE= B*(r*r)*EI
                                    
                                    Energy[x][0]=t1.timestamp
                                    Energy[x][1]=EE
                
                                except:
                                    try:
                                        sta = 'LS01' # STATION 
                                        cha = 'EHZ' # CHANNEL
                                        t1 = UTCDateTime(events[x]) #the format is year:day_of_the_year:month
                                        t2 = t1 + 80
                                        st = Stream()
                                        st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
                                    
                                        st.detrend(type='linear')
                                        st.detrend(type='demean')
                                        st.filter(type='bandpass',freqmin=0.2, freqmax=10)
                                        
                                        r=rs1
                                        tr=st[0]
                                        st_c = calibrate1(tr)
                                
                                        B=2*pi*rhoE*cE*(1/A)
                                        dl=len(st_c[0].data)
                                        p = np.linspace(0,dl/100, num=dl)
                                        
                                        y= st_c[0].data
                                        y2= np.square(st_c[0].data)
                                        
                                        
                                        y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
                                        EI = y_int2[-1]
                                        EE= B*(r*r)*EI
                                        
                                        Energy[x][0]=t1.timestamp
                                        Energy[x][1]=EE
                                        
                                        
                                        
                                    except:
                                        try:
                                            sta = 'LS02' # STATION 
                                            cha = 'EHZ' # CHANNEL
                                            t1 = UTCDateTime(events[x]) #the format is year:day_of_the_year:month
                                            t2 = t1 + 80
                                            st = Stream()
                                            st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
                                        
                                            st.detrend(type='linear')
                                            st.detrend(type='demean')
                                            st.filter(type='bandpass',freqmin=0.2, freqmax=10)
                                            
                                            r=rs2
                                            tr=st[0]
                                            st_c = calibrate1(tr)
                                    
                                            B=2*pi*rhoE*cE*(1/A)
                                            dl=len(st_c[0].data)
                                            p = np.linspace(0,dl/100, num=dl)
                                            
                                            y= st_c[0].data
                                            y2= np.square(st_c[0].data)
                                            
                                            
                                            y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
                                            EI = y_int2[-1]
                                            EE= B*(r*r)*EI
                                            
                                            Energy[x][0]=t1.timestamp
                                            Energy[x][1]=EE
                                        
                                        except:
                                            try:
                                                sta = 'LS03' # STATION 
                                                cha = 'EHZ' # CHANNEL
                                                t1 = UTCDateTime(events[x]) #the format is year:day_of_the_year:month
                                                t2 = t1 + 80
                                                st = Stream()
                                                st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
                                            
                                                st.detrend(type='linear')
                                                st.detrend(type='demean')
                                                st.filter(type='bandpass',freqmin=0.2, freqmax=10)
                                                
                                                r=rs3
                                                tr=st[0]
                                                st_c = calibrate1(tr)
                                        
                                                B=2*pi*rhoE*cE*(1/A)
                                                dl=len(st_c[0].data)
                                                p = np.linspace(0,dl/100, num=dl)
                                                
                                                y= st_c[0].data
                                                y2= np.square(st_c[0].data)
                                                
                                                
                                                y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
                                                EI = y_int2[-1]
                                                EE= B*(r*r)*EI
                                                
                                                Energy[x][0]=t1.timestamp
                                                Energy[x][1]=EE
                                        
                                            except:
                                                try:
                                                    sta = 'LS04' # STATION 
                                                    cha = 'EHZ' # CHANNEL
                                                    t1 = UTCDateTime(events[x]) #the format is year:day_of_the_year:month
                                                    t2 = t1 + 80
                                                    st = Stream()
                                                    st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
                                                
                                                    st.detrend(type='linear')
                                                    st.detrend(type='demean')
                                                    st.filter(type='bandpass',freqmin=0.2, freqmax=10)
                                                    
                                                    r=rs4
                                                    tr=st[0]
                                                    st_c = calibrate1(tr)
                                            
                                                    B=2*pi*rhoE*cE*(1/A)
                                                    dl=len(st_c[0].data)
                                                    p = np.linspace(0,dl/100, num=dl)
                                                    
                                                    y= st_c[0].data
                                                    y2= np.square(st_c[0].data)
                                                    
                                                    
                                                    y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
                                                    EI = y_int2[-1]
                                                    EE= B*(r*r)*EI
                                                    
                                                    Energy[x][0]=t1.timestamp
                                                    Energy[x][1]=EE
                                        
                                                except:
                                                    try:
                                                        sta = 'LS05' # STATION 
                                                        cha = 'EHZ' # CHANNEL
                                                        t1 = UTCDateTime(events[x]) #the format is year:day_of_the_year:month
                                                        t2 = t1 + 80
                                                        st = Stream()
                                                        st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
                                                    
                                                        st.detrend(type='linear')
                                                        st.detrend(type='demean')
                                                        st.filter(type='bandpass',freqmin=0.2, freqmax=10)
                                                        
                                                        r=rs5
                                                        tr=st[0]
                                                        st_c = calibrate1(tr)
                                                
                                                        B=2*pi*rhoE*cE*(1/A)
                                                        dl=len(st_c[0].data)
                                                        p = np.linspace(0,dl/100, num=dl)
                                                        
                                                        y= st_c[0].data
                                                        y2= np.square(st_c[0].data)
                                                        
                                                        
                                                        y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
                                                        EI = y_int2[-1]
                                                        EE= B*(r*r)*EI
                                                        
                                                        Energy[x][0]=t1.timestamp
                                                        Energy[x][1]=EE
                                        
                                                    except:
                                                        try:
                                                            sta = 'LS06' # STATION 
                                                            cha = 'EHZ' # CHANNEL
                                                            t1 = UTCDateTime(events[x]) #the format is year:day_of_the_year:month
                                                            t2 = t1 + 80
                                                            st = Stream()
                                                            st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
                                                        
                                                            st.detrend(type='linear')
                                                            st.detrend(type='demean')
                                                            st.filter(type='bandpass',freqmin=0.2, freqmax=10)
                                                            
                                                            r=rs6
                                                            tr=st[0]
                                                            st_c = calibrate1(tr)
                                                    
                                                            B=2*pi*rhoE*cE*(1/A)
                                                            dl=len(st_c[0].data)
                                                            p = np.linspace(0,dl/100, num=dl)
                                                            
                                                            y= st_c[0].data
                                                            y2= np.square(st_c[0].data)
                                                            
                                                            
                                                            y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
                                                            EI = y_int2[-1]
                                                            EE= B*(r*r)*EI
                                                            
                                                            Energy[x][0]=t1.timestamp
                                                            Energy[x][1]=EE
                                        
                                                        except:
                                                            print('not worked for', UTCDateTime(events[x]) )
                                





np.savetxt("/Users/william/Documents/scanner/Oliver_Energy.csv", Energy,delimiter=",",header="Time,Energy")























