#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 09:26:26 2019

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







#%% Plot catalogue
  
#    
#    
#   # STATION, CHANNEL (DDF --> 400 Hz), NETWWORK AND LOCATION CODES 
#color1='b'
#net = 'Z4'  # 
#loc = ''    # location, it depends mostly of which network you are in. 
#client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
#num_26=0
#for x in range(0,len(events)):
#    if UTCDateTime(events[x]).year == 2017:
#        num_26 +=1
#        if 100 <= num_26  < 200:
#            try:
#                sta = 'LB01' # STATION 
#                cha = 'HHZ' # CHANNEL
#                t1 = UTCDateTime(events[x]-30) #the format is year:day_of_the_year:month
#                t2 = t1 + 90
#                st = Stream()
#                st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
#            
#                st.detrend(type='linear')
#                st.detrend(type='demean')
#                st.filter(type='bandpass',freqmin=0.2, freqmax=10)
#                st.plot(color=color1 ,starttime=t1, endtime=t2)
#            except:
#                try:
#                    sta = 'LB02' # STATION 
#                    cha = 'HHZ' # CHANNEL
#                    t1 = UTCDateTime(events[x]-30) #the format is year:day_of_the_year:month
#                    t2 = t1 + 90
#                    st = Stream()
#                    st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
#                
#                    st.detrend(type='linear')
#                    st.detrend(type='demean')
#                    st.filter(type='bandpass',freqmin=0.2, freqmax=10)
#                    st.plot(color=color1,starttime=t1, endtime=t2)
#                except:
#                    try:
#                        sta = 'LB03' # STATION 
#                        cha = 'HHZ' # CHANNEL
#                        t1 = UTCDateTime(events[x]-30) #the format is year:day_of_the_year:month
#                        t2 = t1 + 90
#                        st = Stream()
#                        st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
#                    
#                        st.detrend(type='linear')
#                        st.detrend(type='demean')
#                        st.filter(type='bandpass',freqmin=0.2, freqmax=10)
#                        st.plot(color=color1,starttime=t1, endtime=t2)
#                    except:
#                        try:
#                            sta = 'LS04' # STATION 
#                            cha = 'EHZ' # CHANNEL
#                            t1 = UTCDateTime(events[x]-30) #the format is year:day_of_the_year:month
#                            t2 = t1 + 90
#                            st = Stream()
#                            st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
#                        
#                            st.detrend(type='linear')
#                            st.detrend(type='demean')
#                            st.filter(type='bandpass',freqmin=0.2, freqmax=10)
#                            st.plot(color=color1,starttime=t1, endtime=t2)
#                        except:
#                            try:
#                                sta = 'LS05' # STATION 
#                                cha = 'EHZ' # CHANNEL
#                                t1 = UTCDateTime(events[x]-30) #the format is year:day_of_the_year:month
#                                t2 = t1 + 90
#                                st = Stream()
#                                st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
#                            
#                                st.detrend(type='linear')
#                                st.detrend(type='demean')
#                                st.filter(type='bandpass',freqmin=0.2, freqmax=10)
#                                st.plot(color=color1,starttime=t1, endtime=t2)
#                            except:
#                                print(events[x])
#            
                
#%% Find how many of similarities between catalogues                
#
#my_cat = genfromtxt("/Users/william/Documents/scanner/all_stations/Explosion_Catalogue_V4c.csv", delimiter=',',skip_header=1)
#ev_t = my_cat[:,0]
#            
#same_events = np.zeros(shape=(0,1))
#only_oliver = np.zeros(shape=(0,1))
#only_will = np.zeros(shape=(0,1))
#
#num_s=0
#num_o=0
#num_w=0
#
#for x in range(0,6104):
#    OE = events[x]
#    near_event,ind=find_nearest(ev_t[:],OE )  #find how near an event in Olivers catalogue is to mine          
#    if abs(OE-near_event) <  120:
#        same_events = np.lib.pad(same_events, ((0,1),(0,0)), 'constant', constant_values=(0))
#        same_events[num_s][0] = OE
#        num_s += 1
#    else:
#        only_oliver= np.lib.pad(only_oliver, ((0,1),(0,0)), 'constant', constant_values=(0))
#        only_oliver[num_o][0] = OE
#        num_o += 1
#            
#for x in range(0,len(ev_t)):
#    WE=ev_t[x]
#    near_same,indx=find_nearest(same_events[:],WE )
#    if abs(WE-near_same) > 120:
#        only_will = np.lib.pad(only_will, ((0,1),(0,0)), 'constant', constant_values=(0))
#        only_will[num_w][0]=WE
#        num_w +=1
#            
        
#np.savetxt("/Users/william/Documents/scanner/both_catalogues.csv", same_events,delimiter=",",header="event_time")
#np.savetxt("/Users/william/Documents/scanner/only_olivers_catalogue.csv",  only_oliver,delimiter=",",header="event_time")
#np.savetxt("/Users/william/Documents/scanner/only_wills_catalogue.csv", only_will,delimiter=",",header="event_time")

#%%    Check same detections
        
        
#color1='k'
#net = 'Z4'  # 
#loc = ''    # location, it depends mostly of which network you are in. 
#client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23       
#for x in range(1100,1200):
#    
#    try:
#        sta = 'LB01' # STATION 
#        cha = 'HHZ' # CHANNEL
#        t1 = UTCDateTime(same_events[x]-30) #the format is year:day_of_the_year:month
#        t2 = t1 + 90
#        st = Stream()
#        st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
#    
#        st.detrend(type='linear')
#        st.detrend(type='demean')
#        st.filter(type='bandpass',freqmin=0.2, freqmax=10)
#        st.plot(color=color1 ,starttime=t1, endtime=t2)
#    except:
#        try:
#            sta = 'LB02' # STATION 
#            cha = 'HHZ' # CHANNEL
#            t1 = UTCDateTime(same_events[x]-30) #the format is year:day_of_the_year:month
#            t2 = t1 + 90
#            st = Stream()
#            st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
#        
#            st.detrend(type='linear')
#            st.detrend(type='demean')
#            st.filter(type='bandpass',freqmin=0.2, freqmax=10)
#            st.plot(color=color1,starttime=t1, endtime=t2)
#        except:
#            try:
#                sta = 'LB03' # STATION 
#                cha = 'HHZ' # CHANNEL
#                t1 = UTCDateTime(same_events[x]-30) #the format is year:day_of_the_year:month
#                t2 = t1 + 90
#                st = Stream()
#                st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
#            
#                st.detrend(type='linear')
#                st.detrend(type='demean')
#                st.filter(type='bandpass',freqmin=0.2, freqmax=10)
#                st.plot(color=color1,starttime=t1, endtime=t2)
#            except:
#                try:
#                    sta = 'LS04' # STATION 
#                    cha = 'EHZ' # CHANNEL
#                    t1 = UTCDateTime(same_events[x]-30) #the format is year:day_of_the_year:month
#                    t2 = t1 + 90
#                    st = Stream()
#                    st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
#                
#                    st.detrend(type='linear')
#                    st.detrend(type='demean')
#                    st.filter(type='bandpass',freqmin=0.2, freqmax=10)
#                    st.plot(color=color1,starttime=t1, endtime=t2)
#                except:
#                    try:
#                        sta = 'LS05' # STATION 
#                        cha = 'EHZ' # CHANNEL
#                        t1 = UTCDateTime(same_events[x]-30) #the format is year:day_of_the_year:month
#                        t2 = t1 + 90
#                        st = Stream()
#                        st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
#                    
#                        st.detrend(type='linear')
#                        st.detrend(type='demean')
#                        st.filter(type='bandpass',freqmin=0.2, freqmax=10)
#                        st.plot(color=color1,starttime=t1, endtime=t2)
#                    except:
#                        print(same_events[x])


#%%    Check only_olivers detections
        
only_oliver= genfromtxt("/Users/william/Documents/scanner/only_olivers_catalogue.csv", delimiter=',',skip_header=1)
  

    
color1='k'
net = 'Z4'  # 
loc = ''    # location, it depends mostly of which network you are in. 
client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23       

num2=0

for x in range(0,len(only_oliver)):
#    if UTCDateTime(only_oliver[x]).year == 2017:
        num2 +=1
        if 100 <= num2  < 200:
            try:
                sta = 'LB01' # STATION 
                cha = 'HHZ' # CHANNEL
                t1 = UTCDateTime(only_oliver[x]-30) #the format is year:day_of_the_year:month
                t2 = t1 + 90
                st = Stream()
                st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
            
                st.detrend(type='linear')
                st.detrend(type='demean')
                st.filter(type='bandpass',freqmin=1, freqmax=10)
                st.plot(color=color1 ,starttime=t1, endtime=t2)
            except:
                try:
                    sta = 'LB02' # STATION 
                    cha = 'HHZ' # CHANNEL
                    t1 = UTCDateTime(only_oliver[x]-30) #the format is year:day_of_the_year:month
                    t2 = t1 + 90
                    st = Stream()
                    st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
                
                    st.detrend(type='linear')
                    st.detrend(type='demean')
                    st.filter(type='bandpass',freqmin=1, freqmax=10)
                    st.plot(color=color1,starttime=t1, endtime=t2)
                except:
                    try:
                        sta = 'LB03' # STATION 
                        cha = 'HHZ' # CHANNEL
                        t1 = UTCDateTime(only_oliver[x]-30) #the format is year:day_of_the_year:month
                        t2 = t1 + 90
                        st = Stream()
                        st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
                    
                        st.detrend(type='linear')
                        st.detrend(type='demean')
                        st.filter(type='bandpass',freqmin=1, freqmax=10)
                        st.plot(color=color1,starttime=t1, endtime=t2)
                    except:
                        try:
                            sta = 'LS04' # STATION 
                            cha = 'EHZ' # CHANNEL
                            t1 = UTCDateTime(only_oliver[x]-30) #the format is year:day_of_the_year:month
                            t2 = t1 + 90
                            st = Stream()
                            st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
                        
                            st.detrend(type='linear')
                            st.detrend(type='demean')
                            st.filter(type='bandpass',freqmin=1, freqmax=10)
                            st.plot(color=color1,starttime=t1, endtime=t2)
                        except:
                            try:
                                sta = 'LS05' # STATION 
                                cha = 'EHZ' # CHANNEL
                                t1 = UTCDateTime(only_oliver[x]-30) #the format is year:day_of_the_year:month
                                t2 = t1 + 90
                                st = Stream()
                                st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
                            
                                st.detrend(type='linear')
                                st.detrend(type='demean')
                                st.filter(type='bandpass',freqmin=1, freqmax=10)
                                st.plot(color=color1,starttime=t1, endtime=t2)
                            except:
                                try:
                                    sta = 'LB06' # STATION 
                                    cha = 'HHZ' # CHANNEL
                                    t1 = UTCDateTime(only_oliver[x]-30) #the format is year:day_of_the_year:month
                                    t2 = t1 + 90
                                    st = Stream()
                                    st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
                                
                                    st.detrend(type='linear')
                                    st.detrend(type='demean')
                                    st.filter(type='bandpass',freqmin=1, freqmax=10)
                                    st.plot(color=color1,starttime=t1, endtime=t2)
                                except:
                                    print(only_oliver[x])
        
        
        
#%%
                                    
                                    
                                    
        
        
