#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  3 15:23:28 2018

@author: william
"""
import obspy
from obspy import read
import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
from obspy import Stream
from obspy import UTCDateTime
from obspy.clients.earthworm import Client

cha = 'HHZ' # CHANNEL
net = 'Z4'  # 
loc = ''    # location, it depends mostly of which network you are in. 
client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23

percentage=0
count=0
all_count =0
all_percentage =0
larger_s_tot=0
                
cat= genfromtxt("/Users/william/Documents/scanner/all_stations/Explosion_Catalogue_V4c.csv", delimiter=',',skip_header=1)

pands=np.zeros(shape=(1,3))
pands[0,0]=cat[0,0]
pands[0,1]=1
pands[0,2]=0
pands = np.lib.pad(pands, ((0,1),(0,0)), 'constant', constant_values=(0))

last_p = cat[0,0]

p_tot=1
s_tot=0

event_num=1

sec=np.zeros(shape=(0,3))
num=0
#%%
#for x in range(1,5000):
for x in range(1,len(cat)):
#    if 1459468800 < cat[x,0] < 1475884800:
        pands[event_num,0]=cat[x,0]
        
        if cat[x,0] - last_p < 10*60:
            sec = np.lib.pad(sec, ((0,1),(0,0)), 'constant', constant_values=(0))
            sec[num][0]= last_p
            sec[num][1]= cat[x,0]
            
#%%
            if cat[x,10]==1:
                sta = 'LB01' # STATION 
                
                arrival1=last_p
                arrival2=cat[x,0]
                
                t1 = UTCDateTime(arrival1)-30 #the format is year:day_of_the_year:month
                t2 = t1 + 120
                st = Stream()
#                try:
                st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                st.detrend(type='linear')
                st.detrend(type='demean')
                peak1=max(st[0].data)
#                    g=0
#                except:
#                    g=1
            
                t1 = UTCDateTime(arrival2)-30 #the format is year:day_of_the_year:month
                t2 = t1 + 120
                st = Stream()
#                try:
                st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                st.detrend(type='linear')
                st.detrend(type='demean')
                peak2=max(st[0].data)
#                    f=0
#                except:
#                    f=1
            
            
                peak_percent = (peak2/peak1)*100
                all_count += 1
                all_percentage += peak_percent
                
#                if  f==0 and g==0:
                if peak1 > peak2:
                    pands[event_num,1]=2
                    s_tot += 1
#                    print(peak_percent)
                    count += 1
                    percentage += peak_percent
                
                    arrival=last_p
                    t1 = UTCDateTime(arrival)-60 #the format is year:day_of_the_year:month
                    t2 = t1 + (cat[x,0] - last_p) + 120
                    st = Stream()
                    st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                    st.detrend(type='linear')
                    st.detrend(type='demean')
                    st.filter(type='bandpass',freqmin=0.5, freqmax=6)
                    st.plot(color='b',starttime=t1, endtime=t2)
                    print(cat[x,0] - last_p)
#                        except:
                else:
#                        pands[event_num,1]=1
#                        last_p=cat[x,0]
#                        p_tot +=1
                    larger_s_tot += 1
#%%              
            if cat[x-1,10]==0 and cat[x-1,11]==1:
                sta = 'LB02' # STATION 
                arrival1=last_p
                arrival2=cat[x,0]
                
                t1 = UTCDateTime(arrival1)-30 #the format is year:day_of_the_year:month
                t2 = t1 + 120
                st = Stream()
                try:
                    st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                    st.detrend(type='linear')
                    st.detrend(type='demean')
                    peak1=max(st[0].data)
                    g=0
                except:
                    g=1
            
                t1 = UTCDateTime(arrival2)-30 #the format is year:day_of_the_year:month
                t2 = t1 + 120
                st = Stream()
                try:
                    st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                    st.detrend(type='linear')
                    st.detrend(type='demean')
                    peak2=max(st[0].data)
                    f=0
                except:
                    f=1
            
            
                peak_percent = (peak2/peak1)*100
                all_count += 1
                all_percentage += peak_percent
                
                if  f==0 and g==0:
                    if peak1 > peak2 :
                        pands[event_num,1]=2
                        s_tot += 1
    #                    print(peak_percent)
                        count += 1
                        percentage += peak_percent
                    
    #                    try:
    #                        arrival=last_p
    #                        t1 = UTCDateTime(arrival)-120 #the format is year:day_of_the_year:month
    #                        t2 = t1 + 780
    #                        st = Stream()
    #                        st = client.get_waveforms(net, sta, '', cha, t1 , t2)
    #                        st.detrend(type='linear')
    #                        st.detrend(type='demean')
    #                        st.filter(type='bandpass',freqmin=0.1, freqmax=6)
    #                        st.plot(color='b',starttime=t1, endtime=t2)
    #                    except:
    #                        h=1
                    else:
#                        pands[event_num,1]=1
#                        last_p=cat[x,0]
#                        p_tot +=1
                        larger_s_tot += 1
#%%              
            if cat[x-1,10] + cat[x-1,11]==0 and cat[x-1,12]==1:
                sta = 'LB03' # STATION 
                arrival1=last_p
                arrival2=cat[x,0]
                
                t1 = UTCDateTime(arrival1)-30 #the format is year:day_of_the_year:month
                t2 = t1 + 120
                st = Stream()
                try:
                    st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                    st.detrend(type='linear')
                    st.detrend(type='demean')
                    peak1=max(st[0].data)
                    g=0
                except:
                    g=1
            
                t1 = UTCDateTime(arrival2)-30 #the format is year:day_of_the_year:month
                t2 = t1 + 120
                st = Stream()
                try:
                    st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                    st.detrend(type='linear')
                    st.detrend(type='demean')
                    peak2=max(st[0].data)
                    f=0
                except:
                    f=1
            
            
                peak_percent = (peak2/peak1)*100
                all_count += 1
                all_percentage += peak_percent
                
                if  f==0 and g==0:
                    if peak1 > peak2 :
                        pands[event_num,1]=2
                        s_tot += 1
    #                    print(peak_percent)
                        count += 1
                        percentage += peak_percent
                    
    #                    try:
    #                        arrival=last_p
    #                        t1 = UTCDateTime(arrival)-120 #the format is year:day_of_the_year:month
    #                        t2 = t1 + 600
    #                        st = Stream()
    #                        st = client.get_waveforms(net, sta, '', cha, t1 , t2)
    #                        st.detrend(type='linear')
    #                        st.detrend(type='demean')
    #                        st.plot(color='b',starttime=t1, endtime=t2)
    #                    except:
    #                        h=1
                    else:
#                        pands[event_num,1]=1
#                        last_p=cat[x,0]
#                        p_tot +=1
                        larger_s_tot += 1
 #%%             
            if cat[x-1,10] + cat[x-1,11] + cat[x-1,12]==0 and cat[x-1,13]==1:
                sta = 'LB04' # STATION 
                arrival1=last_p
                arrival2=cat[x,0]
                
                t1 = UTCDateTime(arrival1)-30 #the format is year:day_of_the_year:month
                t2 = t1 + 120
                st = Stream()
                try:
                    st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                    st.detrend(type='linear')
                    st.detrend(type='demean')
                    peak1=max(st[0].data)
                    g=0
                except:
                    g=1
            
                t1 = UTCDateTime(arrival2)-30 #the format is year:day_of_the_year:month
                t2 = t1 + 120
                st = Stream()
                try:
                    st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                    st.detrend(type='linear')
                    st.detrend(type='demean')
                    peak2=max(st[0].data)
                    f=0
                except:
                    f=1
            
                peak_percent = (peak2/peak1)*100
                all_count += 1
                all_percentage += peak_percent
                
                if  f==0 and g==0:
                    if peak1 > peak2 :
                        pands[event_num,1]=2
                        s_tot += 1
    #                    print(peak_percent)
                        count += 1
                        percentage += peak_percent
                   
    #                    try:
    #                        arrival=last_p
    #                        t1 = UTCDateTime(arrival)-120 #the format is year:day_of_the_year:month
    #                        t2 = t1 + 600
    #                        st = Stream()
    #                        st = client.get_waveforms(net, sta, '', cha, t1 , t2)
    #                        st.detrend(type='linear')
    #                        st.detrend(type='demean')
    #                        st.plot(color='b',starttime=t1, endtime=t2)
    #                    except:
    #                        h=1
                    else:
#                        pands[event_num,1]=1
#                        last_p=cat[x,0]
#                        p_tot +=1
                        larger_s_tot += 1
#%%                  
            if cat[x-1,10] + cat[x-1,11] + cat[x-1,12] + cat[x-1,13]==0 and cat[x-1,14]==1:
                sta = 'LB05' # STATION 
                arrival1=last_p
                arrival2=cat[x,0]
                
                t1 = UTCDateTime(arrival1)-30 #the format is year:day_of_the_year:month
                t2 = t1 + 120
                st = Stream()
                try:
                    st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                    st.detrend(type='linear')
                    st.detrend(type='demean')
                    peak1=max(st[0].data)
                    g=0
                except:
                    g=1
            
                t1 = UTCDateTime(arrival2)-30 #the format is year:day_of_the_year:month
                t2 = t1 + 120
                st = Stream()
                try:
                    st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                    st.detrend(type='linear')
                    st.detrend(type='demean')
                    peak2=max(st[0].data)
                    f=0
                except:
                    f=1
            
            
                peak_percent = (peak2/peak1)*100
                all_count += 1
                all_percentage += peak_percent
                
                
                if  f==0 and g==0:
                    if peak1 > peak2:
                        pands[event_num,1]=2
                        s_tot += 1
    #                    print(peak_percent)
                        count += 1
                        percentage += peak_percent
                    
                        
    #                    try:
    #                        arrival=last_p
    #                        t1 = UTCDateTime(arrival)-120 #the format is year:day_of_the_year:month
    #                        t2 = t1 + 600
    #                        st = Stream()
    #                        st = client.get_waveforms(net, sta, '', cha, t1 , t2)
    #                        st.detrend(type='linear')
    #                        st.detrend(type='demean')
    #                        st.plot(color='b',starttime=t1, endtime=t2)
    #                    except:
    #                        h=1
                    else:
#                        pands[event_num,1]=1
#                        last_p=cat[x,0]
#                        p_tot +=1
                        larger_s_tot += 1
#%%                  
            if cat[x-1,10] + cat[x-1,11] + cat[x-1,12] + cat[x-1,13] + cat[x-1,14]==0 and cat[x-1,15]==1:
                sta = 'LB06' # STATION 
                arrival1=last_p
                arrival2=cat[x,0]
                
                t1 = UTCDateTime(arrival1)-30 #the format is year:day_of_the_year:month
                t2 = t1 + 120
                st = Stream()
                try:
                    st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                    st.detrend(type='linear')
                    st.detrend(type='demean')
                    peak1=max(st[0].data)
                    g=0
                except:
                    g=1
            
                t1 = UTCDateTime(arrival2)-30 #the format is year:day_of_the_year:month
                t2 = t1 + 120
                st = Stream()
                try:
                    st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                    st.detrend(type='linear')
                    st.detrend(type='demean')
                    peak2=max(st[0].data)
                    f=0
                except:
                    f=1
            
            
                peak_percent = (peak2/peak1)*100
                all_count += 1
                all_percentage += peak_percent
                
                if f==0 and g==0:
                    if peak1 > peak2 :
                        pands[event_num,1]=2
                        s_tot += 1
    #                    print(peak_percent)
                        count += 1
                        percentage += peak_percent
                    
    #                    try:
    #                        arrival=last_p
    #                        t1 = UTCDateTime(arrival)-120 #the format is year:day_of_the_year:month
    #                        t2 = t1 + 600
    #                        st = Stream()
    #                        st = client.get_waveforms(net, sta, '', cha, t1 , t2)
    #                        st.detrend(type='linear')
    #                        st.detrend(type='demean')
    #                        st.plot(color='b',starttime=t1, endtime=t2)
    #                    except:
    #                        h=1
                    else:
#                        pands[event_num,1]=1
#                        last_p=cat[x,0]
#                        p_tot +=1
                        larger_s_tot += 1
#%%                  
            if cat[x-1,10] + cat[x-1,11] + cat[x-1,12] + cat[x-1,13] + cat[x-1,14] + cat[x-1,15]==0 and cat[x-1,16]==1:
                sta = 'LS01' # STATION 
                arrival1=last_p
                arrival2=cat[x,0]
                
                t1 = UTCDateTime(arrival1)-30 #the format is year:day_of_the_year:month
                t2 = t1 + 120
                st = Stream()
                try:
                    st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                    st.detrend(type='linear')
                    st.detrend(type='demean')
                    peak1=max(st[0].data)
                    g=0
                except:
                    g=1
            
                t1 = UTCDateTime(arrival2)-30 #the format is year:day_of_the_year:month
                t2 = t1 + 120
                st = Stream()
                try:
                    st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                    st.detrend(type='linear')
                    st.detrend(type='demean')
                    peak2=max(st[0].data)
                    f=0
                except:
                    f=1
            
                peak_percent = (peak2/peak1)*100
                all_count += 1
                all_percentage += peak_percent
                
                if  f==0 and g==0:
                    if peak1 > peak2:
                        pands[event_num,1]=2
                        s_tot += 1
    #                    print(peak_percent)
                        count += 1
                        percentage += peak_percent
                    
    #                    try:
    #                        arrival=last_p
    #                        t1 = UTCDateTime(arrival)-120 #the format is year:day_of_the_year:month
    #                        t2 = t1 + 600
    #                        st = Stream()
    #                        st = client.get_waveforms(net, sta, '', cha, t1 , t2)
    #                        st.detrend(type='linear')
    #                        st.detrend(type='demean')
    #                        st.plot(color='b',starttime=t1, endtime=t2)
    #                    except:
    #                        h=1
                    else:
#                        pands[event_num,1]=1
#                        last_p=cat[x,0]
#                        p_tot +=1
                        larger_s_tot += 1
#%%                  
            if cat[x-1,10] + cat[x-1,11] + cat[x-1,12] + cat[x-1,13] + cat[x-1,14] + cat[x-1,15] + cat[x-1,16]==0 and cat[x-1,17]==1:
                sta = 'LS02' # STATION 
                arrival1=last_p
                arrival2=cat[x,0]
                
                t1 = UTCDateTime(arrival1)-30 #the format is year:day_of_the_year:month
                t2 = t1 + 120
                st = Stream()
                try:
                    st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                    st.detrend(type='linear')
                    st.detrend(type='demean')
                    peak1=max(st[0].data)
                    g=0
                except:
                    g=1
                
            
                t1 = UTCDateTime(arrival2)-30 #the format is year:day_of_the_year:month
                t2 = t1 + 120
                st = Stream()
                try:
                    st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                    st.detrend(type='linear')
                    st.detrend(type='demean')
                    peak2=max(st[0].data)
                    f=0
                except:
                    f=1
            
            
                peak_percent = (peak2/peak1)*100
                all_count += 1
                all_percentage += peak_percent
                
                if f==0 and g==0: 
                    if peak1 > peak2:
                        pands[event_num,1]=2
                        s_tot += 1
    #                    print(peak_percent)
                        count += 1
                        percentage += peak_percent
                
    #                    try:
    #                        arrival=last_p
    #                        t1 = UTCDateTime(arrival)-120 #the format is year:day_of_the_year:month
    #                        t2 = t1 + 600
    #                        st = Stream()
    #                        st = client.get_waveforms(net, sta, '', cha, t1 , t2)
    #                        st.detrend(type='linear')
    #                        st.detrend(type='demean')
    #                        st.plot(color='b',starttime=t1, endtime=t2)
    #                    except:
    #                        h=1
                    else:
#                        pands[event_num,1]=1
#                        last_p=cat[x,0]
#                        p_tot +=1
                        larger_s_tot += 1
#%%                
                
            if cat[x-1,10] + cat[x-1,11] + cat[x-1,12] + cat[x-1,13] + cat[x-1,14] + cat[x-1,15] + cat[x-1,16] + cat[x-1,17]==0 and cat[x-1,18]==1:
                sta = 'LS03' # STATION 
                arrival1=last_p
                arrival2=cat[x,0]
                
                t1 = UTCDateTime(arrival1)-30 #the format is year:day_of_the_year:month
                t2 = t1 + 120
                st = Stream()
                try:
                    st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                    st.detrend(type='linear')
                    st.detrend(type='demean')
                    peak1=max(st[0].data)
                    g=0
                except:
                    g=1
            
                t1 = UTCDateTime(arrival2)-30 #the format is year:day_of_the_year:month
                t2 = t1 + 120
                st = Stream()
                try:
                    st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                    st.detrend(type='linear')
                    st.detrend(type='demean')
                    peak2=max(st[0].data)
                    f=0
                except:
                    f=1
            
            
                peak_percent = (peak2/peak1)*100
                all_count += 1
                all_percentage += peak_percent
                
                if  f==0 and g==0:
                    if peak1 > peak2:
                        pands[event_num,1]=2
                        s_tot += 1
    #                    print(peak_percent)
                        count += 1
                        percentage += peak_percent
                    
    #                    try:
    #                        arrival=last_p
    #                        t1 = UTCDateTime(arrival)-120 #the format is year:day_of_the_year:month
    #                        t2 = t1 + 600
    #                        st = Stream()
    #                        st = client.get_waveforms(net, sta, '', cha, t1 , t2)
    #                        st.detrend(type='linear')
    #                        st.detrend(type='demean')
    #                        st.plot(color='b',starttime=t1, endtime=t2)
    #                    except:
    #                        h=1
                    else:
#                        pands[event_num,1]=1
#                        last_p=cat[x,0]
#                        p_tot +=1
                        larger_s_tot += 1
 #%%           
            if cat[x-1,10] + cat[x-1,11] + cat[x-1,12] + cat[x-1,13] + cat[x-1,14] + cat[x-1,15] + cat[x-1,16] + cat[x-1,17] + cat[x-1,18]==0 and cat[x-1,19]==1:                                                        
                sta = 'LS04' # STATION 
                arrival1=last_p
                arrival2=cat[x,0]
                
                t1 = UTCDateTime(arrival1)-30 #the format is year:day_of_the_year:month
                t2 = t1 + 120
                st = Stream()
                try:
                    st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                    st.detrend(type='linear')
                    st.detrend(type='demean')
                    peak1=max(st[0].data)
                    g=0
                except:
                    g=1
                
            
                t1 = UTCDateTime(arrival2)-30 #the format is year:day_of_the_year:month
                t2 = t1 + 120
                st = Stream()
                try:
                    st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                    st.detrend(type='linear')
                    st.detrend(type='demean')
                    peak2=max(st[0].data)
                    f=0
                except:
                    f=1
            
                peak_percent = (peak2/peak1)*100
                all_count += 1
                all_percentage += peak_percent
                
                if  f==0 and g==0:
                    if peak1 > peak2 :
                        pands[event_num,1]=2
                        s_tot += 1
    #                    print(peak_percent)
                        count += 1
                        percentage += peak_percent
                    
    #                    try:
    #                        arrival=last_p
    #                        t1 = UTCDateTime(arrival)-120 #the format is year:day_of_the_year:month
    #                        t2 = t1 + 600
    #                        st = Stream()
    #                        st = client.get_waveforms(net, sta, '', cha, t1 , t2)
    #                        st.detrend(type='linear')
    #                        st.detrend(type='demean')
    #                        st.plot(color='b',starttime=t1, endtime=t2)
    #                    except:
    #                        h=1
                    else:
#                        pands[event_num,1]=1
#                        last_p=cat[x,0]
#                        p_tot +=1
                        larger_s_tot += 1
#%%                
            if cat[x-1,10] + cat[x-1,11] + cat[x-1,12] + cat[x-1,13] + cat[x-1,14] + cat[x-1,15] + cat[x-1,16] + cat[x-1,17] + cat[x-1,18] + cat[x-1,19]==0 and cat[x-1,20]==1:     
                sta = 'LS05' # STATION 
                arrival1=last_p
                arrival2=cat[x,0]
                
                t1 = UTCDateTime(arrival1)-30 #the format is year:day_of_the_year:month
                t2 = t1 + 120
                st = Stream()
                try:
                    st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                    st.detrend(type='linear')
                    st.detrend(type='demean')
                    peak1=max(st[0].data)
                    g=0
                except:
                    g=1
            
                t1 = UTCDateTime(arrival2)-30 #the format is year:day_of_the_year:month
                t2 = t1 + 120
                st = Stream()
                try:
                    st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                    st.detrend(type='linear')
                    st.detrend(type='demean')
                    peak2=max(st[0].data)
                    f=0
                except:
                    f=1
            
            
                peak_percent = (peak2/peak1)*100
                all_count += 1
                all_percentage += peak_percent
                
                if f==0 and g==0:
                    if peak1 > peak2 :
                        pands[event_num,1]=2
                        s_tot += 1
    #                    print(peak_percent)
                        count += 1
                        percentage += peak_percent
                    
    #                    try:
    #                        arrival=last_p
    #                        t1 = UTCDateTime(arrival)-120 #the format is year:day_of_the_year:month
    #                        t2 = t1 + 600
    #                        st = Stream()
    #                        st = client.get_waveforms(net, sta, '', cha, t1 , t2)
    #                        st.detrend(type='linear')
    #                        st.detrend(type='demean')
    #                        st.plot(color='b',starttime=t1, endtime=t2)
    #                    except:
    #                        h=1
                    else:
#                        pands[event_num,1]=1
#                        last_p=cat[x,0]
#                        p_tot +=1
                        larger_s_tot += 1
 #%%               
            if cat[x-1,10] + cat[x-1,11] + cat[x-1,12] + cat[x-1,13] + cat[x-1,14] + cat[x-1,15] + cat[x-1,16] + cat[x-1,17] + cat[x-1,18] + cat[x-1,19] + cat[x-1,20] ==0 and cat[x-1,21]==1:     
                sta = 'LS06' # STATION 
                arrival1=last_p
                arrival2=cat[x,0]
                
                t1 = UTCDateTime(arrival1)-30 #the format is year:day_of_the_year:month
                t2 = t1 + 120
                st = Stream()
                try:
                    st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                    st.detrend(type='linear')
                    st.detrend(type='demean')
                    peak1=max(st[0].data)
                    g=0
                except:
                    g=1
            
                t1 = UTCDateTime(arrival2)-30 #the format is year:day_of_the_year:month
                t2 = t1 + 120
                st = Stream()
                try:
                    st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                    st.detrend(type='linear')
                    st.detrend(type='demean')
                    peak2=max(st[0].data)
                    f=0
                except:
                    f=1
            
            
                peak_percent = (peak2/peak1)*100
                all_count += 1
                all_percentage += peak_percent
                
                if  f==0 and g==0:
                    if peak1 > peak2:
                        pands[event_num,1]=2
                        s_tot += 1
    #                    print(peak_percent)
                        count += 1
                        percentage += peak_percent
                    
    #                    try:
    #                        arrival=last_p
    #                        t1 = UTCDateTime(arrival)-120 #the format is year:day_of_the_year:month
    #                        t2 = t1 + 600
    #                        st = Stream()
    #                        st = client.get_waveforms(net, sta, '', cha, t1 , t2)
    #                        st.detrend(type='linear')
    #                        st.detrend(type='demean')
    #                        st.plot(color='b',starttime=t1, endtime=t2)
    #                    except:
    #                        h=1
                    else:
#                        pands[event_num,1]=1
#                        last_p=cat[x,0]
#                        p_tot +=1
                        larger_s_tot += 1
            num +=1
#%%
            
        else:
            pands[event_num,1]=1
            last_p=cat[x,0]
            p_tot +=1
        pands[event_num,2]=(cat[x,0]-cat[x-1,0] )/60 
        if x < len(cat)-1:   
            pands = np.lib.pad(pands, ((0,1),(0,0)), 'constant', constant_values=(0))
        event_num += 1
        
#np.savetxt("/Users/william/Documents/scanner/all_stations/primary_secondary_v1.csv", pands,delimiter=",",header="time_stanp,p_or_s,repose time")


print('number of Events =',s_tot + p_tot + larger_s_tot)
print('number of Primary EXPs =', p_tot)
print("Number of smaller secondary EXP's =",s_tot)
print("number of larger secondaries =" , larger_s_tot )

av_perc = percentage/count
print("average secondary %age of primary amplitude = ", av_perc)

av_perc_all = all_percentage/all_count
print("average secondary %age of primary amplitude for all secondaries = ", av_perc_all)


#pony = 0
#for x in range(0,len(pands)-1):
#    if pands[x,1] == 1 and pands[x+1,1] == 1:
#        pony +=1
#        
#if pands[-1,1]==1:
#    pony+=1
#
#print('number of P only Events =', pony)
#print('number of Events with secondaries =', p_tot - pony)
    
    
    
    
    
    
    
    
