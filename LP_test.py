#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 13 08:57:55 2018

@author: william
"""


from obspy.clients.earthworm import Client
from obspy import UTCDateTime
from obspy import Stream
from numpy import genfromtxt
from obspy.signal.trigger import recursive_sta_lta
from obspy.signal.trigger import plot_trigger, trigger_onset
import numpy as np

cat= genfromtxt("/Users/william/Documents/scanner/all_stations/Explosion_Catalogue_V4c.csv", delimiter=',',skip_header=1)

sta = 'LB05' # STATION 
cha = 'HHZ' # CHANNEL
net = 'Z4'  # 
loc = ''    # location, it depends mostly of which network you are in. 
client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23

sr = 100
nsta=int(1*sr)                                      
nlta=int(20*sr)                                   
trig_on=7.5
trig_off=0.2 
event = 0
on_off_save=np.zeros(shape=(0,4))
no_data = 0
shift = 15000
crit=0.5


#%%
t1_ref = UTCDateTime(1459324380) 
t2_ref = t1_ref +120
st_ref = Stream()
st_ref = client.get_waveforms(net, 'LB03', '', 'HHZ', t1_ref , t2_ref)

st_ref.detrend(type='linear')
st_ref.detrend(type='demean')
tr_ref=st_ref[0].slice(starttime=t1_ref, endtime=t2_ref)

tr_ref.detrend(type='linear')
tr_ref.detrend(type='demean')
tr_ref.filter(type='bandpass',freqmin=0.001, freqmax=0.1)

st_c_ref = calibrate(tr_ref)
st_c_ref.plot(color='g',starttime=t1_ref, endtime=t1_ref+120)


#%%
data_stream_unfilt = Stream()
data_stream_filt = Stream()

for x in range(0,len(cat)):
#    if 1459042900 < cat[x,0] < 1477791497 : # if  1416791640 < cat[x,0] < 1418430749 :  #  for getting LP's in mid time
        tim = UTCDateTime(cat[x,0])
        
        try:
            t1 = tim - 120 #the format is year:day_of_the_year:month
            t2 = t1 + 5*60
            st = Stream()
            st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)

            st.detrend(type='linear')
            st.detrend(type='demean')
            
            tr=st[0].slice(starttime=t1, endtime=t2)

            tr.filter(type='bandpass',freqmin=0.001, freqmax=0.1)
            tr.detrend(type='linear')
            tr.detrend(type='demean')
            tr_copy = tr.copy()
            st_c = calibrate(tr)
            st1 = calibrate(st[0])
                        
            cft=recursive_sta_lta(tr, nsta, nlta)
            on_off = trigger_onset(cft,trig_on,trig_off)
            
            st1_cut=st1.trim(starttime=t1+(on_off[0,0]/100)-30, endtime=t1+(on_off[0,0]/100)+90)
            st_c_cut=st_c.trim(starttime=t1+(on_off[0,0]/100)-30, endtime=t1+(on_off[0,0]/100)+90)

            if on_off[0,0] > 30 :
#                plot_trigger(tr_copy, cft, trig_on, trig_off)
                top_v,top,corell = corel(st_c_ref[0],st_c[0],shift)
                 
                if abs(top_v) > crit:
                
                    on_off_save = np.lib.pad(on_off_save, ((0,1),(0,0)), 'constant', constant_values=(0))
                    on_off_save[event][0] = on_off[0,0]
                    on_off_save[event][1] = tim
                    on_off_save[event][2] = t1+(on_off[0,0]/100)-30
                    on_off_save[event][3] = t1+(on_off[0,0]/100)+90
                    event += 1
                                        
                    data_stream_unfilt.append(st1_cut[0])
                    data_stream_filt.append(st_c_cut[0])
                    
        except:
            no_data +=1
            
            
for x in range(0,len(data_stream_filt)):  
    data_stream_unfilt[x].plot(color='b')
    data_stream_filt[x].plot(color='k')
            
            
            
            
            
            
            
            
            
            
            
            