#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 12:38:32 2018

@author: william
"""
import obspy
from obspy import read
import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
from obspy import Stream
from obspy.signal.trigger import classic_sta_lta, recursive_sta_lta
from obspy.signal.trigger import plot_trigger, trigger_onset
from obspy import UTCDateTime
from obspy.clients.earthworm import Client

year1=2014
month1=11
day1=24
hour1=0
minute1=0
second1=0

r1=4630
r2=3370
r3=2310
r4=1300
r5=810
r6=7660


rs1=4430
rs2=3190
rs3=2270
rs4=950
rs5=1100
rs6=7510


sso=340
sse = 450
sott=r1/sso
sett=rs1/sse
extd = sott-sett

td=np.zeros(shape=(0,4))
num=0



sta = 'LB01' # STATION 
net = 'Z4'  # 
loc = ''    # location, it depends mostly of which network you are in. 

client = Client('138.253.112.23', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23

def s_inf(last_p,s_arrival,station,next_p):
    
        
    cha = 'HDF' # CHANNEL
    t1 = UTCDateTime(last_p) 
    t2 = t1 + 40
    st = Stream()
    st = client.get_waveforms(net, sta, '', cha, t1-10 , t2)
    st.detrend(type='linear')
    st.detrend(type='demean')
    tr_inf=st[0]
    tr_inf.filter("bandpass", freqmin=0.2,freqmax=2)
    sr = tr_inf.stats.sampling_rate
    nsta=int(0.5*sr)                                      #2
    nlta=int(15*sr)                                     #20
    stream=tr_inf.data
    cft=recursive_sta_lta(stream, nsta, nlta)
    trig_on=10                                           #8
    trig_off=1                                        #0.2
#    plot_trigger(trs, cft, trig_on, trig_off) 

    on_off = trigger_onset(cft,trig_on,trig_off)
    
    if len(on_off) > 0:
#        tr_inf.plot(color='r',starttime=t1, endtime=t2)
   
    
        cha = 'HHZ' # CHANNEL
        t1 = UTCDateTime(last_p) 
        t2 = t1 + 120
        st = Stream()
        st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
        st.detrend(type='linear')
        st.detrend(type='demean')
        trs=st[0]
        trs.filter("bandpass", freqmin=0.2,freqmax=6)
        
        sr = trs.stats.sampling_rate
        nsta=int(1*sr)                                      #2
        nlta=int(20*sr)                                     #20
        stream=trs.data
        cft=recursive_sta_lta(stream, nsta, nlta)
        trig_on=5                                            #8
        trig_off=1                                        #0.2
#        plot_trigger(trs, cft, trig_on, trig_off) 
    
        on_off2 = trigger_onset(cft,trig_on,trig_off)
        
#        if len(on_off2) > 0:
#            trs.plot(color='r',starttime=t1, endtime=t2)
#            print(arrival_s)                
#                diff= arrival_s - arrival
#                print(diff)

#%%
    
    cha = 'HDF' # CHANNEL
    t1 = UTCDateTime(s_arrival) 
    t2 = t1 + 40
    st = Stream()
    st = client.get_waveforms(net, sta, '', cha, t1-10 , t2)
    st.detrend(type='linear')
    st.detrend(type='demean')
    tr_inf=st[0]
    tr_inf.filter("bandpass", freqmin=0.2,freqmax=2)
    sr = tr_inf.stats.sampling_rate
    nsta=int(0.5*sr)                                      #2
    nlta=int(15*sr)                                     #20
    stream=tr_inf.data
    cft=recursive_sta_lta(stream, nsta, nlta)
    trig_on=10                                           #8
    trig_off=1                                        #0.2
#    plot_trigger(trs, cft, trig_on, trig_off) 

    on_off3 = trigger_onset(cft,trig_on,trig_off)
    
    if len(on_off3) > 0:
#        tr_inf.plot(color='b',starttime=t1, endtime=t2)
   
    
        cha = 'HHZ' # CHANNEL
        t1 = UTCDateTime(s_arrival) 
        t2 = t1 + 120
        st = Stream()
        st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
        st.detrend(type='linear')
        st.detrend(type='demean')
        trs=st[0]
        trs.filter("bandpass", freqmin=0.2,freqmax=6)
        
        sr = trs.stats.sampling_rate
        nsta=int(1*sr)                                      #2
        nlta=int(20*sr)                                     #20
        stream=trs.data
        cft=recursive_sta_lta(stream, nsta, nlta)
        trig_on=5                                            #8
        trig_off=1                                        #0.2
#        plot_trigger(trs, cft, trig_on, trig_off) 
    
        on_off4 = trigger_onset(cft,trig_on,trig_off)
        
#        if len(on_off4) > 0:
#            trs.plot(color='b',starttime=t1, endtime=t2)
#            print(arrival_s)                
#                diff= arrival_s - arrival
#                print(diff)


#%%
    cha = 'HDF' # CHANNEL
    t5 = UTCDateTime(next_p) 
    t6 = t5 + 40
    st = Stream()
    st = client.get_waveforms(net, sta, '', cha, t5-10 , t6)
    st.detrend(type='linear')
    st.detrend(type='demean')
    tr_inf=st[0]
    tr_inf.filter("bandpass", freqmin=0.2,freqmax=2)
    sr = tr_inf.stats.sampling_rate
    nsta=int(0.5*sr)                                      #2
    nlta=int(15*sr)                                     #20
    stream=tr_inf.data
    cft=recursive_sta_lta(stream, nsta, nlta)
    trig_on=10                                           #8
    trig_off=1                                        #0.2
#    plot_trigger(trs, cft, trig_on, trig_off) 

    on_off5 = trigger_onset(cft,trig_on,trig_off)
    
    if len(on_off) > 0:
#        tr_inf.plot(color='r',starttime=t1, endtime=t2)
   
    
        cha = 'HHZ' # CHANNEL
        t5 = UTCDateTime(last_p) 
        t6 = t5 + 120
        st = Stream()
        st = client.get_waveforms(net, sta, '', cha, t5-20 , t6)
        st.detrend(type='linear')
        st.detrend(type='demean')
        trs=st[0]
        trs.filter("bandpass", freqmin=0.2,freqmax=6)
        
        sr = trs.stats.sampling_rate
        nsta=int(1*sr)                                      #2
        nlta=int(20*sr)                                     #20
        stream=trs.data
        cft=recursive_sta_lta(stream, nsta, nlta)
        trig_on=5                                            #8
        trig_off=1                                        #0.2
#        plot_trigger(trs, cft, trig_on, trig_off) 
    
        on_off6 = trigger_onset(cft,trig_on,trig_off)


#%%            
    all_4=np.zeros(shape=(1,6))
    
    if len(on_off) > 0 and len(on_off2) > 0 and len(on_off3) > 0 and len(on_off4) > 0:
        all_4 = np.lib.pad(all_4, ((0,1),(0,0)), 'constant', constant_values=(0))
        all_4[0][0]=last_p + on_off[0,0]/100 #Primary infrasound
        all_4[0][1]=last_p #+ on_off2[0,0]/100 #Primary seismic
        all_4[0][2]=s_arrival + on_off3[0,0]/100 #secondary infrasound
        all_4[0][3]=s_arrival #+ on_off4[0,0]/100 #secondary seismic
        all_4[0][4]=next_p + on_off5[0,0]/100 #secondary infrasound
        all_4[0][5]=next_p #+ on_off4[0,0]/100 #secondary seismic
        
        

    return(tr_inf,all_4)