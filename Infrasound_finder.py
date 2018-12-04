#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 11:38:50 2018

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


sso=330
sse = 500
sott=r1/sso
sett=rs1/sse
extd = sott-sett

td=np.zeros(shape=(0,4))
num=0

sta = 'LB01' # STATION 
net = 'Z4'  # 
loc = ''    # location, it depends mostly of which network you are in. 

client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23


cat= genfromtxt("/Users/william/Documents/scanner/all_stations/Explosion_Catalogue_V4c.csv", delimiter=',',skip_header=1)


for x in range(0,len(cat)):
#for x in range(0,500):    
    try:
        event=cat[x,0]
        
    
        cha = 'HDF' # CHANNEL
        t1 = UTCDateTime(event) 
        t2 = t1 + 40
        st = Stream()
        st = client.get_waveforms(net, sta, '', cha, t1-10 , t2)
        st.detrend(type='linear')
        st.detrend(type='demean')
        trs=st[0]
        trs.filter("bandpass", freqmin=0.2,freqmax=2)
        sr = trs.stats.sampling_rate
        nsta=int(0.5*sr)                                      #2
        nlta=int(15*sr)                                     #20
        stream=trs.data
        cft=recursive_sta_lta(stream, nsta, nlta)
        trig_on=10                                           #8
        trig_off=1                                        #0.2
    #    plot_trigger(trs, cft, trig_on, trig_off) 
    
        on_off = trigger_onset(cft,trig_on,trig_off)
        
        if len(on_off) > 0:
            start1=(on_off[0])/100
            start=start1[1]
            arrival = t1 + start
    #        print(t1)
    #        print(arrival)
    #    
        
            cha = 'HHZ' # CHANNEL
            t1 = UTCDateTime(event) 
            t2 = t1 + 120
            st = Stream()
            st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
            st.detrend(type='linear')
            st.detrend(type='demean')
            trs=st[0]
            
            sr = trs.stats.sampling_rate
            nsta=int(1*sr)                                      #2
            nlta=int(20*sr)                                     #20
            stream=trs.data
            cft=recursive_sta_lta(stream, nsta, nlta)
            trig_on=5                                            #8
            trig_off=1                                        #0.2
    #        plot_trigger(trs, cft, trig_on, trig_off) 
        
            on_off2 = trigger_onset(cft,trig_on,trig_off)
            
            if len(on_off) > 0:
                start2=(on_off2[0])/100
                start_s=start2[1]
                arrival_s = t1 + start_s
                
    #            print(arrival_s)
                
                diff= arrival_s - arrival
                
#                print(diff)
                
                if 2 < diff < 15:
                    td = np.lib.pad(td, ((0,1),(0,0)), 'constant', constant_values=(0))
                            
                    td[num][0]=arrival_s
                    td[num][1]=arrival
                    td[num][2]=diff
                    td[num][3]=diff-extd
                    num+=1
                
            
    except:
        p=1
        



plt.figure(1)
plt.plot(td[:,2])
plt.title('arrival difference')

plt.figure(2)
plt.plot(td[:,3])
plt.title('source difference')
#%%
numb=0
tds=np.zeros(shape=(0,2))
for q in range(0,len(td)-20):
    tds = np.lib.pad(tds, ((0,1),(0,0)), 'constant', constant_values=(0))
    
    smoothed=sum(td[q:q+20,2])/20  
    tds[numb][1]=td[q+10,0]          
    tds[numb][0]=smoothed
    numb+=1
    
    
numb2=0   
tds2=np.zeros(shape=(0,2))
for q in range(0,len(td)-20):
    tds2 = np.lib.pad(tds2, ((0,1),(0,0)), 'constant', constant_values=(0))
    
    smoothed2=sum(td[q:q+20,3])/20  
    tds2[numb2][1]=td[q+10,0]          
    tds2[numb2][0]=smoothed2
    numb2+=1
    
plt.figure(3)
plt.plot(tds[:,0])
plt.title('smoothed arrival difference')
#plt.ylim([0,10])

plt.figure(4)
plt.plot(tds2[:,0])
plt.title('smoothed source difference')
#plt.ylim([0,10])

#np.savetxt("/Users/william/Documents/scanner/analysis/acoustic_delay_LB01.csv", td,delimiter=",",header="seismic_time,Acoustic_time,time_delay,adjusted_delay")
#np.savetxt("/Users/william/Documents/scanner/analysis/acoustic_delay_smoothed_LB01.csv", tds,delimiter=",",header="smoothed_time,seismic_arrival")




















    