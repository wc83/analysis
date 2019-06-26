#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 15:02:56 2019

@author: william
"""


import numpy as np
import obspy
import matplotlib.pyplot as plt 
from obspy.clients.earthworm import Client
from obspy import UTCDateTime
from obspy.signal.trigger import classic_sta_lta, recursive_sta_lta
from obspy.signal.trigger import plot_trigger, trigger_onset
from obspy import Stream
from numpy import genfromtxt

cat = genfromtxt("/Users/william/Documents/scanner/all_stations/Final_Catalogue_2014_2018.csv", delimiter=',',skip_header=1)


sta = 'LB01' # STATION 
cha1 = 'HDF' # Acoustic CHANNEL
cha2 = 'HHZ' # Vertical Seismic Channel
net = 'Z4'  # 
loc = ''    # location, it depends mostly of which network you are in. 
client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23

t_diff=np.zeros(shape=(0,5))
num=0
for x in range(0,len(cat)) : 
    
    if cat[x,10] == 1 and cat[x,3] == 2016 and 4 < cat[x,4] < 10:
        try:
            # Seismics
            # read in data
            t1 = UTCDateTime(cat[x,0]-30) #the format is year:day_of_the_year:month
            t2 = t1 + 90           
            st2 = Stream()
            st2 = client.get_waveforms(net, sta, '', cha2, t1 , t2)
#            print(t1)
            # Filter data
            st2.detrend(type='linear')
            st2.detrend(type='demean')
            st2.filter(type='bandpass',freqmin=0.5, freqmax=10)
            trs=st2[0]
            
            # sta_lta for Seismics
            sr = trs.stats.sampling_rate
            ssta=int(3.5 *sr)                                      #2
            slta=int(30*sr)
            
            streams=trs.data
            cfts=recursive_sta_lta(streams, ssta, slta)
            trig_on=2                                          #8
            trig_off=0.5                                        #0.2
#            plot_trigger(trs, cfts, trig_on, trig_off) 
            
            on_off_s = trigger_onset(cfts,trig_on,trig_off)
            
            #Seismic energy
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
            
            event_ss = trs.stats.starttime + (on_off_s[0,0]/100)
            trs2 = trs.slice(starttime = event_ss - 10, endtime = event_ss + 40 )
            
            st_c = calibrate1(trs2)
            B=2*pi*rhoE*cE*(1/A)
            
            EI = sum(np.square(st_c[0].data))
            EE= B*(r1*r1)*EI
#            print(EE)
            # Infrasound
            
            # Read in data
            
            t1a = UTCDateTime(trs.stats.starttime + (on_off_s[0,0]/100) -10) #the format is year:day_of_the_year:month
            t2a = t1a + 40
            
            st1 = Stream()
            st1 = client.get_waveforms(net, sta, '', cha1, t1a , t2a)
            
            # filter data
            st1.detrend(type='linear')
            st1.detrend(type='demean')
            st1.filter(type='bandpass',freqmin=0.2, freqmax=2)
            tra=st1[0]
            
            # sta_lta for acoustic
            sr = tra.stats.sampling_rate
            asta=int(0.2*sr)                                      #2
            alta=int(12*sr)                                     #20
            
            streama=tra.data
            cfta=recursive_sta_lta(streama, asta, alta)
            trig_on=15                                           #8
            trig_off=1                                        #0.2
#            plot_trigger(tra, cfta, trig_on, trig_off) 
            
            on_off_a = trigger_onset(cfta,trig_on,trig_off)
#            print(on_off_a[1]) 
                
            
#            if len(on_off_a) > 0 and len(on_off_s) > 0 :
                
            event_as = tra.stats.starttime + (on_off_a[0,0]/100) 
            print(event_as)
            
            time_diff = event_as - event_ss
            print(time_diff)
            print(EE)
            if 18 > time_diff > 2:

                t_diff = np.lib.pad(t_diff, ((0,1),(0,0)), 'constant', constant_values=(0))
                t_diff[num][0]=trs.stats.starttime
                t_diff[num][1]=time_diff
                t_diff[num][2]=EE
                t_diff[num][3]=event_ss
                t_diff[num][4]=event_as
                
                num +=1

        except:
            p=1