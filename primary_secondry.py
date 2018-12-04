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
from obspy.signal.trigger import classic_sta_lta, recursive_sta_lta
from obspy.signal.trigger import plot_trigger, trigger_onset
from obspy import UTCDateTime
from matplotlib import pyplot
import scipy
from obspy import UTCDateTime

sta = 'LB01' # STATION 
cha = 'HHZ' # CHANNEL
net = 'Z4'  # 
loc = ''    # location, it depends mostly of which network you are in. 
client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23


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

#for x in range(1,200):
for x in range(1,len(cat)):
#    if 1459468800 < cat[x,0] < 1475884800:
        pands[event_num,0]=cat[x,0]
        
        if cat[x,0] - last_p < 10*60:
            pands[event_num,1]=2
            s_tot += 1
            arrival=cat[x,0]
            t1 = last_p - 30
            t2=t1+600
            st = Stream()
            try:
                st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
                st.detrend(type='linear')
                st.detrend(type='demean')
                st.plot(color='b',starttime=t1, endtime=t2)
            except:
                f=1
        else:
            pands[event_num,1]=1
            last_p=cat[x,0]
            p_tot +=1
        pands[event_num,2]=(cat[x,0]-cat[x-1,0] )/60 
        if x < len(cat)-1:   
            pands = np.lib.pad(pands, ((0,1),(0,0)), 'constant', constant_values=(0))
        event_num += 1
        
#np.savetxt("/Users/william/Documents/scanner/all_stations/primary_secondary_v1.csv", pands,delimiter=",",header="time_stanp,p_or_s,repose time")


print('number of Events =',s_tot + p_tot)
print('number of Primary EXPs =', p_tot)
print("Number of secondary EXP's =",s_tot)


pony = 0
for x in range(0,len(pands)-1):
    if pands[x,1] == 1 and pands[x+1,1] == 1:
        pony +=1
        
if pands[-1,1]==1:
    pony+=1

print('number of P only Events =', pony)
print('number of Events with secondaries =', p_tot - pony)
    
    
    
    
    
    
    
    
