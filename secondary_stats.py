#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 18:05:51 2019

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


#%%

secondary= genfromtxt("//Users/william/Documents/scanner/analysis/secondaries_list.csv", delimiter=',',skip_header=1)
corr_list= genfromtxt("/Users/william/Documents/scanner/analysis/LB01_inf_correltion_list_6s_v2.csv", delimiter=',',skip_header=1)
event_times= genfromtxt("/Users/william/Documents/scanner/analysis/LB01_inf_correltion_event_details_6s_v2.csv", delimiter=',',skip_header=1)
corr_m= genfromtxt("/Users/william/Documents/scanner/analysis/LB01_inf_correltion_matrix_6s_v2.csv", delimiter=',',skip_header=1)


#%%

# STATION, CHANNEL (DDF --> 400 Hz), NETWWORK AND LOCATION CODES 
sta = 'LB01' # STATION 
cha = 'HDF' # CHANNEL
net = 'Z4'  # 
loc = ''    # location, it depends mostly of which network you are in.      
client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23

for x in range(0,20):
        
    t1 = UTCDateTime(secondary[x][1]) #the format is year:day_of_the_year:month
    t2 = t1 + 6
    st = Stream()
    st = client.get_waveforms(net, sta, '', cha, t1-20 , t2+60)
    #print(st)
    
    
    st.detrend(type='linear')
    st.detrend(type='demean')
    st.filter(type='bandpass',freqmin=0.2, freqmax=2)
    st.plot(color='b',starttime=t1-15, endtime=t2+30)


#%%

count = 0

#for x in range(0,100):
for x in range(0,len(secondary)):   
    primary = secondary[x][0]
    second = secondary[x][1]
    
    nearp,indp=find_nearest(event_times[:], primary)
    nears,inds=find_nearest(event_times[:], second)
    
    if abs(primary - nearp) < 60:
        if abs(second-nears) < 60:
            print(nearp)
            count += 1
        else:
            print("no secondary in list")
            
    else:   
        print("no primary in list")
        
        
print(len(secondary))
print(count)            
            
            
 #%%

num =0

for y in range(1,len(event_times)):
    if event_times[y] - event_times[y-1] < 600:
        num += 1
        
          
print(num)          



#%%


last_p = event_times[0]
prim_i = 0

numb_s=0
numb_a=0
ps_corr=np.zeros(shape=(0,1))
a_corr=np.zeros(shape=(0,1))

cot = 0
for x in range(1,len(event_times)):
#for x in range(1,250):
    print('x=',x)
    for y in range(0,x):
        cot += 1
        if x - y < 3:
            
            if 0 < event_times[x]-event_times[y] < 600:
                if event_times[x]-event_times[y-1] > 600:
                    ps_corr = np.lib.pad(ps_corr, ((0,1),(0,0)), 'constant', constant_values=(0))
                    ps_corr[numb_s] = corr_m[x][y]
                    numb_s +=1
                else:
                    if event_times[x]-event_times[y-2] > 600:
                        ps_corr = np.lib.pad(ps_corr, ((0,1),(0,0)), 'constant', constant_values=(0))
                        ps_corr[numb_s] = corr_m[x][y-1]
                        numb_s +=1
            else:
                a_corr = np.lib.pad(a_corr, ((0,1),(0,0)), 'constant', constant_values=(0))
                a_corr[numb_a] = corr_m[x][y]
                numb_a +=1
        else:
            a_corr = np.lib.pad(a_corr, ((0,1),(0,0)), 'constant', constant_values=(0))
            a_corr[numb_a] = corr_m[x][y]
            numb_a +=1




#%%





s_mean = np.mean(ps_corr)
s_var = np.var(ps_corr)
s_n = len(ps_corr)

a_mean = np.mean(a_corr)
a_var = np.var(a_corr)
a_n = len(a_corr)



t_calc = (s_mean - a_mean)/np.sqrt((s_var/s_n)+(a_var/a_n))

print("t_calc (any pair to secondary) = ", t_calc) # = 6.800 (3dp) - significant to 0.1%

s_std = np.std(ps_corr)
a_std = np.std(a_corr)

print('secondary mean = ', s_mean)
print('normal mean = ', a_mean)

print('secondary variance = ', s_var)
print('normal variance = ', a_var)

print('secondary standard deviation = ', s_std)
print('normal standard deviation = ', a_std)


#%%


np.savetxt("/Users/william/Documents/scanner/analysis/primary_secondary_correlation_v2.csv", ps_corr,delimiter=",",header="")  
np.savetxt("/Users/william/Documents/scanner/analysis/all_explosions_correlation_v2.csv", a_corr,delimiter=",",header="")  

#%%


n_corr=np.zeros(shape=(0,1))
numb_n =0

for x in range(1,len(event_times)):
    if  event_times[x]-event_times[x-1] > 600:
        n_corr = np.lib.pad(n_corr, ((0,1),(0,0)), 'constant', constant_values=(0))
        n_corr[numb_n] = corr_m[x][x-1]
        numb_n +=1



n_mean = np.mean(n_corr)
n_var = np.var(n_corr)
n_n = len(n_corr)

t_calc = (s_mean - n_mean)/np.sqrt((s_var/s_n)+(n_var/n_n))

print("t_calc (next to secondary) = ", t_calc) # = 4.067 (3dp) - significant to 0.1%




#%%


t_calc = (a_mean - n_mean)/np.sqrt((a_var/a_n)+(n_var/n_n))

print("t_calc (next to any) = ", abs(t_calc)) # = 6.215 (3dp) - significant to 0.1%

        
        
