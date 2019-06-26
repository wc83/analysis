#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 09:49:17 2019

@author: william
"""

import io
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
from scipy.stats import norm

cat = genfromtxt("/Users/william/Documents/scanner/all_stations/Final_Catalogue_2014_2018.csv", delimiter=',',skip_header=1)


Duration = cat[:,24]
Energy = cat[:,23]
evt = cat[:,0]
year = cat[:,3]


plt.figure(200)
plt.hist(Duration,80)
plt.xlim([0,80])
plt.xlabel('Event Duration [s]')
plt.ylabel('Occurance [#]')

        
plt.figure(300)
plt.semilogy(Duration,Energy,'kx')
plt.xlabel('Event Duration [s]')
plt.ylabel('Energy [J]')
plt.ylim([0,1e11])

plt.figure(400)
plt.semilogy(evt,Energy,'bx')
plt.ylim([0,1e10])
plt.xlabel('Event time [timestamp]')
plt.ylabel('Energy [J]')

plt.figure(500)
plt.plot(evt,Duration,'rx')
plt.xlabel('Event time [timestamp]')
plt.ylabel('Duration [s]')

#%% Duration and Energy

#Duration1=np.zeros(shape=(0,1))
#num1=0
Duration2=np.zeros(shape=(0,1))
num2=0
Duration3=np.zeros(shape=(0,1))
num3=0
Duration4=np.zeros(shape=(0,1))
num4=0
Duration5=np.zeros(shape=(0,1))
num5=0
Duration6=np.zeros(shape=(0,1))
num6=0
#Duration7=np.zeros(shape=(0,1))
#num7=0

for x in range(0,len(cat)):
#    if cat[x,23] < 1e4:
#        Duration1 = np.lib.pad(Duration1, ((0,1),(0,0)), 'constant', constant_values=(0))
#        Duration1[num1][0] = cat[x,24]
#        num1+=1
    if 1e4 < cat[x,23] < 1e5:
        Duration2 = np.lib.pad(Duration2, ((0,1),(0,0)), 'constant', constant_values=(0))
        Duration2[num2][0] = cat[x,24]
        num2+=1
    if 1e5 < cat[x,23] < 1e6:
        Duration3 = np.lib.pad(Duration3, ((0,1),(0,0)), 'constant', constant_values=(0))
        Duration3[num3][0] = cat[x,24]
        num3+=1
    if 1e6 < cat[x,23] < 1e7:
        Duration4 = np.lib.pad(Duration4, ((0,1),(0,0)), 'constant', constant_values=(0))
        Duration4[num4][0] = cat[x,24]
        num4+=1
    if 1e7 < cat[x,23] < 1e8:
        Duration5 = np.lib.pad(Duration5, ((0,1),(0,0)), 'constant', constant_values=(0))
        Duration5[num5][0] = cat[x,24]
        num5+=1
    if 1e8 < cat[x,23] < 1e9:
        Duration6 = np.lib.pad(Duration6, ((0,1),(0,0)), 'constant', constant_values=(0))
        Duration6[num6][0] = cat[x,24]
        num6+=1
#    if 1e9 < cat[x,23] < 1e10:
#        Duration7 = np.lib.pad(Duration7, ((0,1),(0,0)), 'constant', constant_values=(0))
#        Duration7[num7][0] = cat[x,24]
#        num7+=1





#plt.figure(1000)
#plt.hist(Duration1,bins=40)
#plt.xlim([0,80])
#plt.xlabel('Event_Duration [s]')
#plt.ylabel('Occurance [#]')
#plt.title('Energy < 1e4')

plt.figure(2000)
plt.hist(Duration2,bins=40)
plt.xlim([0,80])
plt.xlabel('Event_Duration [s]')
plt.ylabel('Occurance [#]')
plt.title('1e4 < Energy < 1e5')

plt.figure(3000)
plt.hist(Duration3,bins=40)
plt.xlim([0,80])
plt.xlabel('Event_Duration [s]')
plt.ylabel('Occurance [#]')
plt.title('1e5 < Energy < 1e6')

plt.figure(4000)
plt.hist(Duration4,bins=40)
plt.xlim([0,80])
plt.xlabel('Event_Duration [s]')
plt.ylabel('Occurance [#]')
plt.title('1e6 < Energy < 1e7')

plt.figure(5000)
plt.hist(Duration5,bins=40)
plt.xlim([0,80])
plt.xlabel('Event_Duration [s]')
plt.ylabel('Occurance [#]')
plt.title('1e7 < Energy < 1e8')

plt.figure(6000)
plt.hist(Duration6,bins=40)
plt.xlim([0,80])
plt.xlabel('Event_Duration [s]')
plt.ylabel('Occurance [#]')
plt.title('1e8 < Energy < 1e9')

#plt.figure(7000)
#plt.hist(Duration7,bins=40)
#plt.xlim([0,80])
#plt.xlabel('Event_Duration [s]')
#plt.ylabel('Occurance [#]')
#plt.title('1e9 < Energy < 1e10')

     
#%% Duration Through Time

Duration10=np.zeros(shape=(0,1))
num1=0
Duration20=np.zeros(shape=(0,1))
num2=0
Duration30=np.zeros(shape=(0,1))
num3=0
Duration40=np.zeros(shape=(0,1))
num4=0
Duration50=np.zeros(shape=(0,1))
num5=0


for x in range(0,len(cat)):
    if year[x] == 2014:
        Duration10 = np.lib.pad(Duration10, ((0,1),(0,0)), 'constant', constant_values=(0))
        Duration10[num1][0] = cat[x,24]
        num1+=1
    if year[x] == 2015:
        Duration20 = np.lib.pad(Duration20, ((0,1),(0,0)), 'constant', constant_values=(0))
        Duration20[num2][0] = cat[x,24]
        num2+=1
    if year[x] == 2016:
        Duration30 = np.lib.pad(Duration30, ((0,1),(0,0)), 'constant', constant_values=(0))
        Duration30[num3][0] = cat[x,24]
        num3+=1
    if year[x] == 2017:
        Duration40 = np.lib.pad(Duration40, ((0,1),(0,0)), 'constant', constant_values=(0))
        Duration40[num4][0] = cat[x,24]
        num4+=1
    if year[x] == 2018:
        Duration50 = np.lib.pad(Duration50, ((0,1),(0,0)), 'constant', constant_values=(0))
        Duration50[num5][0] = cat[x,24]
        num5+=1
    





plt.figure(10000)
plt.hist(Duration10,bins=40)
plt.xlim([0,80])
plt.xlabel('Event_Duration [s]')
plt.ylabel('Occurance [#]')
plt.title('2014')
print('2014 median = ', np.median(Duration10))  

plt.figure(20000)
plt.hist(Duration20,bins=40)
plt.xlim([0,80])
plt.xlabel('Event_Duration [s]')
plt.ylabel('Occurance [#]')
plt.title('2015')
print('2015 median = ', np.median(Duration20))  

plt.figure(30000)
plt.hist(Duration30,bins=40)
plt.xlim([0,80])
plt.xlabel('Event_Duration [s]')
plt.ylabel('Occurance [#]')
plt.title('2016')
print('2016 median = ', np.median(Duration30))  

plt.figure(40000)
plt.hist(Duration40,bins=40)
plt.xlim([0,80])
plt.xlabel('Event_Duration [s]')
plt.ylabel('Occurance [#]')
plt.title('2017')
print('2017 median = ', np.median(Duration40))  

plt.figure(50000)
plt.hist(Duration50,bins=40)
plt.xlim([0,80])
plt.xlabel('Event_Duration [s]')
plt.ylabel('Occurance [#]')
plt.title('2018')
print('2018 median = ', np.median(Duration50))  

           
        
 #%%

E10 = np.zeros(shape=(0,1))
E20 = np.zeros(shape=(0,1))
E30 = np.zeros(shape=(0,1))
E40 = np.zeros(shape=(0,1))
E50 = np.zeros(shape=(0,1))
E60 = np.zeros(shape=(0,1))
E70 = np.zeros(shape=(0,1))
E80 = np.zeros(shape=(0,1))

ne10 = 0
ne20 = 0
ne30 = 0
ne40 = 0
ne50 = 0
ne60 = 0
ne70 = 0
ne80 = 0

for x in range(0,len(cat)):
    if Duration[x] < 10:
        E10 = np.lib.pad(E10, ((0,1),(0,0)), 'constant', constant_values=(0))
        E10[ne10][0]=Energy[x]
        ne10 += 1
    if 10 < Duration[x] < 20:
        E20 = np.lib.pad(E20, ((0,1),(0,0)), 'constant', constant_values=(0))
        E20[ne20][0]=Energy[x]
        ne20 += 1
    if 20 < Duration[x] < 30:
        E30 = np.lib.pad(E30, ((0,1),(0,0)), 'constant', constant_values=(0))
        E30[ne30][0]=Energy[x]
        ne30 += 1
    if 30 < Duration[x] < 40:
        E40 = np.lib.pad(E40, ((0,1),(0,0)), 'constant', constant_values=(0))
        E40[ne40][0]=Energy[x]
        ne40 += 1
    if 40 < Duration[x] < 50:
        E50 = np.lib.pad(E50, ((0,1),(0,0)), 'constant', constant_values=(0))
        E50[ne50][0]=Energy[x]
        ne50 += 1
    if 50 < Duration[x] < 60:
        E60 = np.lib.pad(E60, ((0,1),(0,0)), 'constant', constant_values=(0))
        E60[ne60][0]=Energy[x]
        ne60 += 1
    if 60 < Duration[x] < 70:
        E70 = np.lib.pad(E70, ((0,1),(0,0)), 'constant', constant_values=(0))
        E70[ne70][0]=Energy[x]
        ne70 += 1
    if 70 < Duration[x] < 80:
        E80 = np.lib.pad(E80, ((0,1),(0,0)), 'constant', constant_values=(0))
        E80[ne80][0]=Energy[x]
        ne80 += 1

plt.figure(34231)        
plt.semilogy(Duration,Energy,'kx')
#plt.semilogy(5,np.mean(E10),'bx') 
plt.semilogy(5,np.median(E10),'rx') 
  
#plt.semilogy(15,np.mean(E20),'bx')   
#plt.semilogy(25,np.mean(E30),'bx')   
#plt.semilogy(35,np.mean(E40),'bx')   
#plt.semilogy(45,np.mean(E50),'bx')  
#plt.semilogy(55,np.mean(E60),'bx')      
#plt.semilogy(65,np.mean(E70),'bx')   
#plt.semilogy(75,np.mean(E80),'bx')        
  
plt.semilogy(15,np.median(E20),'rx')   
plt.semilogy(25,np.median(E30),'rx')   
plt.semilogy(35,np.median(E40),'rx')   
plt.semilogy(45,np.median(E50),'rx')  
plt.semilogy(55,np.median(E60),'rx')      
plt.semilogy(65,np.median(E70),'rx')   
plt.semilogy(75,np.median(E80),'rx')     
plt.legend(['raw','median'])
        
        
        
        
        
        
        
        
        
        
        
#Energy2=np.zeros(shape=(0,1))
#num=0
#for x in range(0,len(cat)):
#    if cat[x,23] < 1e7:
#        Energy2 = np.lib.pad(Energy2, ((0,1),(0,0)), 'constant', constant_values=(0))                            
#        Energy2[num][0]=cat[x,23]
#        num+=1



#plt.figure(4001)
#plt.hist(Energy2,bins=20, histtype='step')
#plt.yscale('log')
#plt.xscale('log')
#plt.xlim([1e4,5e7])
#plt.xlabel('Energy [J]')
#plt.ylabel('Occurance [#]')

#%%
#Repose = np.zeros(shape=(0,1))
#num=0
#for x in range(1,len(evt)):
#    rep = evt[x]-evt[x-1]
#    
#    if rep < 6*60*60:
#        Repose = np.lib.pad(Repose, ((0,1),(0,0)), 'constant', constant_values=(0))             
#        Repose[num][0]=rep/60
#        num+=1
#        
#        
#        
#plt.figure(1)
#plt.hist(Repose,bins=180)
#plt.xlabel('Repose Time [mins]')
#plt.ylabel('Occurance [#]')
#plt.title('Repose Times - My Catalogue')
 


































       
        