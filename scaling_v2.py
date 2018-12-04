#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 10 13:43:32 2018

@author: william
"""


from numpy import argmax
import numpy as np
from obspy.clients.earthworm import Client
from obspy import UTCDateTime
from obspy.signal.trigger import trigger_onset
from obspy import Stream
import matplotlib.pyplot as plt

import io
import obspy 

reclen = 512
chunksize = 100000 * reclen # Around 50 MB

scale=np.zeros(shape=(0,9))

event_num = 0

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

first=0

with io.open("/Users/william/Documents/scanner/output_data/m30.mseed", "rb") as fh:
        # just month 2
#with io.open("/Users/william/Documents/scanner/output_data/EXP_all_data_stream_2_month_1.mseed", "rb") as fh:
    while True:
        with io.BytesIO() as buf:
            c = fh.read(chunksize)
            if not c:
                break
            buf.write(c)
            buf.seek(0, 0)
            st = obspy.read(buf)
           
            
        for x in range(3,len(st)):
            if first == 0:
                if st[x].stats.station == "LB01":# or "LB02" or "LB03" or "LB04" or "LB05" or "LB06" :
                    tr = st[x]
                    t1 = st[x].stats.starttime
                    t2 = st[x].stats.endtime
                    window=t2-t1
                    tr_data=st[x].data
                    m=np.mean(tr_data)
                    tr_data = tr_data-m
                    famp = abs(np.fft.rfft(tr_data))
                    famp2=famp*famp
                    
                    # Dominant, Central, Bamdwidth50
                    dom1,cent1, bw, bwid25 = freq_info(tr,t1 ,t2)
                    
                    r=r1
            
                    st_c1 = calibrate(tr)
                    B=2*pi*rhoE*cE*(1/A)
                    EI = sum(np.square(st_c1[0].data))
                    Energy1 = B*(r*r)*EI
            
                    Amp1 = max(st_c1[0].data)
                    st_c1.plot()
                    print("\n")
                        
                    scale = np.lib.pad(scale, ((0,1),(0,0)), 'constant', constant_values=(0))
                    scale[event_num][0]=Energy1
                    scale[event_num][1]=Amp1*1000
                    scale[event_num][2]=Amp1*1000
                    scale[event_num][3]=1
                    scale[event_num][4]=t1
                    scale[event_num][5]=1
                    scale[event_num][6]=Amp1*1000
                    scale[event_num][7]=dom1
                    scale[event_num][8]=dom1
                    event_num += 1
                    first=1
                        
            else:
                if st[x].stats.station == "LB01":# or "LB02" or "LB03" or "LB04" or "LB05" or "LB06" :
                    tr = st[x]
                    t1 = st[x].stats.starttime
                    t2 = st[x].stats.endtime
                    window=t2-t1
                    tr_data=st[x].data
                    m=np.mean(tr_data)
                    tr_data = tr_data-m
                    famp = abs(np.fft.rfft(tr_data))
                    famp2=famp*famp
                    
                    # Dominant, Central, Bamdwidth50
                    dom,cent, bw, bwid25 = freq_info(tr,t1 ,t2)
                    
                    r=r1
            
                    st_c = calibrate(tr)
                    B=2*pi*rhoE*cE*(1/A)
                    EI = sum(np.square(st_c[0].data))
                    Energy = B*(r*r)*EI
            
                    Amp = max(st_c[0].data)
                    
                    
                    fd=Energy/Energy1
                    alpha=fd**(1/3)
                    alph2=alpha**2
                    
                    exp_dom = alph2*dom1
                    
                    expA= Amp1*(fd**(1/3))
                    sqA=Amp1*(fd**(1/2))
                    
                    scale = np.lib.pad(scale, ((0,1),(0,0)), 'constant', constant_values=(0))
                    scale[event_num][0]=Energy
                    scale[event_num][1]=Amp*1000
                    scale[event_num][2]=expA*1000
                    scale[event_num][3]=expA/Amp
                    scale[event_num][4]=t1
                    scale[event_num][5]=fd
                    scale[event_num][6]=sqA*1000
                    scale[event_num][7]=dom
                    scale[event_num][8]=exp_dom
                    
                    event_num +=1
            
                
                    
                    
print('Minimum A scale ratio = ',min(scale[:,3]))
print('Maximum A scale ratio = ',max(scale[:,3]))
print('Largest A scale ratio = ',max(scale[:,3])/min(scale[:,3]))

av_scale= sum(scale[:,3])/len(scale)
print('Average A scale ratio = ',av_scale)

print('Maximum Energy = ', max(scale[:,0]))
print('Minumum Energy = ', min(scale[:,0]))                   
                    

print('number of data points = ', len(scale))
 
#u=np.linspace(0,1.5e14,100001)
#i=u**(1/3)
#e=u**(1/2)
#o=0.0000005
#l=0.000004
#
#plt.figure(1)                   
#plt.plot(scale[:,0],scale[:,1],'bx') 
##plt.plot(u,o*i,'k-')
#plt.plot(u,l*e,'k--')
#plt.xlim((0 ,7e13))            
#plt.ylim((0,50))       
#plt.ylabel('Max Ground Velocity Amplitude [km/s]')
#plt.xlabel('Explosion Energy [J]')                    
#plt.title('Energy vs Amplitude') 
#plt.legend(['Explosions','Cube law','Square law'])
#
#                  
#
#p=np.linspace(0,0.3,1001)
#q=p**(0.5)
#k=0.0025
#
#plt.figure(2)
#plt.xlabel('True Ground Velocity [km/s]')
#plt.ylabel('Predicted Ground Velocity [km/s]') 
#plt.title('Cube Law Scaling') 
#sc = plt.scatter(x=scale[0:1500,1], y=scale[0:1500,2], c=scale[0:1500,0], cmap="rainbow")
#plt.colorbar(sc)
#plt.plot(p,p,'k-')    
##plt.plot(p,k*q,'k--')              
#plt.xlim((0 ,0.2))            
#plt.ylim((0,0.2))
##plt.legend(['Explosions','Cube law','Square law'])
#plt.show()  
#
#v=np.linspace(0,0.3,1001)
#        
#                    
#plt.figure(3)
#plt.xlabel('True Ground Velocity [km/s]')
#plt.ylabel('Predicted Ground Velocity [km/s]') 
#plt.title('Square Law Scaling') 
#sc = plt.scatter(x=scale[0:1500,1], y=scale[0:1500,6], c=scale[0:1500,0], cmap="rainbow")
#plt.colorbar(sc)
#plt.plot(v,v,'k-')                
#plt.xlim((0 ,0.2))            
#plt.ylim((0,0.2))
##plt.legend(['Explosions'])
#plt.show()      
#
#
#
plt.figure(4)                   
plt.plot(scale[0:200,0],scale[0:200,7],'bx') 
#plt.plot(u,o*i,'k-')
#plt.plot(u,l*e,'k--')
#plt.xlim((0 ,7e13))            
#plt.ylim((0,50))       
plt.ylabel('Dominant F [Hz]')
plt.xlabel('Explosion Energy [J]')                    
plt.title('Energy vs Frequency') 
#plt.legend(['Explosions','Cube law','Square law'])
   
        
                    
r=np.linspace(0,5,1001)
        
                    
plt.figure(5)
plt.xlabel('True dominant f [Hz]')
plt.ylabel('Predicted dominant f [Hz]') 
plt.title('Frequency Scaling') 
sc = plt.scatter(x=scale[0:1000,7], y=scale[0:1000,8], c=scale[0:1000,0], cmap="rainbow")
plt.colorbar(sc)
plt.plot(r,r,'k-')                
#plt.xlim((0 ,0.2))            
#plt.ylim((0,20))
#plt.legend(['Explosions'])
plt.show()      
      

  
                    
                    