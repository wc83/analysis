#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 10:37:38 2018

@author: william
"""


from numpy import argmax
import numpy as np
from obspy.clients.earthworm import Client
from obspy import UTCDateTime
from obspy.signal.trigger import trigger_onset
from obspy import Stream
from obspy import Trace
import matplotlib.pyplot as plt
from obspy import read, read_inventory
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

#%% Get large exposion to create source function with
with io.open("/Users/william/Documents/scanner/output_data/m30.mseed", "rb") as fh:
    while True:
        with io.BytesIO() as buf:
            c = fh.read(chunksize)
            if not c:
                break
            buf.write(c)
            buf.seek(0, 0)
            st = obspy.read(buf)
            
            for x in range(0,len(st)):
                if st[x].stats.station == "LB01":
                    if  1438145500 > st[x].stats.starttime.timestamp > 1438145480:
                        tr = st[x]
                        tr.detrend(type='linear')
                        tr.detrend(type='demean')
                        t1 = st[x].stats.starttime
                        t2 = st[x].stats.endtime
                        window=t2-t1
                        tr_data=st[x].data
                        m=np.mean(tr_data)
                        tr_data = tr_data-m
                        
                        B=2*pi*rhoE*cE*(1/A)                    
                        EI = sum(np.square(tr.data))
                        EE1= B*(r1*r1)*EI
                        
                        st_c1 = calibrate(tr)
                        st_c1.plot(color='r')
                        
                        famp = abs(np.fft.rfft(tr_data))
  
#%% Get all explosions to create their source functions using big explosion
scale=np.zeros(shape=(0,2))       
number=0            
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
            
            for x in range(0,len(st)):
#            for x in range(0,20):  
#                print(x)
                if st[x].stats.station == "LB01":
                    if 1451607000 > st[x].stats.starttime.timestamp:
                        tr2 = st[x]
                        tr2.detrend(type='linear')
                        tr2.detrend(type='demean')
                        t1 = st[x].stats.starttime
                        t2 = st[x].stats.endtime
                        window=t2-t1
                        tr2_data=st[x].data
                        m=np.mean(tr2_data)
                        tr2_data = tr2_data-m
                        
                        st_c2 = calibrate(tr2)
#                        st_c2.plot(color='g')
                        
                    
                        famp2 = abs(np.fft.rfft(tr2_data))
                     
                        
                        DC = famp2[2:-2]/famp[2:-2]
                        st_dc = np.fft.ifft(DC)
                        st_r = st_dc[10:-10].real
#                        plt.plot(st_r)
                        
                        trace_obs = Trace()
                        trace_obs.data = st_r
                        st_obs = Stream(trace_obs)
                        st_obs.filter(type='lowpass',freq=0.5/100)
#                        st_obs.plot(color='k')
                        
                        
                        B=2*pi*rhoE*cE*(1/A)                    
                        EI = sum(np.square(tr2.data))
                        EE= B*(r1*r1)*EI
#                        if EE < 1e8:
                        scale = np.lib.pad(scale, ((0,1),(0,0)), 'constant', constant_values=(0))
                        scale[number][0]=max(st_obs[0].data) - (min(st_obs[0].data))
                        scale[number][1]=EE
                        number +=1
                        
#                        print('greatest amplitude = ', abs(min(st_obs[0].data)))
                        
            break                  
                        
plt.figure()
plt.semilogy(scale[:,0],scale[:,1],'rx')   
plt.xlim([0,0.002])
#plt.ylim([1e1,2e8])
plt.title("Source Function Amplitude (P2P) vs Energy")   
plt.xlabel("Peak to Peak Amplitude")                  
plt.ylabel("Explosion Energy")

#np.savetxt("/Users/william/Documents/scanner/analysis/Scaling.csv", scale,delimiter=",",header="Amplitude ,Energy ")                    
#%% predict amplitude of explosion p from energy of both explosion p and x energy + explosion x amplitude
slen = len(scale)

alphaA=np.zeros(shape=(slen,slen)) 
s2 = np.zeros(shape=(0,6))
num=0
 
for x in range(0,slen):
    for p in range(x+1,slen):
        if scale[p,0] < 0.02:
            alphaA[x][p]=(scale[p,1]/scale[x,1])**(1/3)
            s2 = np.lib.pad(s2, ((0,1),(0,0)), 'constant', constant_values=(0))
            s2[num][0]=scale[x,0]
            s2[num][1]=scale[p,0]
            s2[num][2]=(scale[x,0])* (alphaA[x][p])
            s2[num][3]=scale[x,1]
            s2[num][4]=scale[p,1]
            
            E_ratio = scale[x,1]/scale[p,1]
            if E_ratio >= 1:
                s2[num][5]=E_ratio
            else:
                E_ratio2 = scale[p,1]/scale[1,1]
                s2[num][5]=E_ratio2
            
            
            num += 1

r=np.linspace(0,0.1,1001)

plt.figure()
plt.plot(s2[:,1],s2[:,2],'bx')
plt.title("Calculated vs Observed Amplitude")   
plt.xlabel("Observed Peak to Peak Amplitude")                  
plt.ylabel("Calculated Peak to Peak Amplitude")
plt.plot(r,r,'k-') 
plt.xlim([0,0.003])
plt.ylim([0,0.006])

#np.savetxt("/Users/william/Documents/scanner/analysis/Scaling_prediction.csv", s2,delimiter=",",header="S1 ,S2_obs ,S2_est ")                    



plt.figure()
sc = plt.scatter(x=s2[:,1], y=s2[:,2], c=s2[:,3], cmap="rainbow")
plt.colorbar(sc)
plt.plot(r,r,'k-')    
plt.xlim([0,0.003])
plt.ylim([0,0.01])
plt.title("Calculated vs Observed Amplitude of S2")   
plt.xlabel("Observed Peak to Peak Amplitude of S2")                  
plt.ylabel("Calculated Peak to Peak Amplitude of S2")
plt.legend(['y=x','energy of S1'])

plt.figure()
sc = plt.scatter(x=s2[:,1], y=s2[:,2], c=s2[:,0], cmap="rainbow")
plt.colorbar(sc)
plt.plot(r,r,'k-')    
plt.xlim([0,0.0075])
plt.ylim([0,0.015])
plt.title("Calculated vs Observed Amplitude")   
plt.xlabel("Observed Peak to Peak Amplitude")                  
plt.ylabel("Calculated Peak to Peak Amplitude")
plt.legend(['y=x','Amplitude of S1'])

plt.figure()
sc = plt.scatter(x=s2[:,1], y=s2[:,2], c=s2[:,4], cmap="rainbow")
plt.colorbar(sc)
plt.plot(r,r,'k-')    
plt.xlim([0,0.003])
plt.ylim([0,0.01])
plt.title("Calculated vs Observed Amplitude of S2")   
plt.xlabel("Observed Peak to Peak Amplitude of S2")                  
plt.ylabel("Calculated Peak to Peak Amplitude of S2 ")
plt.legend(['y=x','energy of S2'])




       