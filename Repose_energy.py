#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 10:43:54 2018

@author: william
"""



import numpy as np
from numpy import genfromtxt
from obspy import UTCDateTime
import numpy as np
import matplotlib.pyplot as plt
import obspy
import io
from obspy import Stream
from obspy.clients.earthworm import Client

#%%


cat= genfromtxt("/Users/william/Documents/scanner/all_stations/Explosion_Catalogue_V4c.csv", delimiter=',',skip_header=1)


pands=np.zeros(shape=(1,3))
pands[0,0]=cat[0,0]
pands[0,1]=1
pands[0,2]=0
pands = np.lib.pad(pands, ((0,1),(0,0)), 'constant', constant_values=(0))

last_p = cat[0,0]

p_tot=1
s_tot=0


sec=np.zeros(shape=(0,3))
num=0

for x in range(1,len(cat)):
    pands[x,0]=cat[x,0]
    
    if cat[x,0] - last_p < 10*60:
        pands[x,1]=2
        s_tot += 1
        
  
        
        sec = np.lib.pad(sec, ((0,1),(0,0)), 'constant', constant_values=(0))
        sec[num][0]= last_p
        sec[num][1]= cat[x,0]
        num += 1
        
    else:
        pands[x,1]=1
        last_p=cat[x,0]
        p_tot +=1
    pands[x,2]=(cat[x,0]-cat[x-1,0] )/60 
    if x < len(cat)-1:   
        pands = np.lib.pad(pands, ((0,1),(0,0)), 'constant', constant_values=(0))
        



#print(len(sec))

#%%
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




net = 'Z4'  # 
loc = ''    # location, it depends mostly of which network you are in. 
client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23

for x in range(0,len(sec)):
    try:
        sta = 'LB01' # STATION 
        cha = 'HHZ' # Vertical Seismic Channel
        # PRIMARY
        t1 = UTCDateTime(sec[x,0]-10) #the format is year:day_of_the_year:month
        t2 = t1 + 60           
        st = Stream()
        st = client.get_waveforms(net, sta, '', cha, t1 , t2)
        
        # Filter data
        st.detrend(type='linear')
        st.detrend(type='demean')
        st.filter(type='bandpass',freqmin=0.2, freqmax=10)
        trs=st[0]
    
        st_c = calibrate1(trs)
        B=2*pi*rhoE*cE*(1/A)
        EI = sum(np.square(st_c[0].data))
        EEP= B*(r1*r1)*EI
    
        # SECONDARY
        t1 = UTCDateTime(sec[x,1]-10) #the format is year:day_of_the_year:month
        t2 = t1 + 60           
        st = Stream()
        st = client.get_waveforms(net, sta, '', cha, t1 , t2)
        
        # Filter data
        st.detrend(type='linear')
        st.detrend(type='demean')
        st.filter(type='bandpass',freqmin=0.2, freqmax=10)
        trs=st[0]
    
        st_c = calibrate1(trs)
        B=2*pi*rhoE*cE*(1/A)
        EI = sum(np.square(st_c[0].data))
        EES= B*(r1*r1)*EI
    
        # Percentage
        
        
        sec[x][2]= (EES/EEP)*100
    except: 
         try:
            sta = 'LB02' # STATION 
            cha = 'HHZ' # Vertical Seismic Channel
            # PRIMARY
            t1 = UTCDateTime(sec[x,0]-10) #the format is year:day_of_the_year:month
            t2 = t1 + 60           
            st = Stream()
            st = client.get_waveforms(net, sta, '', cha, t1 , t2)
            
            # Filter data
            st.detrend(type='linear')
            st.detrend(type='demean')
            st.filter(type='bandpass',freqmin=0.2, freqmax=10)
            trs=st[0]
        
            st_c = calibrate1(trs)
            B=2*pi*rhoE*cE*(1/A)
            EI = sum(np.square(st_c[0].data))
            EEP= B*(r2*r2)*EI
        
            # SECONDARY
            t1 = UTCDateTime(sec[x,1]-10) #the format is year:day_of_the_year:month
            t2 = t1 + 60           
            st = Stream()
            st = client.get_waveforms(net, sta, '', cha, t1 , t2)
            
            # Filter data
            st.detrend(type='linear')
            st.detrend(type='demean')
            st.filter(type='bandpass',freqmin=0.2, freqmax=10)
            trs=st[0]
        
            st_c = calibrate1(trs)
            B=2*pi*rhoE*cE*(1/A)
            EI = sum(np.square(st_c[0].data))
            EES= B*(r2*r2)*EI
        
            # Percentage
            
            
            sec[x][2]= (EES/EEP)*100
         except: 
             try:
                sta = 'LB03' # STATION 
                cha = 'HHZ' # Vertical Seismic Channel
                # PRIMARY
                t1 = UTCDateTime(sec[x,0]-10) #the format is year:day_of_the_year:month
                t2 = t1 + 60           
                st = Stream()
                st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                
                # Filter data
                st.detrend(type='linear')
                st.detrend(type='demean')
                st.filter(type='bandpass',freqmin=0.2, freqmax=10)
                trs=st[0]
            
                st_c = calibrate1(trs)
                B=2*pi*rhoE*cE*(1/A)
                EI = sum(np.square(st_c[0].data))
                EEP= B*(r3*r3)*EI
            
                # SECONDARY
                t1 = UTCDateTime(sec[x,1]-10) #the format is year:day_of_the_year:month
                t2 = t1 + 60           
                st = Stream()
                st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                
                # Filter data
                st.detrend(type='linear')
                st.detrend(type='demean')
                st.filter(type='bandpass',freqmin=0.2, freqmax=10)
                trs=st[0]
            
                st_c = calibrate1(trs)
                B=2*pi*rhoE*cE*(1/A)
                EI = sum(np.square(st_c[0].data))
                EES= B*(r3*r3)*EI
            
                # Percentage
                
                
                sec[x][2]= (EES/EEP)*100
             except: 
                 try:
                    sta = 'LB04' # STATION 
                    cha = 'HHZ' # Vertical Seismic Channel
                    # PRIMARY
                    t1 = UTCDateTime(sec[x,0]-10) #the format is year:day_of_the_year:month
                    t2 = t1 + 60           
                    st = Stream()
                    st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                    
                    # Filter data
                    st.detrend(type='linear')
                    st.detrend(type='demean')
                    st.filter(type='bandpass',freqmin=0.2, freqmax=10)
                    trs=st[0]
                
                    st_c = calibrate1(trs)
                    B=2*pi*rhoE*cE*(1/A)
                    EI = sum(np.square(st_c[0].data))
                    EEP= B*(r4*r4)*EI
                
                    # SECONDARY
                    t1 = UTCDateTime(sec[x,1]-10) #the format is year:day_of_the_year:month
                    t2 = t1 + 60           
                    st = Stream()
                    st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                    
                    # Filter data
                    st.detrend(type='linear')
                    st.detrend(type='demean')
                    st.filter(type='bandpass',freqmin=0.2, freqmax=10)
                    trs=st[0]
                
                    st_c = calibrate1(trs)
                    B=2*pi*rhoE*cE*(1/A)
                    EI = sum(np.square(st_c[0].data))
                    EES= B*(r4*r4)*EI
                
                    # Percentage
                    
                    
                    sec[x][2]= (EES/EEP)*100
                 except: 
                     try:
                        sta = 'LB05' # STATION 
                        cha = 'HHZ' # Vertical Seismic Channel
                        # PRIMARY
                        t1 = UTCDateTime(sec[x,0]-10) #the format is year:day_of_the_year:month
                        t2 = t1 + 60           
                        st = Stream()
                        st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                        
                        # Filter data
                        st.detrend(type='linear')
                        st.detrend(type='demean')
                        st.filter(type='bandpass',freqmin=0.2, freqmax=10)
                        trs=st[0]
                    
                        st_c = calibrate1(trs)
                        B=2*pi*rhoE*cE*(1/A)
                        EI = sum(np.square(st_c[0].data))
                        EEP= B*(r5*r5)*EI
                    
                        # SECONDARY
                        t1 = UTCDateTime(sec[x,1]-10) #the format is year:day_of_the_year:month
                        t2 = t1 + 60           
                        st = Stream()
                        st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                        
                        # Filter data
                        st.detrend(type='linear')
                        st.detrend(type='demean')
                        st.filter(type='bandpass',freqmin=0.2, freqmax=10)
                        trs=st[0]
                    
                        st_c = calibrate1(trs)
                        B=2*pi*rhoE*cE*(1/A)
                        EI = sum(np.square(st_c[0].data))
                        EES= B*(r5*r5)*EI
                    
                        # Percentage
                        
                        
                        sec[x][2]= (EES/EEP)*100
                     except: 
                        try:
                            sta = 'LB06' # STATION 
                            cha = 'HHZ' # Vertical Seismic Channel
                            # PRIMARY
                            t1 = UTCDateTime(sec[x,0]-10) #the format is year:day_of_the_year:month
                            t2 = t1 + 60           
                            st = Stream()
                            st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                            
                            # Filter data
                            st.detrend(type='linear')
                            st.detrend(type='demean')
                            st.filter(type='bandpass',freqmin=0.2, freqmax=10)
                            trs=st[0]
                        
                            st_c = calibrate1(trs)
                            B=2*pi*rhoE*cE*(1/A)
                            EI = sum(np.square(st_c[0].data))
                            EEP= B*(r6*r6)*EI
                        
                            # SECONDARY
                            t1 = UTCDateTime(sec[x,1]-10) #the format is year:day_of_the_year:month
                            t2 = t1 + 60           
                            st = Stream()
                            st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                            
                            # Filter data
                            st.detrend(type='linear')
                            st.detrend(type='demean')
                            st.filter(type='bandpass',freqmin=0.2, freqmax=10)
                            trs=st[0]
                        
                            st_c = calibrate1(trs)
                            B=2*pi*rhoE*cE*(1/A)
                            EI = sum(np.square(st_c[0].data))
                            EES= B*(r6*r6)*EI
                        
                            # Percentage
                            
                            
                            sec[x][2]= (EES/EEP)*100
                        except: 
                        
                            
                            
                                try:
                                    sta = 'LS01' # STATION 
                                    cha = 'EHZ' 
                                    # PRIMARY
                                    t1 = UTCDateTime(sec[x,0]-10) #the format is year:day_of_the_year:month
                                    t2 = t1 + 60           
                                    st = Stream()
                                    st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                                    
                                    # Filter data
                                    st.detrend(type='linear')
                                    st.detrend(type='demean')
                                    st.filter(type='bandpass',freqmin=0.2, freqmax=10)
                                    trs=st[0]
                                
                                    st_c = calibrate1(trs)
                                    B=2*pi*rhoE*cE*(1/A)
                                    EI = sum(np.square(st_c[0].data))
                                    EEP= B*(rs1*rs1)*EI
                                
                                    # SECONDARY
                                    t1 = UTCDateTime(sec[x,1]-10) #the format is year:day_of_the_year:month
                                    t2 = t1 + 60           
                                    st = Stream()
                                    st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                                    
                                    # Filter data
                                    st.detrend(type='linear')
                                    st.detrend(type='demean')
                                    st.filter(type='bandpass',freqmin=0.2, freqmax=10)
                                    trs=st[0]
                                
                                    st_c = calibrate1(trs)
                                    B=2*pi*rhoE*cE*(1/A)
                                    EI = sum(np.square(st_c[0].data))
                                    EES= B*(rs1*rs1)*EI
                                
                                    # Percentage
                                    
                                    
                                    sec[x][2]= (EES/EEP)*100
                                except: 
                                     try:
                                        sta = 'LS02' # STATION 
                                        cha = 'EHZ' 
                                        # PRIMARY
                                        t1 = UTCDateTime(sec[x,0]-10) #the format is year:day_of_the_year:month
                                        t2 = t1 + 60           
                                        st = Stream()
                                        st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                                        
                                        # Filter data
                                        st.detrend(type='linear')
                                        st.detrend(type='demean')
                                        st.filter(type='bandpass',freqmin=0.2, freqmax=10)
                                        trs=st[0]
                                    
                                        st_c = calibrate1(trs)
                                        B=2*pi*rhoE*cE*(1/A)
                                        EI = sum(np.square(st_c[0].data))
                                        EEP= B*(rs2*rs2)*EI
                                    
                                        # SECONDARY
                                        t1 = UTCDateTime(sec[x,1]-10) #the format is year:day_of_the_year:month
                                        t2 = t1 + 60           
                                        st = Stream()
                                        st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                                        
                                        # Filter data
                                        st.detrend(type='linear')
                                        st.detrend(type='demean')
                                        st.filter(type='bandpass',freqmin=0.2, freqmax=10)
                                        trs=st[0]
                                    
                                        st_c = calibrate1(trs)
                                        B=2*pi*rhoE*cE*(1/A)
                                        EI = sum(np.square(st_c[0].data))
                                        EES= B*(rs2*rs2)*EI
                                    
                                        # Percentage
                                        
                                        
                                        sec[x][2]= (EES/EEP)*100
                                     except: 
                                         try:
                                            sta = 'LS03' # STATION 
                                            cha = 'EHZ' 
                                            # PRIMARY
                                            t1 = UTCDateTime(sec[x,0]-10) #the format is year:day_of_the_year:month
                                            t2 = t1 + 60           
                                            st = Stream()
                                            st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                                            
                                            # Filter data
                                            st.detrend(type='linear')
                                            st.detrend(type='demean')
                                            st.filter(type='bandpass',freqmin=0.2, freqmax=10)
                                            trs=st[0]
                                        
                                            st_c = calibrate1(trs)
                                            B=2*pi*rhoE*cE*(1/A)
                                            EI = sum(np.square(st_c[0].data))
                                            EEP= B*(rs3*rs3)*EI
                                        
                                            # SECONDARY
                                            t1 = UTCDateTime(sec[x,1]-10) #the format is year:day_of_the_year:month
                                            t2 = t1 + 60           
                                            st = Stream()
                                            st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                                            
                                            # Filter data
                                            st.detrend(type='linear')
                                            st.detrend(type='demean')
                                            st.filter(type='bandpass',freqmin=0.2, freqmax=10)
                                            trs=st[0]
                                        
                                            st_c = calibrate1(trs)
                                            B=2*pi*rhoE*cE*(1/A)
                                            EI = sum(np.square(st_c[0].data))
                                            EES= B*(rs3*rs3)*EI
                                        
                                            # Percentage
                                            
                                            
                                            sec[x][2]= (EES/EEP)*100
                                         except: 
                                            try:
                                                sta = 'LS04' # STATION 
                                                cha = 'EHZ' 
                                                # PRIMARY
                                                t1 = UTCDateTime(sec[x,0]-10) #the format is year:day_of_the_year:month
                                                t2 = t1 + 60           
                                                st = Stream()
                                                st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                                                
                                                # Filter data
                                                st.detrend(type='linear')
                                                st.detrend(type='demean')
                                                st.filter(type='bandpass',freqmin=0.2, freqmax=10)
                                                trs=st[0]
                                            
                                                st_c = calibrate1(trs)
                                                B=2*pi*rhoE*cE*(1/A)
                                                EI = sum(np.square(st_c[0].data))
                                                EEP= B*(rs4*rs4)*EI
                                            
                                                # SECONDARY
                                                t1 = UTCDateTime(sec[x,1]-10) #the format is year:day_of_the_year:month
                                                t2 = t1 + 60           
                                                st = Stream()
                                                st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                                                
                                                # Filter data
                                                st.detrend(type='linear')
                                                st.detrend(type='demean')
                                                st.filter(type='bandpass',freqmin=0.2, freqmax=10)
                                                trs=st[0]
                                            
                                                st_c = calibrate1(trs)
                                                B=2*pi*rhoE*cE*(1/A)
                                                EI = sum(np.square(st_c[0].data))
                                                EES= B*(rs4*rs4)*EI
                                            
                                                # Percentage
                                                
                                                
                                                sec[x][2]= (EES/EEP)*100
                                            except: 
                                                try:
                                                    sta = 'LS05' # STATION 
                                                    cha = 'EHZ' 
                                                    # PRIMARY
                                                    t1 = UTCDateTime(sec[x,0]-10) #the format is year:day_of_the_year:month
                                                    t2 = t1 + 60           
                                                    st = Stream()
                                                    st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                                                    
                                                    # Filter data
                                                    st.detrend(type='linear')
                                                    st.detrend(type='demean')
                                                    st.filter(type='bandpass',freqmin=0.2, freqmax=10)
                                                    trs=st[0]
                                                
                                                    st_c = calibrate1(trs)
                                                    B=2*pi*rhoE*cE*(1/A)
                                                    EI = sum(np.square(st_c[0].data))
                                                    EEP= B*(rs5*rs5)*EI
                                                
                                                    # SECONDARY
                                                    t1 = UTCDateTime(sec[x,1]-10) #the format is year:day_of_the_year:month
                                                    t2 = t1 + 60           
                                                    st = Stream()
                                                    st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                                                    
                                                    # Filter data
                                                    st.detrend(type='linear')
                                                    st.detrend(type='demean')
                                                    st.filter(type='bandpass',freqmin=0.2, freqmax=10)
                                                    trs=st[0]
                                                
                                                    st_c = calibrate1(trs)
                                                    B=2*pi*rhoE*cE*(1/A)
                                                    EI = sum(np.square(st_c[0].data))
                                                    EES= B*(rs5*rs5)*EI
                                                
                                                    # Percentage
                                                    
                                                    
                                                    sec[x][2]= (EES/EEP)*100
                                                except: 
                                                    try:
                                                        sta = 'LS06' # STATION 
                                                        cha = 'EHZ' 
                                                        # PRIMARY
                                                        t1 = UTCDateTime(sec[x,0]-10) #the format is year:day_of_the_year:month
                                                        t2 = t1 + 60           
                                                        st = Stream()
                                                        st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                                                        
                                                        # Filter data
                                                        st.detrend(type='linear')
                                                        st.detrend(type='demean')
                                                        st.filter(type='bandpass',freqmin=0.2, freqmax=10)
                                                        trs=st[0]
                                                    
                                                        st_c = calibrate1(trs)
                                                        B=2*pi*rhoE*cE*(1/A)
                                                        EI = sum(np.square(st_c[0].data))
                                                        EEP= B*(rs6*rs6)*EI
                                                    
                                                        # SECONDARY
                                                        t1 = UTCDateTime(sec[x,1]-10) #the format is year:day_of_the_year:month
                                                        t2 = t1 + 60           
                                                        st = Stream()
                                                        st = client.get_waveforms(net, sta, '', cha, t1 , t2)
                                                        
                                                        # Filter data
                                                        st.detrend(type='linear')
                                                        st.detrend(type='demean')
                                                        st.filter(type='bandpass',freqmin=0.2, freqmax=10)
                                                        trs=st[0]
                                                    
                                                        st_c = calibrate1(trs)
                                                        B=2*pi*rhoE*cE*(1/A)
                                                        EI = sum(np.square(st_c[0].data))
                                                        EES= B*(rs6*rs6)*EI
                                                    
                                                        # Percentage
                                                        
                                                        
                                                        sec[x][2]= (EES/EEP)*100
                                                    except: 
                                                        g=0
#%%

for y in range(0,len(sec)): 
    if sec[y][2]==0: 
        print(UTCDateTime(sec[y][0]))
        
        
        
  #%%        
count =0

secondary=np.zeros(shape=(0,3))
for y in range(0,len(sec)): 
    if 0 < sec[y][2] < 25: 
        
        secondary = np.lib.pad(secondary, ((0,1),(0,0)), 'constant', constant_values=(0))
        secondary[count][0]= sec[y][0]
        secondary[count][1]= sec[y][1]
        secondary[count][2]= sec[y][2]
        count +=1
        
        
plt.hist(secondary[:,2],50)        
print('number of secondaries =', count)  
#np.savetxt("/Users/william/Documents/scanner/analysis/secondaries_list.csv",secondary ,delimiter=",",header="")  
    
 #%%      
count =0

secondary2=np.zeros(shape=(0,3))
for y in range(0,len(sec)): 
    if 0 < sec[y][2] < 200: 
        
        secondary2 = np.lib.pad(secondary2, ((0,1),(0,0)), 'constant', constant_values=(0))
        secondary2[count][0]= sec[y][0]
        secondary2[count][1]= sec[y][1]
        secondary2[count][2]= sec[y][2]
        count +=1
 
plt.hist(secondary2[:,2],100)
plt.title('Energy of Second Explosions')
plt.xlabel('Secondary Energy Compared to Primary [%]')
plt.ylabel('Count [#]')
print('number of secondaries =', count)   
  #%%
#for y in range(0,len(secondary)): 
#    print(UTCDateTime(secondary[y,0]))
#  
  
  #%%
        
#
#reclen = 512
#chunksize = 100000 * reclen # Around 50 MB
#ad=0
#sum_events=0
#
#R_Energy=np.zeros(shape=(0,5))
#
#A=1
#rhoE=2500
#cE=2000
#pi=3.14159
#r1=4630
#r2=3370
#r3=2310
#r4=1300
#r5=810
#r6=7660
#
#E_max=1e16
#
#with io.open("/Users/william/Documents/scanner/output_data/m30.mseed", "rb") as fh:
#        # just month 2
##with io.open("/Users/william/Documents/scanner/output_data/EXP_all_data_stream_2_month_2.mseed", "rb") as fh:
#    while True:
#        with io.BytesIO() as buf:
#            c = fh.read(chunksize)
#            if not c:
#                break
#            buf.write(c)
#            buf.seek(0, 0)
#            st = obspy.read(buf)
#        
#
#        for x in range (0,int(len(st))):
#            if st[x].stats.station == "LB01":
#                r=r1
#                
#                tr = st[x]
#                st_c = calibrate1(tr)
##                    print(st_c[0])
#                B=2*pi*rhoE*cE*(1/A)
#                
#                EI = sum(np.square(st_c[0].data))
#                EE= B*(r*r)*EI
#                if EE < E_max:
#                    R_Energy = np.lib.pad(R_Energy, ((0,1),(0,0)), 'constant', constant_values=(0))
#                    R_Energy[sum_events][0]=sum_events + 1
#                    R_Energy[sum_events][1]=EE
#                    R_Energy[sum_events][2]=st[x].stats.starttime
#                    
#                    near,ind=find_nearest(pands[:,0],float(st[x].stats.starttime))
#                    R_Energy[sum_events][3]=pands[ind,1]
#                    R_Energy[sum_events][4]=pands[ind,2]
#                    
#                    sum_events += 1
#
#                
#                    
#plt.figure(10)
#plt.semilogy(R_Energy[:,4],R_Energy[:,1],'rx')
#plt.xlim([0,400])
#plt.ylim([1000,100000000000000])
#plt.xlabel("Repose time [mins]")
#plt.ylabel("Explosion Energy [J]")
#plt.title("Repose Time Vs Energy")
#
##np.savetxt("/Users/william/Documents/scanner/analysis/repose_energy_lb01.csv", R_Energy,delimiter=",",header="event_number,energy,time,p_or_s,repose time")
#
##%%
#reclen = 512
#chunksize = 100000 * reclen # Around 50 MB
#ad=0
#sum_events=0
#
#R2_Energy=np.zeros(shape=(0,5))
#
#A=1
#rhoE=2500
#cE=2000
#pi=3.14159
#r1=4630
#r2=3370
#r3=2310
#r4=1300
#r5=810
#r6=7660
#
#E_max=1e16
#
#with io.open("/Users/william/Documents/scanner/output_data/m30.mseed", "rb") as fh:
#        # just month 2
##with io.open("/Users/william/Documents/scanner/output_data/EXP_all_data_stream_2_month_2.mseed", "rb") as fh:
#    while True:
#        with io.BytesIO() as buf:
#            c = fh.read(chunksize)
#            if not c:
#                break
#            buf.write(c)
#            buf.seek(0, 0)
#            st = obspy.read(buf)
#        
#
#        for x in range (0,int(len(st))):
#            if st[x].stats.station == "LB02":
#                r=r1
#                
#                tr = st[x]
#                st_c = calibrate1(tr)
##                    print(st_c[0])
#                B=2*pi*rhoE*cE*(1/A)
#                
#                EI = sum(np.square(st_c[0].data))
#                EE= B*(r*r)*EI
#                if EE < E_max:
#                    R2_Energy = np.lib.pad(R2_Energy, ((0,1),(0,0)), 'constant', constant_values=(0))
#                    R2_Energy[sum_events][0]=sum_events + 1
#                    R2_Energy[sum_events][1]=EE
#                    R2_Energy[sum_events][2]=st[x].stats.starttime
#                    
#                    near,ind=find_nearest(pands[:,0],float(st[x].stats.starttime))
#                    R2_Energy[sum_events][3]=pands[ind,1]
#                    R2_Energy[sum_events][4]=pands[ind,2]
#                    
#                    sum_events += 1
#
#                
#                    
#plt.figure(11)
#plt.semilogy(R2_Energy[:,4],R2_Energy[:,1],'bx')
#plt.xlim([0,300])
#plt.ylim([1000,100000000000000])
#plt.xlabel("Repose time [mins]")
#plt.ylabel("Explosion Energy [J]")
#plt.title("Repose Time Vs Energy")
#
##np.savetxt("/Users/william/Documents/scanner/analysis/repose_energy_lb02.csv", R2_Energy,delimiter=",",header="event_number,energy,time,p_or_s,repose time")
#
##%%
#reclen = 512
#chunksize = 100000 * reclen # Around 50 MB
#ad=0
#sum_events=0
#
#R3_Energy=np.zeros(shape=(0,5))
#
#A=1
#rhoE=2500
#cE=2000
#pi=3.14159
#r1=4630
#r2=3370
#r3=2310
#r4=1300
#r5=810
#r6=7660
#
#E_max=1e16
#
#with io.open("/Users/william/Documents/scanner/output_data/m30.mseed", "rb") as fh:
#        # just month 2
##with io.open("/Users/william/Documents/scanner/output_data/EXP_all_data_stream_2_month_2.mseed", "rb") as fh:
#    while True:
#        with io.BytesIO() as buf:
#            c = fh.read(chunksize)
#            if not c:
#                break
#            buf.write(c)
#            buf.seek(0, 0)
#            st = obspy.read(buf)
#        
#
#        for x in range (0,int(len(st))):
#            if st[x].stats.station == "LB03":
#                r=r1
#                
#                tr = st[x]
#                st_c = calibrate1(tr)
##                    print(st_c[0])
#                B=2*pi*rhoE*cE*(1/A)
#                
#                EI = sum(np.square(st_c[0].data))
#                EE= B*(r*r)*EI
#                if EE < E_max:
#                    R3_Energy = np.lib.pad(R3_Energy, ((0,1),(0,0)), 'constant', constant_values=(0))
#                    R3_Energy[sum_events][0]=sum_events + 1
#                    R3_Energy[sum_events][1]=EE
#                    R3_Energy[sum_events][2]=st[x].stats.starttime
#                    
#                    near,ind=find_nearest(pands[:,0],float(st[x].stats.starttime))
#                    R3_Energy[sum_events][3]=pands[ind,1]
#                    R3_Energy[sum_events][4]=pands[ind,2]
#                    
#                    sum_events += 1
#
#               
#                    
#plt.figure(12)
#plt.semilogy(R3_Energy[:,4],R3_Energy[:,1],'gx')
#plt.xlim([0,300])
#plt.ylim([1000,10000000000])
#plt.xlabel("Repose time [mins]")
#plt.ylabel("Explosion Energy [J]")
#plt.title("Repose Time Vs Energy")
#
##np.savetxt("/Users/william/Documents/scanner/analysis/repose_energy_lb03.csv", R3_Energy,delimiter=",",header="event_number,energy,time,p_or_s,repose time")
