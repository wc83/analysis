#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 15:31:30 2019

@author: william
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 16:42:44 2019

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

#%% constants

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


#%% read catalogue of events

f=open("/Users/william/Documents/scanner/Olivers_Final_exp_2014_2017.txt", "r")

cat =f.read()

events=np.zeros(shape=(6104,1))
year1=np.zeros(shape=(6104,1))
month1=np.zeros(shape=(6104,1))
day1=np.zeros(shape=(6104,1))
hour1=np.zeros(shape=(6104,1))
minute1=np.zeros(shape=(6104,1))
second1=np.zeros(shape=(6104,1))
num=0

for x in range(0,len(cat),20):  
    
    year = int(float(cat[x])*1000 + float(cat[x+1])*100 + float(cat[x+2])*10 + float(cat[x+3]))
    month = int(float(cat[x+5])*10+float(cat[x+6]))
    day = int(float(cat[x+8])*10 + float(cat[x+9])) 

    
    hour = int(float(cat[x+11])*10 + float(cat[x+12])) 
    minute = int(float(cat[x+14])*10 + float(cat[x+15]))
    second = int(float(cat[x+17])*10 + float(cat[x+18]))


    events[num] = UTCDateTime(year, month, day, hour, minute, second)
    year1[num] = year
    month1[num] = month
    day1[num]= day
    hour1[num]=hour
    minute1[num]=minute
    second1[num]=second
    
    
    
    
    num+=1









#%%
Energy_trip = np.zeros(shape=(len(events),2))
Energy_trip[:,0]=events[:,0]


#%% get all 2014-2017 data
net = 'Z4' 
st_n = Stream()
st_e = Stream() 
client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23

for x in range(0,len(events)):

        if  Energy_trip[x,1] ==0 :

            t1=UTCDateTime(events[x] -10)
            t2=t1+80
            try:
                r=r1
                sta='LB01'
                cha1='HHZ'
                cha2='HHE'
                cha3='HHN'
                
                st_z = client.get_waveforms(net, sta, '', cha1, t1 , t2)
                trsz= st_z[0]
                trsz.detrend(type='linear')
                trsz.detrend(type='demean')
                trsz.filter(type='bandpass',freqmin=0.1, freqmax=10)
                st_acz = calibrate1(trsz)
                    
                st_e = client.get_waveforms(net, sta, '', cha2, t1 , t2)
                trse= st_e[0]
                trse.detrend(type='linear')
                trse.detrend(type='demean')
                trse.filter(type='bandpass',freqmin=0.1, freqmax=10)
                st_ace = calibrate1(trse)
                
                st_n = client.get_waveforms(net, sta, '', cha3, t1 , t2)
                trsn= st_n[0]
                trsn.detrend(type='linear')
                trsn.detrend(type='demean')
                trsn.filter(type='bandpass',freqmin=0.1, freqmax=10)
                st_acn = calibrate1(trsn)
                

                B=2*pi*rhoE*cE*(1/A)
                dl=len(st_acz[0].data)
                p = np.linspace(0,dl/100, num=dl)
                
                y= np.sqrt(st_acz[0].data**2 + st_ace[0].data**2 + st_acn[0].data**2)
                y2= np.square(y)
                y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
                EI = y_int2[-1]
                EE1= B*(r*r)*EI
                
                Energy_trip[x,1]= EE1
                
                
            except:
                try:
                    r=r2
                    sta='LB02'
                    cha1='HHZ'
                    cha2='HHE'
                    cha3='HHN'
                    
                    st_z = client.get_waveforms(net, sta, '', cha1, t1 , t2)
                    trsz= st_z[0]
                    trsz.detrend(type='linear')
                    trsz.detrend(type='demean')
                    trsz.filter(type='bandpass',freqmin=0.1, freqmax=10)
                    st_acz = calibrate1(trsz)
                        
                    st_e = client.get_waveforms(net, sta, '', cha2, t1 , t2)
                    trse= st_e[0]
                    trse.detrend(type='linear')
                    trse.detrend(type='demean')
                    trse.filter(type='bandpass',freqmin=0.1, freqmax=10)
                    st_ace = calibrate1(trse)
                    
                    st_n = client.get_waveforms(net, sta, '', cha3, t1 , t2)
                    trsn= st_n[0]
                    trsn.detrend(type='linear')
                    trsn.detrend(type='demean')
                    trsn.filter(type='bandpass',freqmin=0.1, freqmax=10)
                    st_acn = calibrate1(trsn)
                    
    
                    B=2*pi*rhoE*cE*(1/A)
                    dl=len(st_acz[0].data)
                    p = np.linspace(0,dl/100, num=dl)
                    
                    y= np.sqrt(st_acz[0].data**2 + st_ace[0].data**2 + st_acn[0].data**2)
                    y2= np.square(y)
                    y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
                    EI = y_int2[-1]
                    EE1= B*(r*r)*EI
                    
                    Energy_trip[x,1]= EE1
                    
                    
                except:
                    try:
                        r=r3
                        sta='LB03'
                        cha1='HHZ'
                        cha2='HHE'
                        cha3='HHN'
                        
                        st_z = client.get_waveforms(net, sta, '', cha1, t1 , t2)
                        trsz= st_z[0]
                        trsz.detrend(type='linear')
                        trsz.detrend(type='demean')
                        trsz.filter(type='bandpass',freqmin=0.1, freqmax=10)
                        st_acz = calibrate1(trsz)
                            
                        st_e = client.get_waveforms(net, sta, '', cha2, t1 , t2)
                        trse= st_e[0]
                        trse.detrend(type='linear')
                        trse.detrend(type='demean')
                        trse.filter(type='bandpass',freqmin=0.1, freqmax=10)
                        st_ace = calibrate1(trse)
                        
                        st_n = client.get_waveforms(net, sta, '', cha3, t1 , t2)
                        trsn= st_n[0]
                        trsn.detrend(type='linear')
                        trsn.detrend(type='demean')
                        trsn.filter(type='bandpass',freqmin=0.1, freqmax=10)
                        st_acn = calibrate1(trsn)
                        
        
                        B=2*pi*rhoE*cE*(1/A)
                        dl=len(st_acz[0].data)
                        p = np.linspace(0,dl/100, num=dl)
                        
                        y= np.sqrt(st_acz[0].data**2 + st_ace[0].data**2 + st_acn[0].data**2)
                        y2= np.square(y)
                        y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
                        EI = y_int2[-1]
                        EE1= B*(r*r)*EI
                        
                        Energy_trip[x,1]= EE1
                    except:
                        try:
                            r=r6
                            sta='LB06'
                            cha1='HHZ'
                            cha2='HHE'
                            cha3='HHN'
                            
                            st_z = client.get_waveforms(net, sta, '', cha1, t1 , t2)
                            trsz= st_z[0]
                            trsz.detrend(type='linear')
                            trsz.detrend(type='demean')
                            trsz.filter(type='bandpass',freqmin=0.1, freqmax=10)
                            st_acz = calibrate1(trsz)
                                
                            st_e = client.get_waveforms(net, sta, '', cha2, t1 , t2)
                            trse= st_e[0]
                            trse.detrend(type='linear')
                            trse.detrend(type='demean')
                            trse.filter(type='bandpass',freqmin=0.1, freqmax=10)
                            st_ace = calibrate1(trse)
                            
                            st_n = client.get_waveforms(net, sta, '', cha3, t1 , t2)
                            trsn= st_n[0]
                            trsn.detrend(type='linear')
                            trsn.detrend(type='demean')
                            trsn.filter(type='bandpass',freqmin=0.1, freqmax=10)
                            st_acn = calibrate1(trsn)
                            
            
                            B=2*pi*rhoE*cE*(1/A)
                            dl=len(st_acz[0].data)
                            p = np.linspace(0,dl/100, num=dl)
                            
                            y= np.sqrt(st_acz[0].data**2 + st_ace[0].data**2 + st_acn[0].data**2)
                            y2= np.square(y)
                            y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
                            EI = y_int2[-1]
                            EE1= B*(r*r)*EI
                            
                            Energy_trip[x,1]= EE1   
                        except:
                            try:
                                r=rs1
                                sta='LS01'
                                cha1='EHZ'
                                cha2='EHE'
                                cha3='EHN'
                                
                                st_z = client.get_waveforms(net, sta, '', cha1, t1 , t2)
                                trsz= st_z[0]
                                trsz.detrend(type='linear')
                                trsz.detrend(type='demean')
                                trsz.filter(type='bandpass',freqmin=0.1, freqmax=10)
                                st_acz = calibrate1(trsz)
                                    
                                st_e = client.get_waveforms(net, sta, '', cha2, t1 , t2)
                                trse= st_e[0]
                                trse.detrend(type='linear')
                                trse.detrend(type='demean')
                                trse.filter(type='bandpass',freqmin=0.1, freqmax=10)
                                st_ace = calibrate1(trse)
                                
                                st_n = client.get_waveforms(net, sta, '', cha3, t1 , t2)
                                trsn= st_n[0]
                                trsn.detrend(type='linear')
                                trsn.detrend(type='demean')
                                trsn.filter(type='bandpass',freqmin=0.1, freqmax=10)
                                st_acn = calibrate1(trsn)
                                
                
                                B=2*pi*rhoE*cE*(1/A)
                                dl=len(st_acz[0].data)
                                p = np.linspace(0,dl/100, num=dl)
                                
                                y= np.sqrt(st_acz[0].data**2 + st_ace[0].data**2 + st_acn[0].data**2)
                                y2= np.square(y)
                                y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
                                EI = y_int2[-1]
                                EE1= B*(r*r)*EI
                                
                                Energy_trip[x,1]= EE1
                            except:
                                try:
                                    r=rs2
                                    sta='LS02'
                                    cha1='EHZ'
                                    cha2='EHE'
                                    cha3='EHN'
                                    
                                    st_z = client.get_waveforms(net, sta, '', cha1, t1 , t2)
                                    trsz= st_z[0]
                                    trsz.detrend(type='linear')
                                    trsz.detrend(type='demean')
                                    trsz.filter(type='bandpass',freqmin=0.1, freqmax=10)
                                    st_acz = calibrate1(trsz)
                                        
                                    st_e = client.get_waveforms(net, sta, '', cha2, t1 , t2)
                                    trse= st_e[0]
                                    trse.detrend(type='linear')
                                    trse.detrend(type='demean')
                                    trse.filter(type='bandpass',freqmin=0.1, freqmax=10)
                                    st_ace = calibrate1(trse)
                                    
                                    st_n = client.get_waveforms(net, sta, '', cha3, t1 , t2)
                                    trsn= st_n[0]
                                    trsn.detrend(type='linear')
                                    trsn.detrend(type='demean')
                                    trsn.filter(type='bandpass',freqmin=0.1, freqmax=10)
                                    st_acn = calibrate1(trsn)
                                    
                    
                                    B=2*pi*rhoE*cE*(1/A)
                                    dl=len(st_acz[0].data)
                                    p = np.linspace(0,dl/100, num=dl)
                                    
                                    y= np.sqrt(st_acz[0].data**2 + st_ace[0].data**2 + st_acn[0].data**2)
                                    y2= np.square(y)
                                    y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
                                    EI = y_int2[-1]
                                    EE1= B*(r*r)*EI
                                    
                                    Energy_trip[x,1]= EE1
                                except:
                                    try:
                                        r=rs3
                                        sta='LS03'
                                        cha1='EHZ'
                                        cha2='EHE'
                                        cha3='EHN'
                                        
                                        st_z = client.get_waveforms(net, sta, '', cha1, t1 , t2)
                                        trsz= st_z[0]
                                        trsz.detrend(type='linear')
                                        trsz.detrend(type='demean')
                                        trsz.filter(type='bandpass',freqmin=0.1, freqmax=10)
                                        st_acz = calibrate1(trsz)
                                            
                                        st_e = client.get_waveforms(net, sta, '', cha2, t1 , t2)
                                        trse= st_e[0]
                                        trse.detrend(type='linear')
                                        trse.detrend(type='demean')
                                        trse.filter(type='bandpass',freqmin=0.1, freqmax=10)
                                        st_ace = calibrate1(trse)
                                        
                                        st_n = client.get_waveforms(net, sta, '', cha3, t1 , t2)
                                        trsn= st_n[0]
                                        trsn.detrend(type='linear')
                                        trsn.detrend(type='demean')
                                        trsn.filter(type='bandpass',freqmin=0.1, freqmax=10)
                                        st_acn = calibrate1(trsn)
                                        
                        
                                        B=2*pi*rhoE*cE*(1/A)
                                        dl=len(st_acz[0].data)
                                        p = np.linspace(0,dl/100, num=dl)
                                        
                                        y= np.sqrt(st_acz[0].data**2 + st_ace[0].data**2 + st_acn[0].data**2)
                                        y2= np.square(y)
                                        y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
                                        EI = y_int2[-1]
                                        EE1= B*(r*r)*EI
                                        
                                        Energy_trip[x,1]= EE1
                                    except:
                                        try:
                                            r=rs4
                                            sta='LS04'
                                            cha1='EHZ'
                                            cha2='EHE'
                                            cha3='EHN'
                                            
                                            st_z = client.get_waveforms(net, sta, '', cha1, t1 , t2)
                                            trsz= st_z[0]
                                            trsz.detrend(type='linear')
                                            trsz.detrend(type='demean')
                                            trsz.filter(type='bandpass',freqmin=0.1, freqmax=10)
                                            st_acz = calibrate1(trsz)
                                                
                                            st_e = client.get_waveforms(net, sta, '', cha2, t1 , t2)
                                            trse= st_e[0]
                                            trse.detrend(type='linear')
                                            trse.detrend(type='demean')
                                            trse.filter(type='bandpass',freqmin=0.1, freqmax=10)
                                            st_ace = calibrate1(trse)
                                            
                                            st_n = client.get_waveforms(net, sta, '', cha3, t1 , t2)
                                            trsn= st_n[0]
                                            trsn.detrend(type='linear')
                                            trsn.detrend(type='demean')
                                            trsn.filter(type='bandpass',freqmin=0.1, freqmax=10)
                                            st_acn = calibrate1(trsn)
                                            
                            
                                            B=2*pi*rhoE*cE*(1/A)
                                            dl=len(st_acz[0].data)
                                            p = np.linspace(0,dl/100, num=dl)
                                            
                                            y= np.sqrt(st_acz[0].data**2 + st_ace[0].data**2 + st_acn[0].data**2)
                                            y2= np.square(y)
                                            y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
                                            EI = y_int2[-1]
                                            EE1= B*(r*r)*EI
                                            
                                            Energy_trip[x,1]= EE1
                                        except:
                                            try:
                                                r=rs5
                                                sta='LS05'
                                                cha1='EHZ'
                                                cha2='EHE'
                                                cha3='EHN'
                                                
                                                st_z = client.get_waveforms(net, sta, '', cha1, t1 , t2)
                                                trsz= st_z[0]
                                                trsz.detrend(type='linear')
                                                trsz.detrend(type='demean')
                                                trsz.filter(type='bandpass',freqmin=0.1, freqmax=10)
                                                st_acz = calibrate1(trsz)
                                                    
                                                st_e = client.get_waveforms(net, sta, '', cha2, t1 , t2)
                                                trse= st_e[0]
                                                trse.detrend(type='linear')
                                                trse.detrend(type='demean')
                                                trse.filter(type='bandpass',freqmin=0.1, freqmax=10)
                                                st_ace = calibrate1(trse)
                                                
                                                st_n = client.get_waveforms(net, sta, '', cha3, t1 , t2)
                                                trsn= st_n[0]
                                                trsn.detrend(type='linear')
                                                trsn.detrend(type='demean')
                                                trsn.filter(type='bandpass',freqmin=0.1, freqmax=10)
                                                st_acn = calibrate1(trsn)
                                                
                                
                                                B=2*pi*rhoE*cE*(1/A)
                                                dl=len(st_acz[0].data)
                                                p = np.linspace(0,dl/100, num=dl)
                                                
                                                y= np.sqrt(st_acz[0].data**2 + st_ace[0].data**2 + st_acn[0].data**2)
                                                y2= np.square(y)
                                                y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
                                                EI = y_int2[-1]
                                                EE1= B*(r*r)*EI
                                                
                                                Energy_trip[x,1]= EE1
                                            except:
                                                try:
                                                    r=rs6
                                                    sta='LS06'
                                                    cha1='EHZ'
                                                    cha2='EHE'
                                                    cha3='EHN'
                                                    
                                                    st_z = client.get_waveforms(net, sta, '', cha1, t1 , t2)
                                                    trsz= st_z[0]
                                                    trsz.detrend(type='linear')
                                                    trsz.detrend(type='demean')
                                                    trsz.filter(type='bandpass',freqmin=0.1, freqmax=10)
                                                    st_acz = calibrate1(trsz)
                                                        
                                                    st_e = client.get_waveforms(net, sta, '', cha2, t1 , t2)
                                                    trse= st_e[0]
                                                    trse.detrend(type='linear')
                                                    trse.detrend(type='demean')
                                                    trse.filter(type='bandpass',freqmin=0.1, freqmax=10)
                                                    st_ace = calibrate1(trse)
                                                    
                                                    st_n = client.get_waveforms(net, sta, '', cha3, t1 , t2)
                                                    trsn= st_n[0]
                                                    trsn.detrend(type='linear')
                                                    trsn.detrend(type='demean')
                                                    trsn.filter(type='bandpass',freqmin=0.1, freqmax=10)
                                                    st_acn = calibrate1(trsn)
                                                    
                                    
                                                    B=2*pi*rhoE*cE*(1/A)
                                                    dl=len(st_acz[0].data)
                                                    p = np.linspace(0,dl/100, num=dl)
                                                    
                                                    y= np.sqrt(st_acz[0].data**2 + st_ace[0].data**2 + st_acn[0].data**2)
                                                    y2= np.square(y)
                                                    y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
                                                    EI = y_int2[-1]
                                                    EE1= B*(r*r)*EI
                                                    
                                                    Energy_trip[x,1]= EE1
                                                       
                                                except:
                                                        print('still not found', UTCDateTime(events[x]))
                            
                                
                                
                                
        

#%%  
#E1 = sum(Energy_trip[:,1])
#EE = sum(Energy_trip[:,4])
#Er=EE/E1
#print(EE/E1) 
#for x in range(0,len(Energy_trip)):
#     if Energy_trip[x,2] == 0:
#         Energy_trip[x,1] = cat[x,23]
#         Energy_trip[x,4] = cat[x,23]*Er


#%%

for x in range(0,len(Energy_trip)):
    if Energy_trip[x,1] ==0:
        print(Energy_trip[x,0])

        
#%%
        
Energy3 = Energy_trip[:,1]
        
new_cat=np.zeros(shape=(len(events),8))
new_cat[:,0]=events[:,0]
new_cat[:,1]=year1[:,0]
new_cat[:,2]=month1[:,0]
new_cat[:,3]=day1[:,0]
new_cat[:,4]=hour1[:,0]
new_cat[:,5]=minute1[:,0]
new_cat[:,6]=second1[:,0]
new_cat[:,7]=Energy3

        
 



#%%
        
#np.savetxt("/Users/william/Documents/scanner/Olivers_Catalogue_2014_2017_Energy.csv", new_cat,delimiter=",",header="Time_stamp,year,month,day,hour,min,sec,Energy")

#%%  
#for x in range(0,len(new_cat)):
#    if new_cat[x,23] ==0 or new_cat[x,24] ==0:
#        print(new_cat[x,0])
#  
  
  
  
  
  
  
  
  
  
  
  
  