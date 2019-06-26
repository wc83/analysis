#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 10:20:29 2019

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
from numpy.polynomial.polynomial import polyfit

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

cat = genfromtxt("/Users/william/Documents/scanner/all_stations/Final_Catalogue_2014_2018.csv", delimiter=',',skip_header=1)
ev_t = cat[:,0]
day_num= cat[:,1]
doy=cat[:,2]
year= cat[:,3]
month = cat[:,4]
day = cat[:,5]
hour = cat[:,6]
minute = cat[:,7]
second = cat[:,8]
milisec = cat[:,9]
lb01 = cat[:,10]
lb02 = cat[:,11]
lb03 = cat[:,12]
lb04 = cat[:,13]
lb05 = cat[:,14]
lb06 = cat[:,15]
ls01 = cat[:,16]
ls02 = cat[:,17]
ls03 = cat[:,18]
ls04 = cat[:,19]
ls05 = cat[:,20]
ls06 = cat[:,21]
trust = cat[:,22]
Energy = cat[:,23]
dur = cat[:,24]


#%%
ac = 2
el = 2*ac
ol = (2*ac) +1

ran = len(cat)
nume_test=0
for x in range(0,ran):
    if lb01[x] == 1 and  lb04[x] == 1: # and lb03[x] ==1 :
#        print(UTCDateTime(ev_t[x]))
        
        
#        active = lb01[x] + lb02[x] + lb03[x] + lb04[x] + lb05[x] + lb06[x] + ls01[x] + ls02[x] + ls03[x] + ls04[x] + ls05[x] + ls06[x]
#        if  active == ac :
#            print(UTCDateTime(ev_t[x]))
        nume_test += 1

print(nume_test)

#%%
ed = np.zeros(shape=(nume_test,el))
enl = np.zeros(shape=(0,1))

al = np.zeros(shape=(0,1))
ad = np.zeros(shape=(nume_test,el))

numb=0
numl=-1
net = 'Z4' 

lb1=0
lb2=0
lb3=0
lb4=0
lb5=0
lb6=0
ls1=0
ls2=0
ls3=0
ls4=0
ls5=0
ls6=0

lba1=0
lba2=0
lba3=0
lba4=0
lba5=0
lba6=0
lsa1=0
lsa2=0
lsa3=0
lsa4=0
lsa5=0
lsa6=0

lb1n=0
lb2n=0
lb3n=0
lb4n=0
lb5n=0
lb6n=0
ls1n=0
ls2n=0
ls3n=0
ls4n=0
ls5n=0
ls6n=0



client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23

nume=0

for x in range(0,ran):
#    active = lb01[x] + lb02[x] + lb03[x] + lb04[x] + lb05[x] + lb06[x] + ls01[x] + ls02[x] + ls03[x] + ls04[x] + ls05[x] + ls06[x]
   if lb01[x] == 1 and  lb04[x] == 1:
       
       try:
#            print(UTCDateTime(ev_t[x]))
            numb=0
            nume += 1
            t1=UTCDateTime(ev_t[x] -10)
            t2=UTCDateTime(ev_t[x] +70)
            numl+=1
            if lb01[x] == 1:
                
                sta='LB01'
                cha1='HHZ'
                cha2='HHN'
                cha3='HHE'
                r=r1
                
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
                EE= B*(r*r)*EI
                
                En1= EE
                A1 = EI
                
                enl = np.lib.pad(enl, ((0,1),(0,0)), 'constant', constant_values=(0))
                enl[numb][0] = r
                al = np.lib.pad(al, ((0,1),(0,0)), 'constant', constant_values=(0))
                al[numb][0] = r
                numb+=1
                enl = np.lib.pad(enl, ((0,1),(0,0)), 'constant', constant_values=(0))
                enl[numb][0] = En1
                al = np.lib.pad(al, ((0,1),(0,0)), 'constant', constant_values=(0))
                al[numb][0] = A1
                numb+=1
                
                lb1 += En1
                lba1 += A1
                lb1n += 1
            
#            if lb02[x] == 1:
#                sta='LB02'
#                cha1='HHZ'
#                cha2='HHN'
#                cha3='HHE'
#                r=r2
#                
#                st_z = client.get_waveforms(net, sta, '', cha1, t1 , t2)
#                trsz= st_z[0]
#                trsz.detrend(type='linear')
#                trsz.detrend(type='demean')
#                trsz.filter(type='bandpass',freqmin=0.1, freqmax=10)
#            
#                st_acz = calibrate1(trsz)
#                
#                st_e = client.get_waveforms(net, sta, '', cha2, t1 , t2)
#                trse= st_e[0]
#                trse.detrend(type='linear')
#                trse.detrend(type='demean')
#                trse.filter(type='bandpass',freqmin=0.1, freqmax=10)
#            
#                st_ace = calibrate1(trse)
#                
#                st_n = client.get_waveforms(net, sta, '', cha3, t1 , t2)
#                trsn= st_n[0]
#                trsn.detrend(type='linear')
#                trsn.detrend(type='demean')
#                trsn.filter(type='bandpass',freqmin=0.1, freqmax=10)
#            
#                st_acn = calibrate1(trsn)
#                
#                
#                B=2*pi*rhoE*cE*(1/A)
#                dl=len(st_acz[0].data)
#                p = np.linspace(0,dl/100, num=dl)
#                
#                y= np.sqrt(st_acz[0].data**2 + st_ace[0].data**2 + st_acn[0].data**2)
#                y2= np.square(y)
#                y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
#                           
#                EI = y_int2[-1]
#                EE= B*(r*r)*EI
#                
#                En2= EE
#                A2 = EI
#                
#                enl = np.lib.pad(enl, ((0,1),(0,0)), 'constant', constant_values=(0))
#                enl[numb][0] = r
#                al = np.lib.pad(al, ((0,1),(0,0)), 'constant', constant_values=(0))
#                al[numb][0] = r
#                numb+=1
#                enl = np.lib.pad(enl, ((0,1),(0,0)), 'constant', constant_values=(0))
#                enl[numb][0] = En2
#                al = np.lib.pad(al, ((0,1),(0,0)), 'constant', constant_values=(0))
#                al[numb][0] = A2
#                numb+=1
#                
#                lb2 += En2
#                lba2 += A2
#                lb2n += 1
#                
#            if lb03[x] == 1:
#                sta='LB03'
#                cha1='HHZ'
#                cha2='HHN'
#                cha3='HHE'
#                r=r3
#                
#                st_z = client.get_waveforms(net, sta, '', cha1, t1 , t2)
#                trsz= st_z[0]
#                trsz.detrend(type='linear')
#                trsz.detrend(type='demean')
#                trsz.filter(type='bandpass',freqmin=0.1, freqmax=10)
#            
#                st_acz = calibrate1(trsz)
#                
#                st_e = client.get_waveforms(net, sta, '', cha2, t1 , t2)
#                trse= st_e[0]
#                trse.detrend(type='linear')
#                trse.detrend(type='demean')
#                trse.filter(type='bandpass',freqmin=0.1, freqmax=10)
#            
#                st_ace = calibrate1(trse)
#                
#                st_n = client.get_waveforms(net, sta, '', cha3, t1 , t2)
#                trsn= st_n[0]
#                trsn.detrend(type='linear')
#                trsn.detrend(type='demean')
#                trsn.filter(type='bandpass',freqmin=0.1, freqmax=10)
#            
#                st_acn = calibrate1(trsn)
#                
#                
#                B=2*pi*rhoE*cE*(1/A)
#                dl=len(st_acz[0].data)
#                p = np.linspace(0,dl/100, num=dl)
#                
#                y= np.sqrt(st_acz[0].data**2 + st_ace[0].data**2 + st_acn[0].data**2)
#                y2= np.square(y)
#                y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
#                           
#                EI = y_int2[-1]
#                EE= B*(r*r)*EI
#                
#                En3= EE
#                A3 = EI
#                
#                enl = np.lib.pad(enl, ((0,1),(0,0)), 'constant', constant_values=(0))
#                enl[numb][0] = r
#                al = np.lib.pad(al, ((0,1),(0,0)), 'constant', constant_values=(0))
#                al[numb][0] = r
#                numb+=1
#                enl = np.lib.pad(enl, ((0,1),(0,0)), 'constant', constant_values=(0))
#                enl[numb][0] = En3
#                al = np.lib.pad(al, ((0,1),(0,0)), 'constant', constant_values=(0))
#                al[numb][0] = A3
#                numb+=1
#                
#                lb3 += En3
#                lba3 += A3
#                lb3n += 1
#                
            if lb04[x] == 1:
                sta='LB04'
                cha1='HHZ'
                cha2='HHN'
                cha3='HHE'
                r=r4
                
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
                EE= B*(r*r)*EI
                
                En4= EE
                A4 = EI
                
                enl = np.lib.pad(enl, ((0,1),(0,0)), 'constant', constant_values=(0))
                enl[numb][0] = r
                al = np.lib.pad(al, ((0,1),(0,0)), 'constant', constant_values=(0))
                al[numb][0] = r
                numb+=1
                enl = np.lib.pad(enl, ((0,1),(0,0)), 'constant', constant_values=(0))
                enl[numb][0] = En4
                al = np.lib.pad(al, ((0,1),(0,0)), 'constant', constant_values=(0))
                al[numb][0] = A4
                numb+=1
                
                lb4 += En4
                lba4 += A4
                lb4n += 1
#                
#            if lb05[x] == 1:
#                sta='LB05'
#                cha1='HHZ'
#                cha2='HHN'
#                cha3='HHE'
#                r=r5
#                
#                st_z = client.get_waveforms(net, sta, '', cha1, t1 , t2)
#                trsz= st_z[0]
#                trsz.detrend(type='linear')
#                trsz.detrend(type='demean')
#                trsz.filter(type='bandpass',freqmin=0.1, freqmax=10)
#            
#                st_acz = calibrate1(trsz)
#                
#                st_e = client.get_waveforms(net, sta, '', cha2, t1 , t2)
#                trse= st_e[0]
#                trse.detrend(type='linear')
#                trse.detrend(type='demean')
#                trse.filter(type='bandpass',freqmin=0.1, freqmax=10)
#            
#                st_ace = calibrate1(trse)
#                
#                st_n = client.get_waveforms(net, sta, '', cha3, t1 , t2)
#                trsn= st_n[0]
#                trsn.detrend(type='linear')
#                trsn.detrend(type='demean')
#                trsn.filter(type='bandpass',freqmin=0.1, freqmax=10)
#            
#                st_acn = calibrate1(trsn)
#                
#                
#                B=2*pi*rhoE*cE*(1/A)
#                dl=len(st_acz[0].data)
#                p = np.linspace(0,dl/100, num=dl)
#                
#                y= np.sqrt(st_acz[0].data**2 + st_ace[0].data**2 + st_acn[0].data**2)
#                y2= np.square(y)
#                y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
#                           
#                EI = y_int2[-1]
#                EE= B*(r*r)*EI
#                
#                En5= EE
#                A5 = EI
#                
#                enl = np.lib.pad(enl, ((0,1),(0,0)), 'constant', constant_values=(0))
#                enl[numb][0] = r
#                al = np.lib.pad(al, ((0,1),(0,0)), 'constant', constant_values=(0))
#                al[numb][0] = r
#                numb+=1
#                enl = np.lib.pad(enl, ((0,1),(0,0)), 'constant', constant_values=(0))
#                enl[numb][0] = En5
#                al = np.lib.pad(al, ((0,1),(0,0)), 'constant', constant_values=(0))
#                al[numb][0] = A5
#                numb+=1
#                
#                
#                lb5 += En5
#                lba5 += A5
#                lb5n += 1
#                
#            if lb06[x] == 1:
#                sta='LB06'
#                cha1='HHZ'
#                cha2='HHN'
#                cha3='HHE'
#                r=r6
#                
#                st_z = client.get_waveforms(net, sta, '', cha1, t1 , t2)
#                trsz= st_z[0]
#                trsz.detrend(type='linear')
#                trsz.detrend(type='demean')
#                trsz.filter(type='bandpass',freqmin=0.1, freqmax=10)
#            
#                st_acz = calibrate1(trsz)
#                
#                st_e = client.get_waveforms(net, sta, '', cha2, t1 , t2)
#                trse= st_e[0]
#                trse.detrend(type='linear')
#                trse.detrend(type='demean')
#                trse.filter(type='bandpass',freqmin=0.1, freqmax=10)
#            
#                st_ace = calibrate1(trse)
#                
#                st_n = client.get_waveforms(net, sta, '', cha3, t1 , t2)
#                trsn= st_n[0]
#                trsn.detrend(type='linear')
#                trsn.detrend(type='demean')
#                trsn.filter(type='bandpass',freqmin=0.1, freqmax=10)
#            
#                st_acn = calibrate1(trsn)
#                
#                
#                B=2*pi*rhoE*cE*(1/A)
#                dl=len(st_acz[0].data)
#                p = np.linspace(0,dl/100, num=dl)
#                
#                y= np.sqrt(st_acz[0].data**2 + st_ace[0].data**2 + st_acn[0].data**2)
#                y2= np.square(y)
#                y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
#                           
#                EI = y_int2[-1]
#                EE= B*(r*r)*EI
#                
#                En6= EE
#                A6 = EI
#                
#                enl = np.lib.pad(enl, ((0,1),(0,0)), 'constant', constant_values=(0))
#                enl[numb][0] = r6
#                al = np.lib.pad(al, ((0,1),(0,0)), 'constant', constant_values=(0))
#                al[numb][0] = r6
#                numb+=1
#                enl = np.lib.pad(enl, ((0,1),(0,0)), 'constant', constant_values=(0))
#                enl[numb][0] = En6
#                al = np.lib.pad(al, ((0,1),(0,0)), 'constant', constant_values=(0))
#                al[numb][0] = A6
#                numb+=1
#                
#                lb6 += En6
#                lba6 += A6
#                lb6n += 1
#                
#            if ls01[x] == 1:
#                sta='LS01'
#                cha1='EHZ'
#                cha2='EHN'
#                cha3='EHE'
#                r=rs1
#                
#                st_z = client.get_waveforms(net, sta, '', cha1, t1 , t2)
#                trsz= st_z[0]
#                trsz.detrend(type='linear')
#                trsz.detrend(type='demean')
#                trsz.filter(type='bandpass',freqmin=0.1, freqmax=10)
#            
#                st_acz = calibrate1(trsz)
#                
#                st_e = client.get_waveforms(net, sta, '', cha2, t1 , t2)
#                trse= st_e[0]
#                trse.detrend(type='linear')
#                trse.detrend(type='demean')
#                trse.filter(type='bandpass',freqmin=0.1, freqmax=10)
#            
#                st_ace = calibrate1(trse)
#                
#                st_n = client.get_waveforms(net, sta, '', cha3, t1 , t2)
#                trsn= st_n[0]
#                trsn.detrend(type='linear')
#                trsn.detrend(type='demean')
#                trsn.filter(type='bandpass',freqmin=0.1, freqmax=10)
#            
#                st_acn = calibrate1(trsn)
#                
#                
#                B=2*pi*rhoE*cE*(1/A)
#                dl=len(st_acz[0].data)
#                p = np.linspace(0,dl/100, num=dl)
#                
#                y= np.sqrt(st_acz[0].data**2 + st_ace[0].data**2 + st_acn[0].data**2)
#                y2= np.square(y)
#                y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
#                           
#                EI = y_int2[-1]
#                EE= B*(r*r)*EI
#                
#                Ens1= EE
#                AS1 = EI
#                
#                enl = np.lib.pad(enl, ((0,1),(0,0)), 'constant', constant_values=(0))
#                enl[numb][0] = r
#                al = np.lib.pad(al, ((0,1),(0,0)), 'constant', constant_values=(0))
#                al[numb][0] = r
#                numb+=1
#                enl = np.lib.pad(enl, ((0,1),(0,0)), 'constant', constant_values=(0))
#                enl[numb][0] = Ens1
#                al = np.lib.pad(al, ((0,1),(0,0)), 'constant', constant_values=(0))
#                al[numb][0] = AS1
#                numb+=1
#                
#                ls1 += Ens1
#                lsa1 += AS1
#                ls1n += 1
#                
#            if ls02[x] == 1:
#                sta='LS02'
#                cha1='EHZ'
#                cha2='EHN'
#                cha3='EHE'
#                r=rs2
#                
#                st_z = client.get_waveforms(net, sta, '', cha1, t1 , t2)
#                trsz= st_z[0]
#                trsz.detrend(type='linear')
#                trsz.detrend(type='demean')
#                trsz.filter(type='bandpass',freqmin=0.1, freqmax=10)
#            
#                st_acz = calibrate1(trsz)
#                
#                st_e = client.get_waveforms(net, sta, '', cha2, t1 , t2)
#                trse= st_e[0]
#                trse.detrend(type='linear')
#                trse.detrend(type='demean')
#                trse.filter(type='bandpass',freqmin=0.1, freqmax=10)
#            
#                st_ace = calibrate1(trse)
#                
#                st_n = client.get_waveforms(net, sta, '', cha3, t1 , t2)
#                trsn= st_n[0]
#                trsn.detrend(type='linear')
#                trsn.detrend(type='demean')
#                trsn.filter(type='bandpass',freqmin=0.1, freqmax=10)
#            
#                st_acn = calibrate1(trsn)
#                
#                
#                B=2*pi*rhoE*cE*(1/A)
#                dl=len(st_acz[0].data)
#                p = np.linspace(0,dl/100, num=dl)
#                
#                y= np.sqrt(st_acz[0].data**2 + st_ace[0].data**2 + st_acn[0].data**2)
#                y2= np.square(y)
#                y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
#                           
#                EI = y_int2[-1]
#                EE= B*(r*r)*EI
#                
#                Ens2= EE
#                AS2 = EI
#                
#                enl = np.lib.pad(enl, ((0,1),(0,0)), 'constant', constant_values=(0))
#                enl[numb][0] = r
#                al = np.lib.pad(al, ((0,1),(0,0)), 'constant', constant_values=(0))
#                al[numb][0] = r
#                numb+=1
#                enl = np.lib.pad(enl, ((0,1),(0,0)), 'constant', constant_values=(0))
#                enl[numb][0] = Ens2
#                al = np.lib.pad(al, ((0,1),(0,0)), 'constant', constant_values=(0))
#                al[numb][0] = AS2
#                numb+=1
#                
#                ls2 += Ens2
#                lsa2 += AS2
#                ls2n += 1
#                
#            if ls03[x] == 1:
#                sta='LS03'
#                cha1='EHZ'
#                cha2='EHN'
#                cha3='EHE'
#                r=rs3
#                
#                st_z = client.get_waveforms(net, sta, '', cha1, t1 , t2)
#                trsz= st_z[0]
#                trsz.detrend(type='linear')
#                trsz.detrend(type='demean')
#                trsz.filter(type='bandpass',freqmin=0.1, freqmax=10)
#            
#                st_acz = calibrate1(trsz)
#                
#                st_e = client.get_waveforms(net, sta, '', cha2, t1 , t2)
#                trse= st_e[0]
#                trse.detrend(type='linear')
#                trse.detrend(type='demean')
#                trse.filter(type='bandpass',freqmin=0.1, freqmax=10)
#            
#                st_ace = calibrate1(trse)
#                
#                st_n = client.get_waveforms(net, sta, '', cha3, t1 , t2)
#                trsn= st_n[0]
#                trsn.detrend(type='linear')
#                trsn.detrend(type='demean')
#                trsn.filter(type='bandpass',freqmin=0.1, freqmax=10)
#            
#                st_acn = calibrate1(trsn)
#                
#                
#                B=2*pi*rhoE*cE*(1/A)
#                dl=len(st_acz[0].data)
#                p = np.linspace(0,dl/100, num=dl)
#                
#                y= np.sqrt(st_acz[0].data**2 + st_ace[0].data**2 + st_acn[0].data**2)
#                y2= np.square(y)
#                y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
#                           
#                EI = y_int2[-1]
#                EE= B*(r*r)*EI
#                
#                Ens3= EE
#                AS3 = EI
#                
#                
#                enl = np.lib.pad(enl, ((0,1),(0,0)), 'constant', constant_values=(0))
#                enl[numb][0] = r
#                al = np.lib.pad(al, ((0,1),(0,0)), 'constant', constant_values=(0))
#                al[numb][0] = r
#                numb+=1
#                enl = np.lib.pad(enl, ((0,1),(0,0)), 'constant', constant_values=(0))
#                enl[numb][0] = Ens3
#                al = np.lib.pad(al, ((0,1),(0,0)), 'constant', constant_values=(0))
#                al[numb][0] = AS3
#                numb+=1
#                
#                ls3 += Ens3
#                lsa3 += AS3
#                ls3n += 1
#                
#            if ls04[x] == 1:
#                sta='LS04'
#                cha1='EHZ'
#                cha2='EHN'
#                cha3='EHE'
#                r=rs4
#                
#                st_z = client.get_waveforms(net, sta, '', cha1, t1 , t2)
#                trsz= st_z[0]
#                trsz.detrend(type='linear')
#                trsz.detrend(type='demean')
#                trsz.filter(type='bandpass',freqmin=0.1, freqmax=10)
#            
#                st_acz = calibrate1(trsz)
#                
#                st_e = client.get_waveforms(net, sta, '', cha2, t1 , t2)
#                trse= st_e[0]
#                trse.detrend(type='linear')
#                trse.detrend(type='demean')
#                trse.filter(type='bandpass',freqmin=0.1, freqmax=10)
#            
#                st_ace = calibrate1(trse)
#                
#                st_n = client.get_waveforms(net, sta, '', cha3, t1 , t2)
#                trsn= st_n[0]
#                trsn.detrend(type='linear')
#                trsn.detrend(type='demean')
#                trsn.filter(type='bandpass',freqmin=0.1, freqmax=10)
#            
#                st_acn = calibrate1(trsn)
#                
#                
#                B=2*pi*rhoE*cE*(1/A)
#                dl=len(st_acz[0].data)
#                p = np.linspace(0,dl/100, num=dl)
#                
#                y= np.sqrt(st_acz[0].data**2 + st_ace[0].data**2 + st_acn[0].data**2)
#                y2= np.square(y)
#                y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
#                           
#                EI = y_int2[-1]
#                EE= B*(r*r)*EI
#                
#                Ens4= EE
#                AS4 = EI
#                
#                enl = np.lib.pad(enl, ((0,1),(0,0)), 'constant', constant_values=(0))
#                enl[numb][0] = r
#                al = np.lib.pad(al, ((0,1),(0,0)), 'constant', constant_values=(0))
#                al[numb][0] = r
#                numb+=1
#                enl = np.lib.pad(enl, ((0,1),(0,0)), 'constant', constant_values=(0))
#                enl[numb][0] = Ens4
#                al = np.lib.pad(al, ((0,1),(0,0)), 'constant', constant_values=(0))
#                al[numb][0] = AS4
#                numb+=1
#                
#                ls4 += Ens4
#                lsa4 += AS4
#                ls4n += 1
#                
#            if ls05[x] == 1:
#                sta='LS05'
#                cha1='EHZ'
#                cha2='EHN'
#                cha3='EHE'
#                r=rs5
#                
#                st_z = client.get_waveforms(net, sta, '', cha1, t1 , t2)
#                trsz= st_z[0]
#                trsz.detrend(type='linear')
#                trsz.detrend(type='demean')
#                trsz.filter(type='bandpass',freqmin=0.1, freqmax=10)
#            
#                st_acz = calibrate1(trsz)
#                
#                st_e = client.get_waveforms(net, sta, '', cha2, t1 , t2)
#                trse= st_e[0]
#                trse.detrend(type='linear')
#                trse.detrend(type='demean')
#                trse.filter(type='bandpass',freqmin=0.1, freqmax=10)
#            
#                st_ace = calibrate1(trse)
#                
#                st_n = client.get_waveforms(net, sta, '', cha3, t1 , t2)
#                trsn= st_n[0]
#                trsn.detrend(type='linear')
#                trsn.detrend(type='demean')
#                trsn.filter(type='bandpass',freqmin=0.1, freqmax=10)
#            
#                st_acn = calibrate1(trsn)
#                
#                
#                B=2*pi*rhoE*cE*(1/A)
#                dl=len(st_acz[0].data)
#                p = np.linspace(0,dl/100, num=dl)
#                
#                y= np.sqrt(st_acz[0].data**2 + st_ace[0].data**2 + st_acn[0].data**2)
#                y2= np.square(y)
#                y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
#                           
#                EI = y_int2[-1]
#                EE= B*(r*r)*EI
#                
#                Ens5= EE
#                AS5 = EI
#                
#                enl = np.lib.pad(enl, ((0,1),(0,0)), 'constant', constant_values=(0))
#                enl[numb][0] = r
#                al = np.lib.pad(al, ((0,1),(0,0)), 'constant', constant_values=(0))
#                al[numb][0] = r
#                numb+=1
#                enl = np.lib.pad(enl, ((0,1),(0,0)), 'constant', constant_values=(0))
#                enl[numb][0] = Ens5
#                al = np.lib.pad(al, ((0,1),(0,0)), 'constant', constant_values=(0))
#                al[numb][0] = AS5
#                numb+=1
#                
#                ls5 += Ens5
#                lsa5 += AS5
#                ls5n += 1
#                
#            if ls06[x] == 1:
#                sta='LS06'
#                cha1='EHZ'
#                cha2='EHN'
#                cha3='EHE'
#                r=rs6
#                
#                st_z = client.get_waveforms(net, sta, '', cha1, t1 , t2)
#                trsz= st_z[0]
#                trsz.detrend(type='linear')
#                trsz.detrend(type='demean')
#                trsz.filter(type='bandpass',freqmin=0.1, freqmax=10)
#            
#                st_acz = calibrate1(trsz)
#                
#                st_e = client.get_waveforms(net, sta, '', cha2, t1 , t2)
#                trse= st_e[0]
#                trse.detrend(type='linear')
#                trse.detrend(type='demean')
#                trse.filter(type='bandpass',freqmin=0.1, freqmax=10)
#            
#                st_ace = calibrate1(trse)
#                
#                st_n = client.get_waveforms(net, sta, '', cha3, t1 , t2)
#                trsn= st_n[0]
#                trsn.detrend(type='linear')
#                trsn.detrend(type='demean')
#                trsn.filter(type='bandpass',freqmin=0.1, freqmax=10)
#            
#                st_acn = calibrate1(trsn)
#                
#                
#                B=2*pi*rhoE*cE*(1/A)
#                dl=len(st_acz[0].data)
#                p = np.linspace(0,dl/100, num=dl)
#                
#                y= np.sqrt(st_acz[0].data**2 + st_ace[0].data**2 + st_acn[0].data**2)
#                y2= np.square(y)
#                y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
#                           
#                EI = y_int2[-1]
#                EE= B*(r*r)*EI
#                
#                Ens6= EE
#                AS6 = EI
#                
#                enl = np.lib.pad(enl, ((0,1),(0,0)), 'constant', constant_values=(0))
#                enl[numb][0] = r
#                al = np.lib.pad(al, ((0,1),(0,0)), 'constant', constant_values=(0))
#                al[numb][0] = r
#                numb+=1
#                enl = np.lib.pad(enl, ((0,1),(0,0)), 'constant', constant_values=(0))
#                enl[numb][0] = Ens6
#                al = np.lib.pad(al, ((0,1),(0,0)), 'constant', constant_values=(0))
#                al[numb][0] = AS6
#                numb+=1
#                
#                ls6 += Ens6
#                lsa6 += AS6
#                ls6n += 1
#            
            for y in range(0,el):
                
                ad[numl][y]=al[y][0]
                ed[numl][y]=enl[y][0]
            enl = np.zeros(shape=(0,1))
            al = np.zeros(shape=(0,1))
            numb=0
            
            
       except:
            print('event lengths messed up')

#%% remove events that didnt work
            
for y in range(0,10):            
    for x in range(0,len(ad)):
        try:
            if ad[x,0] == 0:
                ad = np.delete(ad, (x), axis=0)
        except:
            fwjnwjnvfnskpjwebfgbjad = 1
    for x in range(0,len(ed)):
        try:
            if ed[x,0] == 0:
                ed = np.delete(ed, (x), axis=0) 
        except:
            fwjnwjnvfnskpjwebfgbjad = 11
        
nume = len(ad)
#%%
        

plt.figure()
plt.xlabel('Energy')
plt.ylabel('station distance (m)')
plt.title('24 events')



plt.semilogx(ed[0][1:ol:2],ed[0][0:el:2],'bx')
plt.semilogx(ed[1][1:ol:2],ed[1][0:el:2],'kx')
plt.semilogx(ed[2][1:ol:2],ed[2][0:el:2],'gx')
plt.semilogx(ed[3][1:ol:2],ed[3][0:el:2],'rx')
plt.semilogx(ed[4][1:ol:2],ed[4][0:el:2],'cx')
plt.semilogx(ed[5][1:ol:2],ed[5][0:el:2],'mx')
plt.semilogx(ed[6][1:ol:2],ed[6][0:el:2],'b^')
plt.semilogx(ed[7][1:ol:2],ed[7][0:el:2],'b+')
plt.semilogx(ed[8][1:ol:2],ed[8][0:el:2],'k+')
plt.semilogx(ed[9][1:ol:2],ed[9][0:el:2],'g+')
plt.semilogx(ed[10][1:ol:2],ed[10][0:el:2],'c+')
plt.semilogx(ed[11][1:ol:2],ed[11][0:el:2],'m+')
plt.semilogx(ed[12][1:ol:2],ed[12][0:el:2],'k^')
plt.semilogx(ed[13][1:ol:2],ed[13][0:el:2],'r+')
plt.semilogx(ed[14][1:ol:2],ed[14][0:el:2],'g^')
plt.semilogx(ed[15][1:ol:2],ed[15][0:el:2],'r^')
plt.semilogx(ed[16][1:ol:2],ed[16][0:el:2],'m^')
plt.semilogx(ed[17][1:ol:2],ed[17][0:el:2],'c^')
plt.semilogx(ed[18][1:ol:2],ed[18][0:el:2],'b*')
plt.semilogx(ed[19][1:ol:2],ed[19][0:el:2],'k*')
plt.semilogx(ed[20][1:ol:2],ed[20][0:el:2],'g*')
plt.semilogx(ed[21][1:ol:2],ed[21][0:el:2],'r*')
plt.semilogx(ed[22][1:ol:2],ed[22][0:el:2],'m*')
plt.semilogx(ed[23][1:ol:2],ed[23][0:el:2],'c*')


#%% Average Energy

lb1a = lb1/lb1n
#lb2a = lb2/lb2n
#lb3a = lb3/lb3n
lb4a = lb4/lb4n
#lb5a = lb5/lb5n
#lb6a = lb6/lb6n

#ls1a = ls1/ls1n
#ls2a = ls2/ls2n
#ls3a = ls3/ls3n
#ls4a = ls4/ls4n
#ls5a = ls5/ls5n
#ls6a = ls6/ls6n



plt.figure()
plt.title('Number of Events = {}'.format(nume))
plt.xlabel('Station Average Energy')
plt.ylabel('Distance from Source [m]')
plt.semilogx(lb1a,r1,'bx')
#plt.semilogx(lb2a,r2,'kx')
#plt.semilogx(lb3a,r3,'gx')
plt.semilogx(lb4a,r4,'rx')
#plt.semilogx(lb5a,r5,'cx')
#plt.semilogx(lb6a,r6,'mx')
#plt.semilogx(ls1a,rs1,'b+')
#plt.semilogx(ls2a,rs2,'k+')
#plt.semilogx(ls3a,rs3,'g+')
#plt.semilogx(ls4a,rs4,'r+')
#plt.semilogx(ls5a,rs5,'c+')
#plt.semilogx(ls6a,rs6,'m+')
plt.xlim([5E5, 5E8])
#plt.legend(('LB01','LB02','LB03','LB04','LB05','LS01','LS02','LS03','LS04','LS05','LS06'))



#%% Correction using average Energy compared to LB01


lb1A = lb1a/lb1a
#lb2A = lb1a/lb2a
#lb3A = lb1a/lb3a
lb4A = lb1a/lb4a
#lb5A = lb1a/lb5a
#lb6A = lb1a/lb6a
#ls1A = lb1a/ls1a
#ls2A = lb1a/ls2a
#ls3A = lb1a/ls3a
#ls4A = lb1a/ls4a
#ls5A = lb1a/ls5a
#ls6A = lb1a/ls6a


print('LB01 amplification to LB01 = ', lb1A)
#print('LB02 amplification to LB01 = ', lb2A)
#print('LB03 amplification to LB01 = ', lb3A)
print('LB04 amplification to LB01 = ', lb4A)
#print('LB05 amplification to LB01 = ', lb5A)
#print('LB06 amplification to LB01 = ', lb6A)
#print('LS01 amplification to LB01 = ', ls1A)
#print('LS02 amplification to LB01 = ', ls2A)
#print('LS03 amplification to LB01 = ', ls3A)
#print('LS04 amplification to LB01 = ', ls4A)
#print('LS05 amplification to LB01 = ', ls5A)
#print('LS06 amplification to LB01 = ', ls6A)


plt.figure()
plt.title('Number of Events = {}'.format(nume))
plt.xlabel('Corrected Station Average log10|Energy|')
plt.ylabel('Distance from Source [m]')
plt.plot(lb1a*lb1A,r1,'bx')
#plt.plot(lb2a*lb2A,r2,'kx')
#plt.plot(lb3a*lb3A,r3,'gx')
plt.plot(lb4a*lb4A,r4,'rx')
#plt.plot(lb5a*lb5A,r5,'cx')
#plt.plot(lb6a*lb6A,r6,'mx')
#plt.plot(ls1a*ls1A,rs1,'b+')
#plt.plot(ls2a*ls2A,rs2,'k+')
#plt.plot(ls3a*ls3A,rs3,'g+')
#plt.plot(ls4a*ls4A,rs4,'r+')
#plt.plot(ls5a*ls5A,rs5,'c+')
#plt.plot(ls6a*ls6A,rs6,'m+')
#plt.xlim([5.5, 7.5])
#plt.legend(('LB01','LB02','LB03','LB04','LB05','LS01','LS02','LS03','LS04','LS05','LS06'))


#%%

#%%
        

plt.figure()
plt.xlabel('Sum of Squared Amplitude')
plt.ylabel('station distance (m)')
plt.title('24 events')

plt.loglog(ad[0][1:ol:2],ad[0][0:el:2],'bx')
plt.loglog(ad[1][1:ol:2],ad[1][0:el:2],'kx')
plt.loglog(ad[2][1:ol:2],ad[2][0:el:2],'gx')
plt.loglog(ad[3][1:ol:2],ad[3][0:el:2],'rx')
plt.loglog(ad[4][1:ol:2],ad[4][0:el:2],'cx')
plt.loglog(ad[5][1:ol:2],ad[5][0:el:2],'mx')
plt.loglog(ad[6][1:ol:2],ad[6][0:el:2],'b^')
plt.loglog(ad[7][1:ol:2],ad[7][0:el:2],'b+')
plt.loglog(ad[8][1:ol:2],ad[8][0:el:2],'k+')
plt.loglog(ad[9][1:ol:2],ad[9][0:el:2],'g+')
plt.loglog(ad[10][1:ol:2],ad[10][0:el:2],'c+')
plt.loglog(ad[11][1:ol:2],ad[11][0:el:2],'m+')
plt.loglog(ad[12][1:ol:2],ad[12][0:el:2],'k^')
plt.loglog(ad[13][1:ol:2],ad[13][0:el:2],'r+')
plt.loglog(ad[14][1:ol:2],ad[14][0:el:2],'g^')
plt.loglog(ad[15][1:ol:2],ad[15][0:el:2],'r^')
plt.loglog(ad[16][1:ol:2],ad[16][0:el:2],'m^')
plt.loglog(ad[17][1:ol:2],ad[17][0:el:2],'c^')
plt.loglog(ad[18][1:ol:2],ad[18][0:el:2],'b*')
plt.loglog(ad[19][1:ol:2],ad[19][0:el:2],'k*')
plt.loglog(ad[20][1:ol:2],ad[20][0:el:2],'g*')
plt.loglog(ad[21][1:ol:2],ad[21][0:el:2],'r*')
plt.loglog(ad[22][1:ol:2],ad[22][0:el:2],'m*')
plt.loglog(ad[23][1:ol:2],ad[23][0:el:2],'c*')


#%%

lba1a = lba1/lb1n
#lba2a = lba2/lb2n
#lba3a = lba3/lb3n
lba4a = lba4/lb4n
#lba5a = lba5/lb5n
#lba6a = lba6/lb6n
#lsa1a = lsa1/ls1n
#lsa2a = lsa2/ls2n
#lsa3a = lsa3/ls3n
#lsa4a = lsa4/ls4n
#lsa5a = lsa5/ls5n
#lsa6a = lsa6/ls6n



plt.figure()
plt.title('Number of Events = {}'.format(nume))
plt.xlabel('Station Average Amplidude^2')
plt.ylabel('Distance from Source [m]')
plt.loglog(lba1a,r1,'bx')
#plt.loglog(lba2a,r2,'kx')
#plt.loglog(lba3a,r3,'gx')
plt.loglog(lba4a,r4,'rx')
#plt.loglog(lba5a,r5,'cx')
#plt.loglog(lba6a,r6,'mx')
#plt.loglog(lsa1a,rs1,'b+')
#plt.loglog(lsa2a,rs2,'k+')
#plt.loglog(lsa3a,rs3,'g+')
#plt.loglog(lsa4a,rs4,'r+')
#plt.loglog(lsa5a,rs5,'c+')
#plt.loglog(lsa6a,rs6,'m+')
#plt.xlim([5E5, 5E8])
#plt.legend(('LB01','LB02','LB03','LB04','LB05','LS01','LS02','LS03','LS04','LS05','LS06'))





plt.figure()
plt.title('Number of Events = {}'.format(nume))
plt.xlabel('Corrected Station Average Amplitudes')
plt.ylabel('Distance from Source [m]')
plt.loglog(lba1a*lb1A,r1,'bx')
#plt.loglog(lba2a*lb2A,r2,'kx')
#plt.loglog(lba3a*lb3A,r3,'gx')
plt.loglog(lba4a*lb4A,r4,'rx')
#plt.loglog(lba5a*lb5A,r5,'cx')
#plt.loglog(lba6a*lb6A,r6,'mx')
#plt.loglog(lsa1a*ls1A,rs1,'b+')
#plt.loglog(lsa2a*ls2A,rs2,'k+')
#plt.loglog(lsa3a*ls3A,rs3,'g+')
#plt.loglog(lsa4a*ls4A,rs4,'r+')
#plt.loglog(lsa5a*ls5A,rs5,'c+')
#plt.loglog(lsa6a*ls6A,rs6,'m+')
#plt.xlim([5.5, 7.5])
#plt.legend(('LB01','LB02','LB03','LB04','LB05','LS01','LS02','LS03','LS04','LS05','LS06'))

#%%

dist=np.asarray([r1,r4])
amps=np.asarray([lba1a*lb1A,lba4a*lb4A])#,lba3a*lb3A,lba4a*lb4A,lba5a*lb5A]) #

lx1=np.zeros(shape=(2,1))
ly1=np.zeros(shape=(2,1))

for p in range(0,2):
    lx=np.log10(amps[p])
    ly=np.log10(dist[p]**2)
    
    lx1[p][0]=lx
    ly1[p][0]=ly
    
#log(y) = log(c) + m*log(x)   


lx1=np.concatenate(lx1)
ly1=np.concatenate(ly1)
#
c, m = polyfit(lx1, ly1, 1)

plt.figure()
plt.plot(lx1,ly1,'kx')
plt.plot(lx1, c + m * lx1, 'r-')
plt.xlabel("Log10|Amp^2|")
plt.ylabel("Log10|Dist^2|")
print("y= d*x^m;  m= ",m , "d= ", 10**c)
print("log(y)= c + m*log(x);  m= ",m , "c= ", c)


#%%


#for x in range(0,400):
#    active = lb01[x] + lb02[x] + lb03[x] + lb04[x] + lb05[x] + lb06[x] + ls01[x] + ls02[x] + ls03[x] + ls04[x] + ls05[x] + ls06[x]
#    if  active == 10:        
#        t1=UTCDateTime(ev_t[x] -10)
#        t2=t1+80
#        
#        if lb01[x] == 1:
#            sta='LB01'
#            cha1='HHZ'
#            cha2='HHN'
#            cha3='HHE'
#            r=r1
#            
#            st_z = client.get_waveforms(net, sta, '', cha1, t1 , t2)
#            trsz= st_z[0]
#            trsz.detrend(type='linear')
#            trsz.detrend(type='demean')
#            trsz.filter(type='bandpass',freqmin=0.1, freqmax=10)
#        
#            st_acz = calibrate1(trsz)
#                
#            B=2*pi*rhoE*cE*(1/A)
#            dl=len(st_acz[0].data)
#            p = np.linspace(0,dl/100, num=dl)
#            
#            y= st_acz[0].data
#            y2= np.square(st_acz[0].data)
#            
#            
#            y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
#            EI = y_int2[-1]
#            EE1= B*(r*r)*EI
#            
#            
#            st_e = client.get_waveforms(net, sta, '', cha2, t1 , t2)
#            trse= st_e[0]
#            trse.detrend(type='linear')
#            trse.detrend(type='demean')
#            trse.filter(type='bandpass',freqmin=0.1, freqmax=10)
#        
#            st_ace = calibrate1(trse)
#                
#            B=2*pi*rhoE*cE*(1/A)
#            dl=len(st_ace[0].data)
#            p = np.linspace(0,dl/100, num=dl)
#            
#            y= st_ace[0].data
#            y2= np.square(st_ace[0].data)
#            
#            
#            y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
#            EI = y_int2[-1]
#            EE2= B*(r*r)*EI        
#            
#            st_n = client.get_waveforms(net, sta, '', cha3, t1 , t2)
#            trsn= st_n[0]
#            trsn.detrend(type='linear')
#            trsn.detrend(type='demean')
#            trsn.filter(type='bandpass',freqmin=0.1, freqmax=10)
#        
#            st_acn = calibrate1(trsn)
#                
#            B=2*pi*rhoE*cE*(1/A)
#            dl=len(st_acn[0].data)
#            p = np.linspace(0,dl/100, num=dl)
#            
#            y= st_acn[0].data
#            y2= np.square(st_acn[0].data)
#            
#            
#            y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
#            EI = y_int2[-1]
#            EE3= B*(r*r)*EI
#            
#            Energy1 = (EE1 + EE2 + EE3)
#            Energy1c = (EE1 + EE2 + EE3) * lb1A
#            
#        if lb02[x] == 1:
#            sta='LB02'
#            cha1='HHZ'
#            cha2='HHN'
#            cha3='HHE'
#            r=r2
#            
#            st_z = client.get_waveforms(net, sta, '', cha1, t1 , t2)
#            trsz= st_z[0]
#            trsz.detrend(type='linear')
#            trsz.detrend(type='demean')
#            trsz.filter(type='bandpass',freqmin=0.1, freqmax=10)
#        
#            st_acz = calibrate1(trsz)
#                
#            B=2*pi*rhoE*cE*(1/A)
#            dl=len(st_acz[0].data)
#            p = np.linspace(0,dl/100, num=dl)
#            
#            y= st_acz[0].data
#            y2= np.square(st_acz[0].data)
#            
#            
#            y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
#            EI = y_int2[-1]
#            EE1= B*(r*r)*EI
#            
#            
#            st_e = client.get_waveforms(net, sta, '', cha2, t1 , t2)
#            trse= st_e[0]
#            trse.detrend(type='linear')
#            trse.detrend(type='demean')
#            trse.filter(type='bandpass',freqmin=0.1, freqmax=10)
#        
#            st_ace = calibrate1(trse)
#                
#            B=2*pi*rhoE*cE*(1/A)
#            dl=len(st_ace[0].data)
#            p = np.linspace(0,dl/100, num=dl)
#            
#            y= st_ace[0].data
#            y2= np.square(st_ace[0].data)
#            
#            
#            y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
#            EI = y_int2[-1]
#            EE2= B*(r*r)*EI        
#            
#            st_n = client.get_waveforms(net, sta, '', cha3, t1 , t2)
#            trsn= st_n[0]
#            trsn.detrend(type='linear')
#            trsn.detrend(type='demean')
#            trsn.filter(type='bandpass',freqmin=0.1, freqmax=10)
#        
#            st_acn = calibrate1(trsn)
#                
#            B=2*pi*rhoE*cE*(1/A)
#            dl=len(st_acn[0].data)
#            p = np.linspace(0,dl/100, num=dl)
#            
#            y= st_acn[0].data
#            y2= np.square(st_acn[0].data)
#            
#            
#            y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
#            EI = y_int2[-1]
#            EE3= B*(r*r)*EI
#            
#            Energy2 = (EE1 + EE2 + EE3) 
#            Energy2c = (EE1 + EE2 + EE3) *lb2A
#            
#        if lb03[x] == 1:
#            sta='LB03'
#            cha1='HHZ'
#            cha2='HHN'
#            cha3='HHE'
#            r=r3
#            
#            st_z = client.get_waveforms(net, sta, '', cha1, t1 , t2)
#            trsz= st_z[0]
#            trsz.detrend(type='linear')
#            trsz.detrend(type='demean')
#            trsz.filter(type='bandpass',freqmin=0.1, freqmax=10)
#        
#            st_acz = calibrate1(trsz)
#                
#            B=2*pi*rhoE*cE*(1/A)
#            dl=len(st_acz[0].data)
#            p = np.linspace(0,dl/100, num=dl)
#            
#            y= st_acz[0].data
#            y2= np.square(st_acz[0].data)
#            
#            
#            y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
#            EI = y_int2[-1]
#            EE1= B*(r*r)*EI
#            
#            
#            st_e = client.get_waveforms(net, sta, '', cha2, t1 , t2)
#            trse= st_e[0]
#            trse.detrend(type='linear')
#            trse.detrend(type='demean')
#            trse.filter(type='bandpass',freqmin=0.1, freqmax=10)
#        
#            st_ace = calibrate1(trse)
#                
#            B=2*pi*rhoE*cE*(1/A)
#            dl=len(st_ace[0].data)
#            p = np.linspace(0,dl/100, num=dl)
#            
#            y= st_ace[0].data
#            y2= np.square(st_ace[0].data)
#            
#            
#            y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
#            EI = y_int2[-1]
#            EE2= B*(r*r)*EI        
#            
#            st_n = client.get_waveforms(net, sta, '', cha3, t1 , t2)
#            trsn= st_n[0]
#            trsn.detrend(type='linear')
#            trsn.detrend(type='demean')
#            trsn.filter(type='bandpass',freqmin=0.1, freqmax=10)
#        
#            st_acn = calibrate1(trsn)
#                
#            B=2*pi*rhoE*cE*(1/A)
#            dl=len(st_acn[0].data)
#            p = np.linspace(0,dl/100, num=dl)
#            
#            y= st_acn[0].data
#            y2= np.square(st_acn[0].data)
#            
#            
#            y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
#            EI = y_int2[-1]
#            EE3= B*(r*r)*EI
#            
#            Energy3 = (EE1 + EE2 + EE3) 
#            Energy3c = (EE1 + EE2 + EE3) *lb3A
#      
#        
#        
#print(Energy1 , Energy2 , Energy3)
#print(Energy1c , Energy2c , Energy3c)
#
#











