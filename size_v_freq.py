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

svf=np.zeros(shape=(0,6))

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
            
            
        for x in range (0,int(len(st))):
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
                
                # energy two Hz
                two = int((2*window)+1)
                num2e=0
                for t in range(0,two,1):
                    n2e=famp2[t]                  
                    num2e += n2e
                pe2hz=num2e/sum(famp2)
                
                # Energy of explosion
                
                if tr.stats.station =="LB01":
                    r=r1
                if tr.stats.station =="LB02":
                    r=r2
                if tr.stats.station =="LB03":
                    r=r3
                if tr.stats.station =="LB04":
                    r=r4
                if tr.stats.station =="LB05":
                    r=r5
                if tr.stats.station =="LB06":
                    r=r6
                if tr.stats.station =="LS01":
                    r=rs1
                if tr.stats.station =="LS02":
                    r=rs2
                if tr.stats.station =="LS03":
                    r=rs3
                if tr.stats.station =="LS04":
                    r=rs4
                if tr.stats.station =="LS05":
                    r=rs5
                if tr.stats.station =="LS06":
                    r=rs6
                
                st_c = calibrate1(tr)
                B=2*pi*rhoE*cE*(1/A)
                EI = sum(np.square(st_c[0].data))
                Energy = B*(r*r)*EI

                svf = np.lib.pad(svf, ((0,1),(0,0)), 'constant', constant_values=(0))
  
                svf[event_num][0]= event_num + 1
                svf[event_num][1]=pe2hz
                svf[event_num][2]=dom
                svf[event_num][3]=cent
                svf[event_num][4]=bw
                svf[event_num][5]=Energy
                
                event_num+=1
                    
np.savetxt("/Users/william/Documents/scanner/analysis/size_v_freq_LB01_v2.csv", svf,delimiter=",",header="Event,EB2,dom_f,cent_f,Bwid50,Energy")                   
 
                   
plt.figure(10)
plt.semilogy(svf[:,1],svf[:,5],'rx')
plt.xlim((-0.1,1.1))
plt.ylim([1e1,1e16])
plt.xlabel("Energy below 2Hz [%]")
plt.ylabel("Explosion Energy [J]")
plt.title("Size vs Energy division")                    
                    
                    
plt.figure(11)
plt.semilogy(svf[:,2],svf[:,5],'rx')
plt.xlim((0,5))
plt.ylim([1e1,1e16])
plt.xlabel("Dominant Frequency [Hz]")
plt.ylabel("Explosion Energy [J]")
plt.title("Size vs Dominant Frequency")                      
                    
                    
plt.figure(12)
plt.semilogy(svf[:,3],svf[:,5],'rx')
plt.xlim((0,7))
plt.ylim([1e1,1e16])
plt.xlabel("Central Frequency [Hz]")
plt.ylabel("Explosion Energy [J]")
plt.title("Size vs Central Frequency")                     
                    
                    
plt.figure(13)
plt.semilogy(svf[:,4],svf[:,5],'rx')
plt.xlim((0,7))
plt.ylim([1e1,1e16])
plt.xlabel("50% Bandwidth [Hz]")
plt.ylabel("Explosion Energy [J]")
plt.title("Size vs Bandwidth")                      
                    
                    
                    
reclen = 512
chunksize = 100000 * reclen # Around 50 MB

svf=np.zeros(shape=(0,6))

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
            
            
        for x in range (0,int(len(st))):
            if st[x].stats.station == "LB02":# or "LB02" or "LB03" or "LB04" or "LB05" or "LB06" :
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
                
                # energy two Hz
                two = int((2*window)+1)
                num2e=0
                for t in range(0,two,1):
                    n2e=famp2[t]                  
                    num2e += n2e
                pe2hz=num2e/sum(famp2)
                
                # Energy of explosion
                
                if tr.stats.station =="LB01":
                    r=r1
                if tr.stats.station =="LB02":
                    r=r2
                if tr.stats.station =="LB03":
                    r=r3
                if tr.stats.station =="LB04":
                    r=r4
                if tr.stats.station =="LB05":
                    r=r5
                if tr.stats.station =="LB06":
                    r=r6
                if tr.stats.station =="LS01":
                    r=rs1
                if tr.stats.station =="LS02":
                    r=rs2
                if tr.stats.station =="LS03":
                    r=rs3
                if tr.stats.station =="LS04":
                    r=rs4
                if tr.stats.station =="LS05":
                    r=rs5
                if tr.stats.station =="LS06":
                    r=rs6
                
                st_c = calibrate1(tr)
                B=2*pi*rhoE*cE*(1/A)
                EI = sum(np.square(st_c[0].data))
                Energy = B*(r*r)*EI

                svf = np.lib.pad(svf, ((0,1),(0,0)), 'constant', constant_values=(0))
  
                svf[event_num][0]= event_num + 1
                svf[event_num][1]=pe2hz
                svf[event_num][2]=dom
                svf[event_num][3]=cent
                svf[event_num][4]=bw
                svf[event_num][5]=Energy
                
                event_num+=1
                    
np.savetxt("/Users/william/Documents/scanner/analysis/size_v_freq_LB02_v2.csv", svf,delimiter=",",header="Event,EB2,dom_f,cent_f,Bwid50,Energy")                   
 
                   
plt.figure(10)
plt.semilogy(svf[:,1],svf[:,5],'gx')
plt.figure(11)
plt.semilogy(svf[:,2],svf[:,5],'gx')
plt.figure(12)
plt.semilogy(svf[:,3],svf[:,5],'gx')
plt.figure(13)
plt.semilogy(svf[:,4],svf[:,5],'gx')                    
                    
reclen = 512
chunksize = 100000 * reclen # Around 50 MB

svf=np.zeros(shape=(0,6))

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
            
            
        for x in range (0,int(len(st))):
            if st[x].stats.station == "LB03":# or "LB02" or "LB03" or "LB04" or "LB05" or "LB06" :
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
                
                # energy two Hz
                two = int((2*window)+1)
                num2e=0
                for t in range(0,two,1):
                    n2e=famp2[t]                  
                    num2e += n2e
                pe2hz=num2e/sum(famp2)
                
                # Energy of explosion
                
                if tr.stats.station =="LB01":
                    r=r1
                if tr.stats.station =="LB02":
                    r=r2
                if tr.stats.station =="LB03":
                    r=r3
                if tr.stats.station =="LB04":
                    r=r4
                if tr.stats.station =="LB05":
                    r=r5
                if tr.stats.station =="LB06":
                    r=r6
                if tr.stats.station =="LS01":
                    r=rs1
                if tr.stats.station =="LS02":
                    r=rs2
                if tr.stats.station =="LS03":
                    r=rs3
                if tr.stats.station =="LS04":
                    r=rs4
                if tr.stats.station =="LS05":
                    r=rs5
                if tr.stats.station =="LS06":
                    r=rs6
                
                st_c = calibrate1(tr)
                B=2*pi*rhoE*cE*(1/A)
                EI = sum(np.square(st_c[0].data))
                Energy = B*(r*r)*EI

                svf = np.lib.pad(svf, ((0,1),(0,0)), 'constant', constant_values=(0))
  
                svf[event_num][0]= event_num + 1
                svf[event_num][1]=pe2hz
                svf[event_num][2]=dom
                svf[event_num][3]=cent
                svf[event_num][4]=bw
                svf[event_num][5]=Energy
                
                event_num+=1
                    
np.savetxt("/Users/william/Documents/scanner/analysis/size_v_freq_LB03_v2.csv", svf,delimiter=",",header="Event,EB2,dom_f,cent_f,Bwid50,Energy")                   
 
                   
plt.figure(10)
plt.semilogy(svf[:,1],svf[:,5],'bx')
plt.figure(11)
plt.semilogy(svf[:,2],svf[:,5],'bx')
plt.figure(12)
plt.semilogy(svf[:,3],svf[:,5],'bx')
plt.figure(13)
plt.semilogy(svf[:,4],svf[:,5],'bx')                    
                    
reclen = 512
chunksize = 100000 * reclen # Around 50 MB

svf=np.zeros(shape=(0,6))

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
            
            
        for x in range (0,int(len(st))):
            if st[x].stats.station == "LB04":# or "LB02" or "LB03" or "LB04" or "LB05" or "LB06" :
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
                
                # energy two Hz
                two = int((2*window)+1)
                num2e=0
                for t in range(0,two,1):
                    n2e=famp2[t]                  
                    num2e += n2e
                pe2hz=num2e/sum(famp2)
                
                # Energy of explosion
                
                if tr.stats.station =="LB01":
                    r=r1
                if tr.stats.station =="LB02":
                    r=r2
                if tr.stats.station =="LB03":
                    r=r3
                if tr.stats.station =="LB04":
                    r=r4
                if tr.stats.station =="LB05":
                    r=r5
                if tr.stats.station =="LB06":
                    r=r6
                if tr.stats.station =="LS01":
                    r=rs1
                if tr.stats.station =="LS02":
                    r=rs2
                if tr.stats.station =="LS03":
                    r=rs3
                if tr.stats.station =="LS04":
                    r=rs4
                if tr.stats.station =="LS05":
                    r=rs5
                if tr.stats.station =="LS06":
                    r=rs6
                
                st_c = calibrate1(tr)
                B=2*pi*rhoE*cE*(1/A)
                EI = sum(np.square(st_c[0].data))
                Energy = B*(r*r)*EI

                svf = np.lib.pad(svf, ((0,1),(0,0)), 'constant', constant_values=(0))
  
                svf[event_num][0]= event_num + 1
                svf[event_num][1]=pe2hz
                svf[event_num][2]=dom
                svf[event_num][3]=cent
                svf[event_num][4]=bw
                svf[event_num][5]=Energy
                
                event_num+=1
                    
np.savetxt("/Users/william/Documents/scanner/analysis/size_v_freq_LB04_v2.csv", svf,delimiter=",",header="Event,EB2,dom_f,cent_f,Bwid50,Energy")                   
 
                   
plt.figure(10)
plt.semilogy(svf[:,1],svf[:,5],'cx')
plt.figure(11)
plt.semilogy(svf[:,2],svf[:,5],'cx')
plt.figure(12)
plt.semilogy(svf[:,3],svf[:,5],'cx')
plt.figure(13)
plt.semilogy(svf[:,4],svf[:,5],'cx')                     
                    
reclen = 512
chunksize = 100000 * reclen # Around 50 MB

svf=np.zeros(shape=(0,6))

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
            
            
        for x in range (0,int(len(st))):
            if st[x].stats.station == "LB05":# or "LB02" or "LB03" or "LB04" or "LB05" or "LB06" :
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
                
                # energy two Hz
                two = int((2*window)+1)
                num2e=0
                for t in range(0,two,1):
                    n2e=famp2[t]                  
                    num2e += n2e
                pe2hz=num2e/sum(famp2)
                
                # Energy of explosion
                
                if tr.stats.station =="LB01":
                    r=r1
                if tr.stats.station =="LB02":
                    r=r2
                if tr.stats.station =="LB03":
                    r=r3
                if tr.stats.station =="LB04":
                    r=r4
                if tr.stats.station =="LB05":
                    r=r5
                if tr.stats.station =="LB06":
                    r=r6
                if tr.stats.station =="LS01":
                    r=rs1
                if tr.stats.station =="LS02":
                    r=rs2
                if tr.stats.station =="LS03":
                    r=rs3
                if tr.stats.station =="LS04":
                    r=rs4
                if tr.stats.station =="LS05":
                    r=rs5
                if tr.stats.station =="LS06":
                    r=rs6
                
                st_c = calibrate1(tr)
                B=2*pi*rhoE*cE*(1/A)
                EI = sum(np.square(st_c[0].data))
                Energy = B*(r*r)*EI

                svf = np.lib.pad(svf, ((0,1),(0,0)), 'constant', constant_values=(0))
  
                svf[event_num][0]= event_num + 1
                svf[event_num][1]=pe2hz
                svf[event_num][2]=dom
                svf[event_num][3]=cent
                svf[event_num][4]=bw
                svf[event_num][5]=Energy
                
                event_num+=1
                    
np.savetxt("/Users/william/Documents/scanner/analysis/size_v_freq_LB05_v2.csv", svf,delimiter=",",header="Event,EB2,dom_f,cent_f,Bwid50,Energy")                   
 
                   
plt.figure(10)
plt.semilogy(svf[:,1],svf[:,5],'mx')
plt.figure(11)
plt.semilogy(svf[:,2],svf[:,5],'mx')
plt.figure(12)
plt.semilogy(svf[:,3],svf[:,5],'mx')
plt.figure(13)
plt.semilogy(svf[:,4],svf[:,5],'mx')                  
                    
reclen = 512
chunksize = 100000 * reclen # Around 50 MB

svf=np.zeros(shape=(0,6))

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
            
            
        for x in range (0,int(len(st))):
            if st[x].stats.station == "LB06":# or "LB02" or "LB03" or "LB04" or "LB05" or "LB06" :
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
                
                # energy two Hz
                two = int((2*window)+1)
                num2e=0
                for t in range(0,two,1):
                    n2e=famp2[t]                  
                    num2e += n2e
                pe2hz=num2e/sum(famp2)
                
                # Energy of explosion
                
                if tr.stats.station =="LB01":
                    r=r1
                if tr.stats.station =="LB02":
                    r=r2
                if tr.stats.station =="LB03":
                    r=r3
                if tr.stats.station =="LB04":
                    r=r4
                if tr.stats.station =="LB05":
                    r=r5
                if tr.stats.station =="LB06":
                    r=r6
                if tr.stats.station =="LS01":
                    r=rs1
                if tr.stats.station =="LS02":
                    r=rs2
                if tr.stats.station =="LS03":
                    r=rs3
                if tr.stats.station =="LS04":
                    r=rs4
                if tr.stats.station =="LS05":
                    r=rs5
                if tr.stats.station =="LS06":
                    r=rs6
                
                st_c = calibrate1(tr)
                B=2*pi*rhoE*cE*(1/A)
                EI = sum(np.square(st_c[0].data))
                Energy = B*(r*r)*EI

                svf = np.lib.pad(svf, ((0,1),(0,0)), 'constant', constant_values=(0))
  
                svf[event_num][0]= event_num + 1
                svf[event_num][1]=pe2hz
                svf[event_num][2]=dom
                svf[event_num][3]=cent
                svf[event_num][4]=bw
                svf[event_num][5]=Energy
                
                event_num+=1
                    
np.savetxt("/Users/william/Documents/scanner/analysis/size_v_freq_LB06_v2.csv", svf,delimiter=",",header="Event,EB2,dom_f,cent_f,Bwid50,Energy")                   
 
                   
plt.figure(10)
plt.semilogy(svf[:,1],svf[:,5],'kx')
plt.legend(('LB01','LB02','LB03','LB04','LB05','LB06'))  
plt.figure(11)
plt.semilogy(svf[:,2],svf[:,5],'kx')
plt.legend(('LB01','LB02','LB03','LB04','LB05','LB06'))  
plt.figure(12)
plt.semilogy(svf[:,3],svf[:,5],'kx')
plt.legend(('LB01','LB02','LB03','LB04','LB05','LB06'))   
plt.figure(13)
plt.semilogy(svf[:,4],svf[:,5],'kx')
plt.legend(('LB01','LB02','LB03','LB04','LB05','LB06'))  
                    
                                        
                    
                    
#  plt.legend(('LB01','LB02','LB03','LB04','LB05','LB06'))                   
#plt.legend(('LS01','LS02','LS03','LS04','LS05','LS06'))  
                    
                    
                    
                    
                    
                    
