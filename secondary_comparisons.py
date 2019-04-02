#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 12:35:18 2018

@author: william
"""

# looking for infrasound signals from the secondary events


import obspy
import io
from obspy import read
from obspy.clients.earthworm import Client
from obspy import UTCDateTime
from obspy import Stream
import numpy as np
import matplotlib.pyplot as plt
from obspy.signal.trigger import classic_sta_lta, recursive_sta_lta
from obspy.signal.trigger import plot_trigger, trigger_onset

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
      

saved = np.zeros(shape=(0,4))
num = 0
first =0

sta = 'LB01' # STATION 
#cha = 'HHZ' # CHANNEL
net = 'Z4'  # 
loc = ''    # location, it depends mostly of which network you are in. 
client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23

#%% read in data
reclen = 512
chunksize = 100000 * reclen # 100000 = Around 50 MB 
with io.open("/Users/william/Documents/scanner/output_data/m32.mseed", "rb") as fh:
        # just month 2
#with io.open("/Users/william/Documents/scanner/output_data/EXP_all_data_stream_2_month_2.mseed", "rb") as fh:
    while True:
        with io.BytesIO() as buf:
            c = fh.read(chunksize)
            if not c:
                break
            buf.write(c)
            buf.seek(0, 0)
            st = obspy.read(buf)
            
#%% for first 'load' - save time and energy of each event - save 'lowest ranking' station           
            if first == 0:
                first =1
                if st[0].stats.station == 'LB01':
                    station=1
                    r=r1
                if st[0].stats.station == 'LB02':
                    station=2
                    r=r2
                if st[0].stats.station == 'LB03':
                    station=3
                    r=r3
                if st[0].stats.station == 'LB04':
                    station=4
                    r=r4
                if st[0].stats.station == 'LB05':
                    station=5
                    r=r5
                if st[0].stats.station == 'LB06':
                    station=6
                    r=r6
                if st[0].stats.station == 'LS01':
                    station=7
                    r=rs1
                if st[0].stats.station == 'LS02':
                    station=8
                    r=rs2
                if st[0].stats.station == 'LS03':
                    station=9
                    r=rs3
                if st[0].stats.station == 'LS04':
                    station=10
                    r=rs4
                if st[0].stats.station == 'LS05':
                    station=11
                    r=rs5
                if st[0].stats.station == 'LS06':
                    station=12
                    r=rs6
                
                tr = st[0]
                st_c = calibrate1(tr)
                B=2*pi*rhoE*cE*(1/A)
                EI = sum(np.square(st_c[0].data))
                EE= B*(r*r)*EI
                
                saved = np.lib.pad(saved, ((0,1),(0,0)), 'constant', constant_values=(0))
                saved[num][0]=st[0].stats.starttime.timestamp
                saved[num][1]=max(abs(st[0].data))
                saved[num][2]=station
                saved[num][3]=EE
                num += 1
                
                for x in range (1,int(len(st))):
                    rt = st[x].stats.starttime.timestamp
                    near,ind = find_nearest(saved[:,0], rt)
                    
                    if st[x].stats.station == 'LB01':
                        station=1
                        r=r1
                    if st[x].stats.station == 'LB02':
                        station=2
                        r=r2
                    if st[x].stats.station == 'LB03':
                        station=3
                        r=r3
                    if st[x].stats.station == 'LB04':
                        station=4
                        r=r4
                    if st[x].stats.station == 'LB05':
                        station=5
                        r=r5
                    if st[x].stats.station == 'LB06':
                        station=6
                        r=r6
                    if st[x].stats.station == 'LS01':
                        station=7
                        r=rs1
                    if st[x].stats.station == 'LS02':
                        station=8
                        r=rs2
                    if st[x].stats.station == 'LS03':
                        station=9
                        r=rs3
                    if st[x].stats.station == 'LS04':
                        station=10
                        r=rs4
                    if st[x].stats.station == 'LS05':
                        station=11
                        r=rs5
                    if st[x].stats.station == 'LS06':
                        station=12
                        r=rs6
                        
                    tr = st[x]
                    st_c = calibrate1(tr)
                    B=2*pi*rhoE*cE*(1/A)
                    EI = sum(np.square(st_c[0].data))
                    EE= B*(r*r)*EI
                    
                    if abs(rt-near) > 10:
                        saved = np.lib.pad(saved, ((0,1),(0,0)), 'constant', constant_values=(0))
                        saved[num][0]=st[x].stats.starttime.timestamp
                        saved[num][1]=max((st[x].data))
                        saved[num][2]=station
                        saved[num][3]=EE
                        num += 1
                    else:
                        if station < saved[ind,2]:
                            saved[ind,2] = station
                            saved[ind,1] = max((st[x].data))
                            saved[ind,0] = st[x].stats.starttime.timestamp
                            saved[ind,3]=EE
                            
#%% for subsequent 'loads' - save time and energy of each event - save 'lowest ranking' station              
            else:
                for x in range (0,int(len(st))):
                    rt = st[x].stats.starttime.timestamp
                    near,ind = find_nearest(saved[:,0], rt)
                    
                    if st[x].stats.station == 'LB01':
                        station=1
                        r=r1
                    if st[x].stats.station == 'LB02':
                        station=2
                        r=r2
                    if st[x].stats.station == 'LB03':
                        station=3
                        r=r3
                    if st[x].stats.station == 'LB04':
                        station=4
                        r=r4
                    if st[x].stats.station == 'LB05':
                        station=5
                        r=r5
                    if st[x].stats.station == 'LB06':
                        station=6
                        r=r6
                    if st[x].stats.station == 'LS01':
                        station=7
                        r=rs1
                    if st[x].stats.station == 'LS02':
                        station=8
                        r=rs2
                    if st[x].stats.station == 'LS03':
                        station=9
                        r=rs3
                    if st[x].stats.station == 'LS04':
                        station=10
                        r=rs4
                    if st[x].stats.station == 'LS05':
                        station=11
                        r=rs5
                    if st[x].stats.station == 'LS06':
                        station=12
                        r=rs6
                        
                    tr = st[x]
                    st_c = calibrate1(tr)
                    B=2*pi*rhoE*cE*(1/A)
                    EI = sum(np.square(st_c[0].data))
                    EE= B*(r*r)*EI
                    
                    if abs(rt-near) > 10:
                        saved = np.lib.pad(saved, ((0,1),(0,0)), 'constant', constant_values=(0))
                        saved[num][0]=st[x].stats.starttime.timestamp
                        saved[num][1]=max((st[x].data))
                        saved[num][2]=station
                        saved[num][3]=EE
                        num += 1
                    else:
                        if station < saved[ind,2]:
                            saved[ind,2] = station
                            saved[ind,1] = max((st[x].data))
                            saved[ind,0] = st[x].stats.starttime.timestamp
                            saved[ind,3]=EE
                    
col=0                  
saved=saved[np.argsort(saved[:,col])]          # order in time         
#print(num) 


#%% count secondaries/primaries

last_p = saved[0,0]
peak1 = saved[0,1]
last_station = saved[0,2]
last_E= saved[0,3]
p_tot=1
s_tot=0
sl_tot=0
un_tot =0
count =0

pands = np.zeros(shape=(0,4))
all4 = np.zeros(shape=(0,4))
num=0
for x in range(1,len(saved)):
    
    peak2=saved[x,1]
    new_station = saved[x,2]
    new_E= saved[x,3]
    
    if last_station == new_station:
        if saved[x,0] - last_p < 10*60 :
#            if peak1>peak2:
            if last_E>5*new_E:
                Eper=(new_E/last_E)*100
                s_tot += 1
                pands = np.lib.pad(pands, ((0,1),(0,0)), 'constant', constant_values=(0))
                pands[count][0]=saved[x,0]
                pands[count][1]=last_p
                pands[count][2]=saved[x,2]
                pands[count][3]=Eper
                count +=1
                s_arrival = saved[x,0]
                try:
                    tr_inf,all_4 = s_inf(last_p,s_arrival,new_station)
                    if all_4[0][0] != 0:
                        all4 = np.lib.pad(all4, ((0,1),(0,0)), 'constant', constant_values=(0))
                        all4[num][0]=all_4[0][0]
                        all4[num][1]=all_4[0][1]
                        all4[num][2]=all_4[0][2]
                        all4[num][3]=all_4[0][3]
                        num+=1
                except:
                    z=1
            else:
                sl_tot += 1
               
            
        else:
            p_tot +=1
            last_p = saved[x,0]
            peak1=saved[x,1]
            last_station=saved[x,2]
            last_E= saved[x,3]
    else:
        if saved[x,0] - last_p < 10*60 :
            un_tot += 1
            if last_E>5*new_E:
                Eper=(new_E/last_E)*100
                s_tot += 1
                pands = np.lib.pad(pands, ((0,1),(0,0)), 'constant', constant_values=(0))
                pands[count][0]=saved[x,0]
                pands[count][1]=last_p
                pands[count][2]=saved[x,2]
                pands[count][3]=Eper
                count +=1
                s_arrival = saved[x,0]
                try:
                    tr_inf,all_4 = s_inf(last_p,s_arrival,new_station)
                    if all_4[0][0] != 0:
                        all4 = np.lib.pad(all4, ((0,1),(0,0)), 'constant', constant_values=(0))
                        all4[num][0]=all_4[0][0]
                        all4[num][1]=all_4[0][1]
                        all4[num][2]=all_4[0][2]
                        all4[num][3]=all_4[0][3]
                        num+=1
                except:
                    z=1
            else:
                sl_tot += 1
                
            
        else:
             p_tot +=1
             last_p = saved[x,0]
             peak1=saved[x,1]
             last_station=saved[x,2]
             last_E= saved[x,3]




##%%
#for x in range(0,len(all4)):
#    for y in range(0,4):
#        if y == 0 or y == 2:
#            cha= 'HDF'
#            fm=2
#        else:
#            cha= 'HHZ'
#            fm=6
#        t1=UTCDateTime(all4[x,y]) 
#        t2=t1+120
#        st = Stream()
#        st = client.get_waveforms(net, sta, '', cha, t1-60 , t2) 
#        st.detrend(type='linear')
#        st.detrend(type='demean')
#        st.filter("bandpass", freqmin=0.2,freqmax=fm)
#        if y < 2:
#            #primary = red
#            
#            st.plot(color='r',starttime=t1-20, endtime=t2)  
#        else:
#            #secondry = blue
#            st.plot(color='b',starttime=t1-20, endtime=t2)  

     #%%
shift = 1000
sum_cora=0
sum_cors=0
count=0
for x in range(0,len(all4)):
    try:

        # Secondary Seismic
        t1s=UTCDateTime(all4[x,3]) 
        t2s=t1s+120
        st = Stream()
        st = client.get_waveforms(net, sta, '', 'HHZ', t1s-60 , t2s) 
        st_s2 = st[0]
        st_s2.detrend(type='linear')
        st_s2.detrend(type='demean')
        st_s2.filter("bandpass", freqmin=0.2,freqmax=6)
        
        
        sr = st_s2.stats.sampling_rate
        nsta=int(1*sr)                                      #2
        nlta=int(20*sr) 
        stream=st_s2.data
        cft=recursive_sta_lta(stream, nsta, nlta)
        trig_on=8                                            #8
        trig_off=1                                        #0.2
#        plot_trigger(st_s2, cft, trig_on, trig_off) 
        
        on_off_s2 = trigger_onset(cft,trig_on,trig_off)
        delay2s = on_off_s2[0,0]/100
        
#        st_s2.plot(color='b',starttime=t1-70+delay2s, endtime=t1+55)    
#        st_s2s = st_s2.slice(starttime=t1-70+delay2s, endtime=t1+55) 
#%% seismic frequency - secondary
        
        tr=st_s2s.normalize()
        
        Fs = tr.stats.sampling_rate;  # sampling rate
        Ts = 1/Fs; # sampling interval
        start =t1-70+delay2s
        end=t1+55
        t = np.arange(0,(end-start),Ts) # time vector
        
        y = tr.data
        
        n = len(y) # length of the signal
        k = np.arange(n)
        T = n/Fs
        frq2 = k/T # two sides frequency range
        
        Y2 = np.fft.fft(y)/n # fft computing and normalization
   
        
#%%
             
        #Primary Seismic
        t1=UTCDateTime(all4[x,1]) 
        t2=t1+120
        st = Stream()
        st = client.get_waveforms(net, sta, '', 'HHZ', t1-60 , t2) 
        st_s1 = st[0]
        st_s1.detrend(type='linear')
        st_s1.detrend(type='demean')
        st_s1.filter("bandpass", freqmin=0.2,freqmax=6)
        
                                    #20
        stream=st_s1.data
        cft=recursive_sta_lta(stream, nsta, nlta)
        trig_on=8                                            #8
        trig_off=1                                        #0.2
    #    plot_trigger(st_s1, cft, trig_on, trig_off) 
        
        on_off_s1 = trigger_onset(cft,trig_on,trig_off)
        delay1s = on_off_s1[0,0]/100
        
        st_s1.plot(color='r',starttime=t1-70+delay1s, endtime=t1+55)    
        st_s1s = st_s1.slice(starttime=t1-70+delay1s, endtime=t1+55)
        
        st_s2.plot(color='b',starttime=t1s-70+delay2s, endtime=t1s+55)    
        st_s2s = st_s2.slice(starttime=t1s-70+delay2s, endtime=t1s+55)
  
        top_v2,top2,corell2 = corel(st_s1s,st_s2s,shift)
        print('seismic timeseries correlation = ', top_v2)
        
        
#%% seismic frequency -primary
        tr=st_s1s.normalize()
        
        Fs = tr.stats.sampling_rate;  # sampling rate
        Ts = 1/Fs; # sampling interval
        start =t1-70+delay2s
        end=t1+55
        t = np.arange(0,(end-start),Ts) # time vector
        
        y = tr.data
        
        n = len(y) # length of the signal
        k = np.arange(n)
        T = n/Fs
        frq = k/T # two sides frequency range
        
        Y = np.fft.fft(y)/n # fft computing and normalization


 #%%       
#         Primary Infrasound
        t1=UTCDateTime(all4[x,0]) 
        t2=t1+120
        st = Stream()
        st = client.get_waveforms(net, sta, '', 'HDF', t1-60 , t2)
        st_a1 = st[0]
        st_a1.detrend(type='linear')
        st_a1.detrend(type='demean')
        st_a1.filter("bandpass", freqmin=0.3,freqmax=2.5)
        st_a1.plot(color='r',starttime=t1-13, endtime=t1-5)
        st_a1s = st_a1.slice(starttime = t1-13  , endtime= t1-5)
    
        # Secondary infrasound
        t1=UTCDateTime(all4[x,2]) 
        t2=t1+120
        st = Stream()
        st = client.get_waveforms(net, sta, '', 'HDF', t1-60 , t2) 
        st_a2 = st[0]
        st_a2.detrend(type='linear')
        st_a2.detrend(type='demean')
        st_a2.filter("bandpass", freqmin=0.3,freqmax=2.5)
        st_a2.plot(color='b',starttime=t1-13, endtime=t1-5)
        st_a2s = st_a2.slice(starttime = t1-13  , endtime= t1-5) 

        top_v,top,corell = corel(st_a1s,st_a2s,shift)
        print('acoustic correlation = ', top_v)
        
  #%% plot frequency spectra
        fig, ax = plt.subplots(2, 1)
        ax[0].plot(frq,abs(Y),'r')
        ax[0].set_xlabel('Freq (Hz)')
        ax[0].set_ylabel('|amplitude|')
        ax[0].set_xlim([0, 5])
        ax[1].plot(frq2,abs(Y2),'b')
        ax[1].set_xlabel('Freq (Hz)')
        ax[1].set_ylabel('|amplitude|')
        ax[1].set_xlim([0, 5])
        ax[0].set_title(start)
        print('\n')
        print('\n')
#%%        
        
#        sum_cora +=top_v
#        sum_cors +=top_v2
#        count += 1
        
    except: 
        f=0   

#print('\n')
#print('Average Acoustic cor =', sum_cors/count)
#print('Average Seismic cor =', sum_cora/count)
#%%

                   