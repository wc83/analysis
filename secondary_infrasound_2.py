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

#sta = 'LB01' # STATION 
#cha = 'HHZ' # CHANNEL
net = 'Z4'  # 
loc = ''    # location, it depends mostly of which network you are in. 
client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23

#%% read in data
reclen = 512
chunksize = 100000 * reclen # 100000 = Around 50 MB 
with io.open("/Users/william/Documents/scanner/output_data/m32.mseed", "rb") as fh:
#         just month 2
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

pands = np.zeros(shape=(0,5))
all4 = np.zeros(shape=(0,7))
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
                pands[count][0]=saved[x,0] #secondary time
                pands[count][1]=last_p #last primary time
                pands[count][2]=saved[x,2] #station
                pands[count][3]=Eper #energy percentage of primary
                
                if saved[x+1,0] - saved[x,0] > 600:
                    pands[count][4]=saved[x+1,0] #next exposion time
                    count +=1
                else:
                    if saved[x+2,0] - saved[x,0] > 600:
                        pands[count][4]=saved[x+2,0] #next exposion time
                        count +=1
                    else:
                        pands[count][4]=saved[x+3,0] #next exposion time
                        count +=1 
                        
                            
                    
                    
                s_arrival = saved[x,0]
                next_p = pands[count-1][4]
                try:
                    tr_inf,all_4 = s_inf(last_p,s_arrival,new_station,next_p)
                    if all_4[0][0] != 0:
                        all4 = np.lib.pad(all4, ((0,1),(0,0)), 'constant', constant_values=(0))
                        all4[num][0]=all_4[0][0]
                        all4[num][1]=all_4[0][1]
                        all4[num][2]=all_4[0][2]
                        all4[num][3]=all_4[0][3]
                        all4[num][4]=all_4[0][4]
                        all4[num][5]=all_4[0][5]
                        all4[num][6]=new_station
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
                if saved[x+1,0] - saved[x,0] > 600:
                    pands[count][4]=saved[x+1,0] #next exposion
                    count +=1
                else:
                    if saved[x+2,0] - saved[x,0] > 600:
                        pands[count][4]=saved[x+2,0] #next exposion
                        count +=1
                    else:
                        pands[count][4]=saved[x+3,0] #next exposion
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
                        all4[num][4]=all_4[0][4]
                        all4[num][5]=all_4[0][5]
                        all4[num][6]=new_station
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




 

     #%%
shift = 1400
sum_cora=0
sum_cors=0
count=0
for x in range(0,len(all4)):
    try:

        if all4[x,6] == 1:
            sta = 'LB01'
        if all4[x,6] == 2:
            sta = 'LB02'
        if all4[x,6] == 3:
            sta = 'LB03'
        if all4[x,6] == 4:
            sta = 'LB04'
        if all4[x,6] == 5:
            sta = 'LB05'
        if all4[x,6] == 6:
            sta = 'LB06'
        if all4[x,6] == 7:
            sta = 'LS01'
        if all4[x,6] == 8:
            sta = 'LS702'
        if all4[x,6] == 9:
            sta = 'LS03'
        if all4[x,6] == 10:
            sta = 'LS04'
        if all4[x,6] == 11:
            sta = 'LS05'
        if all4[x,6] == 12:
            sta = 'LS06'
            
#         Primary Infrasound
        t1=UTCDateTime(all4[x,0]) 
        t2=t1+120
        st = Stream()
        st = client.get_waveforms(net, sta, '', 'HDF', t1-60 , t2)
        st_a1 = st[0]
        st_a1.detrend(type='linear')
        st_a1.detrend(type='demean')
        st_a1.filter("bandpass", freqmin=0.3,freqmax=2.5)
        
        st_a1s = st_a1.slice(starttime = t1-13  , endtime= t1-5)
    
        # Secondary infrasound
        t1b=UTCDateTime(all4[x,2]) 
        t2b=t1b+120
        st = Stream()
        st = client.get_waveforms(net, sta, '', 'HDF', t1b-60 , t2b) 
        st_a2 = st[0]
        st_a2.detrend(type='linear')
        st_a2.detrend(type='demean')
        st_a2.filter("bandpass", freqmin=0.3,freqmax=2.5)
        
        st_a2s = st_a2.slice(starttime = t1b-13  , endtime= t1b-5) 
        
        # next primary infrasound
        t1c=UTCDateTime(all4[x,4]) 
        t2c=t1c+120
        st = Stream()
        st = client.get_waveforms(net, sta, '', 'HDF', t1c-60 , t2c) 
        st_a3 = st[0]
        st_a3.detrend(type='linear')
        st_a3.detrend(type='demean')
        st_a3.filter("bandpass", freqmin=0.3,freqmax=2.5)
        
        st_a3s = st_a3.slice(starttime = t1c-13  , endtime= t1c-5) 

        top_v1,top,corell1 = corel(st_a1s,st_a2s,shift)
        top_v2,top,corell2 = corel(st_a2s,st_a3s,shift)
        top_v3,top,corell3 = corel(st_a1s,st_a3s,shift)
        
        
        if abs(top_v1) > 0.6:
            print('acoustic correlation Primary - Secondary = ', top_v1)
            print('acoustic correlation Primary - Next Primary = ', top_v3)
            print('acoustic correlation Secondary - Next Primary = ', top_v2)
            print('primary to primary repose time = ', (t1c-t1)/60 , 'mins')
            st_a1.plot(color='r',starttime=t1-13, endtime=t1-5)
            st_a2.plot(color='b',starttime=t1b-13, endtime=t1b-5)
            st_a3.plot(color='g',starttime=t1c-13, endtime=t1c-5)
            count +=1

#        sum_cora +=top_v
#        count += 1
        
    except: 
        f=0   

#print('number of secondaries = ', count)
#print('Average Acoustic cor =', sum_cora/count)
#print('Average Seismic cor =', sum_cora/count)
#%%
#print("small secondaries (<10%) =", s_tot)
#print("10% + secondaries =", sl_tot)
#print('unmatched stations =' , un_tot)
#print('Primaries =', p_tot)
#print("total_secondaries =", sl_tot + s_tot )
#print("total_events =", sl_tot + s_tot + p_tot )

#%% plot size hiatogram for all small secondaries
#plt.figure(4001)
#plt.hist(pands[:,3],bins=50)
#plt.xlabel("Secondary Energy compared to Primary [%]")
#plt.ylabel("Occurance [#]")
#plt.title("Secondary Energies")

#%% RESULTS

# less energy for small secondary
#small secondaries = 2876
#larger secondaries = 1009
#unmatched stations = 130
#Primaries = 12497
#total_secondaries = 3885
#total_events = 16382

# 1/2 energy or less for small secondary
#small secondaries = 2560
#larger secondaries = 1325
#unmatched stations = 130
#Primaries = 12497
#total_secondaries = 3885
#total_events = 16382

# 1/3 energy or less for small secondary
#small secondaries = 2366
#larger secondaries = 1519
#unmatched stations = 130
#Primaries = 12497
#total_secondaries = 3885
#total_events = 16382

# 1/4 energy or less for small secondary
#small secondaries = 2228
#larger secondaries = 1657
#unmatched stations = 130
#Primaries = 12497
#total_secondaries = 3885
#total_events = 16382


#%% Get waveforms and see which data are not available online
#
#cha = 'HHZ' # CHANNEL
#net = 'Z4'  # 
#loc = ''    # location, it depends mostly of which network you are in. 
#client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
#unavailable=[]
#for x in range(0,len(pands)):
#    if pands[x,2] == 1:
#        sta='LB01'
#    if pands[x,2] == 2:
#        sta='LB02'
#    if pands[x,2] == 3:
#        sta='LB03'
#    if pands[x,2] == 4:
#        sta='LB04'
#    if pands[x,2] == 5:
#        sta='LB05'
#    if pands[x,2] == 6:
#        sta='LB06'
#    if pands[x,2] == 7:
#        sta='LS01'
#    if pands[x,2] == 8:
#        sta='LS02'
#    if pands[x,2] == 9:
#        sta='LS03'
#    if pands[x,2] == 10:
#        sta='LS04'
#    if pands[x,2] == 11:
#        sta='LS05'
#    if pands[x,2] == 12:
#        sta='LS06'
#    
#
#    arrival2=pands[x,0]
#    arrival1=pands[x,1]
#    
#    try:
#        t1 = UTCDateTime(arrival1)- 60 #the format is year:day_of_the_year:month
#        t2 =  UTCDateTime(arrival2) + 120
#        st= Stream()
#        st = client.get_waveforms(net, sta, '', cha, t1 , t2)
#        st.detrend(type='linear')
#        st.detrend(type='demean')
#        st.filter(type='bandpass',freqmin=0.5, freqmax=6)
#        st.plot(color='b',starttime=t1, endtime=t2)
#    except:
#        print("data not available at time: ",UTCDateTime(arrival1))
#        unavailable.append(arrival1)
#    
        
    
    




                   