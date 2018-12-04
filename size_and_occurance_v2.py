#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 18:39:00 2018

@author: william
"""

import obspy
from obspy import read
import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
from obspy import Stream
from obspy.signal.trigger import classic_sta_lta, recursive_sta_lta
from obspy.signal.trigger import plot_trigger, trigger_onset
from obspy import UTCDateTime
from matplotlib import pyplot
import scipy
from obspy import UTCDateTime


#%% import data
#stre = read("/Users/william/Documents/scanner/output_data/EXP_all_data_stream_month_1.mseed")
per_day = genfromtxt("/Users/william/Documents/scanner/all_stations/Explosions_per_day_V3.csv", delimiter=',',skip_header=1,skip_footer=1)
data = genfromtxt("/Users/william/Documents/scanner/all_stations/EXP_all_coincidence_month_1.csv", delimiter=',',skip_header=1)
num_active= genfromtxt("/Users/william/Documents/scanner/all_stations/num_active_stations.csv", delimiter=',',skip_header=1)
cat= genfromtxt("/Users/william/Documents/scanner/all_stations/Explosion_Catalogue_V3.csv", delimiter=',',skip_header=1)
eq_cat = genfromtxt("/Users/william/Documents/scanner/all_stations/Earthquake_Catalogue_v1.csv", delimiter=',',skip_header=1)
#%% Events per day and per week
#plt.figure(10001)
#plt.plot(per_day[:,1],per_day[:,0])
#plt.ylim([0,max(per_day[:,0])+10])
#plt.xlabel('Day Number')
#plt.ylabel('Number of Explosions')
#plt.title('Explosions per Day')

epw=np.zeros(shape=(1,4))
day=0
week=0
nepw=0
stations = 0
for x in range(0,903):
    nepw += per_day[x,0]
    stations += num_active[x,1]
    day += 1
    if day == 7:
        epw[week][0] = nepw
        epw[week][1] = week
        epw[week][2] = stations/7
        epw[week][3] = 1416787200.0 + (week*7*24*60*60)
            
        week +=1
        day=0
        nepw=0
        stations=0
        if len(per_day)+7-x > 7:
            epw = np.lib.pad(epw, ((0,1),(0,0)), 'constant', constant_values=(0))
            
np.savetxt("/Users/william/Documents/scanner/output_data/explosions_per_week_v3.csv", epw ,delimiter=",",header="epw,week,stations,time")
plt.figure(10002)
plt.plot(epw[0:129,1],epw[0:129,0])
plt.ylim([0,max(epw[0:129,0])+50])
plt.xlabel('Week Number')
plt.ylabel('Number of Explosions')
plt.title('Explosions per Week')

#fig, ax1 = plt.subplots()
#ax2 = ax1.twinx()
#ax1.plot(epw[0:129,1],epw[0:129,0],color='red',label='Events')
#ax1.set_xlabel('Week')
#ax1.set_ylabel('Explosions per Week')
#ax2.plot(epw[0:129,1],epw[0:129,2],color='blue',label='Stations')
#ax2.set_ylabel('Average Number of active stations')
#plt.title('Weekly Station Activity and Detection')
#fig.legend()

#
#%% Station Activity
#plt.figure(2001)
#plt.plot(num_active[:,0],num_active[:,1])
#plt.xlabel('Day')
#plt.ylabel('Number of active stations')
#plt.title('All Station Activity')
#plt.ylim([0,12])
#
#plt.figure(2002)
#plt.plot(num_active[:,0],num_active[:,2])
#plt.xlabel('Day')
#plt.ylabel('Number of active LB stations')
#plt.title('Broadband Station Activity')
#plt.ylim([0,7])
#
#plt.figure(2003)
#plt.plot(num_active[:,0],num_active[:,3])
#plt.xlabel('Day')
#plt.ylabel('Number of active LS stations')
#plt.title('Short Period Station Activity')
#plt.ylim([0,7])

#%%   Explosions AND Activity

plt.figure()
#ax2 = ax1.twinx()
plt.plot(per_day[:,1],per_day[:,0],'r')
plt.xlabel('Day')
plt.ylabel('Explosions per day')
#ax2.plot(num_active[:,0],num_active[:,1],color='blue',label='Stations')
#ax2.set_ylabel('Number of active stations')
plt.title('Event Detection')
#fig.legend()
#
#fig, ax1 = plt.subplots()
#ax2 = ax1.twinx()
#ax1.plot(per_day[:,1],per_day[:,0],color='red',label='Events')
#ax1.set_xlabel('Day')
#ax1.set_ylabel('Explosions per day')
#ax2.plot(num_active[:,0],num_active[:,2],color='green',label='Stations')
#ax2.set_ylabel('Number of active LB stations')
#plt.title(' LB Station Activity and Explosion Detection')
#fig.legend()
#
#fig, ax1 = plt.subplots()
#ax2 = ax1.twinx()
#ax1.plot(per_day[:,1],per_day[:,0],color='red',label='Events')
#ax1.set_xlabel('Day')
#ax1.set_ylabel('Explosions per day')
#ax2.plot(num_active[:,0],num_active[:,3],color='black',label='Stations')
#ax2.set_ylabel('Number of active LS stations')
#plt.title(' LS Station Activity and Explosion Detection')
#fig.legend()

#%% Explosions Vs Activity
#plt.figure(3001)
#plt.scatter(num_active[0:904,1],num_active[0:904,16])
#plt.xlim([0,12])
#plt.ylim([0,100])
#plt.xlabel('Number of Active Stations')
#plt.ylabel('Number of Detected Explosions')
#plt.title('Active Stations vs Explosion Detection')

#plt.figure(3002)
#plt.scatter(epw[:,2],epw[:,0])
##plt.xlim([0,12])
##plt.ylim([0,100])
#plt.xlabel('Average Number of Active Stations')
#plt.ylabel('Number of Detected Explosions')
#plt.title('Weekly Active Stations vs Explosion Detection')

#%%

#rho,p=scipy.stats.spearmanr(num_active[0:904,1],per_day[0:904,0])
#print('Daily EXP: rho =',rho)
#rho1,p1=scipy.stats.spearmanr(epw[0:129,2],epw[0:129,0])
#print('Weekly EXP: rho =',rho1)

#%%

repose=[]
for x in range(1,len(cat)):
    dt=(cat[x,0]-cat[x-1,0])/60
    if dt<300:
        repose.append(dt)
    
plt.figure(4001)
plt.hist(repose,bins=100)
plt.xlabel('Repose time [mins]')
plt.ylabel('Occurance [#]')
plt.title('Explosion Repose time')

eq_repose=[]
for x in range(1,len(eq_cat)):
    edt=(eq_cat[x,0]-eq_cat[x-1,0])/60
    if edt<1000:
        eq_repose.append(edt)
    
plt.figure(4002)
plt.hist(eq_repose,bins=100)
plt.xlabel('Repose time [mins]')
plt.ylabel('Occurance [#]')
plt.title('Earthquake Repose time')







##%%
#
#eqpd=np.zeros(shape=(1,3))
#neqpd=0
#start=1416787200 - 24*60*60     #time stamp in seconds of 2014-11-24T00:00:00.000000
#for p in range(0,905):
#    start=start+(24*60*60)
#    end=start+(24*60*60)
##    print(start)
#    for x in range(0,len(eq_cat)):
#        if start < eq_cat[x,0] < end:
#            neqpd += 1
#    eqpd[p][0]=p+1
#    eqpd[p][1]=neqpd
#    eqpd[p][2]=num_active[p,1]
#    neqpd=0
#    eqpd = np.lib.pad(eqpd, ((0,1),(0,0)), 'constant', constant_values=(0))
#    
#plt.figure(50001)
#plt.plot(eqpd[0:904,0],eqpd[0:904,1])
#plt.xlabel('Day')
#plt.ylabel('Number of Earthquakes')
#plt.title('Earthquakes per Day') 
#
#
#eqpw=np.zeros(shape=(1,3))
#neqpw=0
#start=1416787200 - 7*24*60*60     #time stamp in seconds of 2014-11-24T00:00:00.000000
#for p in range(0,129):
#    start=start+(7*24*60*60)
#    end=start+(7*24*60*60)
##    print(start)
#    for x in range(0,len(eq_cat)):
#        if start < eq_cat[x,0] < end:
#            neqpw += 1
#    eqpw[p][0]=p+1
#    eqpw[p][1]=neqpw
#    eqpw[p][2]=epw[p,2]
#    neqpw=0
#    eqpw = np.lib.pad(eqpw, ((0,1),(0,0)), 'constant', constant_values=(0))
#    
#plt.figure(50002)
#plt.plot(eqpw[0:129,0],eqpw[0:129,1])
#plt.xlabel('Week')
#plt.ylabel('Number of Earthquakes')
#plt.title('Earthquakes per Week') 
#
#
#
##fig, ax1 = plt.subplots()
#ax2 = ax1.twinx()
#ax1.plot(eqpd[0:904,0],eqpd[0:904,1],color='red',label='Events')
#ax1.set_xlabel('Day')
#ax1.set_ylim([0,15])
#ax1.set_ylabel('Explosions per day')
#ax2.plot(eqpd[0:904,0],eqpd[0:904,2],color='blue',label='Stations')
#ax2.set_ylabel('Number of active stations')
#ax2.set_ylim([0,12])
#plt.title('Station Activity and Earthquake Detection')
#fig.legend()
#
#
##fig, ax1 = plt.subplots() 
#ax2 = ax1.twinx()
#ax1.plot(eqpw[0:129,0],eqpw[0:129,1],color='red',label='Events')
#ax1.set_xlabel('Week')
#ax1.set_ylim([0,80])
#ax1.set_ylabel('Earthquakes per Week')
#ax2.plot(eqpw[0:129,0],eqpw[0:129,2],color='blue',label='Stations')
#ax2.set_ylabel('Number of active stations')
#ax2.set_ylim([0,12])
#plt.title('Station Activity and Earthquake Detection')
#fig.legend()
#
#
#
#       
##%%    
#
#plt.figure(6001)
#plt.scatter(eqpd[0:905,2],eqpd[0:905,1])
#plt.xlim([0,12])
##plt.ylim([0,35])
#plt.xlabel('Number of Active Stations')
#plt.ylabel('Number of Detected Earthquakes')
#plt.title('Active Stations vs Earthquake Detection')
#
#plt.figure(6002)
#plt.scatter(eqpw[:,2],eqpw[:,1])
##plt.xlim([0,12])
##plt.ylim([0,100])
#plt.xlabel('Average Number of Active Stations')
#plt.ylabel('Number of Detected Earthquakes')
#plt.title('Weekly Active Stations vs Earthquake Detection')
#
##%%
#
#rho2,p2=scipy.stats.spearmanr(eqpd[0:904,2],eqpd[0:904,1])
#print('Daily EQ: rho =',rho2)
#rho3,p3=scipy.stats.spearmanr(eqpw[0:129,2],eqpw[0:129,1])
#print('Weekly EQ: rho =',rho3)
#
##%%
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#














##%%    Seismic calibrations
#
##Broadband
##LB01 is redeployed with different calibration factors
#LB01sc1=0.000001/750            # before 2015-12-05T00:00:01.000000Z
#LB01sc2=(10*0.000001)/(750*256) # 2015-12-05T00:00:01 - 2016-06-15T00:00:01 + maybe till end
#
#LB02sc=0.000001/750
#LB03sc=0.000001/(750*256)
#LB04sc=0.000001/(750*256)
#LB05sc=0.000001/(750*256)
#LB06sc=(10*0.000001)/750
#
## Short period
#LSsc=0.000000122/800
#
##%%    Acoustic calibrations
#
#LB01ac = 0.000001/0.0250
#LB02ac = 0.000001/0.0250
#LB03ac = 0.000001/(0.0250*256)
#LB04ac = 0.000001/(0.0250*256)
#LB05ac = 0.000001/(0.0250*256)
#
##%%   Station distances
#
#
#
##%% constants
#
#rho_atmos = 1.2 # kg/m3 at 15˚C
#c_atmos = 340 # m/s at 15˚C
#
#rho_earth = 2000 # kg/m3
#c_earth = 2500 # m/s 
#
##%%
#
#
#
#

##%% "Energy" of each event
#
##################### Need to use proper equations with calibrations for all below #######################
#
#event_stream = Stream()
#event_list=np.zeros(shape=(1,1))
#event_count=0
#
#event_stream.append(stre[0])
#event_list[event_count]=stre[0].stats.starttime.timestamp
#event_list = np.lib.pad(event_list, ((0,1),(0,0)), 'constant', constant_values=(0))
#event_count +=1                            
#for x in range(1,len(stre)):
#    if stre[x].stats.station == "LB01":
#        event_stream.append(stre[x])
#        event_list[event_count]=stre[x].stats.starttime.timestamp
#        event_list = np.lib.pad(event_list, ((0,1),(0,0)), 'constant', constant_values=(0))
#        event_count +=1 
#    else:
#        rt=stre[x].stats.starttime.timestamp
#        near,ix=find_nearest(event_list[:,0], rt)
#        if abs(near-rt) > 60:
#            event_stream.append(stre[x])
#            event_list[event_count]=stre[x].stats.starttime.timestamp
#            event_list = np.lib.pad(event_list, ((0,1),(0,0)), 'constant', constant_values=(0))
#            event_count +=1 
#
#
##print(len(data_stream))
#
#sr = 100
#nsta=int(1*sr)                                      
#nlta=int(10*sr)                                     
#trig_on=2.5                                          
#trig_off=0.05
#
##for x in range(0,len(data_stream)):
#for x in range(0,10):
#    data_s=event_stream[x].data
#    max_a = data_s.max()
#    min_a = data_s.min()
#    p2p= max_a-min_a
#    cft=recursive_sta_lta(data_s, nsta, nlta)
##    plot_trigger(sq_stream[x], cft, trig_on, trig_off)     
#    on_off = trigger_onset(cft,trig_on,trig_off)                                     
#    start = event_stream[x].stats.starttime
#    tr = event_stream[x].slice(starttime=start+(on_off[0,0]/sr) , endtime=start+(on_off[0,1]/sr)) 
#    print('event:',event_stream[x].stats.starttime ,'from station: ',stre[x].stats.station,', has energy: ',sum(np.square(tr.data)),' and peak to peak:', p2p)
#    plt.figure(x)
#    plt.plot(tr) 
#    plt.figure(x+20)
#    plt.plot(np.square(tr.data))   
#     
##%% energy information
#    
#        
#day_one = 1416787200.0
#day_energy=0
#day_count=0
#e_count=0
#days=0
#energy_each_event = av_energy_list =np.zeros(shape=(1,1))
#av_energy_list =np.zeros(shape=(1,1))
#total_energy_list =np.zeros(shape=(1,1))
#
#for x in range(0,len(per_day)):#len(per_day)
#    day_start = day_one + x*24*60*60
#    day_end = day_start + 24*60*60 - 0.01
##    print(UTCDateTime(day_start),' to', UTCDateTime(day_end))
#    for p in range(0,len(event_stream)):
#        if day_start < event_stream[p].stats.starttime.timestamp < day_end :
#            # ADD IN CALIBRATIONS #
#            day_count += 1
#            data_s=event_stream[p].data
#            max_a = data_s.max()
#            min_a = data_s.min()
#            p2p= max_a-min_a
#            cft=recursive_sta_lta(data_s, nsta, nlta)
#        #    plot_trigger(sq_stream[x], cft, trig_on, trig_off)     
#            on_off = trigger_onset(cft,trig_on,trig_off)                                     
#            start = event_stream[p].stats.starttime
#            tr = event_stream[p].slice(starttime=start+(on_off[0,0]/sr) , endtime=start+(on_off[0,1]/sr)) 
#            event_energy = sum(np.square(tr.data))
#            day_energy += event_energy
##            print('event:',event_stream[p].stats.starttime ,'station: ',event_stream[p].stats.station,', energy: ',event_energy,', peak to peak:', p2p)
#            energy_each_event[e_count] = event_energy
#            energy_each_event = np.lib.pad(energy_each_event, ((0,1),(0,0)), 'constant', constant_values=(0))
#            e_count += 1
#            
##            if event_energy > 1e16:
##                tr.plot()
#            
#    print('total energy in the day:',UTCDateTime(day_start),'=', day_energy)
#    av_day_energy = day_energy/day_count
#    print('average energy in day:',UTCDateTime(day_start),'=', av_day_energy)
#    
#    av_energy_list[days]= av_day_energy
#    total_energy_list[days] = day_energy 
#    av_energy_list = np.lib.pad(av_energy_list, ((0,1),(0,0)), 'constant', constant_values=(0))
#    total_energy_list = np.lib.pad(total_energy_list, ((0,1),(0,0)), 'constant', constant_values=(0))
#    day_count = 0
#    day_energy=0
#    days += 1
#    
#    
#
#
#plt.figure(2001)
#plt.plot(av_energy_list)
#
#plt.figure(2002)
#plt.plot(total_energy_list)
#
#plt.figure(2003)
#plt.plot(energy_each_event)























