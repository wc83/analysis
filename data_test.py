#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 19 14:58:34 2018

@author: william
"""

from obspy import read
from obspy import Stream
from obspy.signal.trigger import classic_sta_lta, recursive_sta_lta
from obspy.signal.trigger import plot_trigger, trigger_onset
from obspy.clients.earthworm import Client
from obspy import UTCDateTime
import matplotlib.pyplot as plt
import obspy


#%%
#
#stre1 = read("/Users/william/Documents/Data/2019_array_seismics/UL/ARR0/EHZ.D/UL.ARR0..EHZ.D.2019.012")
#
#
#tr1 = stre1[0]
##tr1.plot()
#t1=tr1.stats.starttime
#t2=tr1.stats.endtime
##print(t1,t2)
#
#
#for x in range(0,24):
#    t1=tr1.stats.starttime + x*60*60 + 60
#    t2=t1 + 60*60
#    
#    tr1.filter("bandpass", freqmin=0.2,freqmax=6)
#    tr1b = tr1.slice(starttime = t1, endtime = t2)
#    
##    if max(tr2.data) > 1500:
##        tr1.plot(starttime = t1, endtime = t2, color = 'b')
#        
##%%
#        
#et1=tr1.stats.starttime + 13*60*60 + 29*60 + 0
#et2 = et1 + 120 #tr1.stats.starttime + 0*60*60 + 11*60 +30
#tr1.plot(starttime = et1, endtime = et2 , color = 'r')
#
#tr1b = tr1.slice(starttime = et1, endtime = et2)

#sr = tr1b.stats.sampling_rate
#nsta=int(0.5*sr)                                      #2
#nlta=int(15*sr)                                     #20
#stream=tr1b.data
#cft=recursive_sta_lta(stream, nsta, nlta)
#trig_on=10                                           #8
#trig_off=1                                        #0.2
##plot_trigger(tr2, cft, trig_on, trig_off) 
#
#on_off = trigger_onset(cft,trig_on,trig_off)
#

#%%


# STATION, CHANNEL (DDF --> 400 Hz), NETWWORK AND LOCATION CODES 
sta = 'ARR0' # STATION 
cha = 'HDF' # CHANNEL
net = 'Z4'  # 
loc = ''    # location, it depends mostly of which network you are in. 

# Corner frequency for high-pass filter
#hp_corner = 0.05

# t1. and t2 are in hours:minutes:seconds
# Get data from (Liverpool Winston default) wave server between times t1 and t2 for all stations in stalist      
client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
t1 = UTCDateTime(2019, 1, 16, 15, 3, 10) #the format is year:day_of_the_year:month
t2 = t1 + 440
st = Stream()
st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
#print(st)


st.detrend(type='linear')
st.detrend(type='demean')
st.filter(type='bandpass',freqmin=0.2, freqmax=2)
st.plot(color='g',starttime=t1, endtime=t2)

#


#%%
#
#stre2 = read("/Users/william/Documents/Data/2019_data_collection/2019/STG5/HDF.D/UL.STG5..HDF.D.2019.016")
#
#
#tr2 = stre2[0]
#tr2.filter("bandpass", freqmin=0.2,freqmax=2)
#et1=tr2.stats.starttime + 15*60*60 + 3*60 + 0
#et2 = et1 + 150
#tr2.plot(starttime = et1, endtime = et2 , color = 'b')
#
#ett1=tr2.stats.starttime + 15*60*60 + 0*60 + 0
#ett2 = ett1 +600
#
#tr2b = tr2.slice(starttime = ett1 , endtime = ett2)
#
#sr = tr2b.stats.sampling_rate
#nsta=int(0.5*sr)                                      #2
#nlta=int(20*sr)                                     #20
#stream=tr2b.data
#cft=recursive_sta_lta(stream, nsta, nlta)
#trig_on=10                                           #8
#trig_off=1                                        #0.2
#plot_trigger(tr2b, cft, trig_on, trig_off) 
#
#on_off = trigger_onset(cft,trig_on,trig_off)
#
#
#if len(on_off) >0:
#    
#    event_s = ett1 + (on_off[0,0]/100) - 3
#    event_e = ett1 + (on_off[0,1]/100) + 2
#
#    event_inf = tr2b.slice(starttime = event_s, endtime = event_e)
#
#    event_inf.plot()



#%%

stre2 = read("/Users/william/Documents/Data/2019_data_collection/2019/STG8/HDF.D/UL.STG8..HDF.D.2019.016")


tr2 = stre2[0]
tr2.filter("bandpass", freqmin=0.2,freqmax=2)
et1=tr2.stats.starttime + 15*60*60 + 3*60 + 0 - 10
et2 = et1 + 270 + 10
tr2.plot(starttime = et1 + 20 , endtime = et2 , color = 'b')

ett1=tr2.stats.starttime + 15*60*60 + 3*60 + 0
ett2 = ett1 +290

tr2b = tr2.slice(starttime = ett1 , endtime = ett2)

sr = tr2b.stats.sampling_rate
nsta=int(0.5*sr)                                      #2
nlta=int(20*sr)                                     #20
stream=tr2b.data
cft=recursive_sta_lta(stream, nsta, nlta)
trig_on=10                                           #8
trig_off=1                                        #0.2
#plot_trigger(tr2b, cft, trig_on, trig_off) 

on_off = trigger_onset(cft,trig_on,trig_off)


if len(on_off) >0:
    
    event_s = ett1 + (on_off[0,0]/100) - 3
    event_e = ett1 + (on_off[0,1]/100) 

    event_inf = tr2b.slice(starttime = event_s, endtime = event_e)

#    event_inf.plot()
    
    if len(on_off) > 1:
        event2_s = ett1 + (on_off[1,0]/100) - 10
        event2_e = ett1 + (on_off[1,1]/100) + 15
    
        event_inf2 = tr2b.slice(starttime = event2_s, endtime = event2_e)
    
#        event_inf2.plot(color = 'r')

#fig = plt.figure()  # a new figure window
#ax = fig.add_subplot(1, 1, 1)  # specify (nrows, ncols, axnum)
#
#
#tr2b.spectrogram(log=False,axes=ax)
#plt.ylim(0,5)
#plt.title( str(tr2b.stats.starttime))
#plt.ylabel('frequency [Hz]')
#plt.xlabel('time [s]' )
#
#plt.show()
#%%
#
#stre2 = read("/Users/william/Documents/Data/2019_data_collection/2019/STG2/HDF.D/UL.STG2..HDF.D.2019.016")
#
#
#tr2 = stre2[0]
#tr2.filter("bandpass", freqmin=0.2,freqmax=2)
##et1=tr2.stats.starttime + 15*60*60 + 3*60 + 0 - 10
##et2 = et1 + 420 + 10
#tr2.plot(starttime = et1 + 20 , endtime = et2 + 60 , color = 'b')
##
##ett1=tr2.stats.starttime + 15*60*60 + 3*60 + 0
##ett2 = ett1 +300
##
##tr2b = tr2.slice(starttime = ett1 , endtime = ett2)
##
##sr = tr2b.stats.sampling_rate
##nsta=int(0.5*sr)                                      #2
##nlta=int(20*sr)                                     #20
##stream=tr2b.data
##cft=recursive_sta_lta(stream, nsta, nlta)
##trig_on=10                                           #8
##trig_off=1                                        #0.2
##plot_trigger(tr2b, cft, trig_on, trig_off) 
##
##on_off = trigger_onset(cft,trig_on,trig_off)
##
##
##if len(on_off) >0:
##    
##    event_s = ett1 + (on_off[0,0]/100) - 3
##    event_e = ett1 + (on_off[0,1]/100) + 2
##
##    event_inf = tr2b.slice(starttime = event_s, endtime = event_e)
##
##    event_inf.plot()
#


#%%




stre3 = read("/Users/william/Documents/Data/2019_data_collection/2019/STG1/EHZ.D/UL.STG1..EHZ.D.2019.016")


tr3 = stre3[0]
tr3.filter("bandpass", freqmin=0.2,freqmax=5)       
et1=tr3.stats.starttime + 15*60*60 + 3*60 + 0 - 10
et2 = et1 + 270 + 10
tr3.plot(starttime = et1 + 20, endtime = et2 , color = 'r')


fig = plt.figure()  # a new figure window
ax = fig.add_subplot(1, 1, 1)  # specify (nrows, ncols, axnum)

tr3b = tr3.slice(starttime = et1 , endtime = et2)
tr3b.spectrogram(log=False,axes=ax)
plt.ylim(0,10)
plt.title( str(tr3b.stats.starttime))
plt.ylabel('frequency [Hz]')
plt.xlabel('time [s]' )

plt.show()
