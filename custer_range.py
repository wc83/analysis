#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 14:40:22 2019

@author: william
"""

import os
from collections import OrderedDict
import numpy as np
import obspy
import scipy.signal as sgn
import matplotlib.pyplot as plt 
import matplotlib.mlab as mlab
from scipy import integrate
from obspy.clients.earthworm import Client
from obspy import UTCDateTime
from obspy.signal.trigger import trigger_onset
#from scipy.signal import welch
from obspy import Stream

# STATION, CHANNEL (DDF --> 400 Hz), NETWWORK AND LOCATION CODES 
sta = 'LB01' # STATION 
cha = 'HHZ' # CHANNEL
net = 'Z4'  # 
loc = ''    # location, it depends mostly of which network you are in. 
client = Client('138.253.112.23', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23

t15 = UTCDateTime(2015, 7, 1, 10, 56, 30) #the format is year:day_of_the_year:month
t25 = t15 + 65
st5 = Stream()
st5 = client.get_waveforms(net, sta, '', cha, t15 , t25)
st5.detrend(type='linear')
st5.detrend(type='demean')
st5.filter(type='bandpass',freqmin=1, freqmax=10)
tr5=st5[0].slice(starttime=t15, endtime=t25)
st_c5 = calibrate1(tr5)
print('1/10')

t16 = UTCDateTime(2015, 6, 25, 15, 44, 8) #the format is year:day_of_the_year:month
t26 = t16 + 65
st6 = Stream()
st6 = client.get_waveforms(net, sta, '', cha, t16 , t26)
st6.detrend(type='linear')
st6.detrend(type='demean')
st6.filter(type='bandpass',freqmin=1, freqmax=10)
tr6=st6[0].slice(starttime=t16, endtime=t26)
st_c6 = calibrate1(tr6)
print('2/10')

t17 = UTCDateTime(2015, 7, 23, 15, 0, 30) #the format is year:day_of_the_year:month
t27 = t17 + 65
st7 = Stream()
st7 = client.get_waveforms(net, sta, '', cha, t17 , t27)
st7.detrend(type='linear')
st7.detrend(type='demean')
st7.filter(type='bandpass',freqmin=1, freqmax=10)
tr7=st7[0].slice(starttime=t17, endtime=t27)
st_c7 = calibrate1(tr7)
print('3/10')

t18 = UTCDateTime(2015, 7, 14, 6, 55, 42) #the format is year:day_of_the_year:month
t28 = t18 + 65
st8 = Stream()
st8 = client.get_waveforms(net, sta, '', cha, t18 , t28)
st8.detrend(type='linear')
st8.detrend(type='demean')
st8.filter(type='bandpass',freqmin=1, freqmax=10)
tr8=st8[0].slice(starttime=t18, endtime=t28)
st_c8 = calibrate1(tr8)
print('4/10')

t19 = UTCDateTime(2015, 8, 8, 1, 53, 5) #the format is year:day_of_the_year:month
t29 = t19 + 65
st9 = Stream()
st9 = client.get_waveforms(net, sta, '', cha, t19 , t29)
st9.detrend(type='linear')
st9.detrend(type='demean')
st9.filter(type='bandpass',freqmin=1, freqmax=10)
tr9=st9[0].slice(starttime=t19, endtime=t29)
st_c9 = calibrate1(tr9)
print('5/10')

t110 = UTCDateTime(2016, 12, 17, 1, 40, 45) #the format is year:day_of_the_year:month
t210 = t110 + 65
st10 = Stream()
st10 = client.get_waveforms(net, sta, '', cha, t110 , t210)
st10.detrend(type='linear')
st10.detrend(type='demean')
st10.filter(type='bandpass',freqmin=1, freqmax=10)
tr10=st10[0].slice(starttime=t110, endtime=t210)
st_c10 = calibrate1(tr10)
print('6/10')

t111 = UTCDateTime(2016, 12, 22, 11, 57, 25) #the format is year:day_of_the_year:month
t211 = t111 + 65
st11 = Stream()
st11 = client.get_waveforms(net, sta, '', cha, t111 , t211)
st11.detrend(type='linear')
st11.detrend(type='demean')
st11.filter(type='bandpass',freqmin=1, freqmax=10)
tr11=st11[0].slice(starttime=t111, endtime=t211)
st_c11 = calibrate1(tr11)
print('7/10')

t112 = UTCDateTime(2016, 12, 17, 5, 57, 5) #the format is year:day_of_the_year:month
t212 = t112 + 65
st12 = Stream()
st12 = client.get_waveforms(net, sta, '', cha, t112 , t212)
st12.detrend(type='linear')
st12.detrend(type='demean')
st12.filter(type='bandpass',freqmin=1, freqmax=10)
tr12=st12[0].slice(starttime=t112, endtime=t212)
st_c12 = calibrate1(tr12)
print('8/10')

t113 = UTCDateTime(2017, 3, 20, 15, 36, 15) #the format is year:day_of_the_year:month
t213 = t113 + 65
st13 = Stream()
st13 = client.get_waveforms(net, sta, '', cha, t113 , t213)
st13.detrend(type='linear')
st13.detrend(type='demean')
st13.filter(type='bandpass',freqmin=1, freqmax=10)
tr13=st13[0].slice(starttime=t113, endtime=t213)
st_c13 = calibrate1(tr13)
print('9/10')

t114 = UTCDateTime(2017, 4, 9, 17, 13, 30) #the format is year:day_of_the_year:month
t214 = t114 + 65
st14 = Stream()
st14 = client.get_waveforms(net, sta, '', cha, t114 , t214)
st14.detrend(type='linear')
st14.detrend(type='demean')
st14.filter(type='bandpass',freqmin=1, freqmax=10)
tr14=st14[0].slice(starttime=t114, endtime=t214)
st_c14 = calibrate1(tr14)
print('10/10')

#st_c5.plot(color='k',starttime=t15, endtime=t25)
#st_c6.plot(color='k',starttime=t16, endtime=t26)
#st_c7.plot(color='k',starttime=t17, endtime=t27)
#st_c8.plot(color='k',starttime=t18, endtime=t28)
#st_c9.plot(color='k',starttime=t19, endtime=t29)
#st_c10.plot(color='k',starttime=t110, endtime=t210)
#st_c11.plot(color='k',starttime=t111, endtime=t211)
#st_c12.plot(color='k',starttime=t112, endtime=t212)
#st_c13.plot(color='k',starttime=t113, endtime=t213)
#st_c14.plot(color='k',starttime=t114, endtime=t214)




wav1 = st_c14
wav2 = st_c13
wav3 = st_c12
wav4 = st_c11
wav5 = st_c10
wav6 = st_c9
wav7 = st_c8
wav8 = st_c7
wav9 = st_c6
wav10 = st_c5
shift = 1600

top_v1,top1,corell1 = corel_pi(wav1[0],wav2[0],shift)
top_v2,top2,corell2 = corel_pi(wav1[0],wav3[0],shift)
top_v3,top3,corell3 = corel_pi(wav1[0],wav4[0],shift)
top_v4,top4,corell4 = corel_pi(wav1[0],wav5[0],shift)
top_v5,top5,corell5 = corel_pi(wav1[0],wav6[0],shift)
top_v6,top6,corell6 = corel_pi(wav1[0],wav7[0],shift)
top_v7,top7,corell7 = corel_pi(wav1[0],wav8[0],shift)
top_v8,top8,corell8 = corel_pi(wav1[0],wav9[0],shift)
top_v9,top9,corell9 = corel_pi(wav1[0],wav10[0],shift)



t1n14 = t114 +5
t1n13=t113 + top1/100 +5
t1n12=t112 + top2/100 +5 
t1n11=t111 + top3/100 +5
t1n10=t110 + top4/100 +5
t1n9=t19 + top5/100 +5
t1n8=t18 + top6/100 +5
t1n7=t17 + top7/100 +5
t1n6=t16 + top8/100 +5
t1n5=t15 + top9/100 +5

t2n14 = t1n14 + 65
t2n13 = t1n13 + 65
t2n12 = t1n12 + 65
t2n11 = t1n11 + 65
t2n10 = t1n10 + 65
t2n9 = t1n9 + 65
t2n8 = t1n8 + 65
t2n7 = t1n7 + 65
t2n6 = t1n6 + 65
t2n5 = t1n5 + 65


#%%


st5 = Stream()
st5 = client.get_waveforms(net, sta, '', cha, t1n5 , t2n5)
st5.detrend(type='linear')
st5.detrend(type='demean')
st5.filter(type='bandpass',freqmin=1, freqmax=10)
tr5=st5[0].slice(starttime=t1n5, endtime=t2n5)
nst_c5 = calibrate1(tr5)
print('1/10')


st6 = Stream()
st6 = client.get_waveforms(net, sta, '', cha, t1n6 , t2n6)
st6.detrend(type='linear')
st6.detrend(type='demean')
st6.filter(type='bandpass',freqmin=1, freqmax=10)
tr6=st6[0].slice(starttime=t1n6, endtime=t2n6)
nst_c6 = calibrate1(tr6)
print('2/10')


st7 = Stream()
st7 = client.get_waveforms(net, sta, '', cha, t1n7 , t2n7)
st7.detrend(type='linear')
st7.detrend(type='demean')
st7.filter(type='bandpass',freqmin=1, freqmax=10)
tr7=st7[0].slice(starttime=t1n7, endtime=t2n7)
nst_c7 = calibrate1(tr7)
print('3/10')


st8 = Stream()
st8 = client.get_waveforms(net, sta, '', cha, t1n8 , t2n8)
st8.detrend(type='linear')
st8.detrend(type='demean')
st8.filter(type='bandpass',freqmin=1, freqmax=10)
tr8=st8[0].slice(starttime=t1n8, endtime=t2n8)
nst_c8 = calibrate1(tr8)
print('4/10')


st9 = Stream()
st9 = client.get_waveforms(net, sta, '', cha, t1n9 , t2n9)
st9.detrend(type='linear')
st9.detrend(type='demean')
st9.filter(type='bandpass',freqmin=1, freqmax=10)
tr9=st9[0].slice(starttime=t1n9, endtime=t2n9)
nst_c9 = calibrate1(tr9)
print('5/10')


st10 = Stream()
st10 = client.get_waveforms(net, sta, '', cha, t1n10 , t2n10)
st10.detrend(type='linear')
st10.detrend(type='demean')
st10.filter(type='bandpass',freqmin=1, freqmax=10)
tr10=st10[0].slice(starttime=t1n10, endtime=t2n10)
nst_c10 = calibrate1(tr10)
print('6/10')


st11 = Stream()
st11 = client.get_waveforms(net, sta, '', cha, t1n11 , t2n11)
st11.detrend(type='linear')
st11.detrend(type='demean')
st11.filter(type='bandpass',freqmin=1, freqmax=10)
tr11=st11[0].slice(starttime=t1n11, endtime=t2n11)
nst_c11 = calibrate1(tr11)
print('7/10')


st12 = Stream()
st12 = client.get_waveforms(net, sta, '', cha, t1n12 , t2n12)
st12.detrend(type='linear')
st12.detrend(type='demean')
st12.filter(type='bandpass',freqmin=1, freqmax=10)
tr12=st12[0].slice(starttime=t1n12, endtime=t2n12)
nst_c12 = calibrate1(tr12)
print('8/10')


st13 = Stream()
st13 = client.get_waveforms(net, sta, '', cha, t1n13 , t2n13)
st13.detrend(type='linear')
st13.detrend(type='demean')
st13.filter(type='bandpass',freqmin=1, freqmax=10)
tr13=st13[0].slice(starttime=t1n13, endtime=t2n13)
nst_c13 = calibrate1(tr13)
print('9/10')


st14 = Stream()
st14 = client.get_waveforms(net, sta, '', cha, t1n14 , t2n14)
st14.detrend(type='linear')
st14.detrend(type='demean')
st14.filter(type='bandpass',freqmin=1, freqmax=10)
tr14=st14[0].slice(starttime=t1n14, endtime=t2n14)
nst_c14 = calibrate1(tr14)
print('10/10')



#nst_c5.plot(color='r',starttime=t1n5, endtime=t2n5)
#nst_c6.plot(color='r',starttime=t1n6, endtime=t2n6)
#nst_c7.plot(color='r',starttime=t1n7, endtime=t2n7)
#nst_c8.plot(color='r',starttime=t1n8, endtime=t2n8)
#nst_c9.plot(color='r',starttime=t1n9, endtime=t2n9)
#nst_c10.plot(color='r',starttime=t1n10, endtime=t2n10)
#nst_c11.plot(color='r',starttime=t1n11, endtime=t2n11)
#nst_c12.plot(color='r',starttime=t1n12, endtime=t2n12)
#nst_c13.plot(color='r',starttime=t1n13, endtime=t2n13)
#nst_c14.plot(color='r',starttime=t1n14, endtime=t2n14)



#%%

x=np.linspace(0,len(nst_c14[0].data)/100,len(nst_c14[0].data))


plt.figure(figsize=(12,3))
plt.plot(x,nst_c5[0].data)
#plt.xlabel('Time [s]')
plt.ylabel('Ground Velocity [m/s]')
plt.title('3.90E4')
#plt.savefig('w1.png')

plt.figure(figsize=(12,3))
plt.plot(x,nst_c6[0].data)
#plt.xlabel('Time [s]')
plt.ylabel('Ground Velocity [m/s]')
plt.title('9.15E5')
#plt.savefig('w2.png')

plt.figure(figsize=(12,3))
plt.plot(x,nst_c7[0].data)
#plt.xlabel('Time [s]')
plt.ylabel('Ground Velocity [m/s]')
plt.title('7.34E8')
#plt.savefig('w3.png')

plt.figure(figsize=(12,3))
plt.plot(x,nst_c8[0].data)
#plt.xlabel('Time [s]')
plt.ylabel('Ground Velocity [m/s]')
plt.title('4.05E7')
#plt.savefig('w4.png')

plt.figure(figsize=(12,3))
plt.plot(x,nst_c9[0].data)
#plt.xlabel('Time [s]')
plt.ylabel('Ground Velocity [m/s]')
plt.title('1.80E8')
#plt.savefig('w5.png')

plt.figure(figsize=(12,3))
plt.plot(x,nst_c10[0].data)
#plt.xlabel('Time [s]')
plt.ylabel('Ground Velocity [m/s]')
plt.title('2.80E4')
#plt.savefig('w6.png')

plt.figure(figsize=(12,3))
plt.plot(x,nst_c11[0].data)
#plt.xlabel('Time [s]')
plt.ylabel('Ground Velocity [m/s]')
plt.title('9.49E4')
#plt.savefig('w7.png')

plt.figure(figsize=(12,3))
plt.plot(x,nst_c12[0].data)
#plt.xlabel('Time [s]')
plt.ylabel('Ground Velocity [m/s]')
plt.title('6.66E5')
#plt.savefig('w8.png')

plt.figure(figsize=(12,3))
plt.plot(x,nst_c13[0].data)
#plt.xlabel('Time [s]')
plt.ylabel('Ground Velocity [m/s]')
plt.title('3.22E6')
#plt.savefig('w9.png')

plt.figure(figsize=(12,3))
plt.plot(x,nst_c14[0].data)
#plt.xlabel('Time [s]')
plt.ylabel('Ground Velocity [m/s]')
plt.title('3.79E7')
#plt.savefig('w10.png')











 