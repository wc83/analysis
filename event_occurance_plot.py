#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 21 13:34:37 2018

@author: william
"""
from obspy.core import read
from obspy.signal.cross_correlation import correlate
import glob
import time
import obspy
import numpy as np
from numpy import argmax
import matplotlib.pyplot as plt
import csv
from obspy.core import UTCDateTime

ts=[]
jd=[]
yr=[]
mo=[]
da=[]
hr=[]
mi=[]
se=[]
ms=[]
dl=[]
with open('/Users/william/Documents/Python_scripts/EXP_times_test_function.csv') as f:
    reader = csv.reader(f)
    for row in reader:
        ts.append(row[0])   #timestamp
        jd.append(row[1])   #julian day 
        yr.append(row[2])   #year
        mo.append(row[3])   #month
        da.append(row[4])   #day
        hr.append(row[5])   #hour
        mi.append(row[6])   #min
        se.append(row[7])   #sec
        ms.append(row[8])   #milisec
        dl.append(row[9])   #day length
    
#%%

et=[]
for x in range(0,len(ts)-1):
    et.append(UTCDateTime(float(ts[x+1])))
      
diff=[]
for x in range(1,len(et)):
    t_diff=et[x]-et[x-1]
    diff.append(t_diff/60)
    diff.sort()
    diff = [x for x in diff if x <= (6*60)] #longer than 6hrs, considered to have missed an event    
plt.hist(diff,bins=np.arange(min(diff), max(diff) + 10, 10))
plt.xlabel('Time Between Explosions (mins)')
plt.ylabel('Frequency')
plt.title('Repose Times')
plt.xlim(0,300)
#plt.savefig('time_gap_exp_events_v3')
#%%

et=[]
repose=0     
av_diff=[]
det=[]

day=da[1]
day_number=0
week=0
event=0
event2=0
epd=[]
dayn=[]
epw=[]
weekn=[]
week_number=0

for x in range(1,len(ts)):
    if float(dl[x]) < 86399:
        if da[x] != day:
            day=da[x]
            day_number += 1
    else:
        if da[x] == day:
            det.append(float(ts[x]))
        else:
            for p in range(1,len(det)):
                t_diff=det[p]-det[p-1]
                repose += t_diff/60
            av_diff.append(repose/(len(det)-1))
            repose=0

            det.clear()
            det.append(float(ts[x]))
            day=da[x]
            
fig9 = plt.figure() 
plt.hist(av_diff,bins=np.arange(min(av_diff), max(av_diff) + 10, 10))
plt.xlabel('Average Time Between Explosions (mins)')
plt.ylabel('Frequency')
plt.title('Average Repose Times per Day')
plt.xlim(0,300)
#plt.savefig('average_time_gap_exp_events_v3')
#%%

day=da[1]
day_number=0
week=0
event=0
event2=0
epd=[]
dayn=[]
epw=[]
weekn=[]
week_number=0

for x in range(1,len(ts)):
    if float(dl[x]) < 86399:
        if da[x] != day:
            day=da[x]
            day_number += 1
    else:
        if da[x] == day:
            event += 1
            event2 +=1
        else:
            day_number += 1
            dayn.append(day_number)
            epd.append(event)
            
            week += 1
            if week==7:
                week_number+=1
                epw.append(event2)
                weekn.append(week_number)
                week=0
                event2=0
            
            event=1
            event2 +=1
            day=da[x]
day_number += 1
dayn.append(day_number)
epd.append(event)
week += 1
if week==7:
    week_number+=1
    epw.append(event2)
    weekn.append(week_number) 


#
fig2 = plt.figure()
plt.plot(dayn,epd)
plt.ylim([0,max(epd)+10])
plt.xlabel('Day Number')
plt.ylabel('Number of Explosions')
plt.title('Explosions per Day')
#plt.savefig('Daily_explosions_short_removed_v3')
#
fig3 = plt.figure()
plt.plot(weekn,epw)
plt.ylim([0,max(epw)+100])
plt.xlabel('Week Number')
plt.ylabel('Number of Explosions')
plt.title('Explosions per Week')
#plt.savefig('Weekly_explosions_short_removed_v3')