#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 10:18:40 2018

@author: william
"""

import obspy
import io
from obspy import read
from obspy.clients.earthworm import Client
from obspy import UTCDateTime
from obspy import Stream
import numpy as np
import matplotlib.pyplot as plt

#%% #whole dataset %age 1Hz and 2Hz

reclen = 512
chunksize = 100000 * reclen # Around 50 MB

sum_b1 =0
sum_b2 =0
sum_eb1 =0
sum_eb2 =0

sum_events = 0
week_start = 1416787200
week_end = 1416787200 + 7*24*60*60 
week =0
weekly_freq_p=np.zeros(shape=(0,5))
aw=0
       
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
        
            
        # For each chunck of time, do analysis
#        print(st)
        for x in range (0,int(len(st))):
            if st[x].stats.station == "LB06":# or "LB02" or "LB03" or "LB04" or "LB05" or "LB06" :
                if week_start < st[x].stats.starttime.timestamp < week_end:
                    
                    t1 = st[x].stats.starttime
                    t2 = st[x].stats.endtime
                    window=t2-t1
                    tr_data=st[x].data
                    m=np.mean(tr_data)
                    tr_data = tr_data-m
                    famp = abs(np.fft.rfft(tr_data))
                    famp2=famp*famp
                    
                        # one Hz
                    one = int( window +1)
                    num1=0
                    for t in range(0,one,1):
                        n1=famp[t]                 
                        num1 += n1
                    p1hz=num1/sum(famp)
                    
                    # two Hz
                    two = int((2*window)+1)
                    num2=0
                    for t in range(0,two,1):
                        n2=famp[t]                  
                        num2 += n2
                    p2hz=num2/sum(famp)
                    
                         # energy one Hz
                    one = int( window +1)
                    num1e=0
                    for t in range(0,one,1):
                        n1e=famp2[t]                 
                        num1e += n1e
                    pe1hz=num1e/sum(famp2)
                    
                    # energy two Hz
                    two = int((2*window)+1)
                    num2e=0
                    for t in range(0,two,1):
                        n2e=famp2[t]                  
                        num2e += n2e
                    pe2hz=num2e/sum(famp2)
                    

                    sum_events +=1
                    sum_b1 += p1hz
                    sum_b2 += p2hz
                    sum_eb1 += pe1hz
                    sum_eb2 += pe2hz
                    
                if st[x].stats.starttime.timestamp > week_end:
                    if sum_events > 1:
                        
                        av_b1 = sum_b1/sum_events 
                        av_b2 = sum_b2/sum_events
                        av_eb1 = sum_eb1/sum_events 
                        av_eb2 = sum_eb2/sum_events
                        
                        weekly_freq_p = np.lib.pad(weekly_freq_p, ((0,1),(0,0)), 'constant', constant_values=(0))
                        
                        weekly_freq_p[aw][0]=week + 1
                        weekly_freq_p[aw][1]=av_b1
                        weekly_freq_p[aw][2]=av_b2
                        weekly_freq_p[aw][3]=av_eb1
                        weekly_freq_p[aw][4]=av_eb2
                        
                        aw += 1
   
                    else:
                        print("LB inactive in week", week+1)
                    
                    sum_b1 =0
                    sum_b2 =0
                    sum_eb1 =0
                    sum_eb2 = 0
                    sum_events =0
                    for p in range(0,120):

                        week_start=week_end
                        week_end += 7*24*60*60
                        week += 1
                        if week_start < st[x].stats.starttime.timestamp < week_end:

                            t1 = st[x].stats.starttime
                            t2 = st[x].stats.endtime
                            window=t2-t1
                            tr_data=st[x].data
                            m=np.mean(tr_data)
                            tr_data = tr_data-m
                            famp = abs(np.fft.rfft(tr_data))
        #                    famp2=famp*famp
                            
                                # one Hz
                            one = int( window +1)
                            num1=0
                            for t in range(0,one,1):
                                n1=famp[t]                 
                                num1 += n1
                            p1hz=num1/sum(famp)
                            
                            # two Hz
                            two = int((2*window)+1)
                            num2=0
                            for t in range(0,two,1):
                                n2=famp[t]                  
                                num2 += n2
                            p2hz=num2/sum(famp)
                            
                                # energy one Hz
                            one = int( window +1)
                            num1e=0
                            for t in range(0,one,1):
                                n1e=famp2[t]                 
                                num1e += n1e
                            pe1hz=num1e/sum(famp2)
                            
                            # energy two Hz
                            two = int((2*window)+1)
                            num2e=0
                            for t in range(0,two,1):
                                n2e=famp2[t]                  
                                num2e += n2e
                            pe2hz=num2e/sum(famp2)
                            
        
                            sum_events +=1
                            sum_b1 += p1hz
                            sum_b2 += p2hz
                            sum_eb1 += pe1hz
                            sum_eb2 += pe2hz
                            break

#np.savetxt("/Users/william/Documents/scanner/analysis/weekly_LS06_freq.csv", weekly_freq,delimiter=",",header="week,dominant,central,bandwidth_50")

#plt.figure(10)
#plt.plot(weekly_freq_p[:,0],weekly_freq_p[:,1],'bx')
#plt.ylim((0,1))
#plt.xlim((0,140))
#plt.xlabel("Week") 
#plt.ylabel("Average Percentage Frequency below 1Hz ")
#plt.title("Frequency below 1Hz")
#
plt.figure(11)
plt.plot(weekly_freq_p[:,0],weekly_freq_p[:,2],'rx-')
plt.ylim((0,1))
plt.xlim((0,140))
plt.xlabel("Week")
plt.ylabel("Average Percentage Frequency below 2Hz ")
plt.title("Frequency below 2Hz")
#
#plt.figure(12)
#plt.plot(weekly_freq_p[:,0],1 - weekly_freq_p[:,2],'bx')
#plt.ylim((0,1)) 
#plt.xlim((0,140))    
#plt.xlabel("Week")   
#plt.ylabel("Average Percentage Frequency above 2Hz ")
#plt.title("Frequency above 2Hz") 
#
#plt.figure(13)
#plt.plot(weekly_freq_p[:,0],weekly_freq_p[:,2] - weekly_freq_p[:,1],'bx')
#plt.ylim((0,1)) 
#plt.xlim((0,140))    
#plt.xlabel("Week")   
#plt.ylabel("Average Percentage Frequency between 1 and 2Hz ")
#plt.title("Frequency between 1 and 2Hz") 

#plt.figure(14)
#plt.plot(weekly_freq_p[:,0],weekly_freq_p[:,1],'rx-')
#plt.plot(weekly_freq_p[:,0],weekly_freq_p[:,2] - weekly_freq_p[:,1],'bx-')
#plt.plot(weekly_freq_p[:,0],1 - weekly_freq_p[:,2],'kx-')
#plt.ylim((0,1)) 
#plt.xlim((0,140))    
#plt.xlabel("Week")   
#plt.ylabel("Average Frequency Percentages")
#plt.title("Percentage Division") 
#plt.legend(('<1Hz','1-2Hz','2Hz<'))

#plt.figure(15)
#plt.plot(weekly_freq_p[:,0],weekly_freq_p[:,2],'rx-')
#plt.plot(weekly_freq_p[:,0],1 - weekly_freq_p[:,2],'kx-')
#plt.ylim((0,1)) 
#plt.xlim((0,140))    
#plt.xlabel("Week")   
#plt.ylabel("Average Frequency Percentages")
#plt.title("Frequency Amplitude Division") 
#plt.legend(('<2Hz','2Hz<'))




#plt.figure(20)
#plt.plot(weekly_freq_p[:,0],weekly_freq_p[:,3],'bx')
#plt.ylim((0,1))
#plt.xlim((0,140))
#plt.xlabel("Week") 
#plt.ylabel("Average Energy Percentage below 1Hz ")
#plt.title("Energy below 1Hz")

plt.figure(21)
plt.plot(weekly_freq_p[:,0],weekly_freq_p[:,4],'rx-')
plt.ylim((0,1))
plt.xlim((0,140))
plt.xlabel("Week")
plt.ylabel("Average Energy Percentage  below 2Hz ")
plt.title("Energy below 2Hz")
#
#plt.figure(22)
#plt.plot(weekly_freq_p[:,0],1 - weekly_freq_p[:,4],'bx')
#plt.ylim((0,1)) 
#plt.xlim((0,140))    
#plt.xlabel("Week")   
#plt.ylabel("Average Energy Percentage above 2Hz ")
#plt.title("Energy above 2Hz") 
#
#plt.figure(23)
#plt.plot(weekly_freq_p[:,0],weekly_freq_p[:,4] - weekly_freq_p[:,3],'bx')
#plt.ylim((0,1)) 
#plt.xlim((0,140))    
#plt.xlabel("Week")   
#plt.ylabel("Average Energy Percentage between 1 and 2Hz ")
#plt.title("Energy between 1 and 2Hz") 

#plt.figure(24)
#plt.plot(weekly_freq_p[:,0],weekly_freq_p[:,4],'rx-')
#plt.plot(weekly_freq_p[:,0],1 - weekly_freq_p[:,4],'kx-')
#plt.ylim((0,1)) 
#plt.xlim((0,140))    
#plt.xlabel("Week")   
#plt.ylabel("Average Energy Percentages")
#plt.title("Energy Division") 
#plt.legend(('<2Hz','2Hz<'))

#%% whole data set - dominant, central, 50% bandwidth
#                            
#
###stre = read("/Users/william/Documents/scanner/output_data/EXP_all_data_stream_2_month_2.mseed")
###print(stre)
#
#reclen = 512
#chunksize = 100000 * reclen # Around 50 MB
#        
#
#        
#sum_peak=0
#sum_cf=0
#sum_bw =0
#sum_events = 0
#week_start = 1416787200
#week_end = 1416787200 + 7*24*60*60 
#week =0
#weekly_freq=np.zeros(shape=(0,4))
#aw=0
#
#
#        #whole dataset
#with io.open("/Users/william/Documents/scanner/output_data/m32.mseed", "rb") as fh:
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
#        # For each chunck of time, do analysis
##        print(st)
#        for x in range (0,int(len(st))):
##            if st[x].stats.station == "LB01" or "LB02" or "LB03" or "LB04" or "LB05" or "LB06" :
##            if st[x].stats.station == "LS06":
#                if week_start < st[x].stats.starttime.timestamp < week_end:
#                    t1 = st[x].stats.starttime
#                    t2 = st[x].stats.endtime
#                    peak,cf, bwid50, bwid25 = freq_info(st[x],t1 ,t2)
#                    
#                    sum_events +=1
#                    sum_peak += peak
#                    sum_cf += cf
#                    sum_bw += bwid50 
# 
#                if st[x].stats.starttime.timestamp > week_end:
#                    if sum_events > 1:
#                        
#                        av_peak = sum_peak/sum_events
#                        av_cf = sum_cf/sum_events
#                        av_bw = sum_bw/sum_events
#    #                    print('Average Dominant f = ',av_peak,' for ' ,sum_events,' events in week', week)
#    #                    print('Average central f = ',av_cf,' for ' ,sum_events,' events in week', week)
#    #                    print('Average bandwidth f = ',av_bw,' for ' ,sum_events,' events in week', week)
#    #                    print("\n")
#                        
#                        weekly_freq = np.lib.pad(weekly_freq, ((0,1),(0,0)), 'constant', constant_values=(0))
#                        
#                        weekly_freq[aw][0]=week + 1
#                        weekly_freq[aw][1]=av_peak
#                        weekly_freq[aw][2]=av_cf
#                        weekly_freq[aw][3]=av_bw     
#                        aw += 1
#   
#                    else:
#                        print("LS inactive in week", week+1)
#                    
#                    sum_peak=0
#                    sum_cf=0
#                    sum_bw=0
#                    sum_events = 0
#
#                    for p in range(0,120):
#
#                        week_start=week_end
#                        week_end += 7*24*60*60
#                        week += 1
#                        if week_start < st[x].stats.starttime.timestamp < week_end:
#
#                            t1 = st[x].stats.starttime
#                            t2 = st[x].stats.endtime
#                            peak,cf, bwid50, bwid25 = freq_info(st[x],t1 ,t2)
#                            
#                            sum_events += 1
#                            sum_peak += peak
#                            sum_cf += cf
#                            sum_bw += bwid50
#                            break
#                    
##print("weekly_peak=",weekly_freq)
#                 
#                
#np.savetxt("/Users/william/Documents/scanner/analysis/weekly_LS06_freq.csv", weekly_freq,delimiter=",",header="week,dominant,central,bandwidth_50")
#
#plt.figure(31)
#plt.plot(weekly_freq[:,0],weekly_freq[:,1],'rx-')
#ym1 = max(weekly_freq[:,1])
#plt.ylim((0,int(round(ym1 +1 ))))
#plt.xlim((0,140))
#plt.xlabel("Week") 
#plt.ylabel("Average Dominant Frequency")
#plt.title("Dominant Frequency of Explosions through Time")
#
#plt.figure(32)
#plt.plot(weekly_freq[:,0],weekly_freq[:,2],'rx-')
#ym2 = max(weekly_freq[:,2])
#plt.ylim((0,int(round(ym2 +1 )))) 
#plt.xlim((0,140))
#plt.xlabel("Week")
#plt.ylabel("Average Central Frequency") 
#plt.title("Central Frequency of Explosions through Time")
#
#plt.figure(33)
#plt.plot(weekly_freq[:,0],weekly_freq[:,3],'rx-')
#ym3 = max(weekly_freq[:,3])
#plt.ylim((0,int(round(ym3 +1 )))) 
#plt.xlim((0,140))    
#plt.xlabel("Week")   
#plt.ylabel("Average 50% Bandwidth") 
#plt.title("50% Bandwidth of Explosions through Time")    
##        
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
