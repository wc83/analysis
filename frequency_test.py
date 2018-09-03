#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 31 09:45:58 2018

@author: william
"""



from numpy import argmax
import numpy as np
from obspy.clients.earthworm import Client
from obspy import UTCDateTime
from obspy.signal.trigger import trigger_onset
from obspy import Stream
import matplotlib.pyplot as plt

## STATION, CHANNEL (DDF --> 400 Hz), NETWWORK AND LOCATION CODES 
#sta = 'LB01' # STATION 
#cha = 'HHZ' # CHANNEL
#net = 'Z4'  # 
#loc = ''    # location, it depends mostly of which network you are in. 
#
## Corner frequency for high-pass filter
##hp_corner = 0.05
#
## t1. and t2 are in hours:minutes:seconds
## Get data from (Liverpool Winston default) wave server between times t1 and t2 for all stations in stalist      
#client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
#t1 = UTCDateTime(2015, 1, 2, 13, 15, 0) #the format is year:day_of_the_year:month
#t2 = t1 + 1*60
#st = Stream()
#st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
##print(st)
#
##st.filter(type='bandpass',freqmin=0.1, freqmax=10)
#st.detrend(type='linear')
#st.detrend(type='demean')
##st.plot(color='b',starttime=t1, endtime=t2)
#
#tr=st[0].slice(starttime=t1, endtime=t2)
#
##frequency info
#T=1/tr.stats.sampling_rate
#window=t2-t1
#tr_data=tr.data
#m=np.mean(tr_data)
#tr_data = tr_data-m
#famp = abs(np.fft.rfft(tr_data))
#famp2=famp*famp
#X = np.fft.rfftfreq(tr_data.size, d=T)
#
#
#
#
##%% plot frequency amplitude
#
#plt.figure(1)
#plt.plot(X,famp)
#plt.xlabel('frequency [Hz]')
#plt.ylabel('|amplitude|')
#plt.title('Station '+tr.stats.station+', Channel '+tr.stats.channel)
#plt.xlim([0, 5])
##
#x=[1 ,1 ,2,2]
#y=[ 0,  4500000, 0, 4500000]
#def connectpoints(x,y,p1,p2):
#    x1, x2 = x[p1], x[p2]
#    y1, y2 = y[p1], y[p2]
#    plt.plot([x1,x2],[y1,y2],'k-')
#
#connectpoints(x,y,0,1)
#connectpoints(x,y,2,3)
#
##plot frequency amplitude squared
#plt.figure(2)
#plt.plot(X,famp2)
#plt.xlabel('frequency [Hz]')
#plt.ylabel('|energy|')
#plt.title('Station '+tr.stats.station+', Channel '+tr.stats.channel)
#plt.xlim([0, 5])
#
#x=[1 ,1 ,2,2]
#y=[ 0,  20000000000000, 0, 20000000000000]
#def connectpoints(x,y,p1,p2):
#    x1, x2 = x[p1], x[p2]
#    y1, y2 = y[p1], y[p2]
#    plt.plot([x1,x2],[y1,y2],'k-')
#
#connectpoints(x,y,0,1)
#connectpoints(x,y,2,3)
#
##%%
#
#
##peak,cf, bwid50, bwid25= freq_info(tr,t1,t2)
##print(peak,cf,bwid50)
#
##%%     percentage amplitudes
#
## one Hz
#one = int( window +1)
#num1=0
#for t in range(0,one,1):
#    n1=famp[t]                 
#    num1 += n1
#
#p1hz=num1/sum(famp)
#print('percentage of spectrum below 1Hz = ',p1hz)
#
## two Hz
#two = int((2*window)+1)
#num2=0
#for t in range(0,two,1):
#    n2=famp[t]                  
#    num2 += n2
#
#p2hz=num2/sum(famp)
# 
#print('percentage between 1 and 2Hz = ',p2hz-p1hz)
#print('percentage of spectrum above 2Hz = ',1-p2hz)
#print('percentage of spectrum below 2Hz = ',p2hz)
#print("\n")
#
##%%  percentage energies
#
#e_famp= famp2
#
#    # one Hz
#one = int( window +1)
#num1=0
#for t in range(0,one):
#    n1=e_famp[t]                 
#    num1 += n1
#
#p1hz=num1/sum(e_famp)
#print('percentage of energy below 1Hz = ',p1hz)
#
#
#    # two Hz
#two = int((2*window)+1)
#num2=0
#for t in range(0,two):
#    n2=e_famp[t]                  
#    num2 += n2
#
#p2hz=num2/sum(e_famp)
#print('percentage energy between 1 and 2Hz = ',p2hz-p1hz)
#print('percentage of energy above 2Hz = ',1-p2hz)
#print('percentage of energy below 2Hz = ',p2hz)
#
#        
#        
   
#%% #whole dataset %age 1Hz and 2Hz
import io
import obspy 

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
weekly_freq_p=np.zeros(shape=(0,6))
aw=0
       
with io.open("/Users/william/Documents/scanner/output_data/m32.mseed", "rb") as fh:
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
        
            
        # For each chunck of time, do analysis
#        print(st)
        for x in range (0,int(len(st))):
            if st[x].stats.station == "LB01":# or "LB02" or "LB03" or "LB04" or "LB05" or "LB06" :
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
                        weekly_freq_p[aw][5]=week_start
                        
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

#np.savetxt("/Users/william/Documents/scanner/analysis/weekly_LB_freq_amp_energy.csv", weekly_freq_p,delimiter=",",header="Week,AB1,AB2,EB1,EB2,time")


#
plt.figure(11)
plt.plot(weekly_freq_p[:,0],weekly_freq_p[:,2],'rx')
plt.ylim((0,1))
plt.xlim((0,160))
plt.xlabel("Week")
plt.ylabel("Average Percentage Frequency below 2Hz ")
plt.title("Frequency below 2Hz")

plt.figure(12)
plt.plot(weekly_freq_p[:,0],weekly_freq_p[:,4],'rx')
plt.ylim((0,1))
plt.xlim((0,160))
plt.xlabel("Week")
plt.ylabel("Average Percentage Energy below 2Hz ")
plt.title("Energy below 2Hz")




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
weekly_freq_p=np.zeros(shape=(0,6))
aw=0
       
with io.open("/Users/william/Documents/scanner/output_data/m32.mseed", "rb") as fh:
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
        
            
        # For each chunck of time, do analysis
#        print(st)
        for x in range (0,int(len(st))):
            if st[x].stats.station == "LB02": #or "LB02" or "LB03" or "LB04" or "LB05" or "LB06" :
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
                        weekly_freq_p[aw][5]=week_start
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

#np.savetxt("/Users/william/Documents/scanner/analysis/weekly_LB02_freq_amp_energy.csv", weekly_freq_p,delimiter=",",header="Week,AB1,AB2,EB1,EB2,time")
#
plt.figure(11)
plt.plot(weekly_freq_p[:,0],weekly_freq_p[:,2],'bx')
plt.figure(12)
plt.plot(weekly_freq_p[:,0],weekly_freq_p[:,4],'bx')

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
weekly_freq_p=np.zeros(shape=(0,6))
aw=0
       
with io.open("/Users/william/Documents/scanner/output_data/m32.mseed", "rb") as fh:
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
        
            
        # For each chunck of time, do analysis
#        print(st)
        for x in range (0,int(len(st))):
            if st[x].stats.station == "LB03": #or "LB02" or "LB03" or "LB04" or "LB05" or "LB06" :
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
                        weekly_freq_p[aw][5]=week_start
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

#np.savetxt("/Users/william/Documents/scanner/analysis/weekly_LB03_freq_amp_energy.csv", weekly_freq_p,delimiter=",",header="Week,AB1,AB2,EB1,EB2,time")

plt.figure(11)
plt.plot(weekly_freq_p[:,0],weekly_freq_p[:,2],'kx')
plt.figure(12)
plt.plot(weekly_freq_p[:,0],weekly_freq_p[:,4],'kx')



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
weekly_freq_p=np.zeros(shape=(0,6))
aw=0
       
with io.open("/Users/william/Documents/scanner/output_data/m32.mseed", "rb") as fh:
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
        
            
        # For each chunck of time, do analysis
#        print(st)
        for x in range (0,int(len(st))):
            if st[x].stats.station == "LB04": #or "LB02" or "LB03" or "LB04" or "LB05" or "LB06" :
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
                        weekly_freq_p[aw][5]=week_start
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

#np.savetxt("/Users/william/Documents/scanner/analysis/weekly_LB04_freq_amp_energy.csv", weekly_freq_p,delimiter=",",header="Week,AB1,AB2,EB1,EB2,time")
#
plt.figure(11)
plt.plot(weekly_freq_p[:,0],weekly_freq_p[:,2],'cx')
plt.figure(12)
plt.plot(weekly_freq_p[:,0],weekly_freq_p[:,4],'cx')



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
weekly_freq_p=np.zeros(shape=(0,6))
aw=0
       
with io.open("/Users/william/Documents/scanner/output_data/m32.mseed", "rb") as fh:
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
        
            
        # For each chunck of time, do analysis
#        print(st)
        for x in range (0,int(len(st))):
            if st[x].stats.station == "LB05": #or "LB02" or "LB03" or "LB04" or "LB05" or "LB06" :
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
                        weekly_freq_p[aw][5]=week_start
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

#np.savetxt("/Users/william/Documents/scanner/analysis/weekly_LB05_freq_amp_energy.csv", weekly_freq_p,delimiter=",",header="Week,AB1,AB2,EB1,EB2,time")
#
plt.figure(11)
plt.plot(weekly_freq_p[:,0],weekly_freq_p[:,2],'yx')
plt.figure(12)
plt.plot(weekly_freq_p[:,0],weekly_freq_p[:,4],'yx')


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
weekly_freq_p=np.zeros(shape=(0,6))
aw=0
       
with io.open("/Users/william/Documents/scanner/output_data/m32.mseed", "rb") as fh:
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
        
            
        # For each chunck of time, do analysis
#        print(st)
        for x in range (0,int(len(st))):
            if st[x].stats.station == "LB06": #or "LB02" or "LB03" or "LB04" or "LB05" or "LB06" :
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
                        weekly_freq_p[aw][5]=week_start
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

#np.savetxt("/Users/william/Documents/scanner/analysis/weekly_LB06_freq_amp_energy.csv", weekly_freq_p,delimiter=",",header="Week,AB1,AB2,EB1,EB2,time")
#
plt.figure(11)
plt.plot(weekly_freq_p[:,0],weekly_freq_p[:,2],'gx')
plt.legend(('LB01','LB02','LB03','LB04','LB05','LB06'))
plt.figure(12)
plt.plot(weekly_freq_p[:,0],weekly_freq_p[:,4],'gx')
plt.legend(('LB01','LB02','LB03','LB04','LB05','LB06'))












        
        
        