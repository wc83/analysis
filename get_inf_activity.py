#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 12:28:24 2018

@author: william
"""



def get_inf_activity(day):

    from obspy.clients.earthworm import Client
    from obspy import UTCDateTime
    from obspy import Stream
    
    year1=2014
    month1=11
    day1=24
    hour1=0
    minute1=0
    second1=0
    
    num=12
    lb_num=6
    ls_num=6
#%% LB01    
    try:  
        
        sta = 'LB01' # STATION LB01
        cha = 'HDF' # CHANNEL - Vertical
        net = 'Z4'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 

        client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
        t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
        t1 = t0 + day*24*60*60
        t2 = t1 + 23*60*60 + 59*60 +59.999 #UTCDateTime(year2, month2, day2, hour2, minute2, second2) # notice we have here 10 minutes, but we can select our times. 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, '', cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        sorted_data = st1.copy()
        sorted_data= abs(sorted_data[0].data)
        sorted_data.sort()
        mid_dat= sorted_data[int(len(sorted_data)/2)]
        
        LB01=1
        
        if sum(abs(st1[0].data)) < 10 or st1[0].stats.npts < 7920000 or mid_dat < 0.1 :
    
            num = num -1
            lb_num = lb_num - 1
            LB01=0
            
    except: # give 2 seconds of data instead
        
        num = num -1
        lb_num = lb_num - 1
        LB01=0
#%%    LB02

    try:  
     
        sta = 'LB02' # STATION LB01
        cha = 'HDF' # CHANNEL - Vertical
        net = 'Z4'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 
        
        client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
        t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
        t1 = t0 + day*24*60*60
        t2 = t1 + 23*60*60 + 59*60 +59.999 #UTCDateTime(year2, month2, day2, hour2, minute2, second2) # notice we have here 10 minutes, but we can select our times. 
        st2 = Stream()
        st2 = client.get_waveforms(net, sta, '', cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st2.detrend(type='linear')
        st2.detrend(type='demean')
        break_test=st2
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)

        sorted_data = st2.copy()
        sorted_data= abs(sorted_data[0].data)
        sorted_data.sort()
        mid_dat= sorted_data[int(len(sorted_data)/2)]
        
        LB02=1
        
        if sum(abs(st2[0].data)) < 10 or st2[0].stats.npts < 7920000 or mid_dat < 0.1 :

            num = num -1
            lb_num = lb_num - 1
            LB02=0
    except: # give 2 seconds of data instead


        num = num -1
        lb_num = lb_num - 1
        LB02=0
#%% LB03    
    try:  
     
        sta = 'LB03' # STATION LB01
        cha = 'HDF' # CHANNEL - Vertical
        net = 'Z4'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 
               
        client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
        t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
        t1 = t0 + day*24*60*60
        t2 = t1 + 23*60*60 + 59*60 +59.999 #UTCDateTime(year2, month2, day2, hour2, minute2, second2) # notice we have here 10 minutes, but we can select our times. 
        st3 = Stream()
        st3 = client.get_waveforms(net, sta, '', cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st3.detrend(type='linear')
        st3.detrend(type='demean')
        break_test=st3
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
                
        sorted_data = st3.copy()
        sorted_data= abs(sorted_data[0].data)
        sorted_data.sort()
        mid_dat= sorted_data[int(len(sorted_data)/2)]
        
        LB03=1
        
        if sum(abs(st3[0].data)) < 10 or st3[0].stats.npts < 7920000 or mid_dat < 0.1 :
        
            num = num -1
            lb_num = lb_num - 1
            LB03=0
    except: # give 2 seconds of data instead

        num = num -1
        lb_num = lb_num - 1
        LB03=0
#%% LB04    
    try:  
     
        sta = 'LB04' # STATION LB01
        cha = 'HDF' # CHANNEL - Vertical
        net = 'Z4'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 
               
        client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
        t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
        t1 = t0 + day*24*60*60
        t2 = t1 + 23*60*60 + 59*60 +59.999 #UTCDateTime(year2, month2, day2, hour2, minute2, second2) # notice we have here 10 minutes, but we can select our times. 
        st4 = Stream()
        st4 = client.get_waveforms(net, sta, '', cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st4.detrend(type='linear')
        st4.detrend(type='demean')
        break_test=st4
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)

        sorted_data = st4.copy()
        sorted_data= abs(sorted_data[0].data)
        sorted_data.sort()
        mid_dat= sorted_data[int(len(sorted_data)/2)]
        
        LB04=1
        
        if sum(abs(st4[0].data)) < 10 or st4[0].stats.npts < 7920000 or mid_dat < 0.1 :
   
                num = num -1
                lb_num = lb_num - 1
                LB04=0
    except: # give 2 seconds of data instead
  
        num = num -1
        lb_num = lb_num - 1
        LB04=0
#%% LB05    
    try:  
     
        sta = 'LB05' # STATION LB01
        cha = 'HDF' # CHANNEL - Vertical
        net = 'Z4'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 
               
        client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
        t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
        t1 = t0 + day*24*60*60
        t2 = t1 + 23*60*60 + 59*60 +59.999 #UTCDateTime(year2, month2, day2, hour2, minute2, second2) # notice we have here 10 minutes, but we can select our times. 
        st5 = Stream()
        st5 = client.get_waveforms(net, sta, '', cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st5.detrend(type='linear')
        st5.detrend(type='demean')
        break_test=st5
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        sorted_data = st5.copy()
        sorted_data= abs(sorted_data[0].data)
        sorted_data.sort()
        mid_dat= sorted_data[int(len(sorted_data)/2)]
        
        LB05=1
        
        if sum(abs(st5[0].data)) < 10 or st5[0].stats.npts < 7920000 or mid_dat < 0.1 :
  
                num = num -1
                lb_num = lb_num - 1
                LB05=0
    except: # give 2 seconds of data instead
        
        num = num -1
        lb_num = lb_num - 1
        LB05=0
#%% LB06    
    try:  
     
        sta = 'LB06' # STATION LB01
        cha = 'HDF' # CHANNEL - Vertical
        net = 'Z4'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 
               
        client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
        t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
        t1 = t0 + day*24*60*60
        t2 = t1 + 23*60*60 + 59*60 +59.999 #UTCDateTime(year2, month2, day2, hour2, minute2, second2) # notice we have here 10 minutes, but we can select our times. 
        st6 = Stream()
        st6 = client.get_waveforms(net, sta, '', cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st6.detrend(type='linear')
        st6.detrend(type='demean')
        break_test=st6
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        sorted_data = st6.copy()
        sorted_data= abs(sorted_data[0].data)
        sorted_data.sort()
        mid_dat= sorted_data[int(len(sorted_data)/2)]
        LB06=1
        if sum(abs(st6[0].data)) < 10 or st6[0].stats.npts < 7920000 or mid_dat < 0.1 :
            
                num = num -1
                lb_num = lb_num - 1
                LB06=0
    except: # give 2 seconds of blank data instead
       
        num = num -1
        lb_num = lb_num - 1
        LB06=0
#%% LS01    
    try:  
        
        sta = 'LS01' # STATION LS01
        cha = 'HDF' # CHANNEL - Vertical
        net = 'Z4'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 

        client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
        t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
        t1 = t0 + day*24*60*60
        t2 = t1 + 23*60*60 + 59*60 +59.999 #UTCDateTime(year2, month2, day2, hour2, minute2, second2) # notice we have here 10 minutes, but we can select our times. 
        sts1 = Stream()
        sts1 = client.get_waveforms(net, sta, '', cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        sts1.detrend(type='linear')
        sts1.detrend(type='demean')
        break_test=sts1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        sorted_data = sts1.copy()
        sorted_data= abs(sorted_data[0].data)
        sorted_data.sort()
        mid_dat= sorted_data[int(len(sorted_data)/2)]
        LS01=1
        if sum(abs(sts1[0].data)) < 10 or sts1[0].stats.npts < 7920000 or mid_dat < 0.1 :
    
            num = num -1
            ls_num = ls_num-1
            LS01=0
    except: # give 2 seconds of data instead
        
        num = num -1
        ls_num = ls_num-1
        LS01=0
#%%    LS02

    try:  
     
        sta = 'LS02' # STATION LS02
        cha = 'HDF' # CHANNEL - Vertical
        net = 'Z4'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 
        
        client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
        t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
        t1 = t0 + day*24*60*60
        t2 = t1 + 23*60*60 + 59*60 +59.999 #UTCDateTime(year2, month2, day2, hour2, minute2, second2) # notice we have here 10 minutes, but we can select our times. 
        sts2 = Stream()
        sts2 = client.get_waveforms(net, sta, '', cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        sts2.detrend(type='linear')
        sts2.detrend(type='demean')
        break_test=sts2
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)

        sorted_data = sts2.copy()
        sorted_data= abs(sorted_data[0].data)
        sorted_data.sort()
        mid_dat= sorted_data[int(len(sorted_data)/2)]
        LS02=1
        if sum(abs(sts2[0].data)) < 10 or sts2[0].stats.npts < 7920000 or mid_dat < 0.1 :
         
            num = num -1
            ls_num = ls_num-1
            LS02=0
    except: # give 2 seconds of data instead

        num = num -1
        ls_num = ls_num-1
        LS02=0
#%% LS03    
    try:  
     
        sta = 'LS03' # STATION LS03
        cha = 'HDF' # CHANNEL - Vertical
        net = 'Z4'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 
               
        client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
        t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
        t1 = t0 + day*24*60*60
        t2 = t1 + 23*60*60 + 59*60 +59.999 #UTCDateTime(year2, month2, day2, hour2, minute2, second2) # notice we have here 10 minutes, but we can select our times. 
        sts3 = Stream()
        sts3 = client.get_waveforms(net, sta, '', cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        sts3.detrend(type='linear')
        sts3.detrend(type='demean')
        break_test=sts3
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
                
        sorted_data = sts3.copy()
        sorted_data= abs(sorted_data[0].data)
        sorted_data.sort()
        mid_dat= sorted_data[int(len(sorted_data)/2)]
        LS03=1
        
        if sum(abs(sts3[0].data)) < 10 or sts3[0].stats.npts < 7920000 or mid_dat < 0.1 or sts3[0].stats.starttime.timestamp >  1442707190: 
        
            num = num -1
            ls_num = ls_num-1
            LS03=0
    except: # give 2 seconds of data instead

        num = num -1
        ls_num = ls_num-1
        LS03=0
#%% LS04    
    try:  
     
        sta = 'LS04' # STATION LS04
        cha = 'HDF' # CHANNEL - Vertical
        net = 'Z4'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 
               
        client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
        t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
        t1 = t0 + day*24*60*60
        t2 = t1 + 23*60*60 + 59*60 +59.999 #UTCDateTime(year2, month2, day2, hour2, minute2, second2) # notice we have here 10 minutes, but we can select our times. 
        sts4 = Stream()
        sts4 = client.get_waveforms(net, sta, '', cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        sts4.detrend(type='linear')
        sts4.detrend(type='demean')
        break_test=sts4
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)

        sorted_data = sts4.copy()
        sorted_data= abs(sorted_data[0].data)
        sorted_data.sort()
        mid_dat= sorted_data[int(len(sorted_data)/2)]
        LS04=1
        if sum(abs(sts4[0].data)) < 10 or sts4[0].stats.npts < 7920000 or mid_dat < 0.1 :
            
                num = num -1
                ls_num = ls_num-1
                LS04=0
    except: # give 2 seconds of data instead
           
        num = num -1
        ls_num = ls_num-1
        LS04=0
#%% LS05    
    try:  
     
        sta = 'LS05' # STATION LB0S
        cha = 'HDF' # CHANNEL - Vertical
        net = 'Z4'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 
               
        client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
        t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
        t1 = t0 + day*24*60*60
        t2 = t1 + 23*60*60 + 59*60 +59.999 #UTCDateTime(year2, month2, day2, hour2, minute2, second2) # notice we have here 10 minutes, but we can select our times. 
        sts5 = Stream()
        sts5 = client.get_waveforms(net, sta, '', cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        sts5.detrend(type='linear')
        sts5.detrend(type='demean')
        break_test=sts5
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        sorted_data = sts5.copy()
        sorted_data= abs(sorted_data[0].data)
        sorted_data.sort()
        mid_dat= sorted_data[int(len(sorted_data)/2)]
        LS05=1
        if sum(abs(sts5[0].data)) < 10 or sts5[0].stats.npts < 7920000 or mid_dat < 0.1 :
            
                num = num -1
                ls_num = ls_num-1
                LS05=0
    except: # give 2 seconds of data instead
        
        num = num -1
        ls_num = ls_num-1
        LS05=0
#%% LS06    
    try:  
     
        sta = 'LS06' # STATION LS06
        cha = 'HDF' # CHANNEL - Vertical
        net = 'Z4'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 
               
        client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
        t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
        t1 = t0 + day*24*60*60
        t2 = t1 + 23*60*60 + 59*60 +59.999 #UTCDateTime(year2, month2, day2, hour2, minute2, second2) # notice we have here 10 minutes, but we can select our times. 
        sts6 = Stream()
        sts6 = client.get_waveforms(net, sta, '', cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        sts6.detrend(type='linear')
        sts6.detrend(type='demean')
        break_test=sts6
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        sorted_data = sts6.copy()
        sorted_data= abs(sorted_data[0].data)
        sorted_data.sort()
        mid_dat= sorted_data[int(len(sorted_data)/2)]
        LS06=1
        if sum(abs(sts6[0].data)) < 10 or sts6[0].stats.npts < 7920000 or mid_dat < 0.1 :
            
                num = num -1
                ls_num = ls_num-1
                LS06=0
    except: # give 2 seconds of blank data instead
        
        num = num -1
        ls_num = ls_num-1
        LS06=0
    #%% return all stations
    
    return(num,lb_num,ls_num,LB01,LB02,LB03,LB04,LB05,LB06,LS01,LS02,LS03,LS04,LS05,LS06) #st7