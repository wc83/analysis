#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 12:13:35 2019

@author: william
"""

def event_mag(t1,a):
    
    import numpy as np
    from obspy.clients.earthworm import Client
    from obspy import UTCDateTime
    from obspy import Stream
    #%%
    
    magnitudes= np.zeros(shape=(a,1))
    
    #NEED TO BE CAREFUL FOR INACTIVE STATIONS
    sta = 'LB01' # STATION 
    cha = 'HHZ' # CHANNEL
    net = 'Z4'  # 
    
    client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
    #t1 = UTCDateTime(2014, 12, 3, 2, 43, 0) #the format is year:day_of_the_year:month
    t2 = t1 + 60*2
    st = Stream()
    st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
    
    st.detrend(type='linear')
    st.detrend(type='demean')
    st.filter(type='bandpass',freqmin=0.25, freqmax=10)
    tr=st[0].slice(starttime=t1-20, endtime=t2)
    st_c = calibrate1(tr)
            
    peak2peak=max(st_c[0].data)- min(st_c[0].data)
    magnitudes[0]=peak2peak
    
    #%%    
    #NEED TO BE CAREFUL FOR INACTIVE STATIONS
    sta = 'LB02' # STATION 
    cha = 'HHZ' # CHANNEL
    net = 'Z4'  # 
    
    client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
    #    t1 = UTCDateTime(2014, 12, 3, 2, 43, 0) #the format is year:day_of_the_year:month
    t2 = t1 + 60*2
    st = Stream()
    st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
    
    st.detrend(type='linear')
    st.detrend(type='demean')
    st.filter(type='bandpass',freqmin=0.25, freqmax=10)
    tr=st[0].slice(starttime=t1-20, endtime=t2)
    st_c = calibrate1(tr)
            
    peak2peak=max(st_c[0].data)- min(st_c[0].data)
    magnitudes[1]=peak2peak
    
    #%%    
    #NEED TO BE CAREFUL FOR INACTIVE STATIONS
    sta = 'LB03' # STATION 
    cha = 'HHZ' # CHANNEL
    net = 'Z4'  # 
    
    client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
    #    t1 = UTCDateTime(2014, 12, 3, 2, 43, 0) #the format is year:day_of_the_year:month
    t2 = t1 + 60*2
    st = Stream()
    st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
    
    st.detrend(type='linear')
    st.detrend(type='demean')
    st.filter(type='bandpass',freqmin=0.25, freqmax=10)
    tr=st[0].slice(starttime=t1-20, endtime=t2)
    st_c = calibrate1(tr)
            
    peak2peak=max(st_c[0].data)- min(st_c[0].data)
    magnitudes[2]=peak2peak
    
    #%%    
#    NEED TO BE CAREFUL FOR INACTIVE STATIONS
    sta = 'LB04' # STATION 
    cha = 'HHZ' # CHANNEL
    net = 'Z4'  # 
    
    client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
    #    t1 = UTCDateTime(2014, 12, 3, 2, 43, 0) #the format is year:day_of_the_year:month
    t2 = t1 + 60*2
    st = Stream()
    st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
    
    st.detrend(type='linear')
    st.detrend(type='demean')
    st.filter(type='bandpass',freqmin=0.25, freqmax=10)
    tr=st[0].slice(starttime=t1-20, endtime=t2)
    st_c = calibrate1(tr)
            
    peak2peak=max(st_c[0].data)- min(st_c[0].data)
    magnitudes[3]=peak2peak
    
    #%%    
    #NEED TO BE CAREFUL FOR INACTIVE STATIONS
    sta = 'LB05' # STATION 
    cha = 'HHZ' # CHANNEL
    net = 'Z4'  # 
    
    client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
    #    t1 = UTCDateTime(2014, 12, 3, 2, 43, 0) #the format is year:day_of_the_year:month
    t2 = t1 + 60*2
    st = Stream()
    st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
    
    st.detrend(type='linear')
    st.detrend(type='demean')
    st.filter(type='bandpass',freqmin=0.25, freqmax=10)
    tr=st[0].slice(starttime=t1-20, endtime=t2)
    st_c = calibrate1(tr)
            
    peak2peak=max(st_c[0].data)- min(st_c[0].data)
    magnitudes[4]=peak2peak
    
        #%%    
#    #NEED TO BE CAREFUL FOR INACTIVE STATIONS
#    sta = 'LB06' # STATION 
#    cha = 'HHZ' # CHANNEL
#    net = 'Z4'  # 
#    
#    client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
#    #    t1 = UTCDateTime(2014, 12, 3, 2, 43, 0) #the format is year:day_of_the_year:month
#    t2 = t1 + 60*2
#    st = Stream()
#    st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
#    
#    st.detrend(type='linear')
#    st.detrend(type='demean')
#    st.filter(type='bandpass',freqmin=0.25, freqmax=10)
#    tr=st[0].slice(starttime=t1-20, endtime=t2)
#    st_c = calibrate1(tr)
#            
#    peak2peak=max(st_c[0].data)- min(st_c[0].data)
#    magnitudes[5]=peak2peak
    
    #%%    
    #NEED TO BE CAREFUL FOR INACTIVE STATIONS
    sta = 'LS01' # STATION 
    cha = 'EHZ' # CHANNEL
    net = 'Z4'  # 
    
    client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
    #    t1 = UTCDateTime(2014, 12, 3, 2, 43, 0) #the format is year:day_of_the_year:month
    t2 = t1 + 60*2
    st = Stream()
    st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
    
    st.detrend(type='linear')
    st.detrend(type='demean')
    st.filter(type='bandpass',freqmin=0.25, freqmax=10)
    tr=st[0].slice(starttime=t1-20, endtime=t2)
    st_c = calibrate1(tr)
            
    peak2peak=max(st_c[0].data)- min(st_c[0].data)
    magnitudes[5]=peak2peak
    
    #%%    
    #NEED TO BE CAREFUL FOR INACTIVE STATIONS
    sta = 'LS02' # STATION 
    cha = 'EHZ' # CHANNEL
    net = 'Z4'  # 
    
    client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
    #    t1 = UTCDateTime(2014, 12, 3, 2, 43, 0) #the format is year:day_of_the_year:month
    t2 = t1 + 60*2
    st = Stream()
    st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
    
    st.detrend(type='linear')
    st.detrend(type='demean')
    st.filter(type='bandpass',freqmin=0.25, freqmax=10)
    tr=st[0].slice(starttime=t1-20, endtime=t2)
    st_c = calibrate1(tr)
            
    peak2peak=max(st_c[0].data)- min(st_c[0].data)
    magnitudes[6]=peak2peak
    
    #%%    
    #NEED TO BE CAREFUL FOR INACTIVE STATIONS
    sta = 'LS03' # STATION 
    cha = 'EHZ' # CHANNEL
    net = 'Z4'  # 
    
    client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
    #    t1 = UTCDateTime(2014, 12, 3, 2, 43, 0) #the format is year:day_of_the_year:month
    t2 = t1 + 60*2
    st = Stream()
    st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
    
    st.detrend(type='linear')
    st.detrend(type='demean')
    st.filter(type='bandpass',freqmin=0.25, freqmax=10)
    tr=st[0].slice(starttime=t1-20, endtime=t2)
    st_c = calibrate1(tr)
            
    peak2peak=max(st_c[0].data)- min(st_c[0].data)
    magnitudes[7]=peak2peak
    
    #%%    
    #NEED TO BE CAREFUL FOR INACTIVE STATIONS
    sta = 'LS04' # STATION 
    cha = 'EHZ' # CHANNEL
    net = 'Z4'  # 
    
    client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
    #    t1 = UTCDateTime(2014, 12, 3, 2, 43, 0) #the format is year:day_of_the_year:month
    t2 = t1 + 60*2
    st = Stream()
    st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
    
    st.detrend(type='linear')
    st.detrend(type='demean')
    st.filter(type='bandpass',freqmin=0.25, freqmax=10)
    tr=st[0].slice(starttime=t1-20, endtime=t2)
    st_c = calibrate1(tr)
            
    peak2peak=max(st_c[0].data)- min(st_c[0].data)
    magnitudes[8]=peak2peak
    
    #%%
    #NEED TO BE CAREFUL FOR INACTIVE STATIONS
    sta = 'LS05' # STATION 
    cha = 'EHZ' # CHANNEL
    net = 'Z4'  # 
    
    client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
    #    t1 = UTCDateTime(2014, 12, 3, 2, 43, 0) #the format is year:day_of_the_year:month
    t2 = t1 + 60*2
    st = Stream()
    st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
    
    st.detrend(type='linear')
    st.detrend(type='demean')
    st.filter(type='bandpass',freqmin=0.25, freqmax=10)
    tr=st[0].slice(starttime=t1-20, endtime=t2)
    st_c = calibrate1(tr)
            
    peak2peak=max(st_c[0].data)- min(st_c[0].data)
    magnitudes[9]=peak2peak
    
      #%%
    #NEED TO BE CAREFUL FOR INACTIVE STATIONS
    sta = 'LS06' # STATION 
    cha = 'EHZ' # CHANNEL
    net = 'Z4'  # 
    
    client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
    #    t1 = UTCDateTime(2014, 12, 3, 2, 43, 0) #the format is year:day_of_the_year:month
    t2 = t1 + 60*2
    st = Stream()
    st = client.get_waveforms(net, sta, '', cha, t1-20 , t2)
    
    st.detrend(type='linear')
    st.detrend(type='demean')
    st.filter(type='bandpass',freqmin=0.25, freqmax=10)
    tr=st[0].slice(starttime=t1-20, endtime=t2)
    st_c = calibrate1(tr)
            
    peak2peak=max(st_c[0].data)- min(st_c[0].data)
    magnitudes[10]=peak2peak
      
    




    return(magnitudes)