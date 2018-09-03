#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 13 16:24:31 2018

@author: william
"""

# need to produce stream for returned trace
from obspy import Stream

def calibrate(st):

    #%%   LB Seismic calibration factors
    
    # need to check LB01
    LB01c1=0.000001/750            # before 2015-12-05T00:00:01.000000Z
    LB01c2=(10*0.000001)/(750*256) # 2015-12-05T00:00:01 - 2016-06-15T00:00:01 + maybe till end
    LB01c3 = 0.000001/750
    LB02c=0.000001/750
    LB03c=0.000001/(750*256)
    LB04c=0.000001/(750*256)
    LB05c=0.000001/(750*256)
    LB06c=(10*0.000001)/750
    
    #%% LS Seismic calibration fators
    
    LSsc=0.000000122/800
    
    #%% Acoustic Calibration factors
    
    # need to check all - especially LB02,LB06
    LB01ac = 0.000001/0.0250
    LB02ac = 0.000001/(0.0250*256)
    LB03ac = 0.000001/(0.0250*256)
    LB04ac = 0.000001/(0.0250*256)
    LB05ac = 0.000001/(0.0250*256)
    LB06ac = 0.000001/(0.01*256)
    LS01ac = 0.000001/(0.0250*256)    
    
    
    #%% Alternative calibrations
    # LB Seismic calibrations
#    LB01c1=0.000001/750            # before 2015-12-05T00:00:01.000000Z
#    LB01c2=(10*0.0000041)/(750*256) # 2015-12-05T00:00:01 - 2016-06-15T00:00:01 + maybe till end
#    LB01c3 = 0.000001/750 
#    
#    LB02c=0.0000000268/750
#    LB03c=0.0000000268/(750*256)
#    LB04c=0.0000000268/(750*256)
#    LB05c=0.0000000268/(750*256)
#    LB06c=(10*0.0000041)/750
#    
    
    
    
    
    
    
    
    
    
    #%%
    st_c=Stream()
    
    #%% Seismic Calibrations
    
    if st.stats.channel == "HHZ":
        if st.stats.station == "LB01":
            if st.stats.starttime.timestamp < 1457013601.0:
                st.data=st.data * LB01c1
                st_c.append(st)
#                print("seismic lb01 start")
            if 1457013601.1 < st.stats.starttime.timestamp < 1465948799.9:
                st.data=st.data * LB01c2
                st_c.append(st)
#                print("seismic lb01 mid")
            if 1465948800.0 < st.stats.starttime.timestamp:
                st.data=st.data * LB01c3
                st_c.append(st)
#                print("seismic lb01 end")
        if st.stats.station == "LB02":
            st.data=st.data * LB02c
            st_c.append(st)
#            print("seismic lb02")
        if st.stats.station == "LB03":
            st.data=st.data * LB03c
            st_c.append(st)
#            print("seismic lb03")
        if st.stats.station == "LB04":
            st.data=st.data * LB04c
            st_c.append(st)
#            print("seismic lb04")
        if st.stats.station == "LB05":
            st.data=st.data * LB05c
            st_c.append(st)
#            print("seismic lb05")
        if st.stats.station == "LB06":
            st.data=st.data * LB06c
            st_c.append(st)
#            print("seismic lb06")
    if st.stats.channel == "EHZ":
        if st.stats.station == "LS01" or "LS02" or "LS03" or "LS04" or "LS05" or "LS06":
            st.data=st.data * LSsc
            st_c.append(st) 
#            print("seismic lS")
            
       #%% Acoustic Calibrations
    if st.stats.channel == "HDF" :
        if st.stats.station == "LB01":
            st.data=st.data * LB01ac
            st_c.append(st)
#            print("acoustic lb01")
        if st.stats.station == "LB02":
            st.data=st.data * LB02ac
            st_c.append(st)
#            print("acoustic lb02")
        if st.stats.station == "LB03":
            st.data=st.data * LB03ac
            st_c.append(st)
#            print("acoustic lb03")
        if st.stats.station == "LB04":
            st.data=st.data * LB04ac
            st_c.append(st)
#            print("acoustic lb04")
        if st.stats.station == "LB05":
            st.data=st.data * LB05ac
            st_c.append(st)
#            print("acoustic lb05")
        if st.stats.station == "LB06":
            st.data=st.data * LB06ac
            st_c.append(st)
#            print("acoustic lb06")
        if st.stats.station == "LS01":
            st.data=st.data * LS01ac
            st_c.append(st)
#            print("acoustic ls01")
#    #
    return(st_c)