#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 13 16:24:31 2018

@author: william
"""

# need to produce stream for returned trace
from obspy import Stream

def calibrate1(st1):

    #%%   LB Seismic calibration factors
    
    # need to check LB01
    LB01c1=1.33e-9           # before 1449273600.0 2015/12/05
    LB01c2=1.33e-9/256 # 2015-12-05T00:00:01 - 2016-06-16T00:00:01 + maybe till end
    LB01c3 = 3.25e-9         # After 2016-06-15T00:00:01
    LB01c4= 1.33e-9/256
    
    LB02c1=1.33e-9   # 0.000001/(750*256) ????
    LB02c2=1.33e-9/256
    
    LB03c=1.33e-9/256
    LB03c2=1.33e-9
    LB04c=1.33e-9/256
    LB05c=1.33e-9/256
    LB06c=3.25e-9
#    
#LB01, 19.11.14 - 5.12.15: 1.33e-9
#LB01, 5.12.15 - 16.3.16: 1.33e-9/256
#LB01, 16.3.16 - 15.6.16: 3.25e-9 (I assume because of the configuration- there were no suitable EQs to check)
#LB01, 15.6.16 - end: 1.33e-9/256
#
#LB02, 22.11.14 - data gap in March 15: 1.33e-9
#LB02, after data gap in March: 1.33e-9/256
#
#LB03 - LB05: 1.33e-9/256
#LB06 - LB07: 3.25e-9 (only checked for LB06, because again no suitable EQs at LB07)
#
#LS02, LS05 - LS08: 1.526e-9 (this would mean a BOB with factor 2.5)
#LS01, LS03, LS04: should be the same as the other LS stations, but amplitudes seem too large


    
    #%% LS Seismic calibration fators
    
    LSsc=1.526e-9
    
    
    #%% Acoustic Calibration factors
    
    LB01ac = 0.000001/0.0250
    LB02ac = 0.000001/(0.0250*256)
    LB03ac = 0.000001/(0.0250*256)
    LB04ac = 0.000001/(0.0250*256)
    LB05ac = 0.000001/(0.0250*256)
    LB06ac = 0.000001/(0.01*256)
    LS01ac = 0.000001/(0.0250*256)    
    

    #%%
    st_c=Stream()
    
    #%% Seismic Calibrations
#    print(st1.stats.station)
    if st1.stats.channel == "HHZ":
        if st1.stats.station == "LB01":
            if st1.stats.starttime.timestamp < 1449273600.0:
                st1.data=st1.data * LB01c1
                st_c.append(st1)
#                print("seismic lb01 start")
            if 1449273600.01 < st1.stats.starttime.timestamp < 1458086400.0:
                st1.data=st1.data * LB01c2
                st_c.append(st1)
#                print("seismic lb01 mid")
            if 1458086400.01 < st1.stats.starttime.timestamp < 1465948800.0:
                st1.data=st1.data * LB01c3
                st_c.append(st1)
#                print("seismic lb01 mid2")
            if 1465948800.01 < st1.stats.starttime.timestamp:
                st1.data=st1.data * LB01c4
                st_c.append(st1)
#                print("seismic lb01 end")
        if st1.stats.station == "LB02":
            if st1.stats.starttime.timestamp < 1428000000.0:
                st1.data=st1.data * LB02c1
                st_c.append(st1)
#                print("seismic lb02 start")
            if st1.stats.starttime.timestamp > 1428000000.01:
                st1.data=st1.data * LB02c2
                st_c.append(st1)
#                print("seismic lb02 end")
        if st1.stats.station == "LB03" or st1.stats.station == "STG8":
            if st1.stats.starttime.timestamp < 1514764800.0:
                st1.data=st1.data * LB03c
                st_c.append(st1)
#                print("seismic lb03")
            if st1.stats.starttime.timestamp > 1514764800.01:
                st1.data=st1.data * LB03c2
                st_c.append(st1)
#                print("seismic lSTG8")
        if st1.stats.station == "LB04":
            st1.data=st1.data * LB04c
            st_c.append(st1)
#            print("seismic lb04")
        if st1.stats.station == "LB05":
            st1.data=st1.data * LB05c
            st_c.append(st1)
#            print("seismic lb05")
        if st1.stats.station == "LB06":
            st1.data=st1.data * LB06c
            st_c.append(st1)
#            print("seismic lb06")
    if st1.stats.channel == "EHZ":
        if st1.stats.station == "LS01" or st1.stats.station =="LS02" or st1.stats.station =="LS03" or st1.stats.station =="LS04" or st1.stats.station =="LS05" or st1.stats.station =="LS06":
            st1.data=st1.data * LSsc
            st_c.append(st1) 
#            print("seismic lS")
            
       #%% Acoustic Calibrations
    if st1.stats.channel == "HDF" :
        if st1.stats.station == "LB01":
            st1.data=st1.data * LB01ac
            st_c.append(st1)
#            print("acoustic lb01")
        if st1.stats.station == "LB02":
            st1.data=st1.data * LB02ac
            st_c.append(st1)
#            print("acoustic lb02")
        if st1.stats.station == "LB03":
            st1.data=st1.data * LB03ac
            st_c.append(st1)
#            print("acoustic lb03")
        if st1.stats.station == "LB04":
            st1.data=st1.data * LB04ac
            st_c.append(st1)
#            print("acoustic lb04")
        if st1.stats.station == "LB05":
            st1.data=st1.data * LB05ac
            st_c.append(st1)
#            print("acoustic lb05")
        if st1.stats.station == "LB06":
            st1.data=st1.data * LB06ac
            st_c.append(st1)
#            print("acoustic lb06")
        if st1.stats.station == "LS01":
            st1.data=st1.data * LS01ac
            st_c.append(st1)
#            print("acoustic ls01")
#    #
    return(st_c)