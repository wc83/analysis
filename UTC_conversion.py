#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 13 13:46:19 2018

@author: william
"""

from obspy.clients.earthworm import Client
from obspy import UTCDateTime
from obspy import Stream
from numpy import genfromtxt
from obspy.signal.trigger import recursive_sta_lta
from obspy.signal.trigger import plot_trigger, trigger_onset
import numpy as np




cat= genfromtxt("/Users/william/Documents/scanner/all_stations/Explosion_Catalogue_V3.csv", delimiter=',',skip_header=1)

UTC=np.zeros(shape=(0,1))

for x in range(0,len(cat)):
    ts=cat[x,0]
    utc=UTCDateTime(ts)
    UTC = np.lib.pad(UTC, ((0,1),(0,0)), 'constant', constant_values=(0))
                    
    UTC[x][0]=utc













