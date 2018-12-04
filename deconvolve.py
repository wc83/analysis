#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 16:31:55 2018

@author: william
"""

from numpy import argmax
import numpy as np
from obspy.clients.earthworm import Client
from obspy import UTCDateTime
from obspy.signal.trigger import trigger_onset
from obspy import Stream
import matplotlib.pyplot as plt
from obspy import read, read_inventory


import io
import obspy 

reclen = 512
chunksize = 100000 * reclen # Around 50 MB

#with io.open("/Users/william/Documents/scanner/output_data/m30.mseed", "rb") as fh:
        # just month 2
with io.open("/Users/william/Documents/scanner/output_data/EXP_all_data_stream_2_month_1.mseed", "rb") as fh:
    while True:
        with io.BytesIO() as buf:
            c = fh.read(chunksize)
            if not c:
                break
            buf.write(c)
            buf.seek(0, 0)
            st = obspy.read(buf)
           
for x in range(0,10):
    tr=st[x]
    inv = read_inventory()
    tr.remove_response(inventory=inv, output="VEL",water_level=60, plot=True)     
   

            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            