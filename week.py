#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 16:36:54 2018

@author: william
"""

import obspy
from obspy import UTCDateTime

week_start = 1419379200
week_end = 1419379200 + 7*24*60*60 

for x in range(0,130):
    print('week',x+1, 'from', UTCDateTime(week_start),'to',UTCDateTime(week_end))
    week_start += 7*24*60*60
    week_end += 7*24*60*60

