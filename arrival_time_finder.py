#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 14:42:40 2018

@author: william
"""

from obspy.geodetics import locations2degrees
from obspy.taup import TauPyModel
import matplotlib.pyplot as plt


#32.718°S 67.123°W
# Location of EQ source
lat1=-32.7          #Latitude positive in north  
lon1=-67.1 #Longitude positive to East
depth=172
# Location of Santiaguito
lat2=14.74                      #Latitude positive in north
lon2=-91.566                    #Longitude positive to East


degree=locations2degrees(lat1, lon1, lat2, lon2)
print('degrees from source to station =',degree)

model = TauPyModel(model="iasp91")

arrivals = model.get_travel_times(source_depth_in_km=depth, distance_in_degree=degree,phase_list=["ttbasic"])
print(arrivals)  
