#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 17:47:54 2019

@author: william
"""

import obspy
from obspy import read
import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
from obspy import Stream
from obspy.signal.trigger import classic_sta_lta, recursive_sta_lta
from obspy.signal.trigger import plot_trigger, trigger_onset
from obspy import UTCDateTime
from obspy.clients.earthworm import Client
from mpl_toolkits.mplot3d import Axes3D



corr_m = genfromtxt("/Users/william/Documents/scanner/analysis/LB01_inf_correltion_matrix_6s.csv", delimiter=',')
Z = genfromtxt("/Users/william/Documents/scanner/analysis/LB01_inf_correltion_list_6s.csv", delimiter=',')

events=genfromtxt("/Users/william/Documents/scanner/analysis/LB01_inf_correltion_event_details_6s.csv", delimiter=',',skip_header=1)  
num_events = len(corr_m)

number_c = int(((num_events-1)**2)/2)+ int(num_events/2)

cor_list=np.zeros(shape=(number_c,1))

num=0
for x in range(0,num_events-1):
    for y in range(x+1,num_events):
        cor_list[num][0]=corr_m[x][y]
        num+=1


        
mean_cor =np.mean(cor_list)

std_cor = np.std(cor_list)


std_1_min = mean_cor - std_cor
std_1_max = mean_cor + std_cor

print("mean correlation = ", mean_cor)
print("1std from the mean = ", std_1_min, " to ", std_1_max)

plt.figure(1)
plt.hist(cor_list,bins=100, histtype='step')




#%% create correlation matrix
        
X=np.linspace(0,num_events-1,num_events)
Y=np.linspace(0,num_events-1,num_events)

X2, Y2 = np.meshgrid(X, Y)
Z2 = Z.reshape(len(X), len(Y))


#%% plot correlation matrix
plt.pcolor(X2, Y2, Z2)
plt.colorbar()
plt.show()






