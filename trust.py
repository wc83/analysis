#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 12:06:45 2018

@author: william
"""

import obspy
from obspy import read
from obspy import Stream
from numpy import genfromtxt
import numpy as np

cat= genfromtxt("/Users/william/Documents/scanner/all_stations/Explosion_Catalogue_V3.csv", delimiter=',',skip_header=1)


LB01=cat[:,8]
np.transpose(LB01)
LB02=cat[:,9]
np.transpose(LB02)
LB03=cat[:,10]
np.transpose(LB03)
LB04=cat[:,11]
np.transpose(LB04)
LB05=cat[:,12]
np.transpose(LB05)
LB06=cat[:,13]
np.transpose(LB06)
LS01=cat[:,14]
np.transpose(LS01)
LS02=cat[:,15]
np.transpose(LS02)
LS03=cat[:,16]
np.transpose(LS03)
LS04=cat[:,17]
np.transpose(LS04)
LS05=cat[:,18]
np.transpose(LS05)
LS06=cat[:,19]
np.transpose(LS06)

T=np.zeros(shape=(len(LB01),1))
T2=np.zeros(shape=(len(LB01),1))
for x in range(0,len(LB01)):
    
    if LB01[x]+LB02[x]+LB03[x]+LB04[x]+LB05[x]+LB06[x] > 1:
        T2[x]=4
    if (LB01[x]+LB02[x]+LB03[x]+LB04[x]+LB05[x]+LB06[x]) == 1 and (LS01[x]+LS02[x]+LS03[x]+LS04[x]+LS05[x]+LS06[x] > 0):
        T2[x]=3
    if (LB01[x]+LB02[x]+LB03[x]+LB04[x]+LB05[x]+LB06[x] == 0) and (LS01[x]+LS02[x]+LS03[x]+LS04[x]+LS05[x]+LS06[x] > 1):
        T2[x]=2    
    if (LB01[x]+LB02[x]+LB03[x]+LB04[x]+LB05[x]+LB06[x] == 1) and (LS01[x]+LS02[x]+LS03[x]+LS04[x]+LS05[x]+LS06[x] == 0):
        T2[x]=2  
    if (LB01[x]+LB02[x]+LB03[x]+LB04[x]+LB05[x]+LB06[x] == 0) and (LS01[x]+LS02[x]+LS03[x]+LS04[x]+LS05[x]+LS06[x] == 1):
        T2[x]=1  
    

    
  #    if LB01[x]+LB02[x]+LB03[x]+LB04[x]+LB05[x] > 1:
#        T[x]=5
#    if (LB01[x]+LB02[x]+LB03[x]+LB04[x]+LB05[x] == 1) and (LB06[x] == 1):
#        T[x]=4   
#    if (LB01[x]+LB02[x]+LB03[x]+LB04[x]+LB05[x] == 1) and (LB06[x] == 0) and (LS01[x]+LS02[x]+LS03[x]+LS04[x]+LS05[x]+LS06[x] > 0):
#        T[x]=4
#    if (LB01[x]+LB02[x]+LB03[x]+LB04[x]+LB05[x] == 1) and (LB06[x] == 0) and (LS01[x]+LS02[x]+LS03[x]+LS04[x]+LS05[x]+LS06[x] == 0):
#        T[x]=3
#    if (LB01[x]+LB02[x]+LB03[x]+LB04[x]+LB05[x] == 0) and (LB06[x] == 1) and (LS01[x]+LS02[x]+LS03[x]+LS04[x]+LS05[x]+LS06[x] > 0):
#        T[x]=3    
#    if (LB01[x]+LB02[x]+LB03[x]+LB04[x]+LB05[x] == 0) and (LB06[x] == 1) and (LS01[x]+LS02[x]+LS03[x]+LS04[x]+LS05[x]+LS06[x] == 0):
#        T[x]=2    
#    if (LB01[x]+LB02[x]+LB03[x]+LB04[x]+LB05[x]+LB06[x] == 0) and (LS01[x]+LS02[x]+LS03[x]+LS04[x]+LS05[x]+LS06[x] > 1):
#        T[x]=2    
#    if (LB01[x]+LB02[x]+LB03[x]+LB04[x]+LB05[x]+LB06[x] == 0) and (LS01[x]+LS02[x]+LS03[x]+LS04[x]+LS05[x]+LS06[x] == 1):
#        T[x]=1    
#    
      



