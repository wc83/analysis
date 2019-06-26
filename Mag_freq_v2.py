#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 11:16:46 2018

@author: william
"""

import obspy
import io
from obspy import read
from obspy.clients.earthworm import Client
from obspy import UTCDateTime
from obspy import Stream
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from numpy import genfromtxt
from numpy.polynomial.polynomial import polyfit


            
 #%%

cat = genfromtxt("/Users/william/Documents/scanner/all_stations/Final_Catalogue_2014_2018_v2.csv", delimiter=',',skip_header=1)

#cat = genfromtxt("/Users/william/Documents/scanner/Olivers_Catalogue_2014_2017_Energy.csv", delimiter=',',skip_header=1)


evt = cat[:,0]
#Energy= cat[:,7]
#year = cat[:,1]
Energy= cat[:,23]
year = cat[:,3]


e1=0
e2=0
e3=0
e4=0
e5=0
e6=0
e7=0
e8=0
e9=0
e10=0
e11=0
e12=0
e13=0




for x in range(0,len(cat)):
    if year[x] < 2019 :
        if 1e4 < Energy[x] :
            e1 +=1
        if 1e5 < Energy[x] :
            e2 +=1
        if 1e6 < Energy[x] :
            e3 +=1
        if 1e7 < Energy[x] :
            e4 +=1
        if 1e8 < Energy[x] :
            e5 +=1           
        if 1e9 < Energy[x] :
            e6 +=1
        if 1e10 < Energy[x] :
            e7 +=1            
        if 1e11 < Energy[x] :
            e8 +=1
#    if 1e11 < Energy[x] :
#        e9 +=1            
#    if 9e8 < Energy[x] :
#        e10 +=1
#    if 10e8 < Energy[x] :
#        e11 +=1            

        
        
        
               
E_his=np.zeros(shape=(8,2)) 
E_his[0][0]=e1          
E_his[1][0]=e2  
E_his[2][0]=e3  
E_his[3][0]=e4     
E_his[4][0]=e5  
E_his[5][0]=e6  
E_his[6][0]=e7  
E_his[7][0]=e8  
#E_his[8][0]=e9  
#E_his[9][0]=e10  
#E_his[10][0]=e11  
#E_his[11][0]=e12 



    
E_his[0,1]=1e4
E_his[1,1]=1e5
E_his[2,1]=1e6
E_his[3,1]=1e7
E_his[4,1]=1e8
E_his[5,1]=1e9
E_his[6,1]=1e10
E_his[7,1]=1e11
#E_his[8,1]=1e12
       
x1a=E_his[:,1]
y1a=E_his[:,0]

x1=E_his[2:8,1]
y1=E_his[2:8,0]


lx1=np.zeros(shape=(len(x1),1))
ly1=np.zeros(shape=(len(y1),1))

for p in range(0,len(x1)):
    lx=np.log10(x1[p])
    ly=np.log10(y1[p])
    
    lx1[p][0]=lx
    ly1[p][0]=ly
    
#log(y) = log(c) + m*log(x)   


plt.figure()
plt.plot(lx1[:],ly1[:],'bx-')

lx1=np.concatenate(lx1)
ly1=np.concatenate(ly1)
#
c, m = polyfit(lx1, ly1, 1)

plt.plot(lx1, ly1, 'kx-')
plt.plot(lx1, c + m * lx1, 'r-')
plt.xlabel("Log10|Energy|")
plt.ylabel("Log10|Occurance|")
print("y= d*x^m;  m= ",m , "d= ", 10**c)
print("log(y)= c + m*log(x);  m= ",m , "c= ", c)





plt.figure()
plt.loglog(x1a,y1a,'bx') 
plt.xlim([1e3,1e13])           
plt.ylim([0.5,1e5])  
plt.xlabel("Energy [J]")
plt.ylabel("Occurance [#]")  
plt.title("Will's Catalogue")


x2=np.linspace(1e3,1e15, 1e6)
#y2=10**(8.321919997363567  -0.7140633710717429*np.log10(x2))
y3= (10**c)*(x2**(m))
plt.loglog(x2,y3,'r-')












#E_h = np.zeros(shape=(0,1))
#num=0
#for x in range(0,len(cat)):
#    if cat[x,23] <  150e5:
#        E_h = np.lib.pad(E_h, ((0,1),(0,0)), 'constant', constant_values=(0))  
#        E_h[num][0] = cat[x,23] 
#        num+=1
#               
            
#plt.figure()           
#plt.hist(E_h,bins=30, histtype='step', color='b')          
#plt.yscale('log')
#plt.xscale('log')


#plt.plot([1e9,1e9],[1,10000],'k--')
##plt.plot([1e6,1e6],[1,10000],'k--')
#plt.plot([0,1e9],[78,78],'k--')
#plt.plot([0,5e7],[1110,1110],'k--')


plt.legend(('Energy-Frequency Distribution','Power-Law Fit'))
#%%          
#


#%%
#lE_h = np.zeros(shape=(len(E_h),1))
#for p in range(0,len(E_h)):
#    
#    leh=np.log10(E_h[p])
#    lE_h[p][0]=leh
#
#
#lEh = np.concatenate(lE_h)
#
#plt.figure()
#plt.hist(lE_h,bins=30, histtype='step')     








         
            
            
            
            
            
            