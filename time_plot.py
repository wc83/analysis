#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 12:38:29 2018

@author: william
"""

import matplotlib.pyplot as plt
from numpy import genfromtxt
import numpy as np
import matplotlib.dates as md
from obspy import UTCDateTime


cat = genfromtxt("/Users/william/Documents/scanner/all_stations/num_active_stations.csv", delimiter=",")

day=cat[:,0]
np.transpose(day)
LB01=cat[:,4]
np.transpose(LB01)
LB02=cat[:,5]
np.transpose(LB02)
LB03=cat[:,6]
np.transpose(LB03)
LB04=cat[:,7]
np.transpose(LB04)
LB05=cat[:,8]
np.transpose(LB05)
LB06=cat[:,9]
np.transpose(LB06)
LS01=cat[:,10]
np.transpose(LS01)
LS02=cat[:,11]
np.transpose(LS02)
LS03=cat[:,12]
np.transpose(LS03)
LS04=cat[:,13]
np.transpose(LS04)
LS05=cat[:,14]
np.transpose(LS05)
LS06=cat[:,15]
np.transpose(LS06)
TS=cat[:,16]
np.transpose(TS)
tot_num=cat[:,1]
np.transpose(tot_num)


#%% LB01 
LB01a=np.zeros(shape=(0,1))

switch_on = 0
change=0

for x in range(0,len(LB01)):
    time_change = TS[x]
    
    
    if LB01[x]==1 and switch_on==0:
        switch_on=1
        LB01a=np.lib.pad(LB01a, ((0,1),(0,0)), 'constant', constant_values=(0))
        LB01a[change]=time_change
        change +=1
        
    if LB01[x] == 0 and switch_on ==1:
        switch_on=0
        LB01a=np.lib.pad(LB01a, ((0,1),(0,0)), 'constant', constant_values=(0))
        LB01a[change]=time_change - LB01a[change-1]
        change +=1
        
#%% LB02 
LB02a=np.zeros(shape=(0,1))

switch_on = 0
change=0

for x in range(0,len(LB02)):
    time_change = TS[x]
    
    
    if LB02[x]==1 and switch_on==0:
        switch_on=1
        LB02a=np.lib.pad(LB02a, ((0,1),(0,0)), 'constant', constant_values=(0))
        LB02a[change]=time_change
        change +=1
        
    if LB02[x] == 0 and switch_on ==1:
        switch_on=0
        LB02a=np.lib.pad(LB02a, ((0,1),(0,0)), 'constant', constant_values=(0))
        LB02a[change]=time_change - LB02a[change-1]
        change +=1
 
#%% LB03
LB03a=np.zeros(shape=(0,1))

switch_on = 0
change=0

for x in range(0,len(LB03)):
    time_change = TS[x]
    
    
    if LB03[x]==1 and switch_on==0:
        switch_on=1
        LB03a=np.lib.pad(LB03a, ((0,1),(0,0)), 'constant', constant_values=(0))
        LB03a[change]=time_change
        change +=1
        
    if LB03[x] == 0 and switch_on ==1:
        switch_on=0
        LB03a=np.lib.pad(LB03a, ((0,1),(0,0)), 'constant', constant_values=(0))
        LB03a[change]=time_change - LB03a[change-1]
        change +=1
 
#%% LB04 
LB04a=np.zeros(shape=(0,1))

switch_on = 0
change=0

for x in range(0,len(LB04)):
    time_change = TS[x]
    
    
    if LB04[x]==1 and switch_on==0:
        switch_on=1
        LB04a=np.lib.pad(LB04a, ((0,1),(0,0)), 'constant', constant_values=(0))
        LB04a[change]=time_change
        change +=1
        
    if LB04[x] == 0 and switch_on ==1:
        switch_on=0
        LB04a=np.lib.pad(LB04a, ((0,1),(0,0)), 'constant', constant_values=(0))
        LB04a[change]=time_change - LB04a[change-1]
        change +=1
 
#%% LB05 
LB05a=np.zeros(shape=(0,1))

switch_on = 0
change=0

for x in range(0,len(LB05)):
    time_change = TS[x]
    
    
    if LB05[x]==1 and switch_on==0:
        switch_on=1
        LB05a=np.lib.pad(LB05a, ((0,1),(0,0)), 'constant', constant_values=(0))
        LB05a[change]=time_change
        change +=1
        
    if LB05[x] == 0 and switch_on ==1:
        switch_on=0
        LB05a=np.lib.pad(LB05a, ((0,1),(0,0)), 'constant', constant_values=(0))
        LB05a[change]=time_change - LB05a[change-1]
        change +=1
  
#%% LB06 
LB06a=np.zeros(shape=(0,1))

switch_on = 0
change=0

for x in range(0,len(LB06)):
    time_change = TS[x]
    
    
    if LB06[x]==1 and switch_on==0:
        switch_on=1
        LB06a=np.lib.pad(LB06a, ((0,1),(0,0)), 'constant', constant_values=(0))
        LB06a[change]=time_change
        change +=1
        
    if LB06[x] == 0 and switch_on ==1:
        switch_on=0
        LB06a=np.lib.pad(LB06a, ((0,1),(0,0)), 'constant', constant_values=(0))
        LB06a[change]=time_change - LB06a[change-1]
        change +=1

#%%  

LS01a=np.zeros(shape=(0,1))

switch_on = 0
change=0
 
for x in range(0,len(LS01)):
    time_change = TS[x]
    
    
    if LS01[x]==1 and switch_on==0:
        switch_on=1
        LS01a=np.lib.pad(LS01a, ((0,1),(0,0)), 'constant', constant_values=(0))
        LS01a[change]=time_change
        change +=1
        
    if LS01[x] == 0 and switch_on ==1:
        switch_on=0
        LS01a=np.lib.pad(LS01a, ((0,1),(0,0)), 'constant', constant_values=(0))
        LS01a[change]=time_change - LS01a[change-1]
        change +=1
        
#%% LB02 
LS02a=np.zeros(shape=(0,1))

switch_on = 0
change=0

for x in range(0,len(LS02)):
    time_change = TS[x]
    
    
    if LS02[x]==1 and switch_on==0:
        switch_on=1
        LS02a=np.lib.pad(LS02a, ((0,1),(0,0)), 'constant', constant_values=(0))
        LS02a[change]=time_change
        change +=1
        
    if LS02[x] == 0 and switch_on ==1:
        switch_on=0
        LS02a=np.lib.pad(LS02a, ((0,1),(0,0)), 'constant', constant_values=(0))
        LS02a[change]=time_change - LS02a[change-1]
        change +=1
 
#%% LS03
LS03a=np.zeros(shape=(0,1))

switch_on = 0
change=0
 
for x in range(0,len(LS03)):
    time_change = TS[x]
    
    
    if LS03[x]==1 and switch_on==0:
        switch_on=1
        LS03a=np.lib.pad(LS03a, ((0,1),(0,0)), 'constant', constant_values=(0))
        LS03a[change]=time_change
        change +=1
        
    if LS03[x] == 0 and switch_on ==1:
        switch_on=0
        LS03a=np.lib.pad(LS03a, ((0,1),(0,0)), 'constant', constant_values=(0))
        LS03a[change]=time_change - LS03a[change-1]
        change +=1
 
#%% LS04 
LS04a=np.zeros(shape=(0,1))

switch_on = 0
change=0

for x in range(0,len(LS04)):
    time_change = TS[x]
    
    
    if LS04[x]==1 and switch_on==0:
        switch_on=1
        LS04a=np.lib.pad(LS04a, ((0,1),(0,0)), 'constant', constant_values=(0))
        LS04a[change]=time_change
        change +=1
        
    if LS04[x] == 0 and switch_on ==1:
        switch_on=0
        LS04a=np.lib.pad(LS04a, ((0,1),(0,0)), 'constant', constant_values=(0))
        LS04a[change]=time_change - LS04a[change-1]
        change +=1
 
#%% LS05 
LS05a=np.zeros(shape=(0,1))

switch_on = 0
change=0

for x in range(0,len(LS05)):
    time_change = TS[x]
    
    
    if LS05[x]==1 and switch_on==0:
        switch_on=1
        LS05a=np.lib.pad(LS05a, ((0,1),(0,0)), 'constant', constant_values=(0))
        LS05a[change]=time_change
        change +=1
        
    if LS05[x] == 0 and switch_on ==1:
        switch_on=0
        LS05a=np.lib.pad(LS05a, ((0,1),(0,0)), 'constant', constant_values=(0))
        LS05a[change]=time_change - LS05a[change-1]
        change +=1
  
#%% LS06 
LS06a=np.zeros(shape=(0,1))

switch_on = 0
change=0

for x in range(0,len(LS06)):
    time_change = TS[x]
    
    
    if LS06[x]==1 and switch_on==0:
        switch_on=1
        LS06a=np.lib.pad(LS06a, ((0,1),(0,0)), 'constant', constant_values=(0))
        LS06a[change]=time_change
        change +=1
        
    if LS06[x] == 0 and switch_on ==1:
        switch_on=0
        LS06a=np.lib.pad(LS06a, ((0,1),(0,0)), 'constant', constant_values=(0))
        LS06a[change]=time_change - LS06a[change-1]
        change +=1


#%%

fig, ax = plt.subplots()
for p in range(0,len(LB01a),2):    
    ax.broken_barh([(LB01a[p],LB01a[p+1])],(2.5,5),facecolors='blue')               
for p in range(0,len(LB02a),2):    
    ax.broken_barh([(LB02a[p],LB02a[p+1])],(10,5),facecolors='black') 
for p in range(0,len(LB03a),2):    
    ax.broken_barh([(LB03a[p],LB03a[p+1])],(17.5,5),facecolors='blue') 
for p in range(0,len(LB04a),2):    
    ax.broken_barh([(LB04a[p],LB04a[p+1])],(25,5),facecolors='black') 
for p in range(0,len(LB05a),2):    
    ax.broken_barh([(LB05a[p],LB05a[p+1])],(32.5,5),facecolors='blue') 
for p in range(0,len(LB06a),2):    
    ax.broken_barh([(LB06a[p],LB06a[p+1])],(40,5),facecolors='black') 
for p in range(0,len(LS01a),2):    
    ax.broken_barh([(LS01a[p],LS01a[p+1])],(47.5,5),facecolors='blue')               
for p in range(0,len(LS02a),2):    
    ax.broken_barh([(LS02a[p],LS02a[p+1])],(55,5),facecolors='black') 
for p in range(0,len(LS03a),2):    
    ax.broken_barh([(LS03a[p],LS03a[p+1])],(62.5,5),facecolors='blue') 
for p in range(0,len(LS04a),2):   
    ax.broken_barh([(LS04a[p],LS04a[p+1])],(70,5),facecolors='black') 
for p in range(0,len(LS05a),2):    
    ax.broken_barh([(LS05a[p],LS05a[p+1])],(77.5,5),facecolors='blue') 
for p in range(0,len(LS06a),2):    
    ax.broken_barh([(LS06a[p],LS06a[p+1])],(85,5),facecolors='black') 

ax2 = ax.twinx()
ax2.plot(TS,tot_num,'r')
      
       
      
ax.set_ylim(0, 92.5)
ax.set_xlim(TS[0] , TS[-1])        
        
ax.set_xlabel('TimeStamp Time')
ax.set_ylabel('Station')
ax2.set_ylabel('Number of active stations')
ax.set_yticks([5,12.5,20,27.5,35,42.5,50,57.5,65,72.5,80,87.5])
ax.set_yticklabels(['LB01','LB02','LB03','LB04','LB05','LB06','LS01','LS02','LS03','LS04','LS05','LS06'])        
ax.set_title('Station Activity')        
plt.show()        
        
        
        

#np.savetxt("/Users/william/Documents/scanner/analysis/tp_LB01.csv", LB01a,delimiter=",",header="time/length")
#np.savetxt("/Users/william/Documents/scanner/analysis/tp_LB02.csv", LB02a,delimiter=",",header="time/length")
#np.savetxt("/Users/william/Documents/scanner/analysis/tp_LB03.csv", LB03a,delimiter=",",header="time/length")
#np.savetxt("/Users/william/Documents/scanner/analysis/tp_LB04.csv", LB04a,delimiter=",",header="time/length")        
#np.savetxt("/Users/william/Documents/scanner/analysis/tp_LB05.csv", LB05a,delimiter=",",header="time/length")
#np.savetxt("/Users/william/Documents/scanner/analysis/tp_LB06.csv", LB06a,delimiter=",",header="time/length")
#np.savetxt("/Users/william/Documents/scanner/analysis/tp_LS01.csv", LS01a,delimiter=",",header="time/length")         
#np.savetxt("/Users/william/Documents/scanner/analysis/tp_LS02.csv", LS02a,delimiter=",",header="time/length") 
#np.savetxt("/Users/william/Documents/scanner/analysis/tp_LS03.csv", LS03a,delimiter=",",header="time/length")       
#np.savetxt("/Users/william/Documents/scanner/analysis/tp_LS04.csv", LS04a,delimiter=",",header="time/length")        
#np.savetxt("/Users/william/Documents/scanner/analysis/tp_LS05.csv", LS05a,delimiter=",",header="time/length")        
#np.savetxt("/Users/william/Documents/scanner/analysis/tp_LS06.csv", LS06a,delimiter=",",header="time/length")
        
        
