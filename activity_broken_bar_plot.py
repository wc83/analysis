import os
from collections import OrderedDict
import numpy as np
import obspy
import scipy.signal as sgn
import matplotlib.pyplot as plt 
import matplotlib.mlab as mlab

from obspy.clients.earthworm import Client
from obspy import UTCDateTime
from obspy.signal.trigger import trigger_onset
from numpy import genfromtxt
#from scipy.signal import welch
from obspy import Stream




cat= genfromtxt("/Users/william/Documents/scanner/all_stations/num_active_stations.csv", delimiter=',',skip_header=1)

LB01a = cat[:,4]
LB02a = cat[:,5]
LB03a = cat[:,6]
LB04a = cat[:,7]
LB05a = cat[:,8]
LB06a = cat[:,9]
LS01a = cat[:,10]
LS02a = cat[:,11]
LS03a = cat[:,12]
LS04a = cat[:,13]
LS05a = cat[:,14]
LS06a = cat[:,15]

ts = cat[:,16]


#%% example
#
#
#fig, ax = plt.subplots()
#ax.broken_barh([(110, 30), (150, 10)], (10, 9), facecolors='tab:blue')
#ax.broken_barh([(10, 50), (100, 20), (130, 10)], (20, 9),
#               facecolors=('tab:orange', 'tab:green', 'tab:red'))
#ax.set_ylim(5, 35)
#ax.set_xlim(0, 200)
#ax.set_xlabel('seconds since start')
#ax.set_yticks([15, 25])
#ax.set_yticklabels(['Bill', 'Jim'])
#ax.grid(True)
#ax.annotate('race interrupted', (61, 25),
#            xytext=(0.8, 0.9), textcoords='axes fraction',
#            arrowprops=dict(facecolor='black', shrink=0.05),
#            fontsize=16,
#            horizontalalignment='right', verticalalignment='top')
#
#plt.show()
#

#%%

switch =0
count =0
off_count =0
num=0
turned=0
LB01_a = np.zeros(shape=(0,1))


for x in range(0,len(LB01a)):
    if switch == 0:
        if int(LB01a[x])==0:
            off_count += 1
        if int(LB01a[x])==1:
            LB01_a = np.lib.pad(LB01_a, ((0,1),(0,0)), 'constant', constant_values=(0))
            LB01_a[num]=x
            num+=1
            count = 1 
            off_count = 0
            switch = 1
            turned = 1
            
    if switch ==1 and turned == 0:
        if int(LB01a[x])==1:
            count += 1
        if int(LB01a[x])==0:
            LB01_a = np.lib.pad(LB01_a, ((0,1),(0,0)), 'constant', constant_values=(0))
            LB01_a[num] = count
            count=0
            off_count = 1
            num +=1
            switch = 0
            
    if turned == 1:
        turned =0
        switch =1
        
        
        

numb=0
LB01_b = np.zeros(shape=(2,1))
LB01_b[0]=0
LB01_b[1]=LB01_a[0]
for x in range(2,len(LB01_a),2):
    LB01_b = np.lib.pad(LB01_b, ((0,1),(0,0)), 'constant', constant_values=(0))
    LB01_b[x]=LB01_a[x-2]+LB01_a[x-1]
    LB01_b = np.lib.pad(LB01_b, ((0,1),(0,0)), 'constant', constant_values=(0))
    LB01_b[x+1]= LB01_a[x]-LB01_b[x]
    
    

#%%

switch =0
count =0
num=0
turned=0
LB02_a = np.zeros(shape=(0,1))

for x in range(0,len(LB02a)):
    if switch == 0:
        if int(LB02a[x])==1:
            LB02_a = np.lib.pad(LB02_a, ((0,1),(0,0)), 'constant', constant_values=(0))
            LB02_a[num]=x
            num+=1
            count = 1 
            switch = 1
            turned = 1
            
    if switch ==1 and turned == 0:
        if int(LB02a[x])==1:
            count += 1
        if int(LB02a[x])==0:
            LB02_a = np.lib.pad(LB02_a, ((0,1),(0,0)), 'constant', constant_values=(0))
            LB02_a[num] = count
            count=0
            num +=1
            switch = 0
            
    if turned == 1:
        turned =0
        switch =1


numb=0
LB02_b = np.zeros(shape=(2,1))
LB02_b[0]=0
LB02_b[1]=LB02_a[0]
for x in range(2,len(LB02_a),2):
    LB02_b = np.lib.pad(LB02_b, ((0,1),(0,0)), 'constant', constant_values=(0))
    LB02_b[x]=LB02_a[x-2]+LB02_a[x-1]
    LB02_b = np.lib.pad(LB02_b, ((0,1),(0,0)), 'constant', constant_values=(0))
    LB02_b[x+1]= LB02_a[x]-LB02_b[x]
    
    
#%%

switch =0
count =0
num=0
turned=0
LB03_a = np.zeros(shape=(0,1))

for x in range(0,len(LB03a)):
    if switch == 0:
        if int(LB03a[x])==1:
            LB03_a = np.lib.pad(LB03_a, ((0,1),(0,0)), 'constant', constant_values=(0))
            LB03_a[num]=x
            num+=1
            count = 1 
            switch = 1
            turned = 1
            
    if switch ==1 and turned == 0:
        if int(LB03a[x])==1:
            count += 1
        if int(LB03a[x])==0:
            LB03_a = np.lib.pad(LB03_a, ((0,1),(0,0)), 'constant', constant_values=(0))
            LB03_a[num] = count
            count=0
            num +=1
            switch = 0
            
    if turned == 1:
        turned =0
        switch =1


numb=0
LB03_b = np.zeros(shape=(2,1))
LB03_b[0]=0
LB03_b[1]=LB03_a[0]
for x in range(2,len(LB03_a),2):
    LB03_b = np.lib.pad(LB03_b, ((0,1),(0,0)), 'constant', constant_values=(0))
    LB03_b[x]=LB03_a[x-2]+LB03_a[x-1]
    LB03_b = np.lib.pad(LB03_b, ((0,1),(0,0)), 'constant', constant_values=(0))
    LB03_b[x+1]= LB03_a[x]-LB03_b[x]
#%%

switch =0
count =0
num=0
turned=0
LB04_a = np.zeros(shape=(0,1))

for x in range(0,len(LB04a)):
    if switch == 0:
        if int(LB04a[x])==1:
            LB04_a = np.lib.pad(LB04_a, ((0,1),(0,0)), 'constant', constant_values=(0))
            LB04_a[num]=x
            num+=1
            count = 1 
            switch = 1
            turned = 1
            
    if switch ==1 and turned == 0:
        if int(LB04a[x])==1:
            count += 1
        if int(LB04a[x])==0:
            LB04_a = np.lib.pad(LB04_a, ((0,1),(0,0)), 'constant', constant_values=(0))
            LB04_a[num] = count
            count=0
            num +=1
            switch = 0
            
    if turned == 1:
        turned =0
        switch =1


numb=0
LB04_b = np.zeros(shape=(2,1))
LB04_b[0]=0
LB04_b[1]=LB04_a[0]
for x in range(2,len(LB04_a),2):
    LB04_b = np.lib.pad(LB04_b, ((0,1),(0,0)), 'constant', constant_values=(0))
    LB04_b[x]=LB04_a[x-2]+LB04_a[x-1]
    LB04_b = np.lib.pad(LB04_b, ((0,1),(0,0)), 'constant', constant_values=(0))
    LB04_b[x+1]= LB04_a[x]-LB04_b[x]
#%%

switch =0
count =0
num=0
turned=0
LB05_a = np.zeros(shape=(0,1))

for x in range(0,len(LB05a)):
    if switch == 0:
        if int(LB05a[x])==1:
            LB05_a = np.lib.pad(LB05_a, ((0,1),(0,0)), 'constant', constant_values=(0))
            LB05_a[num]=x
            num+=1
            count = 1 
            switch = 1
            turned = 1
            
    if switch ==1 and turned == 0:
        if int(LB05a[x])==1:
            count += 1
        if int(LB05a[x])==0:
            LB05_a = np.lib.pad(LB05_a, ((0,1),(0,0)), 'constant', constant_values=(0))
            LB05_a[num] = count
            count=0
            num +=1
            switch = 0
            
    if turned == 1:
        turned =0
        switch =1
        
        
numb=0
LB05_b = np.zeros(shape=(2,1))
LB05_b[0]=0
LB05_b[1]=LB05_a[0]
for x in range(2,len(LB05_a),2):
    LB05_b = np.lib.pad(LB05_b, ((0,1),(0,0)), 'constant', constant_values=(0))
    LB05_b[x]=LB05_a[x-2]+LB05_a[x-1]
    LB05_b = np.lib.pad(LB05_b, ((0,1),(0,0)), 'constant', constant_values=(0))
    LB05_b[x+1]= LB05_a[x]-LB05_b[x]
#%%

switch =0
count =0
num=0
turned=0
LB06_a = np.zeros(shape=(0,1))

for x in range(0,len(LB06a)):
    if switch == 0:
        if int(LB06a[x])==1:
            LB06_a = np.lib.pad(LB06_a, ((0,1),(0,0)), 'constant', constant_values=(0))
            LB06_a[num]=x
            num+=1
            count = 1 
            switch = 1
            turned = 1
            
    if switch ==1 and turned == 0:
        if int(LB06a[x])==1:
            count += 1
        if int(LB06a[x])==0:
            LB06_a = np.lib.pad(LB06_a, ((0,1),(0,0)), 'constant', constant_values=(0))
            LB06_a[num] = count
            count=0
            num +=1
            switch = 0
            
    if turned == 1:
        turned =0
        switch =1     


numb=0
LB06_b = np.zeros(shape=(2,1))
LB06_b[0]=0
LB06_b[1]=LB06_a[0]
for x in range(2,len(LB06_a),2):
    LB06_b = np.lib.pad(LB06_b, ((0,1),(0,0)), 'constant', constant_values=(0))
    LB06_b[x]=LB06_a[x-2]+LB06_a[x-1]
    LB06_b = np.lib.pad(LB06_b, ((0,1),(0,0)), 'constant', constant_values=(0))
    LB06_b[x+1]= LB06_a[x]-LB06_b[x]        
#%%

switch =0
count =0
num=0
turned=0
LS01_a = np.zeros(shape=(0,1))

for x in range(0,len(LS01a)):
    if switch == 0:
        if int(LS01a[x])==1:
            LS01_a = np.lib.pad(LS01_a, ((0,1),(0,0)), 'constant', constant_values=(0))
            LS01_a[num]=x
            num+=1
            count = 1 
            switch = 1
            turned = 1
            
    if switch ==1 and turned == 0:
        if int(LS01a[x])==1:
            count += 1
        if int(LS01a[x])==0:
            LS01_a = np.lib.pad(LS01_a, ((0,1),(0,0)), 'constant', constant_values=(0))
            LS01_a[num] = count
            count=0
            num +=1
            switch = 0
            
    if turned == 1:
        turned =0
        switch =1 


numb=0
LS01_b = np.zeros(shape=(2,1))
LS01_b[0]=0
LS01_b[1]=LS01_a[0]
for x in range(2,len(LS01_a),2):
    LS01_b = np.lib.pad(LS01_b, ((0,1),(0,0)), 'constant', constant_values=(0))
    LS01_b[x]=LS01_a[x-2]+LS01_a[x-1]
    LS01_b = np.lib.pad(LS01_b, ((0,1),(0,0)), 'constant', constant_values=(0))
    LS01_b[x+1]= LS01_a[x]-LS01_b[x]        
#%%

switch =0
count =0
num=0
turned=0
LS02_a = np.zeros(shape=(0,1))

for x in range(0,len(LS02a)):
    if switch == 0:
        if int(LS02a[x])==1:
            LS02_a = np.lib.pad(LS02_a, ((0,1),(0,0)), 'constant', constant_values=(0))
            LS02_a[num]=x
            num+=1
            count = 1 
            switch = 1
            turned = 1
            
    if switch ==1 and turned == 0:
        if int(LS02a[x])==1:
            count += 1
        if int(LS02a[x])==0:
            LS02_a = np.lib.pad(LS02_a, ((0,1),(0,0)), 'constant', constant_values=(0))
            LS02_a[num] = count
            count=0
            num +=1
            switch = 0
            
    if turned == 1:
        turned =0
        switch =1          


numb=0
LS02_b = np.zeros(shape=(2,1))
LS02_b[0]=0
LS02_b[1]=LS02_a[0]
for x in range(2,len(LS02_a),2):
    LS02_b = np.lib.pad(LS02_b, ((0,1),(0,0)), 'constant', constant_values=(0))
    LS02_b[x]=LS02_a[x-2]+LS02_a[x-1]
    LS02_b = np.lib.pad(LS02_b, ((0,1),(0,0)), 'constant', constant_values=(0))
    LS02_b[x+1]= LS02_a[x]-LS02_b[x]
#%%

switch =0
count =0
num=0
turned=0
LS03_a = np.zeros(shape=(0,1))

for x in range(0,len(LS03a)):
    if switch == 0:
        if int(LS03a[x])==1:
            LS03_a = np.lib.pad(LS03_a, ((0,1),(0,0)), 'constant', constant_values=(0))
            LS03_a[num]=x
            num+=1
            count = 1 
            switch = 1
            turned = 1
            
    if switch ==1 and turned == 0:
        if int(LS03a[x])==1:
            count += 1
        if int(LS03a[x])==0:
            LS03_a = np.lib.pad(LS03_a, ((0,1),(0,0)), 'constant', constant_values=(0))
            LS03_a[num] = count
            count=0
            num +=1
            switch = 0
            
    if turned == 1:
        turned =0
        switch =1   


numb=0
LS03_b = np.zeros(shape=(2,1))
LS03_b[0]=0
LS03_b[1]=LS03_a[0]
for x in range(2,len(LS03_a),2):
    LS03_b = np.lib.pad(LS03_b, ((0,1),(0,0)), 'constant', constant_values=(0))
    LS03_b[x]=LS03_a[x-2]+LS03_a[x-1]
    LS03_b = np.lib.pad(LS03_b, ((0,1),(0,0)), 'constant', constant_values=(0))
    LS03_b[x+1]= LS03_a[x]-LS03_b[x]        
#%%

switch =0
count =0
num=0
turned=0
LS04_a = np.zeros(shape=(0,1))

for x in range(0,len(LS01a)):
    if switch == 0:
        if int(LS04a[x])==1:
            LS04_a = np.lib.pad(LS04_a, ((0,1),(0,0)), 'constant', constant_values=(0))
            LS04_a[num]=x
            num+=1
            count = 1 
            switch = 1
            turned = 1
            
    if switch ==1 and turned == 0:
        if int(LS04a[x])==1:
            count += 1
        if int(LS04a[x])==0:
            LS04_a = np.lib.pad(LS04_a, ((0,1),(0,0)), 'constant', constant_values=(0))
            LS04_a[num] = count
            count=0
            num +=1
            switch = 0
            
    if turned == 1:
        turned =0
        switch =1          


numb=0
LS04_b = np.zeros(shape=(2,1))
LS04_b[0]=0
LS04_b[1]=LS04_a[0]
for x in range(2,len(LS04_a),2):
    LS04_b = np.lib.pad(LS04_b, ((0,1),(0,0)), 'constant', constant_values=(0))
    LS04_b[x]=LS04_a[x-2]+LS04_a[x-1]
    LS04_b = np.lib.pad(LS04_b, ((0,1),(0,0)), 'constant', constant_values=(0))
    LS04_b[x+1]= LS04_a[x]-LS04_b[x]
#%%

switch =0
count =0
num=0
turned=0
LS05_a = np.zeros(shape=(0,1))

for x in range(0,len(LS05a)):
    if switch == 0:
        if int(LS05a[x])==1:
            LS05_a = np.lib.pad(LS05_a, ((0,1),(0,0)), 'constant', constant_values=(0))
            LS05_a[num]=x
            num+=1
            count = 1 
            switch = 1
            turned = 1
            
    if switch ==1 and turned == 0:
        if int(LS05a[x])==1:
            count += 1
        if int(LS05a[x])==0:
            LS05_a = np.lib.pad(LS05_a, ((0,1),(0,0)), 'constant', constant_values=(0))
            LS05_a[num] = count
            count=0
            num +=1
            switch = 0
            
    if turned == 1:
        turned =0
        switch =1          


numb=0
LS05_b = np.zeros(shape=(2,1))
LS05_b[0]=0
LS05_b[1]=LS05_a[0]
for x in range(2,len(LS05_a),2):
    LS05_b = np.lib.pad(LS05_b, ((0,1),(0,0)), 'constant', constant_values=(0))
    LS05_b[x]=LS05_a[x-2]+LS05_a[x-1]
    LS05_b = np.lib.pad(LS05_b, ((0,1),(0,0)), 'constant', constant_values=(0))
    LS05_b[x+1]= LS05_a[x]-LS05_b[x]        
#%%

switch =0
count =0
num=0
turned=0
LS06_a = np.zeros(shape=(0,1))

for x in range(0,len(LS06a)):
    if switch == 0:
        if int(LS06a[x])==1:
            LS06_a = np.lib.pad(LS06_a, ((0,1),(0,0)), 'constant', constant_values=(0))
            LS06_a[num]=x
            num+=1
            count = 1 
            switch = 1
            turned = 1
            
    if switch ==1 and turned == 0:
        if int(LS06a[x])==1:
            count += 1
        if int(LS06a[x])==0:
            LS06_a = np.lib.pad(LS06_a, ((0,1),(0,0)), 'constant', constant_values=(0))
            LS06_a[num] = count
            count=0
            num +=1
            switch = 0
            
    if turned == 1:
        turned =0
        switch =1  


numb=0
LS06_b = np.zeros(shape=(2,1))
LS06_b[0]=0
LS06_b[1]=LS06_a[0]
for x in range(2,len(LS06_a),2):
    LS06_b = np.lib.pad(LS06_b, ((0,1),(0,0)), 'constant', constant_values=(0))
    LS06_b[x]=LS06_a[x-2]+LS06_a[x-1]
    LS06_b = np.lib.pad(LS06_b, ((0,1),(0,0)), 'constant', constant_values=(0))
    LS06_b[x+1]= LS06_a[x]-LS06_b[x]
    
#%%

cat3= genfromtxt("/Users/william/Documents/scanner/output_data/events_per_day_full_2014_2018.csv",  delimiter=',',skip_header=1)

dayn=cat3[:,0]
detec=cat3[:,1]

switch =0
count =0
num=0
turned=0
STG8_a = np.zeros(shape=(0,1))

for x in range(1000,len(dayn)):
    if switch == 0:
        if detec[x]>0:
            STG8_a = np.lib.pad(STG8_a, ((0,1),(0,0)), 'constant', constant_values=(0))
            STG8_a[num]=x
            num+=1
            count = 1 
            switch = 1
            turned = 1
            
    if switch ==1 and turned == 0:
        if detec[x]>0:
            count += 1
        if detec[x]==0:
            STG8_a = np.lib.pad(STG8_a, ((0,1),(0,0)), 'constant', constant_values=(0))
            STG8_a[num] = count
            count=0
            num +=1
            switch = 0
            
    if turned == 1:
        turned =0
        switch =1  
        

    
#%%
fig, ax = plt.subplots()
ax.set_xlabel('Day')
ax.set_ylabel('Station')
ax.set_title('Station Activity')
ax.set_xlim(0, 910)
ax.set_yticks([3, 8, 13, 18, 23, 28,33,38,43,48,53,58])
ax.set_yticklabels(['LB01', 'LB02','LB03','LB04','LB05','LB06','LS01', 'LS02','LS03','LS04','LS05','LS06'])
for x in range(0,len(LB01_a),2):
#    if LB01_a[x+1]>2:
        ax.broken_barh([(LB01_a[x], LB01_a[x+1])], (1, 4), facecolors='k')
for x in range(0,len(LB02_a),2):
#    if LB02_a[x+1]>2:
        ax.broken_barh([(LB02_a[x], LB02_a[x+1])], (6, 4), facecolors='k')
for x in range(0,len(LB03_a),2):
#    if LB03_a[x+1]>5:
        ax.broken_barh([(LB03_a[x], LB03_a[x+1])], (11,4), facecolors='k')
for x in range(0,len(STG8_a),2):
#    if LB03_a[x+1]>5:
        ax.broken_barh([(STG8_a[x], STG8_a[x+1])], (11,4), facecolors='k')
for x in range(0,len(LB04_a),2):
#    if LB04_a[x+1]>2:
        ax.broken_barh([(LB04_a[x], LB04_a[x+1])], (16,4), facecolors='k')
for x in range(0,len(LB05_a),2):
#    if LB05_a[x+1]>2:
        ax.broken_barh([(LB05_a[x], LB05_a[x+1])], (21,4), facecolors='k')
for x in range(0,len(LB06_a),2):
#    if LB06_a[x+1]>2:
        ax.broken_barh([(LB06_a[x], LB06_a[x+1])], (26,4), facecolors='k')
    
for x in range(0,len(LS01_a),2):
    if LS01_a[x+1]>2:# and LS01_b[x+1] < 2:
        ax.broken_barh([(LS01_a[x], LS01_a[x+1])], (31, 4), facecolors='k')
for x in range(0,len(LS02_a),2):
    if LS02_a[x+1]>2:# and LS02_b[x+1] < 2:
        ax.broken_barh([(LS02_a[x], LS02_a[x+1])], (36, 4), facecolors='k')
for x in range(0,len(LS03_a),2):
#    if LS03_a[x+1]>2 and LS03_b[x+1] < 2:
        ax.broken_barh([(LS03_a[x], LS03_a[x+1])], (41,4), facecolors='k')
for x in range(0,len(LS04_a),2):
#    if LS04_a[x+1]>2 and LS04_b[x+1] < 2:
        ax.broken_barh([(LS04_a[x], LS04_a[x+1])], (46,4), facecolors='k')
for x in range(0,len(LS05_a),2):
#    if LS05_a[x+1]>2 and LS05_b[x+1] < 2:
        ax.broken_barh([(LS05_a[x], LS05_a[x+1])], (51,4), facecolors='k')
for x in range(0,len(LS06_a),2):
#    if LS06_a[x+1]>2 and LS06_b[x+1] < 2:
        ax.broken_barh([(LS06_a[x], LS06_a[x+1])], (56,4), facecolors='k')


for x in range(0,len(LB01_b),2):
    if LB01_b[x+1]<10:
        ax.broken_barh([(LB01_b[x], LB01_b[x+1])], (1, 4), facecolors='k')
for x in range(0,len(LB02_b),2):
    if LB02_b[x+1]<10:
        ax.broken_barh([(LB02_b[x], LB02_b[x+1])], (6, 4), facecolors='k')
for x in range(0,len(LB03_b),2):
    if LB03_b[x+1]<10:
        ax.broken_barh([(LB03_b[x], LB03_b[x+1])], (11, 4), facecolors='k')
for x in range(0,len(LB04_b),2):
    if LB04_b[x+1]<10:
        ax.broken_barh([(LB04_b[x], LB04_b[x+1])], (16, 4), facecolors='k')
for x in range(0,len(LB05_b),2):
    if LB05_b[x+1]<10:
        ax.broken_barh([(LB05_b[x], LB05_b[x+1])], (21, 4), facecolors='k')
for x in range(0,len(LB06_b),2):
    if LB06_b[x+1]<10:
        ax.broken_barh([(LB06_b[x], LB06_b[x+1])], (26, 4), facecolors='k')
for x in range(0,len(LS01_b),2):
    if LS01_b[x+1]<10:
        ax.broken_barh([(LS01_b[x], LS01_b[x+1])], (31, 4), facecolors='k')
for x in range(0,len(LS02_b),2):
    if LS02_b[x+1]<10:
        ax.broken_barh([(LS02_b[x], LS02_b[x+1])], (36, 4), facecolors='k')
for x in range(0,len(LS03_b),2):
    if LS03_b[x+1]<10:
        ax.broken_barh([(LS03_b[x], LS03_b[x+1])], (41, 4), facecolors='k')        
for x in range(0,len(LS04_b),2):
    if LS04_b[x+1]<10:
        ax.broken_barh([(LS04_b[x], LS04_b[x+1])], (46, 4), facecolors='k')  
for x in range(0,len(LS05_b),2):
    if LS05_b[x+1]<10:
        ax.broken_barh([(LS05_b[x], LS05_b[x+1])], (51, 4), facecolors='k')        
for x in range(0,len(LS06_b),2):
    if LS06_b[x+1]<10:
        ax.broken_barh([(LS06_b[x], LS06_b[x+1])], (56, 4), facecolors='k')        
        
        
        
        
        
        


