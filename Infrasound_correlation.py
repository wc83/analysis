#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 16:08:49 2019

@author: william
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 11:38:50 2018

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


year1=2014
month1=11
day1=24
hour1=0
minute1=0
second1=0

r1=4630
r2=3370
r3=2310
r4=1300
r5=810
r6=7660


rs1=4430
rs2=3190
rs3=2270
rs4=950
rs5=1100
rs6=7510


sso=340
sse = 450
sott=r1/sso
sett=rs1/sse
extd = sott-sett

td=np.zeros(shape=(0,4))
num=0

sta = 'LB01' # STATION 
net = 'Z4'  # 
loc = ''    # location, it depends mostly of which network you are in. 
# find acoustic signal
cha = 'HDF' # CHANNEL

client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23


cat= genfromtxt("/Users/william/Documents/scanner/all_stations/Explosion_Catalogue_V4c.csv", delimiter=',',skip_header=1)

#%% get stream wirh all traces

# open stream to save all traces

st_all = Stream()

#print('getting data')
for x in range(0,len(cat)):
#for x in range(0,200):  
#    print('data at', x+1,'of 50')
    try:
        event=cat[x,0]
        
        
        t1 = UTCDateTime(event) 
        t2 = t1 + 40
        st = Stream()
        st = client.get_waveforms(net, sta, '', cha, t1-10 , t2)
        
        st.detrend(type='linear')
        st.detrend(type='demean')
        st.filter("bandpass", freqmin=0.2,freqmax=2)
        
        
        trs=st[0]
        sr = trs.stats.sampling_rate
        nsta=int(0.5*sr)                                      #2
        nlta=int(15*sr)                                     #20
        stream=trs.data
        cft=recursive_sta_lta(stream, nsta, nlta)
        trig_on=10                                           #8
        trig_off=1                                        #0.2
#        plot_trigger(trs, cft, trig_on, trig_off) 
    
        on_off = trigger_onset(cft,trig_on,trig_off)
        
        if len(on_off) > 0:

            start=st[0].stats.starttime + (np.argmin(st[0])/100) - 3
#            start=t1+on_off[0,0]/100
#            st_s = st[0].slice(starttime=start - 12, endtime=start-6)
            st_s = st[0].slice(starttime=start , endtime=start+6)
        
#            st_s.plot(color='b')#,starttime=start - 13, endtime=start-3)
            st_all.append(st_s) #save all traces
                
            
    except:
        p=1

print('Files read in')
#%% get correlations between all traces
#print('correlating')        
shift=200
L=len(st_all)
cor_M = np.zeros(shape=(L,L))

Z = np.zeros(shape=(0,1))

num=0
for x in range(0,len(st_all)):
    print('x=' ,x+1, 'out of: ',len(st_all))
    for y in range(0,len(st_all)):
        
        if x == y:
            cor_M[x,y] = 1 
            Z= np.lib.pad(Z, ((0,1),(0,0)), 'constant', constant_values=(0))
            Z[num][0]=1
            num += 1
        if x < y :
 
            top_v,top,corell = corel_pi(st_all[x],st_all[y],shift)
                
            cor_M[x,y] = top_v 
            cor_M[y,x] = top_v 
    
            Z= np.lib.pad(Z, ((0,1),(0,0)), 'constant', constant_values=(0))
            Z[num][0]=top_v
            
            num += 1
        if x > y :
            Z= np.lib.pad(Z, ((0,1),(0,0)), 'constant', constant_values=(0))
            Z[num][0]= cor_M[y,x]
            num += 1
            
            
print('all correlated')
#%% create correlation matrix
        
X=np.linspace(0,len(st_all)-1,len(st_all))
Y=np.linspace(0,len(st_all)-1,len(st_all))

X2, Y2 = np.meshgrid(X, Y)
Z2 = Z.reshape(len(X), len(Y))


#%% plot correlation matrix
plt.pcolor(X2, Y2, Z2)
plt.colorbar()
plt.show()


#%% secondaries

sec_cor=0
con=0
size_time_cor=np.zeros(shape=(0,3))
for x in range(2,len(st_all)):
    
    if st_all[x].stats.starttime - st_all[x-1].stats.starttime < 600:
        if st_all[x-1].stats.starttime - st_all[x-2].stats.starttime < 600:
            
            sec_cor += cor_M[x-2,x]
            
            size_time_cor = np.lib.pad(size_time_cor, ((0,1),(0,0)), 'constant', constant_values=(0))
            size_time_cor[con][0]= cor_M[x-2,x]
            size_time_cor[con][1]= st_all[x].stats.starttime - st_all[x-2].stats.starttime
            size_time_cor[con][2]= (max(abs(st_all[x].data))/(max(abs(st_all[x-2].data))))*100
            con += 1
            
            print('primary - secondary correlation = ', cor_M[x-1,x])
            if len(st_all)- x >1:
                print('primary - next primary correlation = ', cor_M[x-1,x+1])
            if len(st_all)- x >2:
                print('primary - 2nd following primary correlation = ', cor_M[x-1,x+2])
            
            
            if cor_M[x-2,x] < 1:
                st_all[x-2].plot(color='r')
                st_all[x].plot(color='b')
#                print('acoustic correlation = ', cor_M[x-2,x])
                
            if len(st_all)- x >1:
                st_all[x+1].plot(color='g')
            if len(st_all)- x >2:
                st_all[x+2].plot(color='k')
            
            
        else:
            sec_cor += cor_M[x-1,x]
            size_time_cor = np.lib.pad(size_time_cor, ((0,1),(0,0)), 'constant', constant_values=(0))
            size_time_cor[con][0]= cor_M[x-1,x]
            size_time_cor[con][1]= st_all[x].stats.starttime - st_all[x-1].stats.starttime
            size_time_cor[con][2]= (max(abs(st_all[x].data))/(max(abs(st_all[x-1].data))))*100
            con += 1
            
#            print('primary - secondary chain:')
            print('primary - secondary correlation = ', cor_M[x-1,x])
            if len(st_all)- x >1:
                print('primary - next primary correlation = ', cor_M[x-1,x+1])
                if len(st_all)- x >2:
                    print('primary - 2nd following primary correlation = ', cor_M[x-1,x+2])
            
            if cor_M[x-1,x] < 1:
                st_all[x-1].plot(color='r')
                st_all[x].plot(color='b')
#                print('acoustic correlation = ', cor_M[x-1,x])
            
            if len(st_all)- x >1:
                st_all[x+1].plot(color='g')
                if len(st_all)- x >2:
                    st_all[x+2].plot(color='k')
    
av_sec_cor = sec_cor/con

    
next_cor = 0
number=0
for x in range(0,len(st_all)-1):
    
    next_cor += cor_M[x,x+1]
    number +=1
    
av_n_cor = next_cor/number
    

all_cor = 0
count=0
for x in range(0,len(st_all)-2):
    for y in range(x+2,len(st_all)):
        
        all_cor += cor_M[x,y]
        count +=1

av_cor = all_cor/count

hour_cor =0
cont=0
for x in range(3,len(st_all)):
    
    if st_all[x].stats.starttime - st_all[x-1].stats.starttime < 3600:
        hour_cor += cor_M[x-1,x]
        cont += 1
        
    if st_all[x].stats.starttime - st_all[x-2].stats.starttime < 3600:
        hour_cor += cor_M[x-2,x]
        cont += 1
        
    if st_all[x].stats.starttime - st_all[x-3].stats.starttime < 3600:
        hour_cor += cor_M[x-3,x]
        cont += 1
        
av_hour_cor=hour_cor/cont        
        
        
print('average correlation of following explosion = ', av_n_cor)
print('average correlation of non-neighbouring explosions = ', av_cor)    
print('average correlation of secondary explosion = ', av_sec_cor) 
print('average correlation of explosion within an hour = ', av_hour_cor)
print('number of explosions =', len(st_all))    
    
#%%
np.savetxt("/Users/william/Documents/scanner/analysis/LB01_inf_correltion_matrix_6s_v2.csv", cor_M,delimiter=",",header="")  
#%%

event_details = np.zeros(shape=(len(st_all),1))
for x in range(0,len(st_all)):
    event_details[x][0]=st_all[x].stats.starttime

np.savetxt("/Users/william/Documents/scanner/analysis/LB01_inf_correltion_event_details_6s_v2.csv",event_details ,delimiter=",",header="")  
    
#%%    
np.savetxt("/Users/william/Documents/scanner/analysis/LB01_inf_correltion_list_6s_v2.csv",Z ,delimiter=",",header="")  
    
    

    