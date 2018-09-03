#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 12:33:53 2018

@author: william
"""
import numpy as np
import matplotlib.pyplot as plt

num_active=np.zeros(shape=(1,16))
first = 0
last = 910
for x in range(first,last): 
    print('day',x+1,'of', last)
    num,lb_num,ls_num,LB01,LB02,LB03,LB04,LB05,LB06,LS01,LS02,LS03,LS04,LS05,LS06 =get_activity(x)
    
    num_active[x-first][0]=x+1
    num_active[x-first][1]=num
    num_active[x-first][2]=lb_num
    num_active[x-first][3]=ls_num
    num_active[x-first][4]=LB01
    num_active[x-first][5]=LB02
    num_active[x-first][6]=LB03
    num_active[x-first][7]=LB04
    num_active[x-first][8]=LB05
    num_active[x-first][9]=LB06
    num_active[x-first][10]=LS01
    num_active[x-first][11]=LS02
    num_active[x-first][12]=LS03
    num_active[x-first][13]=LS04
    num_active[x-first][14]=LS05
    num_active[x-first][15]=LS06
    if x != last-1:
        num_active=np.lib.pad(num_active, ((0,1),(0,0)), 'constant', constant_values=(0))

#%% save data
     
np.savetxt("/Users/william/Documents/scanner/all_stations/num_active_stations.csv", num_active,delimiter=",",header="day,num_active,lb_num_active,ls_num_active,LB01,LB02,LB03,LB04,LB05,LB06,LS01,LS02,LS03,LS04,LS05,LS06")

#%% plot data
plt.figure(1)
plt.plot(num_active[:,0],num_active[:,1])
plt.xlabel('Day')
plt.ylabel('Number of active stations')
plt.title('All Station Activity')
plt.ylim([0,12])

plt.figure(2)
plt.plot(num_active[:,0],num_active[:,2])
plt.xlabel('Day')
plt.ylabel('Number of active LB stations')
plt.title('Broadband Station Activity')
plt.ylim([0,7])

plt.figure(3)
plt.plot(num_active[:,0],num_active[:,3])
plt.xlabel('Day')
plt.ylabel('Number of active LS stations')
plt.title('Short Period Station Activity')
plt.ylim([0,7])