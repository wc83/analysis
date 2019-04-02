
import obspy
import io
from obspy import read
from obspy.clients.earthworm import Client
from obspy import UTCDateTime
from obspy import Stream
import numpy as np
import matplotlib.pyplot as plt

#%% constants


A=1
rhoE=2500
cE=2000
pi=3.14159
r1=4630
r2=3370
r3=2310
r4=1300
r5=810
r6=7660
rs1=5710
rs2=5490
rs3=3900
rs4=5520
rs5=4290
rs6=2610
      

saved = np.zeros(shape=(0,4))
num = 0
first =0

#%% read in data
reclen = 512
chunksize = 100000 * reclen # 100000 = Around 50 MB 
with io.open("/Users/william/Documents/scanner/output_data/m32.mseed", "rb") as fh:
        # just month 2
#with io.open("/Users/william/Documents/scanner/output_data/EXP_all_data_stream_2_month_2.mseed", "rb") as fh:
    while True:
        with io.BytesIO() as buf:
            c = fh.read(chunksize)
            if not c:
                break
            buf.write(c)
            buf.seek(0, 0)
            st = obspy.read(buf)
            
#%% for first 'load' - save time and energy of each event - save 'lowest ranking' station           
            if first == 0:
                first =1
                if st[0].stats.station == 'LB01':
                    station=1
                    r=r1
                if st[0].stats.station == 'LB02':
                    station=2
                    r=r2
                if st[0].stats.station == 'LB03':
                    station=3
                    r=r3
                if st[0].stats.station == 'LB04':
                    station=4
                    r=r4
                if st[0].stats.station == 'LB05':
                    station=5
                    r=r5
                if st[0].stats.station == 'LB06':
                    station=6
                    r=r6
                if st[0].stats.station == 'LS01':
                    station=7
                    r=rs1
                if st[0].stats.station == 'LS02':
                    station=8
                    r=rs2
                if st[0].stats.station == 'LS03':
                    station=9
                    r=rs3
                if st[0].stats.station == 'LS04':
                    station=10
                    r=rs4
                if st[0].stats.station == 'LS05':
                    station=11
                    r=rs5
                if st[0].stats.station == 'LS06':
                    station=12
                    r=rs6
                
                tr = st[0]
                st_c = calibrate1(tr)
                B=2*pi*rhoE*cE*(1/A)
                EI = sum(np.square(st_c[0].data))
                EE= B*(r*r)*EI
                
                saved = np.lib.pad(saved, ((0,1),(0,0)), 'constant', constant_values=(0))
                saved[num][0]=st[0].stats.starttime.timestamp
                saved[num][1]=max(abs(st[0].data))
                saved[num][2]=station
                saved[num][3]=EE
                num += 1
                
                for x in range (1,int(len(st))):
                    rt = st[x].stats.starttime.timestamp
                    near,ind = find_nearest(saved[:,0], rt)
                    
                    if st[x].stats.station == 'LB01':
                        station=1
                        r=r1
                    if st[x].stats.station == 'LB02':
                        station=2
                        r=r2
                    if st[x].stats.station == 'LB03':
                        station=3
                        r=r3
                    if st[x].stats.station == 'LB04':
                        station=4
                        r=r4
                    if st[x].stats.station == 'LB05':
                        station=5
                        r=r5
                    if st[x].stats.station == 'LB06':
                        station=6
                        r=r6
                    if st[x].stats.station == 'LS01':
                        station=7
                        r=rs1
                    if st[x].stats.station == 'LS02':
                        station=8
                        r=rs2
                    if st[x].stats.station == 'LS03':
                        station=9
                        r=rs3
                    if st[x].stats.station == 'LS04':
                        station=10
                        r=rs4
                    if st[x].stats.station == 'LS05':
                        station=11
                        r=rs5
                    if st[x].stats.station == 'LS06':
                        station=12
                        r=rs6
                        
                    tr = st[x]
                    st_c = calibrate1(tr)
                    B=2*pi*rhoE*cE*(1/A)
                    EI = sum(np.square(st_c[0].data))
                    EE= B*(r*r)*EI
                    
                    if abs(rt-near) > 10:
                        saved = np.lib.pad(saved, ((0,1),(0,0)), 'constant', constant_values=(0))
                        saved[num][0]=st[x].stats.starttime.timestamp
                        saved[num][1]=max((st[x].data))
                        saved[num][2]=station
                        saved[num][3]=EE
                        num += 1
                    else:
                        if station < saved[ind,2]:
                            saved[ind,2] = station
                            saved[ind,1] = max((st[x].data))
                            saved[ind,0] = st[x].stats.starttime.timestamp
                            saved[ind,3]=EE
                            
#%% for subsequent 'loads' - save time and energy of each event - save 'lowest ranking' station              
            else:
                for x in range (0,int(len(st))):
                    rt = st[x].stats.starttime.timestamp
                    near,ind = find_nearest(saved[:,0], rt)
                    
                    if st[x].stats.station == 'LB01':
                        station=1
                        r=r1
                    if st[x].stats.station == 'LB02':
                        station=2
                        r=r2
                    if st[x].stats.station == 'LB03':
                        station=3
                        r=r3
                    if st[x].stats.station == 'LB04':
                        station=4
                        r=r4
                    if st[x].stats.station == 'LB05':
                        station=5
                        r=r5
                    if st[x].stats.station == 'LB06':
                        station=6
                        r=r6
                    if st[x].stats.station == 'LS01':
                        station=7
                        r=rs1
                    if st[x].stats.station == 'LS02':
                        station=8
                        r=rs2
                    if st[x].stats.station == 'LS03':
                        station=9
                        r=rs3
                    if st[x].stats.station == 'LS04':
                        station=10
                        r=rs4
                    if st[x].stats.station == 'LS05':
                        station=11
                        r=rs5
                    if st[x].stats.station == 'LS06':
                        station=12
                        r=rs6
                        
                    tr = st[x]
                    st_c = calibrate1(tr)
                    B=2*pi*rhoE*cE*(1/A)
                    EI = sum(np.square(st_c[0].data))
                    EE= B*(r*r)*EI
                    
                    if abs(rt-near) > 10:
                        saved = np.lib.pad(saved, ((0,1),(0,0)), 'constant', constant_values=(0))
                        saved[num][0]=st[x].stats.starttime.timestamp
                        saved[num][1]=max((st[x].data))
                        saved[num][2]=station
                        saved[num][3]=EE
                        num += 1
                    else:
                        if station < saved[ind,2]:
                            saved[ind,2] = station
                            saved[ind,1] = max((st[x].data))
                            saved[ind,0] = st[x].stats.starttime.timestamp
                            saved[ind,3]=EE
                    
col=0                  
saved=saved[np.argsort(saved[:,col])]          # order in time         
#print(num) 


#%% count secondaries/primaries

last_p = saved[0,0]
peak1 = saved[0,1]
last_station = saved[0,2]
last_E= saved[0,3]
p_tot=1
s_tot=0
sl_tot=0
un_tot =0
count =0

pands = np.zeros(shape=(0,4))

for x in range(1,len(saved)):
    
    peak2=saved[x,1]
    new_station = saved[x,2]
    new_E= saved[x,3]
    
    if last_station == new_station:
        if saved[x,0] - last_p < 10*60 :
#            if peak1>peak2:
            if last_E>new_E:
                Eper=(new_E/last_E)*100
                s_tot += 1
                pands = np.lib.pad(pands, ((0,1),(0,0)), 'constant', constant_values=(0))
                pands[count][0]=saved[x,0]
                pands[count][1]=last_p
                pands[count][2]=saved[x,2]
                pands[count][3]=Eper
                count +=1
            else:
                sl_tot += 1
               
            
        else:
            p_tot +=1
            last_p = saved[x,0]
            peak1=saved[x,1]
            last_station=saved[x,2]
            last_E= saved[x,3]
    else:
        if saved[x,0] - last_p < 10*60 :
            un_tot += 1
            if last_E>new_E:
                Eper=(new_E/last_E)*100
                s_tot += 1
                pands = np.lib.pad(pands, ((0,1),(0,0)), 'constant', constant_values=(0))
                pands[count][0]=saved[x,0]
                pands[count][1]=last_p
                pands[count][2]=saved[x,2]
                pands[count][3]=Eper
                count +=1
            else:
                sl_tot += 1
                
            
        else:
             p_tot +=1
             last_p = saved[x,0]
             peak1=saved[x,1]
             last_station=saved[x,2]
             last_E= saved[x,3]

print("small secondaries =", s_tot)
print("larger secondaries =", sl_tot)
print('unmatched stations =' , un_tot)
print('Primaries =', p_tot)
print("total_secondaries =", sl_tot + s_tot )
print("total_events =", sl_tot + s_tot + p_tot )

#%% plot size hiatogram for all small secondaries
plt.figure(4001)
plt.hist(pands[:,3],bins=50)
plt.xlabel("Secondary Energy compared to Primary [%]")
plt.ylabel("Occurance [#]")
plt.title("Secondary Energies")

#%% RESULTS

# less energy for small secondary
#small secondaries = 2876
#larger secondaries = 1009
#unmatched stations = 130
#Primaries = 12497
#total_secondaries = 3885
#total_events = 16382

# 1/2 energy or less for small secondary
#small secondaries = 2560
#larger secondaries = 1325
#unmatched stations = 130
#Primaries = 12497
#total_secondaries = 3885
#total_events = 16382

# 1/3 energy or less for small secondary
#small secondaries = 2366
#larger secondaries = 1519
#unmatched stations = 130
#Primaries = 12497
#total_secondaries = 3885
#total_events = 16382

# 1/4 energy or less for small secondary
#small secondaries = 2228
#larger secondaries = 1657
#unmatched stations = 130
#Primaries = 12497
#total_secondaries = 3885
#total_events = 16382


#%% Get waveforms and see which data are not available online
#
#cha = 'HHZ' # CHANNEL
#net = 'Z4'  # 
#loc = ''    # location, it depends mostly of which network you are in. 
#client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
#unavailable=[]
#for x in range(0,len(pands)):
#    if pands[x,2] == 1:
#        sta='LB01'
#    if pands[x,2] == 2:
#        sta='LB02'
#    if pands[x,2] == 3:
#        sta='LB03'
#    if pands[x,2] == 4:
#        sta='LB04'
#    if pands[x,2] == 5:
#        sta='LB05'
#    if pands[x,2] == 6:
#        sta='LB06'
#    if pands[x,2] == 7:
#        sta='LS01'
#    if pands[x,2] == 8:
#        sta='LS02'
#    if pands[x,2] == 9:
#        sta='LS03'
#    if pands[x,2] == 10:
#        sta='LS04'
#    if pands[x,2] == 11:
#        sta='LS05'
#    if pands[x,2] == 12:
#        sta='LS06'
#    
#
#    arrival2=pands[x,0]
#    arrival1=pands[x,1]
#    
#    try:
#        t1 = UTCDateTime(arrival1)- 60 #the format is year:day_of_the_year:month
#        t2 =  UTCDateTime(arrival2) + 120
#        st= Stream()
#        st = client.get_waveforms(net, sta, '', cha, t1 , t2)
#        st.detrend(type='linear')
#        st.detrend(type='demean')
#        st.filter(type='bandpass',freqmin=0.5, freqmax=6)
#        st.plot(color='b',starttime=t1, endtime=t2)
#    except:
#        print("data not available at time: ",UTCDateTime(arrival1))
#        unavailable.append(arrival1)
#    
        
    
    




                   