import numpy as np


print(" ")
print(" ----------------------------- ")
print(" ")
print("SWAVES Report: ")


swaves_data = [] 
with open('swaves_average_20190409_a_lfr.txt') as fobj: 
    for line in fobj: 
        row = line.split() 
        swaves_data.append(row[:]) 


# time array
# swaves_epoch = 

swaves_freq_l = swaves_data[0]
swaves_back_l = swaves_data[1]
data_buffer = []
for i in range(2,len(swaves_data)-1): 
    data_buffer.append(swaves_data[i][1:]) 

data_buffer = np.array(data_buffer)
new_data = []
for i in range(0,len(data_buffer)-1): 
    new_data.append(data_buffer[i].astype('float64')) 

swaves_data_l = np.array(new_data)


swaves_data = [] 
with open('swaves_average_20190409_a_hfr.txt') as fobj: 
    for line in fobj: 
        row = line.split() 
        swaves_data.append(row[:]) 


swaves_freq_h = swaves_data[0]
swaves_back_h = swaves_data[1]
data_buffer = []
for i in range(2,len(swaves_data)-1): 
    data_buffer.append(swaves_data[i][1:]) 

data_buffer = np.array(data_buffer)
new_data = []
for i in range(0,len(data_buffer)-1): 
    new_data.append(data_buffer[i].astype('float64')) 

    
swaves_data_h = np.array(new_data)



import matplotlib.pyplot as plt
from plot_beam_spectro import backsub

data = backsub(swaves_data_h.T)

plt.figure()
plt.imshow(swaves_data_h.T, aspect='auto')  
plt.show()












print("End SWAVES report.")
print(" ")
print(" ----------------------------- ")
print(" ")