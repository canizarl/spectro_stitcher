# Plots dynamic spectra from h5 files with time axes
# Slightly modified to return the data to a different script for plotting

import sys
import h5py
import tkinter
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as ndimage
import matplotlib as mpl
from  pdb import set_trace
from datetime import datetime
from datetime import timedelta
from matplotlib import ticker
from matplotlib import dates
from matplotlib.ticker import MaxNLocator
from pylab import imshow,xlabel,ylabel,title,close

mpl.rcParams.update({'font.size': 10})

def backsub(data, percentile=1.0):

    # Get time slices with standard devs in the bottom nth percentile.
    # Get average spectra from these time slices.
    # Devide through by this average spec.
    # Expects (row, column)

    print('Performing background subtraction.')
    # data = np.log10(data)
    data[np.where(np.isinf(data)==True)] = 0.0
    data_std = np.std(data, axis=0)
    data_std = data_std[np.nonzero(data_std)]
    min_std_indices = np.where( data_std < np.percentile(data_std, percentile) )[0]
    min_std_spec = data[:, min_std_indices]
    min_std_spec = np.mean(min_std_spec, axis=1)
    data = np.transpose(np.divide( np.transpose(data), min_std_spec))
    print('Background subtraction finished.')

    #Alternative: Normalizing frequency channel responses using median of values.
        #for sb in np.arange(data.shape[0]):
        #       data[sb, :] = data[sb, :]/np.mean(data[sb, :])

    return data


def plot_spectro(file, t0plot, t1plot, downsample=1):
    
    #--------------------------------------#
    # Note downsampling allows a bigger time range to be read into RAM. 
    # But it increases the computation time significantly. On an average
    # machine it may be best to keep full sampling (default), but limit
    # your plot to 15 minute blocks.

    runtime0=datetime.now()
    f = h5py.File( file, 'r' )
    t_lines = np.shape(f['SUB_ARRAY_POINTING_000/BEAM_000/STOKES_0'])[0]
    f_lines = np.shape(f['SUB_ARRAY_POINTING_000/BEAM_000/STOKES_0'])[1]
    target = f.attrs['TARGETS'][0]
    print('Target: %s' %(target))


    #---------------------------------#
    #  Sort out the time arrays 
    #
    obs_start = f.attrs['OBSERVATION_START_UTC'].decode("utf-8") 
    obs_stop = f.attrs['OBSERVATION_END_UTC'].decode("utf-8") 
    obs_start = datetime.strptime(obs_start[0:-4], "%Y-%m-%dT%H:%M:%S.%f")
    obs_stop = datetime.strptime(obs_stop[0:-4], "%Y-%m-%dT%H:%M:%S.%f")

    total_duration = f.attrs['TOTAL_INTEGRATION_TIME'] #in seconds
    tres = total_duration/t_lines

    t0sec = t0plot - obs_start
    t1sec = t1plot - obs_start
    t0index = int( t0sec.seconds/tres )
    t1index = int( t1sec.seconds/tres )

    print('Observation start time: %s' %(obs_start))
    print('Observation stop time: %s' %(obs_stop))
    print('Total observation duration %s seconds.' %(total_duration))
    print('Time resolution: %s seconds.' %(tres*downsample))


    tim_mpl = [dates.date2num(obs_start + timedelta(seconds=tres*(i+t0index))) for i in range(t1index-t0index)]
    tim_mpl = tim_mpl[::downsample]

    #----------------------------------#
    #  Sort out the frequency arrays
    #
    start_freq = f.attrs['OBSERVATION_FREQUENCY_MIN']  #in MHz
    end_freq = f.attrs['OBSERVATION_FREQUENCY_MAX'] 
    fres = (end_freq - start_freq)/f_lines #in MHz
    freq = np.linspace(start_freq, end_freq, f_lines)
    print('Frequency resolution: %s MHz.' %(fres*downsample))

    #-------------------------------------#
    #  Downsample and background subtract
    #
    gigabytes = 4*f_lines*(t1index-t0index)/downsample/1e9 # Data is 32 bit (4 byte).
    print('Reading in %s GB of data.' %(gigabytes))
    #import pdb
    #pdb.set_trace()
    data = f['SUB_ARRAY_POINTING_000/BEAM_000/STOKES_0'][t0index:t1index:downsample,::downsample]
    data = np.transpose(data)
    data = backsub(data)

    #----------------------------------#
    #     High pass filter
    #
    #data_lp = ndimage.gaussian_filter(data, sigma=(10, 10))
    #data = data - data_lp

    xyshape = np.shape(data)
    print('Array size: %s Mpixels.' %(xyshape[0]*xyshape[1]/1e6))
    #----------------------------------#
    #    Plot the spectrogram
    #
    
    # plt.figure(0, figsize=(12,6))
    # imshow(data[::-1,::], aspect='auto',
    #         extent=(tim_mpl[0], tim_mpl[-1], start_freq, end_freq),
    #         cmap=plt.get_cmap('plasma'),
    #         vmin=np.percentile(data, 30.0),
    #         vmax=np.percentile(data, 97.9)
    #         )
    
    # yymmdd=f.attrs['OBSERVATION_END_UTC'][0:10]
    # xlabel('Time (UT)')
    # ylabel('Frequency (MHz)')
    # title('LOFAR %s Beam 000 %s' %(yymmdd,target))
    # np.save('dynamic_spec_hba.npy',data)
    # extent=[tim_mpl[0], tim_mpl[-1], end_freq, start_freq]
    # np.save('extent_hba.npy',extent)
    # ax = plt.gca()
    # ax.set_yscale('log')
    # ax.xaxis_date()
    # ax.xaxis.set_major_locator(dates.MinuteLocator(interval=5))
    # ax.xaxis.set_minor_locator(dates.SecondLocator(interval=60))
    # ax.xaxis.set_major_formatter(dates.DateFormatter('%H:%M:%S'))

    # runtime1=datetime.now()
    # runtime = runtime1 - runtime0
    # print('Execution time: %s seconds' %(runtime.seconds))
    # plt.show()
    # #plt.savefig('L701913_TAB_2019-04-13_Beam_000.png')
    # #close()
    return data[::-1,::],tim_mpl[0], tim_mpl[-1], start_freq, end_freq

    # freq = np.arange(start_freq,end_freq,)
    #return data[::-1,::],tim_mpl,start_freq, end_freq

if __name__=='__main__':

    filename = sys.argv[1]
    
    # Define plot time range. Careful here.
    # For 0.01s and 12 kHz resolution, a 15 minutes spectrogram is 1.7 GB.
    # A large time range may eat your RAM!
    t0 = datetime(2019, 4, 9, 12, 30, 0)
    t1 = datetime(2019, 4, 9, 13, 00, 0)
    plot_spectro(filename, t0, t1)
