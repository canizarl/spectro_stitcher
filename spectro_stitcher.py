""" 

    Spectro Stitcher
    Displays dynamic spectra from different instruments

    Supported instruments:
        Fields (PSP) --- Needs psp_dataprep.py
        LOFAR        --- Needs plot_beam_spectro.py

    
    Planned:
        Waves (Solar Orbiter)





    Some notes on how I code.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   Denotes to be changed
   -----------------------------------------   Denotes module heading to be turned into function
   iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii   Denotes information


"""

import numpy as np
from datetime import datetime
from datetime import timedelta
from psp_dataprep import data_spectro
import os
import plot_beam_spectro as pbs
from matplotlib import dates
from matplotlib import lines
import matplotlib.gridspec as gs
import matplotlib.pyplot as plt
import math

from scipy.io import readsav
from spacepy import pycdf
from sunpy.timeseries import TimeSeries
from sunpy.time import TimeRange, parse_time
from sunpy.net import hek, Fido, attrs as a

class mydate:
    def __init__(self,year,month,day):
        self.year = year
        self.month = month
        self.day = day


def load_PSP_data(date,band):
    '''
    Loads data from txt files
    

    Parameters
    ----------
    fnames : list of three strings 
        fname[0]: path and filename to data text file
        fname[1]: path and filename to frequencies text file
        fname[2]: path and filename to time text file
        NOT YET IMPLEMENTED

    Returns
    -------
    data : 2D numpy array of floating points
        Dynamic spectra data
    freq : 1D numpy array of floating points
        Frequencies of the dynamic spectra data
    time : 1D numpy array of floating points
        Time values of the dynamic spectra

    '''

    yearS = date.year
    monthS = date.month
    dayS = date.day 
    cwd = os.getcwd()
    
    directory_extracted =  cwd+"/ExtractedData" 
    directory_year = directory_extracted + "/"+yearS
    directory_month = directory_year + "/"+monthS
    directory = directory_month + "/"+dayS

    fname_data = directory+"/PSP_" + yearS + "_" + monthS + "_" + dayS+"_data_"+band+".txt"
    fname_freq = directory+"/PSP_" + yearS + "_" + monthS + "_" + dayS+"_freq_"+band+".txt"
    fname_time = directory+"/PSP_" + yearS + "_" + monthS + "_" + dayS+"_time_"+band+".txt"

    
    data = np.loadtxt(fname_data, delimiter=",")
    freq = np.loadtxt(fname_freq, delimiter=",")
    time = np.loadtxt(fname_time, delimiter=",")     

    epoch = []
    for i in range(0, len(time)):
        epoch.append(datetime.fromtimestamp(time[i]))

    epoch = np.array(epoch)
    my_spectro = data_spectro(data,epoch,freq)  
    return my_spectro


def freqlabelsmaker(freq):
    """
    No longer used
    """
    # Y axis -----------------------------------------------
    logfreq = np.zeros(len(freq))
    for i in range(0,len(freq)):
        logfreq[i] = math.log10(freq[i])

    ylab0 = str(round(freq[0]/1E6,4))
    # ylab1 = str(round(np.percentile(freq,10)/1E6,4))
    ylab2 = str(round(np.percentile(freq,20)/1E6,4))
    # ylab3 = str(round(np.percentile(freq,30)/1E6,4))
    ylab4 = str(round(np.percentile(freq,40)/1E6,4))
    # ylab5 = str(round(math.pow(10,logfreq[round(len(logfreq)/2)])/1E6,4))  
    ylab6 = str(round(np.percentile(freq,60)/1E6,4))
    # ylab7 = str(round(np.percentile(freq,70)/1E6,4))
    ylab8 = str(round(np.percentile(freq,80)/1E6,4))
    # ylab9 = str(round(np.percentile(freq,90)/1E6,4))
    ylab10 = str(round(freq[-1]/1E6,4))

    # ylabs = [ylab0,ylab1, ylab2, ylab3, ylab4, ylab5, ylab6,ylab7, ylab8, ylab9, ylab10]
    ylabs = [ylab0, ylab2, ylab4, ylab6, ylab8, ylab10]

    return ylabs


def manual_PSP_delay(t0,t1):
    ## align PSP with Earth UTC using PSP_Orbit.txt  file can be obtained at:
    ## https://sppgway.jhuapl.edu/PosCalc

    c = 299792458       # m/s
    au = 149597870700   # m
    dtEarth = au/c


    ## OPEN the file
    data = [] 
    with open('PSP_Orbit.txt') as fobj: 
        for line in fobj: 
            row = line.split() 
            data.append(row[:]) 


    epoch = []
    X = []
    Y = []
    Z = []
    for i in range(1,np.shape(data)[0]): 
        ## store epoch matrix
        epoch.append(datetime(int(data[i][0]),int(data[i][1]), int(data[i][2]),int(data[i][3]),int(data[i][4]), int(data[i][5])))
        ## store x y and z coordinates
        X.append(float(data[i][6]))
        Y.append(float(data[i][7]))
        Z.append(float(data[i][8]))

    ## pick a range within times 
    epoch = np.asarray(epoch)
    X = np.asarray(X)
    Y = np.asarray(Y)
    Z = np.asarray(Z)
    time_range_indices = np.where((epoch>=t0) & (epoch<=t1))
    X = X[time_range_indices[0]]
    Y = Y[time_range_indices[0]]
    Z = Z[time_range_indices[0]]



    ## square x y and z
    Xsquared = np.square(X)
    Ysquared = np.square(Y)
    Zsquared = np.square(Z)

    ## add them  w = x^2 + y^2 + z^2
    w = np.add(Xsquared,np.add(Ysquared,Zsquared))

    ## sqrt
    distance = np.sqrt(w)
    ## change to metres
    d_inmetres = np.multiply(distance,1000)

    d_ave = np.mean(d_inmetres)
    ## divide by c to get the time in seconds
    dtPSP = d_ave/c
   
    delay = dtEarth-dtPSP

    return delay


def get_goes(tr,t0,t1):
    results = Fido.search(a.Time(tr), a.Instrument('XRS'))
    files = Fido.fetch(results)
    goes = TimeSeries(files, source='XRS')

    year = []
    month = []
    day = []
    hour = []
    minute = []
    for items in goes.index:
        year.append(int(str(items)[0:4])) 
        month.append(int(str(items)[5:7])) 
        day.append(int(str(items)[8:10])) 
        hour.append(int(str(items)[11:13])) 
        minute.append(int(str(items)[14:16])) 

    times = []
    for i in range(0,len(year)):
        times.append(datetime(year[i],month[i],day[i],hour[i],minute[i])) 

    times = np.array(times)

    # t0 = datetime(2019, 4, 9, 10, 0)
    # t1 = datetime(2019, 4, 9, 15, 0)

    indices = np.where((times>=t0) & (times<=t1))[0]


    xrsa_data = goes.data['xrsa'][indices]
    xrsb_data = goes.data['xrsb'][indices]

    return xrsa_data, xrsb_data


def get_swaves(filename,t0,t1,percentile):
    data = readsav(filename)
    swaves_freqs = data['frequencies']
    swaves_back = data['back']
    swaves_spec = data['spectrum']
    

    observation_start_time = datetime(t0.year,t0.month,t0.day,0,0)
    time_res = timedelta(seconds=60)
    swaves_epoch = []
    swaves_epoch.append(observation_start_time)
    for i in range(0,(24*60)-1):
        swaves_epoch.append(swaves_epoch[-1]+time_res)
    if len(swaves_epoch) != (24*60):
        print(f"WARNING: swaves_epoch does not have the correct lenght. Length: {len(swaves_epoch)}")


    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    """       Backsub              """
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    # spec_buffer = []
    # swaves_back = swaves_back[::-1]
    # for i in range(0, swaves_spec.shape[0]):
    #     spec_buffer.append(np.multiply(swaves_spec[i,:],swaves_back))

    # swaves_spec = backsub(np.array(spec_buffer))


    swaves_spec = backsub(np.array(swaves_spec.T),percentile)
    swaves_spec = swaves_spec.T


    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    """       Time Range           """
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    swaves_epoch = np.array(swaves_epoch)
    swaves_spec = np.array(swaves_spec)
    swaves_freqs = np.array(swaves_freqs)

    # FIND WHERE IS THE DATA OF INTEREST
    swaves_time_range_indices = np.where((swaves_epoch>=t0) & (swaves_epoch<=t1))
    # KEEP ONLY THE DATA OF INTEREST I.E. TIME RANGES
    swaves_spec = swaves_spec[swaves_time_range_indices[0],:]  
    swaves_epoch = swaves_epoch[swaves_time_range_indices]


    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    """        Clipping            """
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    # clipping low frequency
    min_freq = 4E-1  # MHz    //  0 for no clipping
    min_freq = min_freq*1E3 # converting to kHz. OG data is in kHz.

    ndi = np.where(swaves_freqs >= min_freq)      # new data indices
    swaves_spec = swaves_spec[:,ndi[0]]
    swaves_freqs = swaves_freqs[ndi[0]]    

    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    """  Remove bad freqs          """
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    # # NO NEED
    # # average of frequency channels
    # spec_ave_freqs = np.mean(swaves_spec, axis=0)
    # swaves_spec = np.delete(swaves_spec,np.where(spec_ave_freqs==0),axis=1)
    # swaves_freqs = np.delete(swaves_freqs,np.where(spec_ave_freqs==0)) 
    return swaves_spec, swaves_epoch, swaves_freqs


def get_windwaves(myfile,t0,t1,percentile):
    """ data_from_CDF
     outputs dynamic spectra data from WIND WAVES CDF datafile.
     expects: wind waves h1  file. example: wi_h1_wav_20190409_v01.cdf
     inputs:
        myfile: string with location of file
        t0: start time 
        t1: end time
        percentile: percentile to be used in backsub routine 

    output:
        data: 2D numpy matrix of dynamic spectra
        epoch: 1D numpy array of datetimes
        freqs: 1D numpy array of the frequency channels. low band and high band merged. 

    """


    # OPEN FILE
    cdf = pycdf.CDF(myfile)

    # EXTRACT DATA
    l_data = cdf['E_VOLTAGE_RAD1']
    h_data = cdf['E_VOLTAGE_RAD2']

    # EXTRACT FREQUENCY ARRAYS
    l_freqs = cdf['Frequency_RAD1']
    h_freqs = cdf['Frequency_RAD2']

    # EXTRACT TIME ARRAYS
    epoch = cdf['Epoch']


    # CONVERT TO NUMPY ARRAYS
    l_data = np.array(l_data)
    h_data = np.array(h_data)

    epoch = np.array(epoch)


    data = np.concatenate((l_data,h_data), axis=1)

    freqs = []
    for items in l_freqs[:]:
        freqs.append(items) 
    for items in h_freqs[:]:
        freqs.append(items) 

    freqs = np.array(freqs)

    cdf.close()


    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    """       Backsub              """
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    data = backsub(data.T, percentile)
    data = data.T

    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    """       Time Range           """
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    wwaves_time_range_indices = np.where((epoch>=t0) & (epoch<=t1))
    data = data[wwaves_time_range_indices[0],:]  
    epoch = epoch[wwaves_time_range_indices]



    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    """    Freq Clipping            """
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    # clipping low frequency
    min_freq = 2E-1  # MHz    //  0 for no clipping
    min_freq = min_freq*1E3 # converting to kHz. OG data is in kHz.

    ndi = np.where(freqs >= min_freq)      # new data indices
    data = data[:,ndi[0]]
    freqs = freqs[ndi[0]]


    return data, epoch, freqs    


def backsub(data, percentile=1.0):
    # written by Eoin Carley
    # Get time slices with standard devs in the bottom nth percentile.
    # Get average spectra from these time slices.
    # Devide through by this average spec.
    # Expects (row, column)

    print('Performing background subtraction.')
    data = np.log10(data)
    data[np.where(np.isinf(data)==True)] = 1.0
    data[np.where(np.isnan(data)==True)] = 1.0
    
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






if __name__=='__main__':
    print(f" RUNNING spectro_stitcher")

    """ iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii """ 
    #             General                         #
    """ iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii """ 
    
    day = "09"
    month = "04"
    year = "2019"
    t0 = datetime(2019, 4, 9, 12, 30, 0)
    t1 = datetime(2019, 4, 9, 13, 00, 0)

    
    """ iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii """ 
    #         Wind WAVES                          #
    """ iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii """ 
    print(" ")
    print(" ----------------------------- ")
    print(" ")
    print("Wind WAVES Report: ")
    
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    """       IMPORT DATA          """
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #


    fname  = "combined_117886/wi_h1_wav_20190409_v01.cdf"
    data_w, epoch_w, freqs_w = get_windwaves(fname,t0,t1,percentile=50)





    print("End Wind WAVES report.")
    print(" ")
    print(" ----------------------------- ")
    print(" ")
   




    """ iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii """ 
    #         SWAVES                          #
    """ iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii """ 
    print(" ")
    print(" ----------------------------- ")
    print(" ")
    print("SWAVES Report: ")
    
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    """       IMPORT DATA          """
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    swaves_spec, swaves_epoch, swaves_freqs = get_swaves("swaves_average_20190409_a.sav",t0,t1, percentile = 100)
    


    print("End SWAVES report.")
    print(" ")
    print(" ----------------------------- ")
    print(" ")
   



    """ iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii """ 
    #             PSP                             #
    """ iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii """ 
    psp_delay = manual_PSP_delay(t0,t1)

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
    #  ADD KERNEL TO AUTOMATICALL EXTRACT PSP LOCATION AND CALCULATE THE TIME DELAY
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
    # orbiter_kernel = spicedata.get_kernel('psp')
    #print(psp_delay)
    # for psp 
    date_open = mydate(year,month,day)

    
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    """       IMPORT DATA          """
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    #import PSP h data FOR THE WHOLE DAY
    h_data = load_PSP_data(date_open,"h")
    #import PSP l data
    l_data = load_PSP_data(date_open,"l")

    #dt = timedelta(seconds=(6.5)*60)             # time shift due to speed of light  (manual )
    dt = timedelta(seconds=psp_delay)

    # GET DATA FOR THE TIME OF INTEREST
    #Data range:
    # FIND WHERE IS THE DATA OF INTEREST
    psp_time_range_indices_h = np.where((h_data.epoch>=t0-dt) & (h_data.epoch<=t1-dt))
    psp_time_range_indices_l = np.where((l_data.epoch>=t0-dt) & (l_data.epoch<=t1-dt))
    # KEEP ONLY THE DATA OF INTEREST I.E. TIME RANGES
    h_data.data = h_data.data[psp_time_range_indices_h[0],:]  
    h_data.epoch = h_data.epoch[psp_time_range_indices_h]
    l_data.data = l_data.data[psp_time_range_indices_l[0],:]  
    l_data.epoch = l_data.epoch[psp_time_range_indices_l]

    # PSP REPORT
    print(" ")
    print(" ----------------------------- ")
    print(" ")
    print("PSP Report: ")
    print(f"Observation start time: {h_data.epoch[0]}" )
    print(f"Observation stop time: {h_data.epoch[-1]}" )


    print("End PSP report.")
    print(" ")
    print(" ----------------------------- ")
    print(" ")


    """ iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii """ 
    #             LOFAR                       #
    """ iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii """ 
    file_lofar = "L701167_SAP000_B000_S0_P000_bf.h5"

    #import lofar data
    print(" ")
    print(" ----------------------------- ")
    print(" ")
    print("LOFAR Report: ")
    
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    """       IMPORT DATA          """
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    data_LOFAR, start_time_LOFAR, end_time_LOFAR, start_freq_LOFAR, end_freq_LOFAR = pbs.plot_spectro(file_lofar, t0, t1, downsample=1)

    print("End LOFAR report.")
    print(" ")
    print(" ----------------------------- ")
    print(" ")


    data_LOFAR = data_LOFAR[::-1]



    """ iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii """ 
    #         GOES                            #
    """ iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii """ 
    # GOES REPORT
    print(" ")
    print(" ----------------------------- ")
    print(" ")
    print("GOES Report: ")

    goes_t0 = datetime(2019,4,9,10,0)
    goes_t1 = datetime(2019,4,9,15,0)
    tr = TimeRange([f'{goes_t0.year}-{goes_t0.month:02}-{goes_t0.day:02} {goes_t0.hour:02}:{goes_t0.minute:02}',
             f'{goes_t1.year}-{goes_t1.month:02}-{goes_t1.day:02} {goes_t1.hour:02}:{goes_t1.minute:02}'])
    
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    """       IMPORT DATA          """
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    xrsa_data, xrsb_data = get_goes(tr,datetime(int(year),int(month),int(day),10,0),datetime(int(year),int(month),int(day),15,0))
    

    print(f" {tr}" )

    print("End GOES report.")
    print(" ")
    print(" ----------------------------- ")
    print(" ")
   


    """ ----------------------------------------------- """
                        ## PLOTS
    """ ----------------------------------------------- """

    """ iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii """
    #             General                     #
    """ iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii """

    # General settings 
    # SAVE plot? 1 = yes  //  0 = No
    savefigure = 1
    plt.rcParams.update({'font.size': 12})

    # ORDER OF PLOTS FROM TOP TO BOTTOM
    instruments = ['wwaves','swaves','psp','lofar','goes']

    displays = {}   # create dictionary for displaying instruments
    # populates dictionary with instrument names
    # this is useful for calling axarr later on
    for i in range(0, len(instruments)):
        displays[instruments[i]] = i

    # CONTROLS THE GLOBAL MARGINS LEFT AND RIGHT OF THE PLOTS
    leftmargin = 0.1
    rightmargin = 0.95

    textmargin = 0.02





    # INITIALISING THE SUBPLOTS
    f, axarr = plt.subplots(len(displays),1,gridspec_kw={'height_ratios': np.full((1, len(displays)), 1)[0]},figsize=(9,13)) 
    
    # THIS CONTROLS GLOBAL SUBPLOTS APPEARANCE DOING NOTHING NOW BECAUSE USING LOCAL SUBPLOTS
    # f.subplots_adjust(left=leftmargin, bottom=0.13, right=rightmargin, top=0.95, wspace=0.35, hspace=0.26)


    # OBJECTS FOR SUBPLOTS 
    # object 1 for goes alone
    gs1 = gs.GridSpec(nrows = 1, ncols = 1)
    gs1.update(left= leftmargin, right=rightmargin,top = 0.20 , bottom = 0.05,  hspace=0.0)

    # object 2 for psp and lofar 
    gs2 = gs.GridSpec(nrows = 3, ncols = 1)
    gs2.update(left=leftmargin, right=rightmargin,top = 0.49 , bottom = 0.25,  hspace=0.0)


    # object 3 for swaves
    gs3 = gs.GridSpec(nrows = 1, ncols = 1)
    gs3.update(left=leftmargin, right=rightmargin,top = 0.72 , bottom = 0.54,  hspace=0.0)

    # object 4 for wwaves
    gs4 = gs.GridSpec(nrows = 1, ncols = 1)
    gs4.update(left=leftmargin, right=rightmargin,top = 0.98 , bottom = 0.77,  hspace=0.0)



    colormap = 'plasma'

    """ iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii """
    #         Wind WAVES                          #
    """ iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii """

    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    """         LEVELS             """
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    v_min_wwaves = np.percentile(data_w,50)
    v_max_wwaves = np.percentile(data_w,96)



    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    """       Plotting             """
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    axarr[displays['wwaves']] = plt.subplot(gs4[0, 0:3])
    plt.pcolormesh(dates.date2num(epoch_w), freqs_w/1E3, data_w.T,
        vmin=v_min_wwaves, vmax=v_max_wwaves,
        cmap = colormap)


    # TITLE 
    axarr[displays['wwaves']].text(textmargin, 0.9,'Wind WAVES',
        horizontalalignment='left',
        fontsize = 'large',
        color = 'w',
        fontweight = 'bold',
        transform=axarr[displays['wwaves']].transAxes)
    axarr[displays['wwaves']].text(rightmargin - textmargin, 0.8,'(a)',
        horizontalalignment='left',
        fontsize = 'large',
        color = 'w',
        fontweight = 'bold',
        transform=axarr[displays['wwaves']].transAxes)
    axarr[displays['wwaves']].set_yscale('log')
    axarr[displays['wwaves']].invert_yaxis()

    axarr[displays['wwaves']].set_ylabel("Frequency [MHz] ")
    axarr[displays['wwaves']].set_xlabel(f"TIME / {year}  -  {month}  -  {day}")


    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    """      time axis labels      """
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    axarr[displays['wwaves']].xaxis_date()
    axarr[displays['wwaves']].xaxis.set_major_locator(dates.MinuteLocator(interval=5))
    axarr[displays['wwaves']].xaxis.set_minor_locator(dates.SecondLocator(interval=60))
    axarr[displays['wwaves']].xaxis.set_major_formatter(dates.DateFormatter('%H:%M:%S'))
        





    """ iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii """
    #         SWAVES                              #
    """ iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii """

    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    """         LEVELS             """
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    v_min_swaves = np.percentile(swaves_spec.data, 20)
    v_max_swaves = np.percentile(swaves_spec.data, 99)

    
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    """       Plotting             """
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    axarr[displays['swaves']] = plt.subplot(gs3[0, 0:3])
    axarr[displays['swaves']].set_yscale('log')
    im_swaves = axarr[displays['swaves']].pcolormesh(dates.date2num(swaves_epoch), swaves_freqs/1E3, swaves_spec.T,
        vmin=v_min_swaves, vmax=v_max_swaves,
        cmap = colormap)
    axarr[displays['swaves']].invert_yaxis()
    axarr[displays['swaves']].set_ylabel("Frequency [MHz] ")

    # TITLE
    axarr[displays['swaves']].text(textmargin,.9,'S/WAVES',
        horizontalalignment='left',
        fontsize = 'large',
        color = 'w',
        fontweight = 'bold',
        transform=axarr[displays['swaves']].transAxes)
    
    axarr[displays['swaves']].text(rightmargin - textmargin, 0.8,'(b)',
        horizontalalignment='left',
        fontsize = 'large',
        color = 'w',
        fontweight = 'bold',
        transform=axarr[displays['swaves']].transAxes)

    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    """      time axis labels      """
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #

    axarr[displays['swaves']].xaxis_date()
    axarr[displays['swaves']].xaxis.set_major_locator(dates.MinuteLocator(interval=5))
    axarr[displays['swaves']].xaxis.set_minor_locator(dates.SecondLocator(interval=60))
    axarr[displays['swaves']].xaxis.set_major_formatter(dates.DateFormatter('%H:%M:%S'))
    axarr[displays['swaves']].set_xlabel(f"TIME / {year}  -  {month}  -  {day}")
        



    """ iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii """
    #         PSP                                 #
    """ iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii """
    # PSP 
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    """        Clipping            """
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    # clipping low frequency
    ndi = np.where(l_data.freq > 8E5)      # new data indices
    l_data.data = l_data.data[:,ndi[0]]
    l_data.freq = l_data.freq[ndi[0]]    

    
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    """  Joining HFR and LFR       """
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    psp_data = data_spectro([],[],[])
    psp_data.data = np.hstack((l_data.data,h_data.data))       
    psp_data.epoch = l_data.epoch
    psp_data.freq = np.concatenate((l_data.freq, h_data.freq))

    # frequency limits
    flims = [psp_data.freq[0], 2E6]


    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    """         LEVELS             """
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    v_min_psp = np.percentile(psp_data.data, 1)
    v_max_psp = np.percentile(psp_data.data, 99)
    
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    """         Extend             """
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    start_time_psp = psp_data.epoch[0]
    end_time_psp = psp_data.epoch[-1]
    end_freq_psp = 2E6
    start_freq_psp = psp_data.freq[0]


    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    """       Plotting             """
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    axarr[displays['psp']] = plt.subplot(gs2[0:2, 0:3])
    axarr[displays['psp']].set_yscale('log')
    im_psp = axarr[displays['psp']].pcolormesh(dates.date2num(psp_data.epoch), psp_data.freq/1E6, psp_data.data.T,
        vmin=v_min_psp, vmax=v_max_psp,
        cmap = colormap)
    axarr[displays['psp']].invert_yaxis()
    axarr[displays['psp']].set_xticks([])
    axarr[displays['psp']].set_ylabel("Frequency [MHz] ")

    # TITLE
    axarr[displays['psp']].text(textmargin,.9,'PSP/FIELDS',
        horizontalalignment='left',
        fontsize = 'large',
        color = 'w',
        fontweight = 'bold',
        transform=axarr[displays['psp']].transAxes)

    # Additional info
    axarr[displays['psp']].text(textmargin,.8,f"T+{psp_delay:0.02f}s",
        horizontalalignment='left',
        fontsize = 'small',
        color = 'w',
        fontweight = 'bold',
        transform=axarr[displays['psp']].transAxes)

    axarr[displays['psp']].text(rightmargin - textmargin, 0.8,'(c)',
        horizontalalignment='left',
        fontsize = 'large',
        color = 'w',
        fontweight = 'bold',
        transform=axarr[displays['psp']].transAxes)
    
    # axarr[displays['psp']].yaxis.get_offset_text().set_visible(False)
    # axarr[displays['psp']].axis([start_time_psp, end_time_psp, end_freq_psp,start_freq_psp])


    """ iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii """
    #         LOFAR                               #
    """ iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii """

    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    """       Plotting             """
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
    # CHANGE TO PCOLORMESH to see if solves y axis problems 
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
    axarr[displays['lofar']] = plt.subplot(gs2[2, 0:3])
    imlofar = axarr[displays['lofar']].imshow(data_LOFAR, aspect='auto',
            extent=(start_time_LOFAR, end_time_LOFAR, end_freq_LOFAR,start_freq_LOFAR),
            vmin=np.percentile(data_LOFAR, 30.0),
            vmax=np.percentile(data_LOFAR, 97.9),
            cmap = colormap) 

    # TITLE 
    axarr[displays['lofar']].text(textmargin,.8,'LOFAR',
        horizontalalignment='left',
        fontsize = 'large',
        color = 'w',
        fontweight = 'bold',
        transform=axarr[displays['lofar']].transAxes)

    axarr[displays['lofar']].text(rightmargin-textmargin,.8,'(d)',
        horizontalalignment='left',
        fontsize = 'large',
        color = 'w',
        fontweight = 'bold',
        transform=axarr[displays['lofar']].transAxes)



    # PLOT SETTINGS
    axarr[displays['lofar']].set_yscale('log')
    axarr[displays['lofar']].set_ylim((1E2,2E1))
    axarr[displays['lofar']].set_xlabel(f"TIME / {year}  -  {month}  -  {day}")


    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    """      time axis labels      """
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #

    axarr[displays['lofar']].xaxis_date()
    axarr[displays['lofar']].xaxis.set_major_locator(dates.MinuteLocator(interval=5))
    axarr[displays['lofar']].xaxis.set_minor_locator(dates.SecondLocator(interval=60))
    axarr[displays['lofar']].xaxis.set_major_formatter(dates.DateFormatter('%H:%M:%S'))
    




    """ iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii """
    #         GOES                                #
    """ iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii """
    
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    """       Plotting             """
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    axarr[displays['goes']] = plt.subplot(gs1[0, 0:3])
    goes_y_min = 1E-9
    goes_y_max = 1E-5
    axarr[displays['goes']].set_ylabel('Watts $m^{-2}$')
    axarr[displays['goes']].set_yscale('log')
    axarr[displays['goes']].set_ylim((goes_y_min,goes_y_max))
    axarr[displays['goes']].set_xlim((goes_t0,goes_t1))
    axarr[displays['goes']].plot(xrsb_data,'r-',label = r'GOES 1.0 -- 8.0 $\AA$')
    axarr[displays['goes']].plot(xrsa_data,'b-',label = r'GOES 0.5 -- 4.0 $\AA$')
 
    axarr[displays['goes']].legend(loc='upper left')
    axarr[displays['goes']].axvline(x=t0, color='r', linestyle='--')
    axarr[displays['goes']].axvline(x=t1, color='r', linestyle='--')

    axarr[displays['goes']].grid(axis = 'y')



    ax_right_goes = axarr[displays['goes']].twinx()
    ax_right_goes.set_ylim((goes_y_min,goes_y_max))
    ax_right_goes.plot(xrsb_data,'r-')
    ax_right_goes.plot(xrsa_data,'b-')
    ax_right_goes.set_yscale('log')
    ax_right_goes.set_yticks([1E-8,1E-7,1E-6,1E-5])
    ax_right_goes.tick_params(axis='both', which='both', length=0)
    ax_right_goes.set_yticklabels(['A','B','C','M'])
    

    axarr[displays['goes']].xaxis.set_major_locator(dates.HourLocator(interval=1))
    axarr[displays['goes']].xaxis.set_minor_locator(dates.MinuteLocator(interval=15))
    axarr[displays['goes']].xaxis.set_major_formatter(dates.DateFormatter('%H:%M:%S'))
    axarr[displays['goes']].set_xlabel(f"TIME / {year}  -  {month}  -  {day}")




    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    """       LINES                """
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiii #

    l1 = lines.Line2D([0.61, rightmargin], [0.2, 0.25], transform=f.transFigure, figure=f, color='r',linestyle='--')
    l2 = lines.Line2D([0.53, leftmargin], [0.2, 0.25], transform=f.transFigure, figure=f, color='r',linestyle='--')
    f.lines.extend([l1, l2])




    if savefigure == 1:
        f.savefig('PSP_LOFAR_WAVES', dpi=300, format='png')


    plt.show()


    




    



