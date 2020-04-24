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
   ¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡   Denotes information


"""




import numpy as np
from datetime import datetime
from datetime import timedelta
from psp_dataprep import data_spectro
import os
import plot_beam_spectro as pbs
from matplotlib import dates
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import math



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


if __name__=='__main__':
    print(f" RUNNING spectro_stitcher")
    
    day = "09"
    month = "04"
    year = "2019"

    file_lofar = "L701167_SAP000_B000_S0_P000_bf.h5"
    # times for lofar
    t0 = datetime(2019, 4, 9, 12, 30, 0)
    t1 = datetime(2019, 4, 9, 13, 00, 0)
    
    # for psp 
    date_open = mydate(year,month,day)

    

    ## import data
    #import PSP h data
    h_data = load_PSP_data(date_open,"h")
    #import PSP l data
    l_data = load_PSP_data(date_open,"l")

    dt = timedelta(seconds=(6.5)*60)             # time shift due to speed of light  (manual )


    #Data range:
    psp_time_range_indices_h = np.where((h_data.epoch>=t0-dt) & (h_data.epoch<=t1-dt))
    psp_time_range_indices_l = np.where((l_data.epoch>=t0-dt) & (l_data.epoch<=t1-dt))
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


    #import lofar data
    print(" ")
    print(" ----------------------------- ")
    print(" ")
    print("LOFAR Report: ")
    data_LOFAR, start_time_LOFAR, end_time_LOFAR, start_freq_LOFAR, end_freq_LOFAR = pbs.plot_spectro(file_lofar, t0, t1, downsample=1)
    #data_LOFAR, time_LOFAR, start_freq_LOFAR, end_freq_LOFAR = pbs.plot_spectro(file_lofar, t0, t1, downsample=1)

    print("End LOFAR report.")
    print(" ")
    print(" ----------------------------- ")
    print(" ")


    data_LOFAR = data_LOFAR[::-1]


    """ ----------------------------------------------- """
                        ## PLOT
    """ ----------------------------------------------- """

    
    f, axarr = plt.subplots(4,1,gridspec_kw={'height_ratios': [1.465, 1.166,0.08, 0.8]}) 
    f.subplots_adjust(left=0.15, bottom=0.13, right=0.90, top=0.95, wspace=0.35, hspace=0)
    
    
    axarr[0].set_title(f"PSP and LOFAR {year}  -  {month}  -  {day}")
    axarr[2].set_xlabel("Time (UT)")
    axarr[1].set_ylabel("Frequency [MHz] ")


    # PSP 
    # clipping low frequency
    ndi = np.where(l_data.freq > 1E5)      # new data indices
    l_data.data = l_data.data[:,ndi[0]]
    l_data.freq = l_data.freq[ndi[0]]    


    v_minl = np.percentile(l_data.data, 1)
    v_maxl = np.percentile(l_data.data, 99.9)

    v_minh = np.percentile(h_data.data, 1)
    v_maxh = np.percentile(h_data.data, 99)

    x_lims = [dates.date2num(h_data.epoch[0]), dates.date2num(h_data.epoch[-1])]

    hflim = [h_data.freq[0]/1E6, h_data.freq[-1]/1E6]
    lflim = [l_data.freq[0]/1E6, l_data.freq[-1]/1E6]



    
    axarr[1].set_yscale('log')
    imh = axarr[1].pcolormesh(dates.date2num(h_data.epoch), h_data.freq/1E6 , h_data.data.T,vmin=v_minh, vmax=v_maxh)
    axarr[1].invert_yaxis()
    axarr[1].set_xticks([])
    #plt.gca().invert_yaxis()
    # imh = axarr[1].imshow(h_data.data.T,
	# 	aspect = 'auto',
	# 	extent = [x_lims[0], x_lims[1], hflim[-1], hflim[0]],
	# 	vmin=v_minh, vmax=v_maxh)
    # axarr[1].set_xticks([])


    axarr[0].set_yscale('log')
    iml = axarr[0].pcolormesh(dates.date2num(l_data.epoch), l_data.freq/1E6 , l_data.data.T,vmin=v_minl, vmax=v_maxl)
    axarr[0].invert_yaxis()
    axarr[0].set_xticks([])
    #plt.gca().invert_yaxis()
    # iml = axarr[0].imshow(l_data.data.T,
	# 	aspect = 'auto',
	# 	extent = [x_lims[0], x_lims[1], lflim[-1], lflim[0]],
	# 	vmin=v_minl, vmax=v_maxl)
    # axarr[0].set_xticks([])


    # # LOFAR
    # freq_LOFAR = np.linspace(start_freq_LOFAR, end_freq_LOFAR, num=np.shape(data_LOFAR)[1])
    # axarr[2].set_yscale('log')
    # imlofar = axarr[2].pcolormesh(time_LOFAR, freq_LOFAR , data_LOFAR)
    

    axarr[2].set_xticks([])
    axarr[2].set_yticks([])

    imlofar = axarr[3].imshow(data_LOFAR, aspect='auto',
            extent=(start_time_LOFAR, end_time_LOFAR, end_freq_LOFAR,start_freq_LOFAR),
            vmin=np.percentile(data_LOFAR, 30.0),
            vmax=np.percentile(data_LOFAR, 97.9)
            )
    
    
 

    """  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   """
            # LABELS SHOULD BE IN LOG SCALE
    """  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   """
    #axarr[0].set_yscale('log')
    #axarr[1].set_yscale('log')
    axarr[3].set_yscale('log')
    

    # """ PSP LOW """
    # ylabs = freqlabelsmaker(l_data.freq)
    # labsrange = np.arange(lflim[0], lflim[-1], step=(lflim[-1]+lflim[0])/6)
    # print(f"labsrange :  {labsrange}")
    # axarr[0].set_yticks(labsrange)
    # axarr[0].set_yticklabels(ylabs)
    # # axarr[0].set_ylabel('Freq. [MHz]')

    # """ PSP HIGH """
    # ylabs = freqlabelsmaker(h_data.freq)
    # labsrange = np.arange(hflim[0], hflim[-1], step=(hflim[-1]+hflim[0])/6)
    # print(f"labsrange :  {labsrange}")
    # axarr[1].set_yticks(labsrange)
    # axarr[1].set_yticklabels(ylabs)
    # # axarr[0].set_ylabel('Freq. [MHz]')






    """ ¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡ """
                # time axis labels 
    """ ¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡ """    

    axarr[3].xaxis_date()
    axarr[3].xaxis.set_major_locator(dates.MinuteLocator(interval=5))
    axarr[3].xaxis.set_minor_locator(dates.SecondLocator(interval=60))
    axarr[3].xaxis.set_major_formatter(dates.DateFormatter('%H:%M:%S'))
    plt.show()


    




    



