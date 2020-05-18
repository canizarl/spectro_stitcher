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
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import math


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












if __name__=='__main__':
    print(f" RUNNING spectro_stitcher")

    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    #             General                     #
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    
    day = "09"
    month = "04"
    year = "2019"
    t0 = datetime(2019, 4, 9, 12, 30, 0)
    t1 = datetime(2019, 4, 9, 13, 00, 0)



    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    #             PSP                         #
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    psp_delay = manual_PSP_delay(t0,t1)
    # orbiter_kernel = spicedata.get_kernel('psp')
    #print(psp_delay)
    # for psp 
    date_open = mydate(year,month,day)

    

    ## import data
    #import PSP h data
    h_data = load_PSP_data(date_open,"h")
    #import PSP l data
    l_data = load_PSP_data(date_open,"l")

    #dt = timedelta(seconds=(6.5)*60)             # time shift due to speed of light  (manual )
    dt = timedelta(seconds=psp_delay)

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


    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    #             LOFAR                       #
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    file_lofar = "L701167_SAP000_B000_S0_P000_bf.h5"

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



    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    #         GOES                            #
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    # GOES REPORT
    print(" ")
    print(" ----------------------------- ")
    print(" ")
    print("GOES Report: ")

    tr = TimeRange(['2019-04-09 10:00', '2019-04-09 15:00'])
    xrsa_data, xrsb_data = get_goes(tr,datetime(int(year),int(month),int(day),10,0),datetime(int(year),int(month),int(day),15,0))
    



    print(f" {tr}" )

    print("End GOES report.")
    print(" ")
    print(" ----------------------------- ")
    print(" ")
   



    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    #         SWAVES                          #
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    print(" ")
    print(" ----------------------------- ")
    print(" ")
    print("SWAVES Report: ")

    
    swaves_data = [] 
    with open('swaves_average_20190409_a_lfr.txt') as fobj: 
        for line in fobj: 
            row = line.split() 
            swaves_data.append(row[:]) 
    
    
    print("End SWAVES report.")
    print(" ")
    print(" ----------------------------- ")
    print(" ")
   

















    """ ----------------------------------------------- """
                        ## PLOT
    """ ----------------------------------------------- """

    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    #             General                     #
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    
    displays = {'psp':0, 'lofar':1,'goes':2} 


    # f = plt.figure(figsize=(10,10), constrained_layout=True) 
    # gs = f.add_gridspec(5,5)

    f, axarr = plt.subplots(len(displays),1,gridspec_kw={'height_ratios': [1, 0.5, 0.2]},figsize=(15,10)) 
    f.subplots_adjust(left=0.15, bottom=0.13, right=0.90, top=0.95, wspace=0.35, hspace=0.26)



    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    #         GOES                            #
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    goes_y_min = 1E-9
    goes_y_max = 1E-5
    axarr[displays['goes']].set_ylabel('Watts $m^{-2}$')
    axarr[displays['goes']].set_yscale('log')
    axarr[displays['goes']].set_ylim((goes_y_min,goes_y_max))
    axarr[displays['goes']].plot(xrsb_data,'r-',label = r'GOES 1.0 -- 8.0 $\AA$')
    axarr[displays['goes']].plot(xrsa_data,'b-',label = r'GOES 0.5 -- 4.0 $\AA$')
    axarr[displays['goes']].xaxis.set_major_formatter(dates.DateFormatter("%H-%M"))
    axarr[displays['goes']].legend(loc='upper left')
    axarr[displays['goes']].axvline(x=t0, color='r', linestyle='--')
    axarr[displays['goes']].axvline(x=t1, color='r', linestyle='--')
    axarr[displays['goes']].set_xlabel(f"TIME / {year}  -  {month}  -  {day}")
    # axarr[displays['goes']].set_title('GOES Xray Flux')


    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    #         PSP                             #
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    # # PSP 
    # # clipping low frequency
    # ndi = np.where(l_data.freq > 1E5)      # new data indices
    # l_data.data = l_data.data[:,ndi[0]]
    # l_data.freq = l_data.freq[ndi[0]]    


    # v_minl = np.percentile(l_data.data, 1)
    # v_maxl = np.percentile(l_data.data, 99.9)

    # v_minh = np.percentile(h_data.data, 1)
    # v_maxh = np.percentile(h_data.data, 99)

    # x_lims = [dates.date2num(h_data.epoch[0]), dates.date2num(h_data.epoch[-1])]

    # hflim = [h_data.freq[0]/1E6, h_data.freq[-1]/1E6]
    # lflim = [l_data.freq[0]/1E6, l_data.freq[-1]/1E6]



    
    # axarr[displays['psp_h']].set_yscale('log')
    # imh = axarr[displays['psp_h']].pcolormesh(dates.date2num(h_data.epoch), h_data.freq/1E6 , h_data.data.T,vmin=v_minh, vmax=v_maxh)
    # axarr[displays['psp_h']].invert_yaxis()
    # axarr[displays['psp_h']].set_xticks([])
    # axarr[displays['psp_h']].set_ylabel("Frequency [MHz] ")

    # #plt.gca().invert_yaxis()
    # # imh = axarr[displays['psp_h']].imshow(h_data.data.T,
	# # 	aspect = 'auto',
	# # 	extent = [x_lims[0], x_lims[1], hflim[-1], hflim[0]],
	# # 	vmin=v_minh, vmax=v_maxh)
    # # axarr[displays['psp_h']].set_xticks([])


    # axarr[displays['psp_l']].set_yscale('log')
    # iml = axarr[displays['psp_l']].pcolormesh(dates.date2num(l_data.epoch), l_data.freq/1E6 , l_data.data.T,vmin=v_minl, vmax=v_maxl)
    # axarr[displays['psp_l']].invert_yaxis()
    # axarr[displays['psp_l']].set_xticks([])
    # axarr[displays['psp_l']].set_ylabel("Frequency [MHz] ")

    # #plt.gca().invert_yaxis()
    # # iml = axarr[displays['psp_l']].imshow(l_data.data.T,
	# # 	aspect = 'auto',
	# # 	extent = [x_lims[0], x_lims[1], lflim[-1], lflim[0]],
	# # 	vmin=v_minl, vmax=v_maxl)
    # # axarr[displays['psp_l']].set_xticks([])


    # # # LOFAR
    # # freq_LOFAR = np.linspace(start_freq_LOFAR, end_freq_LOFAR, num=np.shape(data_LOFAR)[1])
    # # axarr[displays['lofar']].set_yscale('log')
    # # imlofar = axarr[displays['lofar']].pcolormesh(time_LOFAR, freq_LOFAR , data_LOFAR)
    




    # PSP 
    # clipping low frequency
    ndi = np.where(l_data.freq > 8E5)      # new data indices
    l_data.data = l_data.data[:,ndi[0]]
    l_data.freq = l_data.freq[ndi[0]]    

    # x_lims = [dates.date2num(h_data.epoch[0]), dates.date2num(h_data.epoch[-1])]

    # hflim = [h_data.freq[0]/1E6, h_data.freq[-1]/1E6]
    # lflim = [l_data.freq[0]/1E6, l_data.freq[-1]/1E6]

    
    
    psp_data = data_spectro([],[],[])
    psp_data.data = np.hstack((l_data.data,h_data.data))       
    psp_data.epoch = l_data.epoch
    psp_data.freq = np.concatenate((l_data.freq, h_data.freq))

    flims = [psp_data.freq[0], 2E6]

    v_min_psp = np.percentile(psp_data.data, 0.5)
    v_max_psp = np.percentile(psp_data.data, 99.9)
    
    start_time_psp = psp_data.epoch[0]
    end_time_psp = psp_data.epoch[-1]
    end_freq_psp = 2E6
    start_freq_psp = psp_data.freq[0]

    axarr[displays['psp']].set_yscale('log')
    iml = axarr[displays['psp']].pcolormesh(dates.date2num(psp_data.epoch), psp_data.freq/1E6, psp_data.data.T,
        vmin=v_min_psp, vmax=v_max_psp)
    axarr[displays['psp']].invert_yaxis()
    axarr[displays['psp']].set_xticks([])
    axarr[displays['psp']].set_ylabel("Frequency [MHz] ")
    #axarr[displays['psp']].set_title('PSP/FIELDS')
    axarr[displays['psp']].text(.1,.9,'PSP/FIELDS',
        horizontalalignment='center',
        fontsize = 'large',
        color = 'w',
        fontweight = 'bold',
        transform=axarr[displays['psp']].transAxes)

    axarr[displays['psp']].text(.1,.1,f"TIME Shifted {psp_delay:0.02f}s",
        horizontalalignment='center',
        fontsize = 'small',
        color = 'w',
        fontweight = 'bold',
        transform=axarr[displays['psp']].transAxes)
    # axarr[displays['psp']].axis([start_time_psp, end_time_psp, end_freq_psp,start_freq_psp])


    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    #         LOFAR                           #
    # iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii #
    imlofar = axarr[displays['lofar']].imshow(data_LOFAR, aspect='auto',
            extent=(start_time_LOFAR, end_time_LOFAR, end_freq_LOFAR,start_freq_LOFAR),
            vmin=np.percentile(data_LOFAR, 30.0),
            vmax=np.percentile(data_LOFAR, 97.9)
            ) 

    # axarr[displays['lofar']].set_title('LOFAR')
    axarr[displays['lofar']].text(.1,.8,'LOFAR',
    horizontalalignment='center',
    fontsize = 'large',
    color = 'w',
    fontweight = 'bold',
    transform=axarr[displays['lofar']].transAxes)
    axarr[displays['lofar']].set_yscale('log')
    axarr[displays['lofar']].set_ylim((1E2,2E1))
    axarr[displays['lofar']].set_ylabel("Frequency [MHz] ")
    axarr[displays['lofar']].set_xlabel(f"TIME / {year}  -  {month}  -  {day}")



    # axarr[displays['lofar']].set_yticks([])






    """ iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii """
                # time axis labels 
    """ iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii """    

    axarr[displays['lofar']].xaxis_date()
    axarr[displays['lofar']].xaxis.set_major_locator(dates.MinuteLocator(interval=5))
    axarr[displays['lofar']].xaxis.set_minor_locator(dates.SecondLocator(interval=60))
    axarr[displays['lofar']].xaxis.set_major_formatter(dates.DateFormatter('%H:%M:%S'))
    
    
    # axarr[displays['psp']].xaxis_date()
    # axarr[displays['psp']].xaxis.set_major_locator(dates.MinuteLocator(interval=5))
    # axarr[displays['psp']].xaxis.set_minor_locator(dates.SecondLocator(interval=60))
    # axarr[displays['psp']].xaxis.set_major_formatter(dates.DateFormatter('%H:%M:%S'))
    
    
    
    
    plt.show()


    




    



