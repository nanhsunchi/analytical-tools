import numpy as np
import matplotlib.pyplot as plt
from scipy import signal, interpolate

def movingavg(x, m, ncrt):
    # movingavg.m from REN-CHEIH LIEN's mymatlab library
    # x : original time series
    # m : number of points of averaging
    # ncrt : the least number of real number is needed to compute the mean
    x = x[:]
    if m%2 == 1:
        nlost = int((m-1)/2)
    else:
        nlost = int(m/2)
    
    # For 1D only (modified from the original)
    begpt = int(nlost)
    nrow = len(x)
    endpt = int(nrow - nlost);
    # print('nlost & begpt & endpt:',nlost,begpt,endpt)
    avgx = np.nan*np.ones(nrow,)
    for i in range(begpt,endpt):
        bin0 = i-nlost
        bin1 = i+nlost+1
        # print(np.sum(~np.isnan(x[bin0:bin1])))
        if np.sum(~np.isnan(x[bin0:bin1]))> ncrt:
            ax = np.nanmean(x[bin0:bin1])
            avgx[i] = ax
            # print(i,bin0,bin1)
    # print(avgx)
    return avgx


def lpass(x, delt, cutt, npole):
    # lpass.m From Ren-Chieh Lien's matlab library
    # NO XOPLOT input
    #             Low pass time series using Butterworth filter          
    #               xlpass = lpass(xoriginal, tinterval, cutofft, npole,[XO])
    #             EX : xlpass = lpass(x0, 25, 120, 2)
    #                  tinterval = 25 units (eg. secs, mins, hours, ...etc.)
    #                  cutoff time = 120 units
    #                  npole = 2 (poles used in Butterworth filter, default 2)
    #------------------------------------------------------------------------------
    cutt = cutt / 2
    ### CHECK THE INPUT ARGUMENT
    if np.sum(np.isnan(x)):
        # echo time series has NaN or Inf. Moving averages will be proceeded
        yl = movingavg(x, np.floor(cutt/delt),2)#,float('-inf'),float('inf'),np.nan)
        ff = []
        hl = []
    else:
        nptslost = np.ceil(cutt/delt)
    
    ### REΜΟVE THE MEAN FROM THE TIME SERIES FIRST
    originalx = x
    bad = np.isnan(x)
    ibad = np.where(bad)[0]
    good = ~np.isnan(x)
    igood = np.where(good)[0]
    npts = len(x)
    if (np.sum(bad) == 1) & (np.isnan(x[0]) | np.isnan(x[-1])):
        x = np.delete(x,bad); add = 'y'
    if np.sum(bad) > 1:
        if (min(ibad) == 0) & (sum(np.diff(np.diff(ibad))^2) == 0):
            x = np.delete(x,bad); add = 'y'
        if (max(ibad) == npts) & (sum(np.diff(np.diff(ibad))^2) == 0):
            x = np.delete(x,bad); add = 'y'
    x0 = x
    x0 = np.mean(x)
    x = x - x0

    ### SET UP THE SAMPLING AND CUTOFF FREQUENCY
    fs = 1/delt         # sampling frequency
    fcut = 1/cutt       # cuttoff frequency
    lband = fcut/fs     # normalized lpass band
    ### CONSTRUCT THE BUTTERWORTH FILTER
    time = np.array([i for i in range(0,len(x))])*delt
    # npole
    # lband
    [cl,ch] = signal.butter(npole,lband)
    ### CHECK THE FREQUECY RESPONSE FUNCTION
    nfreq = 128*2
    ff = fs*np.arange(0,nfreq)/(2*nfreq)
    hl = signal.freqz(cl,ch,nfreq)
    hl = hl*np.conj(hl)
    ### APPLYING THE FILTER FIRST
    y = signal.lfilter(cl,ch,x)
    ### REVERSING THE TIME SERIES AND APPLYING THE FILTER AGAIN
    ### TO REMOVE THE PHASE SHIFT DUE TO THE FILTERING
    z0 = y[::-1]
    z = signal.lfilter(cl,ch,z0)
    ### REVERSE THE TIME SERIES BACK AGAIN ADN ADD THE MEAN
    w = z[::-1]
    w1 = signal.filtfilt(cl,ch,x)
    xl = w + x0
    yl = w1 + x0

    npts = len(xl)
    add = 'y'
    print('len of x, xl, igood:',len(x), len(xl), len(igood))
    print(np.sum(np.isnan(x)),np.sum(np.isnan(xl)))
    if add == 'y':
        # x[igood] = xl
        x = xl
        xl = x
    
    return yl, ff, hl


def lpass_NaN(x, delt, cutt, npole):
    # lpass_NaN.m from REN-CHIEH LIEN'S matlab library
    #          Low pass time series using Butterworth filter
    #            xlpass = lpass(xoriginal, tinterval, cutofft, npole)
    #          EX : xlpass = lpass(x0, 25, 120, 2)
    #                  tinterval = 25 units (eg. secs, mins, hours, ...etc.)
    #                  cutoff time = 120 units
    #                  npole = 2 (poles used in Butterworth filter, default 2)
    #------------------------------------------------------------------------------
    ###             CHECK THE INPUT ARGUMENT
    yl = np.nan*np.ones(len(x))
    n = len(x);
    nbin = np.arange(1,n+1)
    bad = np.where( np.isnan(x) | np.isinf(x))[0]
    good = np.where( ~np.isnan(x) & ~np.isinf(x))[0]
    interpx = interpolate.interp1d(nbin[good],x[good],kind='linear',bounds_error=False,fill_value='extrapolate')(bad)
    x[bad] = interpx
    # bad_beg_end = np.where( np.isnan(x) )[0]
    # ngood = np.where( ~np.isnan(x) )[0]
    # if len(bad_beg_end)> 0:
    #     x = x(ngood)
    ###             REMOVE THE MEAN FROM THE TIME SERIES FIRST
    x0 = np.nanmean(x)
    x = x - x0
    # print(x0)
    # print(x)
    ###             SET UP THE SAMPLING AND CUTOFF FREQUENCY
    fs = 1/delt;                                   # sampling frequency
    fcut = 1/cutt;                                 # cutoff frequency
    lband = fcut/fs;                               # normalized lpass band
    ###             CONSTRUCT THE BUTTERWORTH FILTER
    [cl,ch] = signal.butter(npole,lband)
    w1 = signal.filtfilt(cl,ch,x)
    yl = w1 + x0
    # yl[np.where(~np.isnan(x))] = w1 + x0
    # yl[bad] = np.nan*np.ones(len(bad))
    ###             CHECK THE FREQUENCY RESPONSE FUNCTION
    nfreq = 128*2;
    ff = fs*np.arange(0,nfreq)/(2*nfreq)
    hl = signal.freqz(cl,ch,nfreq)
    hl = hl*np.conj(hl)

    return yl,ff,hl

def bandpass(x,delt,lowcutt,higcutt,npole):
    # From Ren-Chieh's lienbandpass.m
    # mybandpass.m calls lienbandpass.m. 
    # The lienbandpass.m is the main bandpass filter
    # BANDX = BANDPASS(X,DT,LCTIME,HCTIME,NPOLE,plotind,FILTERTYPE)
    #               BANDX : BANDPASSED TIME SERIES
    #                X : ORIGINAL TIME SERIES
    #           LCTIME : LOWBOUND CUTOFF TIME     
    #           HCTIME : HIGHBOUND CUTOFF TIME     
    #            NPOLE : # OF POLES FOR THE NONLINEAR FILTER
    ###------------------------------------------------------------
    ### define the sampling frequency and cutoff frequency band ###
    ###------------------------------------------------------------
    mx = np.mean(x)
    x= x-mx
    fs = 1/delt
    lowband = 2*delt/lowcutt
    higband = 2*delt/higcutt
    lband = [lowband, higband]
    # time = np.arange(0,len(x),1)*delt #(0:(length(x)-1))*delt;
    # print('lowband, highband',lband)
    ###------------------------------------------------------------
    ### design Butterworth filter                               ###
    ###------------------------------------------------------------
    [cl,ch] = signal.butter(npole,lband,btype='bandpass')
    
    nfreq = 100000
    ff = 0.5*fs*np.arange(0,nfreq,1)/nfreq #(0:nfreq-1)/(nfreq);
    _,hl = signal.freqz(cl,ch,nfreq)
    hl = hl*np.conj(hl)
    minind = np.min( np.where(hl>=0.8)[0] )
    maxind = np.max( np.where(hl>=0.8)[0] )
    
    minff = ff[minind]
    maxff = ff[maxind]
    ###------------------------------------------------------------
    ### apply filter twice to avoid phase shift                 ###
    ###------------------------------------------------------------
    i = np.arange(1,len(x)+1,1)
    y = signal.lfilter(cl,ch,x)
    # plt.plot(i,x,i,y)
    z0 = y[::-1]
    z = signal.lfilter(cl,ch,z0)
    w = z[::-1]
    wtest = signal.filtfilt(cl,ch,x)
    xband = [lowband/2/delt, lowband/2/delt, higband/2/delt, higband/2/delt]
    # plt.plot(ff[2:nfreq+1], hl[2:nfreq+1],label='0')
    ### y is the first filtered output. z is the second filtered output.
    # w is the rearranged z. It is close to filtfilt(cl,ch,x)
    # plt.plot(i,x,i,w,i,wtest)
    # plt.legend(['data','2nd filter output','filtfilt'])
    # plt.savefig('data_2ndfilter_filtfilt.png',dpi=400,bbox_inches='tight',transparent='True')

    return w, minff, maxff, ff, hl

def bandpass_sos(x,delt,lowcutt,higcutt,npole): ### "sos" should be used for general purpose filtering
    # From Ren-Chieh's lienbandpass.m
    # mybandpass.m calls lienbandpass.m. 
    # The lienbandpass.m is the main bandpass filter
    # BANDX = BANDPASS(X,DT,LCTIME,HCTIME,NPOLE,plotind,FILTERTYPE)
    #               BANDX : BANDPASSED TIME SERIES
    #                X : ORIGINAL TIME SERIES
    #           LCTIME : LOWBOUND CUTOFF TIME     
    #           HCTIME : HIGHBOUND CUTOFF TIME     
    #            NPOLE : # OF POLES FOR THE NONLINEAR FILTER
    ###------------------------------------------------------------
    ### define the sampling frequency and cutoff frequency band ###
    ###------------------------------------------------------------
    mx = np.mean(x)
    x= x-mx
    fs = 1/delt
    lowband = 2*delt/lowcutt
    higband = 2*delt/higcutt
    lband = [lowband, higband]
    ###------------------------------------------------------------
    ### design Butterworth filter                               ###
    ###------------------------------------------------------------
    sos = signal.butter(npole,lband,btype='bandpass',output='sos')
    
    ###------------------------------------------------------------
    ### apply filter twice to avoid phase shift                 ###
    ###------------------------------------------------------------
    y = signal.sosfilt(sos,x)
    z0 = y[::-1]
    z = signal.sosfilt(sos,z0)
    w = z[::-1]
    
    return w

def compdemodfun(x,npts,dt,f,modperiod):
    ''' Complext Demodulation for the 1D time series 

    amp, pha = compdemodfun(x, npts, dt, f, modperiod) 
    x : original time series
    npts : total number of points
    dt : sampling interval in hours
    f : centering frequency in cph
    modperiod : modulation period
    
    '''
    x = np.array( x )
    # print('input type:',type(x),x.shape)
    if x.ndim > 1:
        print('compdemodfun can only deal with 1D vector at a time')
        return "Input is not a 1D time series"
    # test = np.cos( 2*np.pi*np.arange(0,npts)*dt*f )
    # print('test type:',type(test), test.shape)
    # print(np.multiply(x,test).shape)
    xcos = np.multiply(x, np.cos( 2*np.pi*np.arange(0,npts)*dt*f ))
    xsin = np.multiply(x, np.sin( 2*np.pi*np.arange(0,npts)*dt*f ))
    bxcos, _,_ = lpass_NaN(xcos,dt,modperiod,2)
    bxsin, _,_ = lpass_NaN(xsin,dt,modperiod,2)
    amp = np.sqrt( np.square(bxcos) + np.square(bxsin) )
    pha = np.arctan2(bxsin, bxcos)*180/np.pi
    return amp, pha
