import numpy as np
import math

def mypower(fev, evl):
    """ p = power(fev,evl)
    Compute the power spectrum of a time series using the eigenspectra estimate in fev.  
    fev: The columns of fev are the independent spectral estimates, and 
    evl: is a column vector of eigenvalues corresponding to the eigenvectors (windows)
    used the compute the columns of fev.
    By John Kuehne, June 1990.
    """
    a = np.multiply(fev, np.conjugate(fev))
    p = np.multiply(a,evl)/sum(evl)
    
    return p

# def spectra(x, evc):
#     """ fev = spectra(x,evc) 
#     Compute complex-valued eigenspectra of column vector x using eigenvectors (windows) in evc. 
#     The eigenvectors are in the columns of evc, which can be supplied by a window function.
#     By John Kuehne, June 1990.
#     """
#     n,k = evc.shape

#     for t in range(k):
#         ev[:n,t] = np.multiply(evc[:n,t],x)
    
#     fev=np.fft.fft(ev)
    
#     return fev
