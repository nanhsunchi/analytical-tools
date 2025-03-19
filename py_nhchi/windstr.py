import numpy as np
import math
import netCDF4 as nc 

def ra_windstr_nc(u,v,roh):
    ### This is translated from ra_windstr_nc.m, which was modified from ra_windstr.m written by 
    # # Ramkrushn S. Patel (ramkrushn.scrv89@gmail.com)

    # DESCRIPTION:  Function to compute wind stress from wind field data Based on Gill, 1982 
    # Formula and a non-linear Cd based on Large and Pond (1981), modified for low wind 
    # speeds (Trenberth et al., 1990)
    # 
    # INPUTS: 
    # u = Zonal wind component [m/s], must be 2D
    # v = Meridional wind component [m/s], must be 2D
    # roh = air density in kg/m^3
    # OUTPUT: 
    # Tx = Zonal wind stress [N/m^2]
    # Ty = Meridional wind stress [N/m^2]
    # 
    # DISCLAIMER: 
    # Albeit this function is designed only for academic purpose, it can be implemented in 
    # research. Nonetheless, author does not guarantee the accuracy.
    # 
    # REFERENCE:
    # A.E. Gill, 1982, �Atmosphere-Ocean Dynamics�, Academy Press, Vol. 30.
    # W. G. Large & S. Pond., 1981,�Open Ocean Measurements in Moderate to Strong Winds�, 
    # J. Physical Oceanography, Vol. 11, pp. 324 - 336.
    # K.E. Trenberth, W.G. Large & J.G. Olson, 1990, �The Mean Annual Cycle in Global Ocean 
    # Wind Stress�, J.Physical Oceanography, Vol. 20, pp. 1742 � 1760.

    # Defining Constant 
    # roh=1.2; % kg/m^3, air density
    # Computation of Wind Stresses
     if u.ndim == 1:
          lt = len(u)
          ln = 1
     if u.ndim == 2:
          lt, ln = u.shape
     
     Tx = np.nan*np.ones( (lt, ln) )
     Ty = Tx.copy()

     for ii in range(lt):
          for jj in range(ln):
               if u.ndim == 1:
                    U = np.sqrt( np.square(u[ii]) + np.square(v[ii])) # Wind speed
               if u.ndim == 2:
                    U = np.sqrt( np.square(u[ii,jj]) + np.square(v[ii,jj])) # Wind speed
               ###
               if U <= 1:
                    Cd = 0.00218
               elif (U > 1) & (U <= 3):
                    Cd = (0.62+1.56/U)*0.001
               elif (U > 3) & (U < 10):
                    Cd = 0.00114; 
               else:
                    Cd = (0.49+0.065*U)*0.001
               ###
               if u.ndim == 1:
                    Tx[ii]=Cd*roh*U*u[ii]# kg/m^3*m/s*m/s= N/m^2
                    Ty[ii]=Cd*roh*U*v[ii]
               if u.ndim == 2:
                    Tx[ii,jj]=Cd*roh*U*u[ii, jj]# kg/m^3*m/s*m/s= N/m^2
                    Ty[ii,jj]=Cd*roh*U*v[ii, jj]
               # print(U, Tx, Ty)

     Tx = np.squeeze(Tx)
     Ty = np.squeeze(Ty)
     return Tx, Ty

# ra_windstr_nc(np.array([[0.01,0.1],[0.01,0.1]]),np.array([[0.015,0.15],[0.015,0.15]]),0)

def wind_adj2height_PL(u1, z1, z2):
     ''' Adjust wind speeds from height, z1 (m), to height z2 (m) following a simple method descirbed by S.A. Hsu et al. (1994).
      See https://www.ndbc.noaa.gov/faq/adjust_wind.shtml#:~:text=NDBC%20adjusts%20wind%20speeds%20to,that%20typical%20of%20ship%20anemometers.
      This method was tested and found to compare favorably with the more elaborate method under near-neutral stability. 
      This is the condition most frequently encountered at sea and occurs when air and water temperatures are not too far apart. 
      The method, referred to as Power Law Method, is offered here for those who may want to explore the nature of the marine wind speed profile
      without having to deal with the complexity of the W.T. Liu et al (1979) method that NDBC uses. 
      
      input: u1: float, int or array of wind speeds, or wind components (u, v). It could be negative. 
      output u2: float/int or array of adjusted wind speeds or wind components. '''
     z2_z1 = np.divide( z2, z1 )
     u2 = np.multiply(u1, np.power(z2_z1,0.11))
     return u2

def Cd_e13_CH2020(U10,V10):
     ''' get Cd(10m) values for U10/V(10) 
     input: U10, V10: zonal and meridional wind velocity at 10m in m/s.
     output: Cd value at 10m at the input wind velocity. 
     '''
     ''' Look at input data '''
     if isinstance(U10, np.ndarray) | isinstance(U10, list) \
          | isinstance(U10, float) | isinstance(U10, int):
          U10M = np.sqrt( np.square(U10) + np.square(V10) )
     else:
          print('Input data type has to the following: np.ndarray, list, float, int.')
          return
     ''' READ Cd Data '''
     ''' Read Edson et al. (2013), field data '''
     path = '/Users/chi/Documents/Data/Curcic-Haus-2020-1.0.0/data/'
     e13 = np.loadtxt(path+'edson_etal_2013.txt', skiprows=1)
     U10_e13, CD_e13, CD_std_e13 = e13[:,0], e13[:,1]*1e-3, e13[:,2]*1e-3
     ''' Read Curcic & Haus (2020), lab data '''
     CH2020 = nc.Dataset( path+'asist-windonly-salt.nc' )
     U10_CH2020 = CH2020.variables['U10_momentum_budget'][:]
     CD_CH2020 = CH2020.variables['CD_momentum_budget'][:]
     # remove high-error momentum budget estimates
     CD_CH2020[:6] = np.nan
     # print(U10_e13)
     # print(CD_e13)
     # print(U10_CH2020)
     # print(CD_CH2020)
     ''' Combine Cd Data '''
     iCH2020 = U10_CH2020> U10_e13[-1]
     U10_combine = np.concatenate( (U10_e13, U10_CH2020[iCH2020]),axis=0 )
     CD10_combine = np.concatenate( (CD_e13, CD_CH2020[iCH2020]),axis=0 )
     # print(U10_combine)
     # print(CD10_combine)
     ''' interpolate to get Cd values '''
     CD10_interp = np.interp( U10M, U10_combine, CD10_combine, \
                             left=CD10_combine[0], right=CD10_combine[-1] )
     return CD10_interp

     