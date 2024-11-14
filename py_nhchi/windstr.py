import numpy as np
import math

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
               elif (U > 1) | (U <= 3):
                    Cd = (0.62+1.56/U)*0.001
               elif (U > 3) | (U < 10):
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