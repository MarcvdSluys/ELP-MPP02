#!/bin/env python

from elp_mpp02 import mpp02 as mpp
import numpy as np
import astrotool.date_time as atdt
import astroconst as ac

mpp.dataDir = 'data'  # Set the dir where the ELP_*.S* data files can be found

mode = 1      # Historical mode
jd = 2451545  # 2000-01-01

# Compute the spherical geocentric ecliptical coordinates of the Moon:
# - note: on the first call, the data files must be read, which can take ~10s or so!
lon,lat,dist = mpp.compute_lbr(jd, mode)  
print('jd =',jd, ':   lon =',lon,'rad,  lat =',lat, 'rad,  dist =',dist,'km.')

# Compute the rectangular geocentric ecliptical coordinates and velocities of the Moon:
xyz,vxyz,ierr = mpp.compute_xyz(jd, mode)  
print('jd =',jd, ':   X =',xyz[0],'km,  Y =',xyz[1], 'km,  Z =',xyz[2],'km.')
print('jd =',jd, ':   Vx =',vxyz[0],'km/day,  Vy =',vxyz[1], 'km/day,  Vz =',vxyz[2],'km/day.')


# Compute the spherical geocentric ecliptical coordinates of the Moon for a range of dates:
print()
print(' year,                         JD,      longitude (rad),        latitude (rad),       distance (km)')
for year in range(-3000,2501,500):
    jd = atdt.jd_from_date(year,1,1.5)
    
    # Compute the spherical geocentric ecliptical coordinates of the Moon:
    # - note: on the first call, the data files must be read, which can take ~10s or so!
    lon,lat,dist = mpp.compute_lbr(jd, mode)
    print("%5i,  %25.15f,  %19.15f,  %20.17f,  %18.10f" % (year, jd, np.mod(lon,ac.pi2),lat,dist))


