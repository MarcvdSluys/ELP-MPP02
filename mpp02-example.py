#!/bin/env python3

from elp_mpp02 import mpp02 as mpp

mpp.dataDir = 'data'  # Set the dir where the ELP_*.S* data files can be found

mode = 1      # Historical mode
jd = 2451545  # 2000-01-01

# Compute the spherical geocentric ecliptical coordinates of the Moon:
# - note: on the first call, the data files must be read, which can take ~10s or so!
lon,lat,dist = mpp.compute_lbr(jd, mode)  
print('jd =',jd, ':   lon =',lon,'rad,  lat =',lat, 'rad,  dist =',dist,'km.')

xyz,vxyz,ierr = mpp.compute_xyz(jd, mode)  
print('jd =',jd, ':   X =',xyz[0],'km,  Y =',xyz[1], 'km,  Z =',xyz[2],'km.')
print('jd =',jd, ':   Vx =',vxyz[0],'km/day,  Vy =',vxyz[1], 'km/day,  Vz =',vxyz[2],'km/day.')


