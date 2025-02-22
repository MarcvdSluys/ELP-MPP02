"""Accurate Moon positions using the Lunar solution ELP/MPP02

This module uses the semi-analytical Lunar solution ELP2000/MPP02 to compute the geocentric position of the
Moon in the dynamical mean ecliptic and equinox of J2000.  This Python code has been adapted from the Fortran
version in libTheSky.



Remarks:

The nominal values of some constants have to be corrected.  There are two sets of corrections, one of which
can be chosen using the parameter 'mode' (used in e.g. initialise()):

- mode=0, the constants are fitted to LLR observations provided from 1970 to 2001 (default);

- mode=1, the constants are fitted to DE405 ephemeris over one century (1950-2060); the lunar angles W1, W2,
  W3 receive also additive corrections to the secular coefficients.  This is known as the 'historical mode'.

When the mode is changed, the constants will be reinitialised and the data file reread.


References:
 
- ELPdoc: Lunar solution ELP, version ELP/MPP02, Jean Chapront and Gerard Francou, October 2002

- Data files and Fortran code: ftp://cyrano-se.obspm.fr/pub/2_lunar_solutions/2_elpmpp02/

- Refereed article: Chapront J., Francou G., A&A 404, 735 (2003)

- libTheSky: https://libthesky.sourceforge.net/


Dependencies:

Apart from the standard Python modules math and sys, numpy and fortranformat must be installed.


Copyright:

Copyright (c) 2019-2025 Marc van der Sluys, Department of Physics, Utrecht University, The Netherlands
and Nikhef (Netherlands institute for high-energy physics and gravitational waves), Amsterdam, The Netherlands -
https://www.nikhef.nl/~sluys/  (this Python code)
 
This file is part of the ELP/MPP02 Python package, 
see: https://pypi.org/project/elp_mpp02/ / https://github.com/MarcvdSluys/ELP-MPP02
 
This is free software: you can redistribute it and/or modify it under the terms of the European Union Public
Licence 1.2 (EUPL 1.2).  This software is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
EU Public Licence for more details.  You should have received a copy of the European Union Public Licence
along with this code.  If not, see <https://www.eupl.eu/1.2/en/>.


Acknowledgements:
- Thanks to A. Vähäkangas for bug reports.

"""


import math as m
import numpy.core as np
import sys


# Global constants:
d2r    = m.radians(1)  # Degrees to radians
r2d    = m.degrees(1)  # Radians to degrees
r2as   = r2d*3600      # Radians to arcseconds
pi     = m.pi
pi2    = pi*2
pio2   = pi/2
jd2000 = 2451545


# Global variables:
modeInit = 999  # := uninitialised
dataDir = '.'   # Directory where the data files are located (default '.': the current directory)

nmpb = np.zeros((3,3))
cmpb = np.zeros(2645)
fmpb = np.zeros((5,2645))
nper = np.zeros((3,4,3))
cper = np.zeros(33256)
fper = np.zeros((5,33256))

w    = np.zeros((3,5))
p1=0.;p2=0.;p3=0.;p4=0.;p5=0.;  q1=0.;q2=0.;q3=0.;q4=0.;q5=0.

eart = np.zeros(5)
peri = np.zeros(5)
dela = np.zeros((4,5))
zeta = np.zeros(5)
p    = np.zeros((8,5))
delnu=0.; dele=0.; delg=0.; delnp=0.; delep=0.; dtasm=0.; am=0.


def initialise_and_read_files(mode=0):
    """Initialise ELP/MPP02 constants and read the data files if necessary
    
    Input parameters:
    
    - mode: Index of the corrections to the constants: 0-Fit to LLR observations (default), 1-Fit to DE405
      1950-2060 (historical)

    """
    
    global modeInit
    # print("Initialise and read files:", modeInit)
    
    # Initializing of constants and reading the files:
    ierr = 0
    if mode!=modeInit:
        initialise(mode)
        ierr = read_files()
        if ierr!=0: return ierr
        
        modeInit = mode  # Indicate that the data have been initialised
        # print("modeInit:", modeInit)
        return ierr
        
    return ierr



def initialise(mode=0):
    """Initialization of the constants and parameters used for the evaluation of the ELP/MPP02 series
    
    Input parameters:
    
    - mode: Index of the corrections to the constants: 0-Fit to LLR observations (default), 1-Fit to DE405
      1950-2060 (historical)

    """
    
    
    global modeInit,  w,eart,peri, dela,zeta,  p,delnu,dele,delg,delnp,delep,dtasm,am
    global p1,p2,p3,p4,p5, q1,q2,q3,q4,q5
    # print("Initialise:", modeInit)
    
    Dprec = -0.29965  # Constant for the correction to the constant of precession - source: IAU 2000A
    
    bp = np.array([
        [  0.311079095,  -0.103837907  ],
        [ -0.4482398e-2,  0.668287e-3  ],
        [ -0.1102485e-2, -0.1298072e-2 ],
        [  0.1056062e-2, -0.178028e-3  ],
        [  0.50928e-4,   -0.37342e-4   ]
    ])  # (2,5)
    
    if mode<0 or mode>1: sys.exit('elp_mpp02.mpp02.initialise(): mode must have value 0 or 1, not %i' % mode)
    
    # Constants for the evaluation of the partial derivatives:
    am     =  0.074801329           # Ratio of the mean motions (EMB / Moon)
    alpha  =  0.002571881           # Ratio of the semi-major axis (Moon / EMB)
    dtasm  =  (2*alpha)/(3*am)      # (2*alpha) / (3*am)
    xa     =  (2*alpha)/3
    
    
    # Corrections to constants:
    if mode==0:  # LLR
        # Values of the corrections to the constants fitted to LLR.  Fit 13-05-02 (2 iterations) except Phi
        # and eps w2_1 and w3_1.  See ELPdoc, Table 3 and paper, Table 1
        Dw1_0   = -0.10525
        Dw2_0   =  0.16826
        Dw3_0   = -0.10760
        Deart_0 = -0.04012
        Dperi   = -0.04854
        Dw1_1   = -0.32311
        Dgam    =  0.00069
        De      =  0.00005
        Deart_1 =  0.01442
        Dep     =  0.00226
        Dw2_1   =  0.08017
        Dw3_1   = -0.04317
        Dw1_2   = -0.03794
    else:  # DE 405 ('historical')
        # Values of the corrections to the constants fitted to DE405 over the time interval (1950-2060)
        Dw1_0   = -0.07008
        Dw2_0   =  0.20794
        Dw3_0   = -0.07215
        Deart_0 = -0.00033
        Dperi   = -0.00749
        Dw1_1   = -0.35106
        Dgam    =  0.00085
        De      = -0.00006
        Deart_1 =  0.00732
        Dep     =  0.00224
        Dw2_1   =  0.08017
        Dw3_1   = -0.04317
        Dw1_2   = -0.03743
        
    # Fundamental arguments (Moon and EMB - ELPdoc, Table 1):
    # W1: mean longitude of the Moon:
    w[0,0]  = elp_dms2rad(218,18,59.95571+Dw1_0)      # Source: ELP - mean motion of the Moon (nu)
    w[0,1]  = (1732559343.73604+Dw1_1)/r2as           # Source: ELP
    w[0,2]  = (        -6.8084 +Dw1_2)/r2as           # Source: DE405
    w[0,3]  =          0.66040e-2/r2as                # Source: ELP
    w[0,4]  =         -0.31690e-4/r2as                # Source: ELP
    
    # W2: mean longitude of the lunar perigee:
    w[1,0]  = elp_dms2rad( 83,21,11.67475+Dw2_0)      # Source: ELP
    w[1,1]  = (  14643420.3171 +Dw2_1)/r2as           # Source: DE405
    w[1,2]  = (       -38.2631)/r2as                  # Source: DE405
    w[1,3]  =         -0.45047e-1/r2as                # Source: ELP
    w[1,4]  =          0.21301e-3/r2as                # Source: ELP
    
    # W3: mean longitude of the lunar ascending node:
    w[2,0]  = elp_dms2rad(125, 2,40.39816+Dw3_0)      # Source: ELP
    w[2,1]  = (  -6967919.5383 +Dw3_1)/r2as           # Source: DE405
    w[2,2]  =          6.3590/r2as                    # Source: DE405
    w[2,3]  =          0.76250e-2/r2as                # Source: ELP
    w[2,4]  =         -0.35860e-4/r2as                # Source: ELP
    
    # Earth-Moon (EMB) elements:
    # Te: mean longitude of EMB:
    eart[0] = elp_dms2rad(100,27,59.13885+Deart_0)    # Source: VSOP2000 - mean motion of EMB (n')
    eart[1] = (129597742.29300 +Deart_1)/r2as         # Source: VSOP2000
    eart[2] =         -0.020200/r2as                  # Source: ELP
    eart[3] =          0.90000e-5/r2as                # Source: ELP
    eart[4] =          0.15000e-6/r2as                # Source: ELP
    
    # Pip: mean longitude of the perihelion of EMB:
    peri[0] = elp_dms2rad(102,56,14.45766+Dperi)      # Source: VSOP2000
    peri[1] =       1161.24342/r2as                   # Source: VSOP2000
    peri[2] =          0.529265/r2as                  # Source: VSOP2000
    peri[3] =         -0.11814e-3/r2as                # Source: VSOP2000
    peri[4] =          0.11379e-4/r2as                # Source: VSOP2000
    
    # Corrections to the secular terms of Moon angles.  This gives a better (long-term?) fit
    #   to DE 406.  See ELPdoc, Table 6/paper, Table 4, line 2:
    if mode==1:  # DE 405 / DE 406
        w[0,3] -= 0.00018865/r2as
        w[0,4] -= 0.00001024/r2as
        
        w[1,2] += 0.00470602/r2as
        w[1,3] -= 0.00025213/r2as
        
        w[2,2] -= 0.00261070/r2as
        w[2,3] -= 0.00010712/r2as
    
    
    # Corrections to the mean motions of the Moon angles W2 and W3, infered from the modifications of the
    # constants:
    x2     =   w[1,1] / w[0,1]
    x3     =   w[2,1] / w[0,1]
    y2     =   am*bp[0,0] + xa*bp[4,0]
    y3     =   am*bp[0,1] + xa*bp[4,1]
    
    d21    =   x2 - y2
    d22    =   w[0,1] * bp[1,0]
    d23    =   w[0,1] * bp[2,0]
    d24    =   w[0,1] * bp[3,0]
    d25    =   y2/am
    
    d31    =   x3 - y3
    d32    =   w[0,1] * bp[1,1]
    d33    =   w[0,1] * bp[2,1]
    d34    =   w[0,1] * bp[3,1]
    d35    =   y3/am
    
    Cw2_1  =  d21*Dw1_1+d25*Deart_1+d22*Dgam+d23*De+d24*Dep
    Cw3_1  =  d31*Dw1_1+d35*Deart_1+d32*Dgam+d33*De+d34*Dep
    
    w[1,1] +=  Cw2_1/r2as
    w[2,1] +=  Cw3_1/r2as
    
    # Arguments of Delaunay:
    for iD in range(5):
        dela[0,iD] = w[0,iD]  - eart[iD]                 # D   =  W1 - Te + 180 degrees
        dela[1,iD] = w[0,iD]  - w[2,iD]                  # F   =  W1 - W3
        dela[2,iD] = w[0,iD]  - w[1,iD]                  # l   =  W1 - W2   mean anomaly of the Moon
        dela[3,iD] = eart[iD] - peri[iD]                 # l'  =  Te - Pip  mean anomaly of EMB
    
    dela[0,0] += pi
    
    # Planetary arguments: mean longitudes for J2000 (from VSOP2000):
    p[0,0] = elp_dms2rad(252, 15,  3.216919)         # Mercury
    p[1,0] = elp_dms2rad(181, 58, 44.758419)         # Venus
    p[2,0] = elp_dms2rad(100, 27, 59.138850)         # EMB (eart(0))
    p[3,0] = elp_dms2rad(355, 26,  3.642778)         # Mars
    p[4,0] = elp_dms2rad( 34, 21,  5.379392)         # Jupiter
    p[5,0] = elp_dms2rad( 50,  4, 38.902495)         # Saturn
    p[6,0] = elp_dms2rad(314,  3,  4.354234)         # Uranus
    p[7,0] = elp_dms2rad(304, 20, 56.808371)         # Neptune
    
    # Planetary arguments: mean motions (from VSOP2000):
    p[0,1] = 538101628.66888/r2as                    # Mercury
    p[1,1] = 210664136.45777/r2as                    # Venus
    p[2,1] = 129597742.29300/r2as                    # EMB (eart(1))
    p[3,1] =  68905077.65936/r2as                    # Mars
    p[4,1] =  10925660.57335/r2as                    # Jupiter
    p[5,1] =   4399609.33632/r2as                    # Saturn
    p[6,1] =   1542482.57845/r2as                    # Uranus
    p[7,1] =    786547.89700/r2as                    # Neptune
    
    p[0:8,2:5] = 0
    
    
    # Mean longitude of the Moon W1 + Rate of precession (pt):
    zeta[0] = w[0,0]
    zeta[1] = w[0,1] + (5029.0966+Dprec)/r2as
    zeta[2] = w[0,2]
    zeta[3] = w[0,3]
    zeta[4] = w[0,4]
    
    # Corrections to the parameters: Nu, E, Gamma, n' et e' (Source: ELP):
    delnu  = (+0.55604+Dw1_1)/r2as/w[0,1]         # Correction to the mean motion of the Moon
    dele   = (+0.01789+De)/r2as                   # Correction to the half coefficient of sin(l) in longitude
    delg   = (-0.08066+Dgam)/r2as                 # Correction to the half coefficient of sin(F) in latitude
    delnp  = (-0.06424+Deart_1)/r2as/w[0,1]       # Correction to the mean motion of EMB
    delep  = (-0.12879+Dep)/r2as                  # Correction to the eccentricity of EMB
    
    # Precession of the longitude of the ascending node of the mean ecliptic of date on fixed ecliptic J2000
    # (Laskar, 1986):
    # P: sine coefficients:
    p1 =  0.10180391e-4
    p2 =  0.47020439e-6
    p3 = -0.5417367e-9
    p4 = -0.2507948e-11
    p5 =  0.463486e-14
    
    # Q: cosine coefficients:
    q1 = -0.113469002e-3
    q2 =  0.12372674e-6
    q3 =  0.1265417e-8
    q4 = -0.1371808e-11
    q5 = -0.320334e-14
    
    return
  
  
def read_files():
    """Read the six data files containing the ELP/MPP02 series
    
    Return values:
    
    - ierr: File error index: ierr=0: no error, ierr=1: file error
    
    """
    
    # print("Read files:")
    global nmpb,cmpb,fmpb, nper,cper,fper
    global w,eart,peri, dela,zeta,  p,delnu,dele,delg,delnp,delep,dtasm,am
    
    # Read the Main Problem series:
    ir  = 0
    ilu = np.zeros(4)  # will contain ints
    a   = 0.
    b   = np.zeros(5)
    # ierr=1
    # nerr=0
    
    import fortranformat as ff
    formatMainHeader = ff.FortranRecordReader('(25x,I10)')              # Block header format
    formatMainBody   = ff.FortranRecordReader('(4I3,2x,F13.5,5F12.2)')  # Block body format
    
    
    for iFile in range(3):  # Main-problem files 1-3
        fileName = dataDir+'/ELP_MAIN.S'+str(iFile+1)
        try:
            inFile = open(fileName,'r')
        except Exception as err:
            sys.stderr.write(str(err)+'\n\n')
            sys.stderr.write('Please ensure that:\n')
            sys.stderr.write(' 1) you downloaded the data files ELP_*.S[123] from '+
                             'ftp://cyrano-se.obspm.fr/pub/2_lunar_solutions/2_elpmpp02/\n')
            sys.stderr.write(' 2) you set the variable mpp.dataDir to the correct value\n')
            exit(1)
        
        line = inFile.readline()
        nmpb[iFile,0] = formatMainHeader.read(line)[0]
        # if nerr!=0: return 3
        
        nmpb[iFile,1] = ir+1
        nmpb[iFile,2] = nmpb[iFile,0] + nmpb[iFile,1] - 1

        nLines = int(round(nmpb[iFile,0]))
        for iLine in range(nLines):
            line = inFile.readline()
            ilu[0],ilu[1],ilu[2],ilu[3], a, b[0],b[1],b[2],b[3],b[4] = formatMainBody.read(line)
            # if nerr!=0: return 4
            
            tgv = b[0] + dtasm*b[4]
            if iFile==2:  a -= 2*a*delnu/3
            cmpb[ir] = a + tgv*(delnp-am*delnu) + b[1]*delg + b[2]*dele + b[3]*delep
            
            for k in range(5):
                fmpb[k,ir] = 0
                for i in range(4):
                    fmpb[k,ir] += ilu[i] * dela[i,k]
                    
            if iFile==2: fmpb[0,ir] += pio2
            ir += 1
            
        inFile.close()
        
    
    # Read the Perturbations series:
    ir = 0
    ipt = 0
    icount = 0
    s = 0.0
    c = 0.0
    ifi = np.zeros(16)  # will contain ints
    
    formatPertHeader = ff.FortranRecordReader('(25x,2I10)')         # Perturbation header format
    formatPertBody   = ff.FortranRecordReader('(I5,2D20.13,16I3)')  # Perturbation body format
    
    for iFile in range(3):  # Perturbation files 1-3
        fileName = dataDir+'/ELP_PERT.S'+str(iFile+1)
        try:
            inFile = open(fileName,'r')
        except Exception as err:
            sys.stderr.write(str(err)+'\n')
            exit(1)
            
        for it in range(4):
            # if nerr!=0: return 6
            line = inFile.readline()
            nper[iFile,it,0],ipt = formatPertHeader.read(line)
            
            nper[iFile,it,1] = ir+1
            nper[iFile,it,2] = nper[iFile,it,0] + nper[iFile,it,1] - 1
            if nper[iFile,it,0]==0: continue  # cycle to next it
            
            nLines = int(round(nper[iFile,it,0]))
            for iLine in range(nLines):
                line = inFile.readline()
                ( icount,s,c,ifi[0],ifi[1],ifi[2],ifi[3],ifi[4],ifi[5],ifi[6],ifi[7],ifi[8],ifi[9],ifi[10], 
                  ifi[11],ifi[12],ifi[13],ifi[14],ifi[15] ) = formatPertBody.read(line)
                # if nerr!=0: return 7
                
                cper[ir] = m.sqrt(c**2+s**2)
                pha = m.atan2(c,s)
                if pha<0: pha = pha+pi2
                
                for k in range(5):
                    fper[k,ir] = 0
                    if k==0: fper[k,ir] = pha
                    for i in range(4):
                        fper[k,ir] += ifi[i] * dela[i,k]
                        
                    for i in range(4,12):
                        fper[k,ir] += ifi[i] * p[i-4,k]
                    
                    fper[k,ir] += ifi[12] * zeta[k]
                ir = ir+1
            
        inFile.close()        
        
    return 0
    


def elp_dms2rad(deg,min,sec):
    """Function for the conversion sexagesimal degrees -> radians in initialise()"""
    
    return (deg+min/60+sec/3600) * d2r



def compute_lbr(jd, mode=0):
    """Compute the spherical lunar coordinates using the ELP2000/MPP02 lunar theory in the dynamical mean ecliptic
          and equinox of J2000.
    
    Input parameters:
    
    - jd: Julian day to compute Moon position for
    
    - mode: Index of the corrections to the constants: 0-Fit to LLR observations (default), 1-Fit to DE405
      1950-2060 (historical)
    
    Return values:
    
    - lon: Ecliptic longitude (rad)
    
    - lat: Ecliptic latitude (rad)
    
    - rad: Distance (km)

    """

    # print("Compute lbr:")
    
    xyz,vxyz, ierr = compute_xyz(jd, mode)
    
    # Compute ecliptic l,b,r:
    rad = m.sqrt(sum(xyz**2))
    lon = m.atan2(xyz[1], xyz[0])
    lat = m.asin(xyz[2]/rad)
    
    # rad = rad/1.49597870700e8  # km -> AU
    
    # print('jd, xyz: ', jd, xyz[0:3])
    # print(jd, (lon%pi2)*r2d, lat*r2d, rad)
    
    return lon,lat,rad



def compute_xyz(jd, mode=0):
    """Compute the rectangular lunar coordinates using the ELP/MPP02 lunar theory in the dynamical mean ecliptic
    and equinox of J2000.
    
    Input parameters:
    
    - jd: Julian day to compute Moon position for
    
    - mode: Index of the corrections to the constants: 0-Fit to LLR observations, 1-Fit to DE405 1950-2060
      (historical)
    
    Return value:
    
    - xyz(3): Geocentric rectangular coordinates (km)
    
    - vxyz(3): Geocentric rectangular velocities (km/day):
    
    - ierr: File error index - ierr=0: no error, ierr=1: file error

    """
    
    # print("Compute xyz:")
    global nmpb,cmpb,fmpb, nper,cper,fper
    global w, p1,p2,p3,p4,p5, q1,q2,q3,q4,q5
    
    # Constants:
    a405 = 384747.9613701725  # Moon mean distance for DE405
    aelp = 384747.980674318   # Moon mean distance for ELP
    sc   = 36525              # Julian century in days
    
    # Initialise data and read files if needed:
    ierr = initialise_and_read_files(mode)
    if ierr!=0: sys.exit('Could not read ELP-MPP02 files')
    
    
    # Initialization of time powers:
    rjd  = jd - jd2000  # Reduced JD - JD since 2000
    t = np.zeros(5)
    t[0] = 1
    t[1] = rjd/sc       # t: time since 2000 in Julian centuries
    t[2] = t[1]**2      # t^2
    t[3] = t[2]*t[1]    # t^3
    t[4] = t[2]**2      # t^4
    
    # Evaluation of the series: substitution of time in the series
    v = np.zeros(6)
    for iVar in range(3):  # iVar=0,1,2: Longitude, Latitude, Distance
        
        # Main Problem series:
        nLineStart = int(round(nmpb[iVar,1]))
        nLineEnd   = int(round(nmpb[iVar,2]))
        for iLine in range(nLineStart-1, nLineEnd):
            x = cmpb[iLine]
            y = fmpb[0,iLine]
            yp = 0
            
            for k in range(1,5):  # k=1,4
                y  +=   fmpb[k,iLine] * t[k]
                yp += k*fmpb[k,iLine] * t[k-1]
            
            v[iVar]   += x * m.sin(y)
            v[iVar+3] += x * yp * m.cos(y)
            
        
        # Perturbations series:
        for it in range(4):
            nLineStart = int(round(nper[iVar,it,1]))
            nLineEnd   = int(round(nper[iVar,it,2]))
            for iLine in range(nLineStart-1, nLineEnd):
                x = cper[iLine]
                y = fper[0,iLine]
                xp = 0
                yp = 0
                if it!=0: xp = it * x * t[it-1]
                
                for k in range(1,5):  # k=1,4
                    y  +=     fper[k,iLine] * t[k]
                    yp += k * fper[k,iLine] * t[k-1]
             
                v[iVar]   += x * t[it] * m.sin(y)
                v[iVar+3] += xp * m.sin(y)  +  x * t[it] * yp * m.cos(y)
            
        
    # Compute the spherical coordinates for the mean inertial ecliptic and equinox of date:
    v[0] = v[0]/r2as + w[0,0] + w[0,1]*t[1] + w[0,2]*t[2] + w[0,3]*t[3] + w[0,4]*t[4]  # Longitude + mean lon. (rad)
    v[1] = v[1]/r2as                                                                   # Latitude (rad)
    v[2] = v[2] * a405 / aelp                                                          # Distance (km)
    
    # v[0] = v[0] % pi2
    
    # print('t: ', t[0],t[1],t[2],t[3],t[4])
    # print('v:       ', v[0]*r2d,v[1]*r2d,v[2], v[3],v[4],v[5])
    
    # compute the rectangular coordinates (for the EoD?):
    clamb  = m.cos(v[0])
    slamb  = m.sin(v[0])
    cbeta  = m.cos(v[1])
    sbeta  = m.sin(v[1])
    cw     = v[2]*cbeta
    sw     = v[2]*sbeta
    # print("c/s l/b: ", clamb,slamb, cbeta,sbeta)
    
    x0     = cw*clamb
    x1     = cw*slamb
    x2     = sw
    # print("x1,x2,x3: ", x0,x1,x2)
    # print("p,q: ", p1,p2,p3,p4,p5, q1,q2,q3,q4,q5)
    
    # Is this simply precession in rectangular coordinates to J2000?  From?
    pw     = (p1 + p2*t[1] + p3*t[2] + p4*t[3] + p5*t[4]) * t[1]
    qw     = (q1 + q2*t[1] + q3*t[2] + q4*t[3] + q5*t[4]) * t[1]
    
    ra     = 2*m.sqrt(1 - pw**2 - qw**2)
    pwqw   = 2*pw*qw
    pw2    = 1 - 2*pw**2
    qw2    = 1 - 2*qw**2
    pwra   = pw*ra
    qwra   = qw*ra
    
    xyz = np.zeros(3)
    xyz[0] =  pw2*x0  + pwqw*x1 + pwra*x2
    xyz[1] =  pwqw*x0 + qw2*x1  - qwra*x2
    xyz[2] = -pwra*x0 + qwra*x1 + (pw2+qw2-1)*x2
    
    
    # Compute the rectangular velocities for the equinox J2000:
    v[3]   = v[3]/r2as + w[0,1] + 2*w[0,2]*t[1] + 3*w[0,3]*t[2] + 4*w[0,4]*t[3]
    v[4]   = v[4]/r2as
    
    xp1    = (v[5]*cbeta - v[4]*sw)*clamb - v[3]*x1
    xp2    = (v[5]*cbeta - v[4]*sw)*slamb + v[3]*x0
    xp3    = v[5]*sbeta  + v[4]*cw
    
    ppw    = p1 + (2*p2 + 3*p3*t[1] + 4*p4*t[2] + 5*p5*t[3]) * t[1]
    qpw    = q1 + (2*q2 + 3*q3*t[1] + 4*q4*t[2] + 5*q5*t[3]) * t[1]
    ppw2   = -4*pw*ppw
    qpw2   = -4*qw*qpw
    ppwqpw = 2*(ppw*qw + pw*qpw)
    rap    = (ppw2+qpw2)/ra
    ppwra  = ppw*ra + pw*rap
    qpwra  = qpw*ra + qw*rap
    
    vxyz = np.zeros(3)
    vxyz[0] = (pw2*xp1 + pwqw*xp2 + pwra*xp3  +  ppw2*x0 + ppwqpw*x1 + ppwra*x2) / sc
    vxyz[1] = (pwqw*xp1 + qw2*xp2 - qwra*xp3  +  ppwqpw*x0 + qpw2*x1 - qpwra*x2) / sc
    vxyz[2] = (-pwra*xp1 + qwra*xp2 + (pw2+qw2-1)*xp3  -  ppwra*x0 + qpwra*x1 + (ppw2+qpw2)*x2) / sc
    
    return xyz,vxyz, ierr

