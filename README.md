# ELP_MPP02: accurate Moon positions using the lunar solution ELP/MPP02 in Python #

Compute accurate Moon positions using the semi-analytical lunar solution ELP2000/MPP02 by Chapront & Francou
(2003) in Python.


## Installation ##

This package can be installed using `pip install elp-mpp02`.  You will need to manually download the six data
files `ELP_MAIN/PERT.S1/2/3` from
[ftp://cyrano-se.obspm.fr/pub/2_lunar_solutions/2_elpmpp02](ftp://cyrano-se.obspm.fr/pub/2_lunar_solutions/2_elpmpp02)
and save them in a directory of your choice.


## Using the package ##

You can import the package as follows:
```python
from elp_mpp02 import mpp02 as mpp
```

Then, make sure you define the directory where the data files are located (if not in the current dir). 
For the subdir `data/` of the current directory, do:
```python
mpp.dataDir = 'data'  # Set the dir where the ELP_*.S* data files can be found
```

Choose whether to run the code in LLR (`mode=0`; default) or DE405 (`mode=1`; 'historical') mode, select a
Julian day and compute the Moon position for that instance:
```python
mode = 1  # Historical mode
jd = 2451545
lon,lat,dist = mpp.compute_lbr(jd, mode)  
print('jd =',jd, ':   lon =',lon,'rad,  lat =',lat, 'rad,  dist =',dist,'km.')
```
The result should be
```
jd = 2451545 :   lon = -2.385534575256455 rad,  lat = 0.09024868423130429 rad,  dist = 402448.6385830673 km
```

The ecliptical longitude and latitude are expressed in radians, the distance is in kilometres.  The
coordinates are valid for the mean equinox of J2000.  Note that on the first call, the constants must be
initialised and the data files need to be read, which can take ~10s.  If `mode` is changed between calls, the
data must be reinitialised.


## Copyright and licence ##

* Copyright:  2019-2025 Marc van der Sluys, Department of Physics, Utrecht University and Nikhef (institute for high-energy physics and gravitational waves), Amsterdam (NL)
* Contact:    https://marc.vandersluys.nl
* Website:    [Github](https://github.com/MarcvdSluys/ELP-MPP02), [PyPI](https://pypi.org/project/elp_mpp02/)
* Licence:    [GPLv3+](https://www.gnu.org/licenses/gpl.html)


## References ##

* [Chapront & Francou (2003)](https://ui.adsabs.harvard.edu/abs/2003A%26A...404..735C/abstract)
* [FTP data files](ftp://cyrano-se.obspm.fr/pub/2_lunar_solutions/2_elpmpp02) &mdash; in case [FTP urls don't work in Markdown](https://github.com/gollum/gollum/issues/759): ftp://cyrano-se.obspm.fr/pub/2_lunar_solutions/2_elpmpp02
* This Python code is adapted from the Fortran implementation in [libTheSky](https://libthesky.sourceforge.net/)
