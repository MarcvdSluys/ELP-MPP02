# -*- coding: utf-8 -*-
# SPDX-License-Identifier: EUPL-1.2
#  
#  Copyright (c) 2019-2025  Marc van der Sluys - Nikhef/Utrecht University - marc.vandersluys.nl
#   
#  This file is part of the evTool Python package:
#  Analyse and display output from the binary stellar-evolution code ev (a.k.a. STARS or TWIN).
#  See: https://github.com/MarcvdSluys/evTool
#   
#  This is free software: you can redistribute it and/or modify it under the terms of the European Union
#  Public Licence 1.2 (EUPL 1.2).  This software is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
#  PURPOSE.  See the EU Public Licence for more details.  You should have received a copy of the European
#  Union Public Licence along with this code.  If not, see <https://www.eupl.eu/1.2/en/>.


"""Accurate Moon positions using the Lunar solution ELP/MPP02

This module uses the semi-analytical Lunar solution ELP2000/MPP02 to compute the geocentric position of the
Moon in the dynamical mean ecliptic and equinox of J2000.  This Python code has been adapted from the Fortran
version in libTheSky.

"""


name = "elp-mpp02"
__author__ = 'Marc van der Sluys - Nikhef/Utrecht University - marc.vandersluys.nl'
