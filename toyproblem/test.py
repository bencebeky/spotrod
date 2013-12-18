#!/usr/bin/env python
# 
# Copyright 2013 Bence Béky
#
# This file is part of Spotrod.
# 
# Spotrod is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Spotrod is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with Spotrod.  If not, see <http://www.gnu.org/licenses/>.
# 

from timeit import timeit;

n = 50;

time1 = timeit(stmt="toypython.circleangleloop(r, 0.4, 0.5)", setup="import numpy, toypython; r=numpy.linspace(0.0, 1.0, 10000)", number=n);
time2 = timeit(stmt="toypython.circleanglemask(r, 0.4, 0.5)", setup="import numpy, toypython; r=numpy.linspace(0.0, 1.0, 10000)", number=n);
time3 = timeit(stmt="toycython.circleangleloop(r, 0.4, 0.5)", setup="import numpy, toycython; r=numpy.linspace(0.0, 1.0, 10000)", number=n);
time4 = timeit(stmt="toycython.circleanglemask(r, 0.4, 0.5)", setup="import numpy, toycython; r=numpy.linspace(0.0, 1.0, 10000)", number=n);
time5 = timeit(stmt="toyc.circleangle(r=r, p=0.4, z=0.5)", setup="import numpy, toyc; r=numpy.linspace(0.0, 1.0, 10000)", number=n);

print("Python loop: {0:5.2f} ms.".format(1000*time1/n));
print("Python mask: {0:5.2f} ms.".format(1000*time2/n));
print("Cython loop: {0:5.2f} ms.".format(1000*time3/n));
print("Cython mask: {0:5.2f} ms.".format(1000*time4/n));
print("C:           {0:5.2f} ms.".format(1000*time5/n));

