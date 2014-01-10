#!/usr/bin/env python
# coding=utf8
# 
# Copyright 2013, 2014 Bence Béky
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

r = ", numpy; r = numpy.linspace(0.0, 1.0, 10000)";
arg = "(r, 0.1, 0.5)";

time1 = timeit(stmt="toypython.circleangleloop" + arg, setup="import toypython" + r, number=n);
time2 = timeit(stmt="toypython.circleanglemask" + arg, setup="import toypython" + r, number=n);
time3 = timeit(stmt="toypython.circleanglesorted" + arg, setup="import toypython" + r, number=n);
time4 = timeit(stmt="toycython.circleangleloop" + arg, setup="import toycython" + r, number=n);
time5 = timeit(stmt="toycython.circleanglemask" + arg, setup="import toycython" + r, number=n);
time6 = timeit(stmt="toycython.circleanglesorted" + arg, setup="import toycython" + r, number=n);
time7 = timeit(stmt="toyc.circleangleloop" + arg, setup="import toyc" + r, number=n);
time8 = timeit(stmt="toyc.circleanglesorted" + arg, setup="import toyc" + r, number=n);

print("Python loop:   {0:5.2f} ms.".format(1000*time1/n));
print("Python mask:   {0:5.2f} ms.".format(1000*time2/n));
print("Python sorted: {0:5.2f} ms.".format(1000*time3/n));
print("Cython loop:   {0:5.2f} ms.".format(1000*time4/n));
print("Cython mask:   {0:5.2f} ms.".format(1000*time5/n));
print("Cython sorted: {0:5.2f} ms.".format(1000*time6/n));
print("C      loop:   {0:5.2f} ms.".format(1000*time7/n));
print("C      sorted: {0:5.2f} ms.".format(1000*time8/n));
