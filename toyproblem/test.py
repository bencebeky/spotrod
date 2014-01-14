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

import numpy;
import toypython;
import toycython;
import toyc;
from matplotlib import pyplot;
pyplot.ion();

r = numpy.linspace(0.0, 1.0, 100);

test1 = toypython.circleanglesorted;
test2 = toyc.circleanglesorted;

def test(p, z):
  pyplot.figure();
  testdata1 = test1(r, p, z);
  testdata2 = test2(r, p, z);
  pyplot.plot(r, testdata1, "r-", r, testdata2, "b-");
  print numpy.max(numpy.abs(testdata1 - testdata2));

test(0.1, 0.5);
test(0.3, 0.2);
