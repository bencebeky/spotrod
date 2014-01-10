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

from distutils.core import setup, Extension;
import numpy.distutils.misc_util;

c_ext = Extension("toyc", ["toyc-wrapper.c", "toyc.c"], extra_compile_args=['-Ofast']);

setup(ext_modules=[c_ext], include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs());
