#!/usr/bin/env python
# coding=utf8
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

from distutils.core import setup;
from distutils.extension import Extension;
from Cython.Distutils import build_ext;

setup(cmdclass = {'build_ext': build_ext}, ext_modules = [Extension("toycython", ["toycython.pyx"])]);
