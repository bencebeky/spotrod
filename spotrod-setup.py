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

from distutils.core import setup, Extension
import numpy

name = "spotrod"
version = "1.0"
description = "A semi-analytic model for transits of spotted stars"
author = "Bence Béky"
author_email = "zsebkecske@gmail.com"
maintainer = author
maintainer_email = author_email
url = "https://github.com/bencebeky/spotrod"

module = Extension(name, sources=['spotrod.c', 'spotrod-python-wrapper.c'],
                   include_dirs= [numpy.get_include()],
                   extra_compile_args=['-Ofast'])

setup(name = name, version = version, description = description,
       author = author, author_email = author_email,
         maintainer = maintainer, maintainer_email = maintainer_email,
           url = url, ext_modules = [module])
