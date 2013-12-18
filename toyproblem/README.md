# spotrod/toyproblem

This directory contains five different implementations of a simple function
used in [spotrod](../spotrod), along with a Makefile and a test.py to do benchmarking. The purpose of this is to justify the choice of C for the final product. The five implementations are:

toypython.py: in Python, using a loop or indexing arrays with Boolean masks;
toycython.pyx: in Cython, using a loop or indexing arrays with Boolean masks;
toyc.c: in C, using a loop.

There are \*-setup.py files for the Cython and C version, as well toyc-wrapper.c for interfacing C with numpy.
 
## TODO

For the upcoming Code Coffee presentation, I should create a short slide show and an ipython-notebook.

```
Copyright 2013 Bence Béky

This file is part of Spotrod.

Spotrod is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Spotrod is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Spotrod.  If not, see <http://www.gnu.org/licenses/>.
```
