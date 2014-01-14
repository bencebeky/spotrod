# spotrod/toyproblem
### A semi-analytic model for transits of spotted stars.

```
Copyright 2013, 2014 Bence Béky

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

## Contents

This directory contains eight different implementations of a simple function
used in [spotrod](../spotrod), along with a Makefile and a test.py to do benchmarking. The purpose of this is to justify the choice of C for the final product. The eight implementations are:

- [toypython.py](toypython.py) in Python, using a loop or indexing arrays with Boolean masks;
- [toycython.pyx](toycython.pyx) in Cython, using a loop or indexing arrays with Boolean masks;
- [toyc.c](toyc.c) in C, using a loop.

Each loop method is implemented two times: once executing all comparisons in the loop core, and once saving on that assuming that the input array r is increasing.

Other files in this directory are:

- [Makefile](Makefile) Makefile for generating toycython.so and toyc.so;
- [README.md](README.md) this file;
- [benchmark.py](benchmark.py) a Python script for benchmarking the different implementations;
- [test.py](test.py) sanity check: a Python script for comparing the different implementations;
- [toyc-setup.py](toyc-setup.py) a Python script for compiling toyc.so
- [toyc-wrapper.c](toyc-wrapper.c) C code wrapping toyc.c in the numpy C API;
- [toyc.h](toyc.h) headers for functions in toyc.c;
- [toycython-setup.py](toycython-setup.py) a Python script for compiling toycython.so.

## Compilation

To generate the modules `toycython.so` and `toyc.so`, run `make` without any arguments. To generate `toypython.pyc`, just `import toypython` from Python, it gets automatically compiled.

## TODO

For the upcoming Code Coffee presentation, I should create a short slide show and an ipython-notebook.
