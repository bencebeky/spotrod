# spotrod
### A semi-analytic model for transits of spotted stars.

## Contents

This repository contains `spotrod`, a semi-analytic model for transits of spotted stars. The model is implemented in C and comes with Python API. The following files and directories are included in the root directory:

- [toyproblem](toyproblem) A toy problem benchmarking Python, Cython, and C implementations of circleangle(), the function calculating $\beta$;
- [COPYING](COPYING) the GNU General Public License;
- [Makefile](Makefile) makefile for generating toycython.so and toyc.so;
- [README.md](README.md) this file;
- [spotrod.c](spotrod.c) C implementation of the model;
- [spotrod.h](spotrod.h) headers for functions in spotrod.c;
- [spotrod-python-wrapper.c](spotrod-python-wrapper.c) C code wrapping spotrod.c in the numpy C API;
- [spotrod-setup.py](spotrod-setup.py) a Python script for compiling spotrod.so;
- [test.py](test.py) a minimal script generating a model lightcurve;
- [test.png](test.png) output of minimal script;
- [mcmc.py](mcmc.py) an MCMC simulation using the [emcee](http://dan.iel.fm/emcee/) package, creating an animation.

## Installation / Compilation
To install this package, first ensure you have Python and the necessary build tools installed. Clone the repository to your local machine and navigate to the repository directory. Then, 
run the following command to build and install the package into your current Python environment:  

```
$ pip install .
```

This command will automatically handle dependencies and compile any C extensions. If you prefer an editable install for development, use: 
  
```
$ pip install -e .
```

This setup allows you to make changes to the codebase and see the effects without reinstalling.

To compile the package without installing it into the currently activated Python environment run:
```
$ python setup.py build
```

To clean up (some of) the build process products run:
```
$ python setup.py clean --all
```

## Citation

If you use `spotrod` in a publication, please consider citing [Béky, Kipping, and Holman, 2014, arXiv:1407.4465](http://adsabs.harvard.edu/abs/2014arXiv1407.4465B).

## License

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
