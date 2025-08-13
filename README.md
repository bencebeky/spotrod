# spotrod: A semi-analytic model for transits of spotted stars

## Contents

This repository contains `spotrod`, a semi-analytic model to calculate
lightcurves of planetary transits of spotted stars. The model is implemented in
pure Python.

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

If you use `spotrod` in a publication, please cite the publication [Béky,
Kipping, and Holman, 2014, MNRAS, 442,
3686](https://ui.adsabs.harvard.edu/abs/2014MNRAS.442.3686B/abstract).

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
