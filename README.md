# spotrod: transit lightcurves of spotted stars

This repository contains `spotrod`, a semi-analytic model to calculate
lightcurves of planetary transits of spotted stars. The model is implemented in
pure Python.

## Installation

Install the latest version of `spotrod` using `pip`:

```
pip install spotrod
```

## Usage

```python
import matplotlib.pyplot as plt
import numpy as np
import spotrod

time = np.linspace(-0.1, 0.1, 200)
eta, xi = spotrod.elements(time, period=5, a=14, k=0.0, h=0.0)

planetx = np.zeros(time.size)
planety = -xi

r = np.linspace(0.0005, 0.9995, 1000)
p = 0.05
z = np.sqrt(planetx * planetx + planety * planety)
planetangle = np.array([spotrod.circleangle(r, p, z[i]) for i in range(z.size)])

f = np.ones(r.size)
spotx = np.array([0.05])
spoty = np.array([0.35])
spotradius = np.array([0.2])
spotcontrast = np.array([0.5])
prediction = spotrod.integratetransit(
    planetx, planety, z, p, r, f, spotx, spoty, spotradius, spotcontrast, planetangle
)

plt.plot(time, prediction)
```

See [examples](examples) directory for more examples.

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
