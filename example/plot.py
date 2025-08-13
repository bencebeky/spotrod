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

"""This file shows how to use `spotrod` in a simple case."""

import matplotlib.pyplot as plt
import numpy as np
import spotrod
import kepler_data

# Transit parameters
period = 4.8878026
periodhour = 24.0 * period
midtransit = 726.01298
rp = 0.05866
semimajoraxis = 14.17
impactparam = 0.203

# Values from Bakos et al. 2010, based on RV fit.
k = 0.216
h = 0.133

# Limb darkening parameters
mu1 = 0.652
mu2 = 0.038

timebkjd = kepler_data.timebkjd
measurement = kepler_data.flux

phase = np.mod((timebkjd - midtransit) / period + 0.5, 1.0) - 0.5


# Quadratic limb darkening function, Claret et al. 2000.
# I(mu)/I(1) = 1 - a(1-mu) - b(1-mu)^2
def quadraticlimbdarkening(r, mu1, mu2):
    answer = np.zeros_like(r)
    mask = r <= 1.0
    oneminusmu = 1.0 - np.sqrt(1.0 - np.power(r[mask], 2))
    answer[mask] = 1.0 - mu1 * oneminusmu - mu2 * np.power(oneminusmu, 2)
    return answer


# Integration radii and weights (midpoint rule)
n = 1000
r = np.linspace(1.0 / (2 * n), 1.0 - 1.0 / (2 * n), n)
f = quadraticlimbdarkening(r, mu1, mu2)

# Orbital elements
eta, xi = spotrod.elements(timebkjd - midtransit, period, semimajoraxis, k, h)
planetx = impactparam * eta / semimajoraxis
planety = -xi
# Distance from center, same as `z` in Mandel, Agol 2002.
z = np.sqrt(np.power(planetx, 2) + np.power(planety, 2))

planetangle = np.array([spotrod.circleangle(r, rp, z[i]) for i in range(z.size)])

spotx = 0.204
spoty = 0.376
spotradius = 0.096
spotcontrast = 0.524
prediction = spotrod.integratetransit(
    planetx.size,
    n,
    1,
    planetx,
    planety,
    z,
    rp,
    r,
    f,
    np.array([spotx]),
    np.array([spoty]),
    np.array([spotradius]),
    np.array([spotcontrast]),
    planetangle,
)

plt.figure(figsize=(11, 7))
plt.plot(phase, measurement, "b.", label="Data")
plt.plot(phase, prediction, "k-", label="Model")
plt.legend()
plt.savefig("test.png")
plt.show()
