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

"""Perform parallel tempered mcmc on a single transit of HAT-P-11.
This might take a dozen minutes or so."""

from emcee import PTSampler
from matplotlib import animation
from matplotlib import pyplot
import kepler_data
import numpy as np
import spotrod

# Set transit parameters.
# Fit values based on Kepler data.
period = 4.8878026
periodhour = 24.0*period
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
flux = kepler_data.flux

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

# Prior for spot parameters: isotropic on the surface of the sphere.
logp = lambda p: -0.5 * numpy.sum(numpy.log((1.0 - numpy.power(p[0::4], 2.0) - numpy.power(p[1::4], 2.0))))

# Likelihood for spot parameters.
logl = lambda p: numpy.sum(numpy.power(spotrod.integratetransit(planetx, planety, z, rp, r, f, p[0::4], p[1::4], p[2::4], p[3::4], planetangle) - flux, 2.0))

# We have one spot, therefore phase space is 4D.
ndim = 4
# Number of temperatures.
ntemps = 10
# Number of parallel walkers at each temperature.
nwalkers = 100
# Number of iterations.
niter = 1000
# Of which burn-in is the first
burnin = 500

# Initial spot parameters.
# [spotx, spoty, spotradius, spotcontrast]
spot = numpy.array([0.204, 0.376, 0.096, 0.524])
# Create 3D matrix for initial state for each temperature and walker.
p0 = numpy.repeat(spot[:,numpy.newaxis].T, ntemps*nwalkers, axis=0).reshape(ntemps, nwalkers, ndim)
# Randomize the initial states in a small neighborhood.
p0 += numpy.random.normal(scale=1e-3, size=p0.shape)

# Initialize sampler.
sampler = PTSampler(ntemps, nwalkers, ndim, logl, logp)

# Run sampler.
pos, prob, state = sampler.run_mcmc(p0, niter)

# Take a view of the T=0 chain.
zerotemp = sampler.chain[0]

# We take iterations at T=0 after burn-in as equilibrium
# distribution. With a 100 walkers, this is 1e4 points.
eq = zerotemp[:,burnin:,:].reshape([nwalkers*(niter-burnin), ndim])

# Plot distribution of every possible pairs.
labels = ["spotx", "spoty", "spotradius", "spotcontrast"]
for ploti in range(ndim-1):
  for plotj in range(ploti+1,ndim):
    pyplot.figure()
    pyplot.plot(eq[:,ploti],eq[:,plotj],"b.")
    pyplot.xlabel(labels[ploti])
    pyplot.ylabel(labels[plotj])
    pyplot.savefig("equilibrium-{0:d}-{1:d}.png".format(ploti,plotj))

# Create an animation in anix and aniy indices.
anix = 0
aniy = 1
fig = pyplot.figure()
ax = pyplot.gca()
chainplot, = pyplot.plot(zerotemp[0,0,anix], zerotemp[0,0,aniy], "b.")
ax.set_xlim(numpy.min(zerotemp[:,:,anix]), numpy.max(zerotemp[:,:,anix]))
ax.set_ylim(numpy.min(zerotemp[:,:,aniy]), numpy.max(zerotemp[:,:,aniy]))
ax.set_xlabel(labels[anix])
ax.set_ylabel(labels[aniy])

def animate(i):
  min = numpy.max([0, i-50])
  max = i
  chainplot.set_data(zerotemp[:,min:max,anix], zerotemp[:,min:max,aniy])
  # Align title to left, otherwise it would jitter
  # due to changing width of rendered digits.
  ax.set_title("Iterations {0:d}:{1:d}".format(min, max), horizontalalignment = "left")

ani = animation.FuncAnimation(fig, animate, frames=niter)
ani.save("mcmc.mp4", fps=30)
