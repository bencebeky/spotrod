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

# Perform parallel tempered mcmc on a single transit of HAT-P-11.
# This might take a dozen minutes or so.

import os;
import pyfits;
import numpy;
import spotrod;
import pickle;
from emcee import PTSampler;
from matplotlib import pyplot;
from matplotlib import animation;

# Set transit parameters.
# Fit values based on Kepler data.
period = 4.8878026;
periodhour = 24.0*period;
midtransit = 726.01298;
rp = 0.05866;
semimajoraxis = 14.17;
impactparam = 0.203;
u1 = 0.652;
u2 = 0.038;
# Values from Bakos et al. 2010, based on RV fit.
k = 0.216;
h = 0.133;

filename = "kplr010748390-2011145075126_slc.fits";
url = "http://archive.stsci.edu/missions/kepler/lightcurves/0107/010748390/" + filename;

# If file does not exist, try downloading it.
if not os.path.exists(filename):
  imported_urllib = True;
  imported_urllib2 = True;
  imported_requests = True;
  try:
    import urllib;
  except ImportError:
    imported_urllib = False;
  try:
    import urllib2;
  except ImportError:
    imported_urllib2 = False;
  try:
    import requests;
  except ImportError:
    imported_requests = False;
  if imported_urllib:
    urllib.urlretrieve(url, filename);
  elif imported_urllib2:
    response = urllib2.urlopen(url);
    with open(filename, 'wb') as f:
      f.write(response.read());
  elif imported_requests:
    r = requests.get(url);
    with open(filename, "wb") as f:
      f.write(r.content);
  else:
    print "Please download {0:s} manually.".format(filename);

# Read Kepler short cadence data.
with pyfits.open(filename) as fitsfile:
  timebkjd = fitsfile[1].data.field("TIME");
  flux = fitsfile[1].data.field("SAP_FLUX");
  fluxerr = fitsfile[1].data.field("SAP_FLUX_ERR");

# Only keep transit 28 with some out of transit data.
mask = numpy.abs(timebkjd-midtransit-28*period) < 0.12 * period;
mask &= numpy.logical_not(numpy.isnan(flux));
timebkjd = timebkjd[mask];
flux = flux[mask];
fluxerr = fluxerr[mask];
del mask;

# Identify out of transit data for detrending.
phase = numpy.mod((timebkjd-midtransit)/period+0.5, 1.0)-0.5;
oot = numpy.abs(phase) > 0.01002;

# Detrend.
p = numpy.polyfit(timebkjd[oot], flux[oot], deg=2);
correction = numpy.polyval(p, timebkjd);
flux /= correction;
fluxerr /= correction;
del correction;

# Throw out some more out of transit data, as we will not need them any more.
mask = numpy.abs(timebkjd-midtransit-28*period) < 0.02 * period;
timebkjd = timebkjd[mask];
flux = flux[mask];
fluxerr = fluxerr[mask];
phase = phase[mask];
oot = oot[mask];
del mask;

# Calculate chi prefactor to be used with MCMC.
minusoneovertwofluxerrsquared = - 0.5 * numpy.power(fluxerr, -2.0);

# Quadratic limb darkening function, Claret et al. 2000.
# I(mu)/I(1) = 1 - a(1-mu) - b(1-mu)^2
def quadraticlimbdarkening(r, u1, u2):
  answer = numpy.zeros_like(r);
  mask = (r<=1.0);
  oneminusmu = 1.0 - numpy.sqrt(1.0 - numpy.power(r[mask],2));
  answer[mask] = 1.0 - u1 * oneminusmu - u2 * numpy.power(oneminusmu,2);
  return answer;

# Initialize spotrod.
# Number of intergration rings.
n = 1000;

# Midpoint rule for integration.
# Integration annulii radii.
r = numpy.linspace(1.0/(2*n), 1.0-1.0/(2*n), n);
# Weights: 2.0 times limb darkening times width of integration annulii.
f = 2.0 * quadraticlimbdarkening(r, u1, u2) / n;

# Alternative: trapeziod rule.
#r = numpy.linspace(0.0, 1.0, n);
#f = 2.0 * quadraticlimbdarkening(r, u1, u2) * numpy.append(numpy.append([0.5], numpy.repeat(1.0, n-2)), [0.5]) / (n-1);

# Calculate orbital elements.
eta, xi = spotrod.elements(timebkjd-midtransit, period, semimajoraxis, k, h);
# Planet coordinates in sky plane, in Rstar units.
planetx = impactparam*eta/semimajoraxis;
planety = -xi;
# Distance from center, same as $z$ in Mandel, Agol 2002.
z = numpy.sqrt(numpy.power(planetx,2) + numpy.power(planety,2));
# Calculate planetangle array.
planetangle = numpy.array([spotrod.circleangle(r, rp, z[i]) for i in xrange(z.shape[0])]);

# Prior for spot parameters: isotropic on the surface of the sphere.
logp = lambda p: -0.5 * numpy.sum(numpy.log((1.0 - numpy.power(p[0::4], 2.0) - numpy.power(p[1::4], 2.0))));

# Likelihood for spot parameters.
logl = lambda p: numpy.sum(numpy.power(spotrod.integratetransit(planetx, planety, z, rp, r, f, p[0::4], p[1::4], p[2::4], p[3::4], planetangle) - flux, 2.0) * minusoneovertwofluxerrsquared);

# We have one spot, therefore phase space is 4D.
ndim = 4;
# Number of temperatures.
ntemps = 10;
# Number of parallel walkers at each temperature.
nwalkers = 100;
# Number of iterations.
niter = 1000;
# Of which burn-in is the first
burnin = 500;

# Initial spot parameters.
# [spotx, spoty, spotradius, spotcontrast]
spot = numpy.array([0.204, 0.376, 0.096, 0.524]);
# Create 3D matrix for initial state for each temperature and walker.
p0 = numpy.repeat(spot[:,numpy.newaxis].T, ntemps*nwalkers, axis=0).reshape(ntemps, nwalkers, ndim);
# Randomize the initial states in a small neighborhood.
p0 += numpy.random.normal(scale=1e-3, size=p0.shape);

# Initialize sampler.
sampler = PTSampler(ntemps, nwalkers, ndim, logl, logp);

# Run sampler.
pos, prob, state = sampler.run_mcmc(p0, niter);

# Take a view of the T=0 chain.
zerotemp = sampler.chain[0];

# We take iterations at T=0 after burn-in as equilibrium
# distribution. With a 100 walkers, this is 1e4 points.
eq = zerotemp[:,burnin:,:].reshape([nwalkers*(niter-burnin), ndim]);

# Plot distribution of every possible pairs.
labels = ["spotx", "spoty", "spotradius", "spotcontrast"];
for ploti in range(ndim-1):
  for plotj in range(ploti+1,ndim):
    pyplot.figure()
    pyplot.plot(eq[:,ploti],eq[:,plotj],"b.");
    pyplot.xlabel(labels[ploti]);
    pyplot.ylabel(labels[plotj]);
    pyplot.savefig("equilibrium-{0:d}-{1:d}.png".format(ploti,plotj));

# Create an animation in anix and aniy indices.
anix = 0;
aniy = 1;
fig = pyplot.figure();
ax = pyplot.gca();
chainplot, = pyplot.plot(zerotemp[0,0,anix], zerotemp[0,0,aniy], "b.");
ax.set_xlim(numpy.min(zerotemp[:,:,anix]), numpy.max(zerotemp[:,:,anix]));
ax.set_ylim(numpy.min(zerotemp[:,:,aniy]), numpy.max(zerotemp[:,:,aniy]));
ax.set_xlabel(labels[anix]);
ax.set_ylabel(labels[aniy]);

def animate(i):
  min = numpy.max([0, i-50]);
  max = i;
  chainplot.set_data(zerotemp[:,min:max,anix], zerotemp[:,min:max,aniy]);
  # Align title to left, otherwise it would jitter
  # due to changing width of rendered digits.
  ax.set_title("Iterations {0:d}:{1:d}".format(min, max), horizontalalignment = "left");

ani = animation.FuncAnimation(fig, animate, frames=niter);
ani.save("mcmc.mp4", fps=30);
