#!/usr/bin/env python
# coding=utf8
#
# Copyright 2025 Bence Béky
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

import os
import pyfits
import numpy as np
import matplotlib.pyplot as plt

filename = "kplr010748390-2011145075126_slc.fits"
url = "http://archive.stsci.edu/missions/kepler/lightcurves/0107/010748390/" + filename

# If file does not exist, try downloading it.
if not os.path.exists(filename):
    imported_urllib = True
    imported_urllib2 = True
    imported_requests = True
    try:
        import urllib
    except ImportError:
        imported_urllib = False
    try:
        import urllib2
    except ImportError:
        imported_urllib2 = False
    try:
        import requests
    except ImportError:
        imported_requests = False
    if imported_urllib:
        urllib.urlretrieve(url, filename)
    elif imported_urllib2:
        response = urllib2.urlopen(url)
        with open(filename, "wb") as f:
            f.write(response.read())
    elif imported_requests:
        r = requests.get(url)
        with open(filename, "wb") as f:
            f.write(r.content)
    else:
        print("Please download {0:s} manually.".format(filename))

# Read Kepler short cadence data.
with pyfits.open(filename) as fitsfile:
    timebkjd = fitsfile[1].data.field("TIME")
    flux = fitsfile[1].data.field("SAP_FLUX")
    fluxerr = fitsfile[1].data.field("SAP_FLUX_ERR")

# Only keep transit 28 with some out of transit data.
mask = np.abs(timebkjd - midtransit - 28 * period) < 0.12 * period
mask &= np.logical_not(np.isnan(flux))
timebkjd = timebkjd[mask]
flux = flux[mask]
fluxerr = fluxerr[mask]
del mask

# Identify out of transit data for detrending.
phase = np.mod((timebkjd - midtransit) / period + 0.5, 1.0) - 0.5
oot = np.abs(phase) > 0.01002

# Detrend.
p = np.polyfit(timebkjd[oot], flux[oot], deg=2)
correction = np.polyval(p, timebkjd)
flux /= correction
fluxerr /= correction
del correction

# Throw out some more out of transit data, as we will not need them any more.
mask = np.abs(timebkjd - midtransit - 28 * period) < 0.02 * period
timebkjd = timebkjd[mask]
flux = flux[mask]
fluxerr = fluxerr[mask]
phase = phase[mask]
oot = oot[mask]
del mask

print(timebkjd)
print(flux)

plt.figure()
plt.plot(timebkjd, flux, "b.")
plt.show()
