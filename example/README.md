# Example usage of spotrod

[plot.py](plot.py) contains a full example of how to use the `spotrod` module.

[kepler_data.py](kepler_data.py) contains hardcoded data of flux observations of
a transit that seems to involve a large stellar spot.

[kepler_download.py](kepler_download.py) is provided for reference on how the
data was downloaded and processed, but the `pyfits` package has been merged into
`astropy.io.fits`, and this script has not been updated.

[mcmc.py](mcmc.py) is provided for reference on how the Monte-Carlo simulation
was run by Béky, Kipping, and Holman, 2014, but the `emcee` package has received
breaking API changes since and this script has not been updated.
