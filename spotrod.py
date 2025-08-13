"""
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
"""

import numpy as np
from numpy.typing import NDArray
from typing import Tuple


def integratetransit(
    planetx: NDArray[np.float64],
    planety: NDArray[np.float64],
    z: NDArray[np.float64],
    p: np.float64,
    r: NDArray[np.float64],
    f: NDArray[np.float64],
    spotx: NDArray[np.float64],
    spoty: NDArray[np.float64],
    spotradius: NDArray[np.float64],
    spotcontrast: NDArray[np.float64],
    planetangle: NDArray[np.float64],
) -> NDArray[np.float64]:
    """
    Calculate integrated flux of a star if it is transited by a planet
    of radius p*R_star, at projected position (planetx, planety)
    measured in R_star units.
    Flux is normalized to out-of-transit flux.
    This algorithm works by integrating over concentric rings,
    the number of which is controlled by n. Use n=1000 for fair results.
    Planetx is the coordinate perpendicular to the transit chord
    normalized to stellar radius units, and planety is the one
    parallel to the transit chord, in a fashion such that it increases
    throughout the transit.
    We assume that the one-dimensional arrays spotx, spoty, spotradius
    and spotcontrast have the same length: the number of the spots.

    Input parameters:

    m             length of time series
    n             number of concentric rings
    k             number of spots
    planet[xy]    planetary center coordinates in stellar radii in sky-projected coordinate system [m]
    z             planetary center distance from stellar disk center in stellar radii     (cached) [m]
    p             planetary radius in stellar radii, scalar
    r             radii of integration annuli in stellar radii, non-decreasing            (cached) [n]
    f             2.0 * limb darkening at r[i] * width of annulus i                       (cached) [n]
    spotx, spoty  spot center coordinates in stellar radii in sky-projected coordinate system      [k]
    spotradius    spot radius in stellar radii [k]
    spotcontrast  spot contrast [k]
    planetangle   value of [circleangle(r, p, z[i]) for i in range(m)]                 (cached) [m,n]

    (cached) means the parameter is redundant, and could be calculated from
    other parameters, but storing it and passing it to this routine speeds up
    iterative execution (fit or MCMC). Note that we do not take limb darkening
    coefficients, all we need is f.

    Output parameters:

    answer        model lightcurve, with oot=1.0 [m]
    """

    # Number of instances
    m = planetx.size
    # Number of integration annulii
    n = r.size
    # Number of spots
    k = spotx.size

    assert planetx.shape == (m,)
    assert planety.shape == (m,)
    assert r.shape == (n,)
    assert f.shape == (n,)
    assert spotx.shape == (k,)
    assert spoty.shape == (k,)
    assert planetangle.shape == (m, n)

    if k == 0:
        ootflux = np.pi * np.sum(r * f)

        answer = np.ones_like(planetx)
        transit_mask = z < 1.0 + p
        transit_count = np.sum(transit_mask)
        answer[transit_mask] = (
            np.sum(
                (np.pi - planetangle[transit_mask, :])
                * np.tile(r, (transit_count, 1))
                * np.tile(f, (transit_count, 1)),
                axis=1,
            )
            / ootflux
        )

        return answer

    # `spotx` and `spoty` are the sky-projected coordinates of the point of the
    # surface of the star that is the center of the spot. However, the
    # projected spot is an ellipse with a center that is the center of the 3D
    # circular disk inside the star. `spotcenterdistance` is the distance of
    # this center from the center of the star in sky projection.
    spotcenterdistance = np.sqrt(
        (spotx * spotx + spoty * spoty) * (1.0 - spotradius * spotradius)
    )

    spotangle = np.empty((k, n))
    for spot in range(k):
        spotangle[spot, :] = ellipseangle(r, spotradius[spot], spotcenterdistance[spot])

    # darkening caused by each spot at each radis; shape (k, n)
    darkening = np.tile(np.expand_dims(spotcontrast - 1.0, axis=1), (1, n)) * spotangle
    ootflux = np.sum((np.pi + np.sum(darkening, axis=0)) * r * f)

    answer = np.ones(m)

    for t in range(m):
        if z[t] >= 1.0 + p:
            # no transit
            continue

        # half arc length not covered by planet
        values = np.pi - planetangle[t, :]

        # Distance squared of spot center and planet center, shape (k)
        dsquared = np.power(
            planetx[t] - spotx * np.sqrt(1.0 - spotradius * spotradius), 2.0
        ) + np.power(planety[t] - spoty * np.sqrt(1.0 - spotradius * spotradius), 2.0)

        for spot in range(k):
            # Calculate central angle between planet and spot.
            if spotcenterdistance[spot] == 0.0 or z[t] == 0:
                planetspotangle = 0.0
            else:
                planetspotangle = np.arccos(
                    (
                        np.power(z[t], 2.0)
                        + np.power(spotcenterdistance[spot], 2.0)
                        - dsquared[spot]
                    )
                    / (2.0 * z[t] * spotcenterdistance[spot])
                )

            for i in range(n):
                # Calculate the integrand at r[i]. The geometry is described by
                # planetspotangle, planetangle and spotangle.

                # Case 1: planet and spot arcs are disjoint, contributions add up.
                if planetspotangle > planetangle[t, i] + spotangle[spot, i]:
                    values[i] += (spotcontrast - 1.0) * spotangle[spot, i]

                # Case 2: planet arc inside spot arc.
                elif spotangle[spot, i] > planetspotangle + planetangle[t, i]:
                    values[i] += (spotcontrast[spot] - 1.0) * (
                        spotangle[spot, i] - planetangle[t, i]
                    )

                # Case 3: triangle inequality holds, partial overlap.
                elif planetangle[t, i] <= planetspotangle + spotangle[spot, i]:
                    # Case 3a: partial overlap on one side only.
                    if (
                        2.0 * np.pi - planetspotangle
                        >= planetangle[t, i] + spotangle[spot, i]
                    ):
                        values[i] += (
                            0.5
                            * (spotcontrast[spot] - 1.0)
                            * (spotangle[spot, i] + planetspotangle - planetangle[t, i])
                        )

                    # Case 3b: partial overlap on two sides.
                    else:
                        values[i] += (spotcontrast[spot] - 1.0) * (
                            np.pi - planetangle[t, i]
                        )

                # Case 4: planet arc covers spot arc, spot arc invisible. No need for
                # correction due to spot contrast.

        answer[t] = np.sum(r * f * values) / ootflux

    return answer


def elements(
    deltaT: NDArray[np.float64],
    period: np.float64,
    a: np.float64,
    k: np.float64,
    h: np.float64,
) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Calculate orbital elements eta and xi.

    Input:

    deltaT   time minus midtransit epoch [n]
    period   planetary period
    a        semimajor axis
    k, h     e cos omega, e sin omega respectively (omega is periastron epoch)

    Output:

    eta, xi  eta and xi at times deltaT [2, n]
    """

    n = deltaT.size
    assert deltaT.shape == (n,)

    # Eccentricity and oblateness.
    e = np.sqrt(k * k + h * h)
    l = 1.0 - np.sqrt(1.0 - k * k - h * h)

    if e == 0:
        # In circular case, phase zero is arbitrarily chosen as the point on orbit
        # towards the observer.
        lam = 0.5 * np.pi + 2 * np.pi * deltaT / period
        # In top view, eta is the coordinate towards the observer,
        # and xi is the perpendicular one.
        eta = a * np.sin(lam)
        xi = a * np.cos(lam)
        return eta, xi

    # ke = k cos E - h sqrt(1-e^2) sin E
    # ke / sqrt(k^2+h^2(1-e^2)) = sin a cos E + cos a sin E
    # ke / sqrt(k^2+h^2(1-e^2)) = sin (a+E)
    #
    # omega  k  h   E     a
    #  0     1  0   pi/2  pi/2
    #  pi/2  0  1   0     pi
    #  pi   -1  0   -pi/2 -pi/2
    #  -pi/2 0 -1   pi    0
    #

    Mdot = 2.0 * np.pi / period
    omega = np.arctan2(h, k)
    # Calculate eccentric anomaly and mean anomaly at midtransit.
    Emid = (
        np.pi
        - np.arcsin(k * e / np.sqrt(k * k + h * h * (1.0 - e * e)))
        - np.arctan2(k, -h * np.sqrt(1.0 - k * k - h * h))
    )
    Mmid = Emid - e * np.sin(Emid)
    # Ten second tolerance (assuming time unit is day).
    tol = 10.0 * np.pi / (43200.0 * period)
    # Now calculate eta and xi throughout the orbit.
    # mean anomaly
    M = Mdot * deltaT + Mmid
    # lambda
    lam = M + omega
    # Mean anomaly is the initial guess for eccentric anomaly.
    E = M * np.ones(n)
    for i in range(n):
        Enew = E[i] - (E[i] - e * np.sin(E[i]) - M[i]) / (1.0 - e * np.cos(E[i]))
        while abs(Enew - E[i]) > tol:
            E[i] = Enew
            Enew = E[i] - (E[i] - e * np.sin(E[i]) - M[i]) / (1.0 - e * np.cos(E[i]))
        E[i] = Enew
    p = e * np.sin(E)
    # In top view, eta is the coordinate towards the observer,
    # and xi is the perpendicular one.
    eta = a * (np.sin(lam + p) - k * p / (2.0 - l) - h)
    xi = a * (np.cos(lam + p) + h * p / (2.0 - l) - k)
    return eta, xi


def circleangle(
    r: NDArray[np.float64], p: np.float64, z: np.float64
) -> NDArray[np.float64]:
    """circleangle(r, p, z)

    Calculate half central angle of the arc of circle of radius r
    (which concentrically spans the inside of the star during integration)
    that is inside a circle of radius p (planet)
    with separation of centers z.
    This is a zeroth order homogeneous function, that is,
    circleangle(alpha*r, alpha*p, alpha*z) = circleangle(r, p, z).

    This version uses a binary search on the sorted r.

    Input:
      r  one dimensional numpy array, must be increasing
      p  scalar
      z  scalar
    They should all be non-negative, but there is no other restriction.

    Output:
      circleangle  one dimensional numpy array, same size as r
    """

    n = r.size
    assert r.shape == (n,)

    # If the circle arc of radius r is disjoint from the circular disk
    # of radius p, then the angle is zero.
    answer = np.empty(n)
    if p > z:
        # Planet covers center of star.
        a, b = np.searchsorted(r, [p - z, p + z], side="right")
        answer[:a] = np.pi
        answer[a:b] = np.arccos((r[a:b] * r[a:b] + z * z - p * p) / (2.0 * z * r[a:b]))
        answer[b:] = 0.0
    else:
        # Planet does not cover center of star.
        a, b = np.searchsorted(r, [z - p, z + p], side="right")
        answer[:a] = 0.0
        answer[a:b] = np.arccos((r[a:b] * r[a:b] + z * z - p * p) / (2.0 * z * r[a:b]))
        answer[b:] = 0.0
    return answer


def ellipseangle(
    r: NDArray[np.float64], a: np.float64, z: np.float64
) -> NDArray[np.float64]:
    """Calculate half central angle of the arc of circle of radius r
    (which concentrically spans the inside of the star during integration)
    that is inside an ellipse of semi-major axis a with separation of centers z.
    The orientation of the ellipse is so that the center of the circle lies on
    the continuation of the minor axis. This is the orientation if the ellipse
    is a circle on the surface of a sphere viewed in projection, and the circle
    is concentric with the projection of the sphere.
    b is calculated from a and z, assuming projection of a circle of radius a
    on the surface of a unit sphere. If a and z are not compatible, a is
    clipped. This is not zeroth order homogeneous function, because it
    calculates b based on a circle of radius a living on the surface of the unit
    sphere. r is an array, a, and z are scalars. They should all be
    non-negative. We store the result on the n double positions starting with
    *answer.

    Input:

    r        radius of circle [n]
    a        semi-major axis of ellipse, non-negative
    z        distance between centers of circle and ellipse,
             non-negative and at most 1

    Output:

    answer   half central angle of arc of circle that lies inside ellipes [n].
    """

    n = r.size
    assert r.shape == (n,)

    # Concentric case
    if z <= 0.0:
        z = 0.0
        answer = np.zeros(n)
        bound = np.searchsorted(r, a)
        answer[:bound] = np.pi

        return answer

    if a <= 0.0:
        return np.zeros(n)

    if z >= 1.0:
        z = 1.0

    zsquared = z * z

    # b = 0
    if a * a + zsquared >= 1.0:
        asquared = 1 - zsquared
        a = np.sqrt(asquared)
        answer = np.zeros(n)
        bound = np.searchsorted(r, z, side="right")
        answer[bound:] = np.arccos(z / r[bound:])
        return answer

    asquared = a * a
    b = a * np.sqrt(1.0 - zsquared / (1.0 - asquared))

    # Square of eccentricity is asquared / bsquared - 1
    A = zsquared / (1.0 - asquared - zsquared)

    answer = np.empty(n)
    bound1, bound2, bound3 = np.searchsorted(r, [b - z, z - b, b + z])

    # The ellipse encloses the cirle of radius r
    answer[:bound1] = np.pi

    # The ellipse and the circle are disjoint
    answer[:bound2] = 0.0

    # The ellipse and the circle intersect
    indices = slice(bound1, bound3) if bound1 > bound2 else slice(bound2, bound3)
    answer[indices] = np.arccos(
        (
            z
            - (
                -z
                + np.sqrt(
                    zsquared - A * (np.power(r[indices], 2.0) - zsquared - asquared)
                )
            )
            / A
        )
        / r[indices]
    )

    # The circle encloses the ellipse
    answer[bound3:] = 0.0

    return answer
