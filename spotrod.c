/* Copyright 2013, 2014 Bence Béky

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
along with Spotrod.  If not, see <http://www.gnu.org/licenses/>. */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "spotrod.h"

void integratetransit(int m, int n, int k, double *planetx, double *planety, double *z, double p, double *r, double *f, double *spotx, double *spoty, double *spotradius, double *spotcontrast, double *planetangle, double *answer) {
/* Calculate integrated flux of a star if it is transited by a planet
  of radius p*R_star, at projected position (planetx, planety)
  in R_star units.
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
  spotradius    spot radius in stellar radii                                                     [k]
  spotcontrast  spot contrast                                                                    [k]
  planetangle   value of [circleangle(r, p, z[i]) for i in xrange(m)]                 (cached) [m,n]

  (cached) means the parameter is redundant, and could be calculated from other parameters,
  but storing it and passing it to this routine speeds up iterative execution (fit or MCMC).
  Note that we do not take limb darkening coefficients, all we need is f.

  Output parameters:

  answer        model lightcurve, with oot=1.0                                                [m] */
  // Running indices for m, n, k, respectively.
  int M, N, K;
  // Out of transit flux to normalize lightcurve with.
  double ootflux;
  // Temporary storage for trapeze area to save on multiplications.
  double trapeze;
  // Temporary storage for what it is called.
  double spotcenterdistancesquared;
  // Cache for spot properties.
  double *spotcenterdistance, *spotangle;
  // Cache for trapezoid integration.
  double *values;
  // Cache for planet-spot center distance and central angle.
  double d, planetspotangle;
  /* Projected distance of center of star and center of spot, in stellar radius,
  accounting for the spot boundary plane being closer to the sphere radius
  than the tangent plane at spot center. An array of length k. */
  // If we have no spot:
  if (k == 0) {
    // Evaluate the integral for ootflux using the trapezoid method.
    ootflux = 0.0;
    for (N=0; N<n; N++) {
      ootflux += M_PI * *(r+N) * *(f+N);
    }
    for (M=0; M<m; M++) {
      // Transit?
      if (*(z+M) < 1.0 + p) {
        // Integrate over rings.
        *answer = 0.0;
        for (N=0; N<n; N++) {
          *answer += *(r+N) * (M_PI - *(planetangle + n*M + N)) * *(f+N);
        }
        // Normalize by trapezoid width and by ootflux.
        *answer /= ootflux;
      } else {
        *answer = 1.0;
      }
      // Advance pointer instead of using *(answer+i) all the times, for performance.
      answer++;
    }
    // Restore pointer.
    answer -= m;
  } else {
    //spotcenterdistance = malloc(k * sizeof(double));
    //spotangle = malloc(k * n * sizeof(double));
    //values = malloc(n * sizeof(double));
    // Instead, do a single malloc for speed.
    spotcenterdistance = malloc((k + k*n + n) * sizeof(double));
    spotangle = spotcenterdistance + k;
    values = spotcenterdistance + k * (n+1);
    // Loop over spots: fill up some arrays that do not depend on z.
    for (K=0; K<k; K++) {
      spotcenterdistancesquared = (*(spotx+K) * *(spotx+K) + *(spoty+K) * *(spoty+K)) * (1.0 - *(spotradius+K) * *(spotradius+K));
      *(spotcenterdistance+K) = sqrt(spotcenterdistancesquared);
      /* Calculate the half central angles of the spot, and store it in a single row of the 2D array.
      These values do not depend on z, that's why we cache them for all spots. */
      ellipseangle(r, *(spotradius+K), *(spotcenterdistance+K), n, spotangle + K*n);
    }
    // Evaluate the integral for ootflux using the trapezoid method.
    ootflux = 0.0;
    for (N=0; N<n; N++) {
      trapeze = M_PI;
      for (K=0; K<k; K++) {
        trapeze += (*(spotcontrast+K)-1.0) * *(spotangle + K*n + N);
      }
      ootflux += trapeze * *(r+N) * *(f+N);
    }
    // Loop over observation times and calculate answer.
    for (M=0; M<m; M++) {
      // Transit?
      if (*(z+M) < 1.0 + p) {
        /* The values array first stores the half arc length integrated contrast for no spots
        that we calculate in the next three lines. Later we add the spot contributions below,
        finally multipy by r*f and sum up. */
        for (N=0; N<n; N++) {
          *(values+N) = M_PI - *(planetangle + n*M + N);
        }
        // Cycle through spots and add their contributions.
        for (K=0; K<k; K++) {
          // Calculate distance of spot center and planet center for this moment.
          d = sqrt(pow(*(planetx+M) - *(spotx+K) * sqrt(1.0 - *(spotradius+K) * *(spotradius+K)), 2.0) + pow(*(planety+M) - *(spoty+K) * sqrt(1.0 - *(spotradius+K) * *(spotradius+K)), 2.0));
          // Calculate central angle between planet and spot.
          if ((*(spotcenterdistance+K) == 0) || (*(z+M) == 0)) {
            planetspotangle = 0.0;
          } else {
            planetspotangle = acos((pow(*(z+M),2.0) + pow(*(spotcenterdistance+K),2.0) - d*d)/(2.0 * *(z+M) * *(spotcenterdistance+K)));
          }
          // Cycle through annuli.
          for (N=0; N<n; N++) {
            /* Evaluate the integrand on r[N].
            The geometry is described by planetspotangle (function of z(M))
            planetangle (function of r(N) and z(M), 2D array),
            and spotangle (function of r(N), does not depend on z, calculated above).
            For each value of r, there are five cases. */
            // Case 1: planet and spot arcs are disjoint, contributions add up.
            if (planetspotangle > *(planetangle + M*n + N) + *(spotangle + K*n + N)) {
              *(values+N) += (*(spotcontrast+K)-1.0) * *(spotangle + K*n + N);
            // Case 2: planet arc inside spot arc.
            } else if (*(spotangle + K*n + N) > planetspotangle + *(planetangle + n*M + N)) {
              *(values+N) += (*(spotcontrast+K)-1.0) * (*(spotangle + K*n + N) - *(planetangle + n*M + N));
            // Case 4: triangle inequality holds, partial overlap.
            } else if (*(planetangle + n*M + N) <= planetspotangle + *(spotangle + K*n + N)) {
              // Case 4a: partial overlap on one side only.
              if (2*M_PI - planetspotangle >= *(planetangle + n*M + N) + *(spotangle + K*n + N)) {
                *(values+N) += 0.5 * (*(spotcontrast+K)-1.0) * (*(spotangle + K*n + N) + planetspotangle - *(planetangle + n*M + N));
              // Case 4b: partial overlap on two sides.
              } else {
                *(values+N) += (*(spotcontrast+K)-1.0) * (M_PI - *(planetangle + n*M + N));
              }
            }
            // Case 3: planet arc covers spot arc, spot arc invisible. No need to do anything.
            //else
              //*(values+N) += 0.0;
          }
        }
        /* Now we multiply the half arc length integrated contrast by r*f
        to get the integrand, and sum it up right away. */
        *answer = 0.0;
        for (N=0; N<n; N++) {
          *answer += *(r+N) * *(f+N) * *(values+N);
        }
        // Finally, we normalize with ootflux.
        *answer /= ootflux;
      } else {
        // If not transit:
        *answer = 1.0;
      }
      // Advance pointer instead of using *(answer+i) all the times, for performance.
      answer++;
    }
    // Restore pointer.
    answer -= m;
    // Free spotcenterdistance, spotangle, and values, that live side by side.
    free(spotcenterdistance);
  }
  return;
}

void elements(double *deltaT, double period, double a, double k, double h, int n, double *eta, double *xi) {
/* Calculate orbital elements eta and xi.

  Input:

  deltaT   time minus midtransit epoch [n]
  period   planetary period
  a        semimajor axis
  k, h     e cos omega, e sin omega respectively, (omega is periastron epoch)
  n        lenght of array deltaT

  Output:

  eta, xi  eta and xi at times deltaT, [n] */
  // Eccentricity and oblateness.
  double e = sqrt(k*k+h*h);
  double l = 1 - sqrt(1-k*k-h*h);
  //
  // ke = k cos E - h sqrt(1-e^2) sin E
  // ke / sqrt(k^2+h^2(1-e^2)) = sin a cos E + cos a sin E
  // ke / sqrt(k^2+h^2(1-e^2)) = sin (a+E)
  //
  // omega  k  h   E     a
  //  0     1  0   pi/2  pi/2
  //  pi/2  0  1   0     pi
  //  pi   -1  0   -pi/2 -pi/2
  //  -pi/2 0 -1   pi    0
  //
  // Aux stuff.
  double Mdot = 2*M_PI/period;
  double atan2hk = atan2(h,k);
  // Calculate eccentric anomaly and mean anomaly at midtransit.
  double Emid = M_PI - asin(k*e / sqrt(k*k + h*h * (1-e*e))) - atan2(k, - h * sqrt(1-k*k-h*h));
  double Mmid = Emid - e * sin(Emid);
  // Ten second tolerance (assuming time unit is day).
  double tol = 10.0 * M_PI / (43200.0 * period);
  // Now calculate eta and xi throughout the orbit.
  int i;
  double M, lam, E, Enew, p;
  for(i=0; i<n; i++) {
    // Mean anomaly.
    M = Mdot * *(deltaT+i) + Mmid;
    // Lambda = M + omega.
    lam = M + atan2hk;
    // Mean anomaly is the initial guess for eccentric anomaly.
    E = M;
    Enew = E - (E-e*sin(E)-M)/(1.0-e*cos(E));
    while (fabs(Enew-E) > tol) {
      E = Enew;
      Enew = E - (E-e*sin(E)-M)/(1.0-e*cos(E));
    }
    E = Enew;
    p = e * sin(E);
    // In top view, eta is the coordinate towards the observer,
    // and xi is the perpendicular one.
    *(eta+i) = a * (sin(lam + p) - k * p/(2.0-l) - h);
    *(xi+i)  = a * (cos(lam + p) + h * p/(2.0-l) - k);
  }
  return;
}

void circleangle(double *r, double p, double z, int n, double *answer) {
/* Calculate half central angle of the arc of circle of radius r
  (which concentrically spans the inside of the star during integration)
  that is inside a circle of radius p (planet)
  with separation of centers z.
  This is a zeroth order homogeneous function, that is,
  circleangle(alpha*r, alpha*p, alpha*z) = circleangle(r, p, z).

  This version uses a binary search. It is only marginally faster
  than using direct comparisons: in the loop, we need to compare 
  i to a and b and n (or 0) n times. A direct loop with comparisons
  would be more expensive by one indirect addressing and by a double
  comparison instead of the integer one (assuming, of course, that the
  input array is sorted, and we only compare to the value that delimits
  the next case, not testing for all cases in each iteration). This is
  barely worth the overhead of the binary search, but is so much cooler.

  Input:
    r       array, non-negative, non-decreasing [n]
    p       scalar, non-negative
    z       scalar, non-negative
    n       number of elements

  Output:
    answer  one dimensional array [n] */
  /* For speed, we increment answer and r (possibly faster than
  adding i in each iteration), and restore it at the end. The compiler
  should not actually compile the restore operation if the value of
  the variable is not used any more. */
  /* If the circle arc of radius r is disjoint from the circular disk 
     of radius p, then the angle is zero. */
  int i, a, b;
  double zsquared = z*z;
  double psquared = p*p;
  if (p > z) {
    // Planet covers center of star.
    a = mybsearch(r, p-z, n);
    b = mybsearch(r, p+z, n);
    for(i=0; i<a; i++) {
      *answer = M_PI;
      answer++;
    }
    r += i;
    for(; i<b; i++) {
      *answer = acos(((*r)*(*r)+zsquared-psquared)/(2*z*(*r)));
      answer++;
      r++;
    }
    r -= i;
    for(; i<n; i++) {
      *answer = 0.0;
      answer++;
    }
    answer -= i;
  } else {
    // Planet does not cover center of star.
    a = mybsearch(r, z-p, n);
    b = mybsearch(r, z+p, n);
    for(i=0; i<a; i++) {
      *answer = 0.0;
      answer++;
    }
    r += i;
    for(; i<b; i++) {
      *answer = acos(((*r)*(*r)+zsquared-psquared)/(2*z*(*r)));
      answer++;
      r++;
    }
    r -= i;
    for(; i<n; i++) {
      *answer = 0.0;
      answer++;
    }
    answer -= i;
  }
  return;
}

void ellipseangle(double *r, double a, double z, int n, double *answer) {
/*  Calculate half central angle of the arc of circle of radius r
  (which concentrically spans the inside of the star during integration)
  that is inside an ellipse of semi-major axis a with separation of centers z.
  The orientation of the ellipse is so that the center of the circle lies on 
  the continuation of the minor axis. This is the orientation if the ellipse
  is a circle on the surface of a sphere viewed in projection, and the circle
  is concentric with the projection of the sphere.
  b is calculated from a and z, assuming projection of a circle of radius a
  on the surface of a unit sphere. If a and z are not compatible, a is clipped.
  This is not zeroth order homogeneous function, because it calculates b based
  on a circle of radius a living on the surface of the unit sphere.
  r is an array, a, and z are scalars. They should all be non-negative.
  We store the result on the n double positions starting with *answer.
  
  Input:

  r        radius of circle [n]
  a        semi-major axis of ellipse, non-negative
  z        distance between centers of circle and ellipse,
           non-negative and at most 1
  n        size of array a

  Output:

  answer   half central angle of arc of circle that lies inside ellipes [n]. */
  int i;
  // Degenerate case.
  if ((a<=0.0) || (z>=1.0)) {
    for (i=0; i<n; i++) {
      *answer = 0.0;
      answer++;
    }
    answer -= i;
  } else if (z<=0) {
  // Concentric case.
    int bound = mybsearch(r, a, n);
    for(i=0; i<bound; i++) {
      *answer = M_PI;
      answer++;
    }
    for(; i<n; i++) {
      *answer = 0.0;
      answer++;
    }
    answer -= i;
  } else if (a*a + z*z >= 1.0) {
  // Unphysical case: clip a. Now b=0.
    a = sqrt(1-z*z);
    int bound = mybsearch(r, z, n);
    for(i=0; i<bound; i++) {
      *answer = 0.0;
      answer++;
    }
    r += i;
    for(; i<n; i++) {
      *answer = acos(z/(*r));
      answer++;
      r++;
    }
    answer -= i;
    r -= i;
  } else {
  // Case of ellipse.
    double b = a * sqrt(1.0-z*z/(1-a*a));
    double zsquared = z*z;
    double asquared = a*a;
    // Calculate A based on a and z to mitigate rounding errors.
    //double A = pow(a/b,2.0) - 1.0;
    double A = zsquared / (1.0 - asquared - zsquared);
    /* First, go through all the cases where the ellipse covers C_r.
    If there is no such r, then bound=0, nothing happens. */
    int bound = mybsearch(r, b-z, n);
    for(i=0; i<bound; i++) {
      *answer = M_PI;
      answer++;
    }
    /* Now go through all the cases when C_r does not reach out to the ellipse.
    Again, the loop body might not get executed. */
    bound = mybsearch(r, z-b, n);
    for(; i<bound; i++) {
      *answer = 0.0;
      answer++;
    }
    /* Now take care of the cases where C_r and the ellipse intersect. z_crit = 1-a^2. */
    if (z < 1.0 - a*a)
      bound = mybsearch(r, z+b, n);
    else
      bound = n;
    r += i;
    for(; i<bound; i++) {
      // If not disjoint from the outside.
      // We must have y_+ > -b now.
      //yp = (-z + sqrt(z*z - A*((*r)*(*r) - zsquared - asquared)))/A;
      //*answer = acos((z - yp)/(*r));
      *answer = acos((z - (-z + sqrt(z*z - A*((*r)*(*r) - zsquared - asquared)))/A)/(*r));
      answer++;
      r++;
    }
    r -= i;
    for(; i<n; i++) {
      *answer = 0.0;
      answer++;
    }
    answer -= i;
  }
  return;
}

int mybsearch(double *array, double val, int n) {
/* Return the smallest of i = 0, 1, ..., n such that val < *(array+i),
with the convention that *(array+n) = inf. */
  if (val < *array)
    return 0;
  if (*(array+n-1) <= val)
    return n;
  int a = 0, b = n-1, c;
  /* From this point on, we always have
  0 <= a <= c <= b <= n-1,
  *(array+a) <= val, and val < *(array+b). */
  while (a+1 < b) {
    c = (a+b)/2;
    if (val < *(array+c))
      b = c;
    else
      a = c;
  }
  return b;
}

void ellipseangle_original(double *r, double a, double z, int n, double *answer) {
/*  Calculate half central angle of the arc of circle of radius r
  (which concentrically spans the inside of the star during integration)
  that is inside an ellipse of semi-axes a and b with separation of centers z.
  b is calculated from a and z, assuming projection of a circle of radius a
  on the surface of a unit sphere.
  The orientation of the ellipse is so that the center of the circle lies on 
  the continuation of the minor axis. This is the orientation if the ellipse
  is a circle on the surface of a sphere viewed in projection, and the circle
  is concentric with the projection of the sphere.
  This is a zeroth order homogeneous function, that is,
  ellispeangle(alpha*r, alpha*a, alpha*z) = ellipseangle(r, a, z).
  r is an array, a, and z are scalars. They should all be non-negative.
  We store the result on the n double positions starting with *answer.
  
  Input:

  r        radius of circle [n]
  a        semi-major axis of ellipse
  b        semi-minor axis of ellipse
  z        distance between centers of circle and ellipse
           (center of circle lies on the straight line
           of the minor axis of the ellipse)
  n        size of array a

  Output:

  answer   half central angle of arc of circle that lies inside ellipes [n]. */
  int i;
  // Degenerate case.
  if ((a==0)) {
    for (i=0; i<n; i++) {
      *(answer+i) = 0.0;
    }
  // Concentric case.
  } else if (z==0) {
    for (i=0; i<n; i++) {
      if (*(r+i) < a)
        *(answer+i) = M_PI;
      else
        *(answer+i) = 0.0;
    }
  // Case of ellipse.
  } else {
    double b = a * sqrt(1.0-z*z/(1-a*a));
    double bminusz = b-z;
    double zsquared = z*z;
    double asquared = a*a;
    double A = pow(a/b,2.0) - 1.0;
    double ri, yp, halfD;
    for (i=0; i<n; i++) {
      ri = *(r+i);
      // If the ellipse entirely covers the circle, the half central angle is pi.
      if (ri <= bminusz) {
        *(answer+i) = M_PI;
      // If not disjoint from the outside.
      } else if (-bminusz < ri) {
        // Try to solve for y_+.
        halfD = z*z - A*(ri*ri - zsquared - asquared);
        // Discriminant negative: no intersection points.
        if (halfD < 0.0) {
          *(answer+i) = 0.0;
        } else {
          yp = (-z + sqrt(halfD))/A;
          // If y_+ > -b, then we have real intersection points.
          if (yp > -b) {
            *(answer+i) = acos((z - yp)/ri);
          } else {
          // If y_+ <= -b, then C_r contains the ellipse.
            *(answer+i) = 0.0;
          }
        }
      // Otherwise, they are disjoint from the outside.
      } else {
        *(answer+i) = 0.0;
      }
    }
  }
  return;
}
