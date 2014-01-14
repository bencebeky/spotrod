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

#include "toyc.h"
#include <math.h>

void circleangleloop(double *r, double p, double z, int n, double *answer) {
/* Calculate half central angle of the arc of circle of radius r
   (which concentrically spans the inside of the star during integration)
   that is inside a circle of radius p (planet)
   with separation of centers z.
   This is a zeroth order homogeneous function, that is,
   circleangle(alpha*r, alpha*p, alpha*z) = circleangle(r, p, z).

   This version uses a loop over r.

   Input:
     r[n] array
     p    scalar
     z    scalar
   They should all be non-negative, but there is no other restriction.

     n    number of elements

   Output:
     answer[n]  one dimensional array, same size as r. */
  /* If the circle arc of radius r is disjoint from the circular disk 
     of radius p, then the angle is zero. */
  int i;
  double pminusz = p-z;
  double pplusz = p+z;
  double zsquared = z*z;
  double psquared = p*p;
  double ri;
  for(i=0; i<n; i++) {
    ri = *(r+i);
    // If the planet entirely covers the circle, the half central angle is pi.
    if (ri <= pminusz)
      *(answer+i) = M_PI;
    /* If the triangle inequalities hold between z, r, and p, 
    then we have partial overlap. If alpha is the half central angle
    in the triangle with sides r, p, and z,
    with p opposite the angle, then p^2 = r^2 + z^2 - 2 rz cos(alpha). */
    else if ((ri < pplusz) && (ri > -pminusz))
      *(answer+i) = acos((ri*ri+zsquared-psquared)/(2*z*ri));
    else
      *(answer+i) = 0;
  }
  return;
}

void circleanglesorted(double *r, double p, double z, int n, double *answer) {
/* Calculate half central angle of the arc of circle of radius r
   (which concentrically spans the inside of the star during integration)
   that is inside a circle of radius p (planet)
   with separation of centers z.
   This is a zeroth order homogeneous function, that is,
   circleangle(alpha*r, alpha*p, alpha*z) = circleangle(r, p, z).

   This version uses a binary search. It might, however, not be faster
   than using direct comparisons: in the loop, we need to compare 
   i to a and b and n (or 0) all the times. Is integer comparison
   faster than double? Is it worth the overhead of the binary search?

   Input:
     r[n] array, must be sorted.
     p    scalar
     z    scalar
   They should all be non-negative, but there is no other restriction.

     n    number of elements

   Output:
     answer[n]  one dimensional array, same size as r. */
  /* If the circle arc of radius r is disjoint from the circular disk 
     of radius p, then the angle is zero. */
  int i, a, b;
  double zsquared = z*z;
  double psquared = p*p;
  double ri;
  if (p > z) {
    // Planet covers center of star.
    a = mybsearch(r, p-z, n);
    b = mybsearch(r, p+z, n);
    for(i=0; i<a; i++)
      *(answer+i) = M_PI;
    for(; i<b; i++) {
      ri = *(r+i);
      *(answer+i) = acos((ri*ri+zsquared-psquared)/(2*z*ri));
    }
    for(; i<n; i++)
      *(answer+i) = 0.0;
  } else {
    // Planet does not cover center of star.
    a = mybsearch(r, z-p, n);
    b = mybsearch(r, z+p, n);
    for(i=0; i<a; i++)
      *(answer+i) = 0.0;
    for(; i<b; i++) {
      ri = *(r+i);
      *(answer+i) = acos((ri*ri+zsquared-psquared)/(2*z*ri));
    }
    for(; i<n; i++)
      *(answer+i) = 0.0;
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
