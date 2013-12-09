#include <math.h>

void circleangle(int n, double *r, double p, double z, double *answer) {
/* Calculate half central angle of the arc of circle of radius r
   (which concentrically spans the inside of the star during integration)
   that is inside a circle of radius p (planet)
   with separation of centers z.
   This is a zeroth order homogeneous function, that is,
   circleangle(alpha*r, alpha*p, alpha*z) = circleangle(r, p, z).

   This version uses a loop over r.

   Input:
     n    number of elements
     r[n] array
     p    scalar
     z    scalar
   They should all be non-negative, but there is no other restriction.

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
