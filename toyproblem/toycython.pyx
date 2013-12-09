# cython: profile=False

import numpy;
cimport numpy;

ctypedef numpy.float_t FLOAT_t
ctypedef numpy.uint8_t BOOL_t

def circleangleloop(numpy.ndarray[FLOAT_t, ndim=1] r, double p, double z):
  """circleangleloop(r, p, z)

  Calculate half central angle of the arc of circle of radius r
  (which concentrically spans the inside of the star during integration)
  that is inside a circle of radius p (planet)
  with separation of centers z.
  This is a zeroth order homogeneous function, that is,
  circleangle(alpha*r, alpha*p, alpha*z) = circleangle(r, p, z).

  This version uses a loop over r.

  Input:
    r  one dimensional numpy array
    p  scalar
    z  scalar
  They should all be non-negative, but there is no other restriction.

  Output:
    circleangle  one dimensional numpy array, same size as r
  """
  # If the circle arc of radius r is disjoint from the circular disk 
  # of radius p, then the angle is zero.
  cdef numpy.ndarray[FLOAT_t, ndim=1] answer = numpy.zeros_like(r)
  cdef double pminusz = p-z
  cdef double pplusz = p+z
  cdef double zsquared = z*z
  cdef double psquared = p*p
  cdef double ri
  for i in xrange(r.shape[0]):
    ri = r[i];
    # If the planet entirely covers the circle, the half central angle is pi.
    if (ri <= pminusz):
      answer[i] = numpy.pi;
    # If the triangle inequalities hold between z, r, and p, 
    # then we have partial overlap.
    # If alpha is the half central angle in the triangle with sides r, p, and z,
    # with p opposite the angle, then p^2 = r^2 + z^2 - 2 rz cos(alpha)
    elif (ri < pplusz) & (ri > -pminusz):
      answer[i] = numpy.arccos((ri*ri+zsquared-psquared)/(2*z*ri));
  return answer;

def circleanglemask(numpy.ndarray[FLOAT_t, ndim=1] r, double p, double z):
  """circleanglemask(r, p, z)

  Calculate half central angle of the arc of circle of radius r
  (which concentrically spans the inside of the star during integration)
  that is inside a circle of radius p (planet)
  with separation of centers z.
  This is a zeroth order homogeneous function, that is,
  circleangle(alpha*r, alpha*p, alpha*z) = circleangle(r, p, z).

  This version uses masked array operations.

  Input:
    r  one dimensional numpy array
    p  scalar
    z  scalar
  They should all be non-negative, but there is no other restriction.

  Output:
    circleangle  one dimensional numpy array, same size as r
  """
  cdef numpy.ndarray[BOOL_t, ndim=1, cast=True] inside = (r < p-z)
  cdef numpy.ndarray[BOOL_t, ndim=1, cast=True] intersect = (r < p+z) & (z < r+p) & numpy.logical_not(inside)
  cdef numpy.ndarray[FLOAT_t, ndim=1] answer = numpy.zeros_like(r)
  answer[inside] = numpy.pi;
  answer[intersect] = numpy.arccos((numpy.power(r[intersect],2)+z*z-p*p)/(2*z*r[intersect]));
  return answer;
