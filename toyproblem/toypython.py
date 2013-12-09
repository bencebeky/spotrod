import numpy;

def circleangleloop(r, p, z):
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
  answer = numpy.zeros_like(r);
  pminusz = p-z;
  pplusz = p+z;
  zsquared = z*z;
  psquared = p*p;
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

def circleanglemask(r, p, z):
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
  inside = (r < p-z);
  intersect = (r < p+z) & (z < r+p) & numpy.logical_not(inside);
  answer = numpy.zeros_like(r);
  answer[inside] = numpy.pi;
  answer[intersect] = numpy.arccos((numpy.power(r[intersect],2)+z*z-p*p)/(2*z*r[intersect]));
  return answer;
