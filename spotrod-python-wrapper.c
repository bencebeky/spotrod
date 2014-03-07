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

#include <Python.h>
#include <numpy/arrayobject.h>
#include "spotrod.h"

/* Docstrings */
static char module_docstring[] = 
"  This module is a fast C implementation of\n"
"  the spotrod lightcurve model.\n"
"  Source is available at https://github.com/bencebeky/spotrod";
static char integratetransit_docstring[] = 
" answer = integratetransit(planetx, planety, z, p, r, f, spotx, spoty, spotradius, spotcontrast, planetangle)\n"
"\n"
"  Calculate integrated flux of a star if it is transited by a planet\n"
"  of radius p*R_star, at projected position (planetx, planety)\n"
"  in R_star units.\n"
"  Flux is normalized to out-of-transit flux.\n"
"  This algorithm works by integrating over concentric rings,\n"
"  the number of which is controlled by n. Use n=1000 for fair results.\n"
"  Planetx is the coordinate perpendicular to the transit chord\n"
"  normalized to stellar radius units, and planety is the one\n"
"  parallel to the transit chord, in a fashion such that it increases\n"
"  throughout the transit.\n"
"  We assume that the one-dimensional arrays spotx, spoty, spotradius\n"
"  and spotcontrast have the same length: the number of the spots.\n"
"\n"
"  Input parameters:\n"
"\n"
"  m             length of time series\n"
"  n             number of concentric rings\n"
"  k             number of spots\n"
"  planet[xy]    planetary center coordinates in stellar radii in sky-projected coordinate system [m]\n"
"  z             planetary center distance from stellar disk center in stellar radii     (cached) [m]\n"
"  p             planetary radius in stellar radii, scalar\n"
"  r             radii of integration annuli in stellar radii, non-decreasing            (cached) [n]\n"
"  f             2.0 * limb darkening at r[i] * width of annulus i                       (cached) [n]\n"
"  spotx, spoty  spot center coordinates in stellar radii in sky-projected coordinate system      [k]\n"
"  spotradius    spot radius in stellar radii                                                     [k]\n"
"  spotcontrast  spot contrast                                                                    [k]\n"
"  planetangle   circleangle(r[i], p, z[j])                                            (cached) [m,n]\n"
"\n"
"  (cached) means the parameter is redundant, and could be calculated from other parameters,\n"
"  but storing it and passing it to this routine speeds up iterative execution (fit or MCMC).\n"
"  Note that we do not take limb darkening coefficients, all we need is f.\n"
"\n"
"  Output parameters:\n"
"\n"
"  answer        model lightcurve, with oot=1.0                                                [m]";
static char elements_docstring[] =
" eta, xi = elements(deltaT, period, a, k, h)\n"
"\n"
"  Calculate orbital elements eta and xi.\n"
"\n"
"  Input:\n"
"\n"
"  deltaT   time minus midtransit epoch [n]\n"
"  period   planetary period\n"
"  a        semimajor axis\n"
"  k, h     e cos omega, e sin omega respectively, (omega is periastron epoch)\n"
"  n        lenght of array deltaT\n"
"\n"
"  Output:\n"
"\n"
"  eta, xi  eta and xi at times deltaT, [n]";
static char circleangle_docstring[] =
" answer = circleangle(r, p, z)\n"
"\n"
"  Calculate half central angle of the arc of circle of radius r\n"
"  (which concentrically spans the inside of the star during integration)\n"
"  that is inside a circle of radius p (planet)\n"
"  with separation of centers z.\n"
"  This is a zeroth order homogeneous function, that is,\n"
"  circleangle(alpha*r, alpha*p, alpha*z) = circleangle(r, p, z).\n"
"\n"
"  This version uses a binary search. It is only marginally faster\n"
"  than using direct comparisons: in the loop, we need to compare \n"
"  i to a and b and n (or 0) n times. A direct loop with comparisons\n"
"  would be more expensive by one indirect addressing and by a double\n"
"  comparison instead of the integer one (assuming, of course, that the\n"
"  input array is sorted, and we only compare to the value that delimits\n"
"  the next case, not testing for all cases in each iteration). This is\n"
"  barely worth the overhead of the binary search, but is so much cooler.\n"
"\n"
"  Input:\n"
"    r       array, non-negative, non-decreasing [n]\n"
"    p       scalar, non-negative\n"
"    z       scalar, non-negative\n"
"    n       number of elements\n"
"\n"
"  Output:\n"
"    answer  one dimensional array [n]";
static char ellipseangle_docstring[] =
" answer = ellipseangle(r, a, z)\n"
"\n"
"  Calculate half central angle of the arc of circle of radius r\n"
"  (which concentrically spans the inside of the star during integration)\n"
"  that is inside an ellipse of semi-major axis a with separation of centers z.\n"
"  The orientation of the ellipse is so that the center of the circle lies on \n"
"  the continuation of the minor axis. This is the orientation if the ellipse\n"
"  is a circle on the surface of a sphere viewed in projection, and the circle\n"
"  is concentric with the projection of the sphere.\n"
"  b is calculated from a and z, assuming projection of a circle of radius a\n"
"  on the surface of a unit sphere. If a and z are not compatible, a is clipped.\n"
"  This is not zeroth order homogeneous function, because it calculates b based\n"
"  on a circle of radius a living on the surface of the unit sphere.\n"
"  r is an array, a, and z are scalars. They should all be non-negative.\n"
"  We store the result on the n double positions starting with *answer.\n"
"  \n"
"  Input:\n"
"\n"
"  r        radius of circle [n]\n"
"  a        semi-major axis of ellipse, non-negative\n"
"  z        distance between centers of circle and ellipse,\n"
"           non-negative and at most 1\n"
"  n        size of array a\n"
"\n"
"  Output:\n"
"\n"
"  answer   half central angle of arc of circle that lies inside ellipes [n].";

/* Function wrappers for external use */
static PyObject *integratetransit_wrapper(PyObject*, PyObject*, PyObject*);
static PyObject *elements_wrapper(PyObject*, PyObject*, PyObject*);
static PyObject *circleangle_wrapper(PyObject*, PyObject*, PyObject*);
static PyObject *ellipseangle_wrapper(PyObject*, PyObject*, PyObject*);

/* Module specification */
static PyMethodDef module_methods[] = {
  {"integratetransit", (PyCFunction)integratetransit_wrapper, METH_VARARGS | METH_KEYWORDS, integratetransit_docstring},
  {"elements", (PyCFunction)elements_wrapper, METH_VARARGS | METH_KEYWORDS, elements_docstring},
  {"circleangle", (PyCFunction)circleangle_wrapper, METH_VARARGS | METH_KEYWORDS, circleangle_docstring},
  {"ellipseangle", (PyCFunction)ellipseangle_wrapper, METH_VARARGS | METH_KEYWORDS, ellipseangle_docstring},
  {NULL, NULL, 0, NULL}
};

/* Initialize the module */
PyMODINIT_FUNC initspotrod(void) {
  PyObject *m = Py_InitModule3("spotrod", module_methods, module_docstring);
  if (m == NULL)
    return;
  /* Load numpy functionality. */
  import_array();
}

/* Wrapper function for integratetransit. */
static PyObject *integratetransit_wrapper(PyObject *self, PyObject *args, PyObject *kwds) {
  /* Input arguments. */
  PyObject *planetx_obj, *planety_obj, *z_obj, *r_obj, *f_obj, *spotx_obj, *spoty_obj, *spotradius_obj, *spotcontrast_obj, *planetangle_obj;
  double p;

  // Keywords.
  static char *kwlist[] = {"planetx", "planety", "z", "p", "r", "f", "spotx", "spoty", "spotradius", "spotcontrast", "planetangle", NULL};

  /* Parse the input tuple */
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOOdOOOOOOO", kwlist, &planetx_obj, &planety_obj, &z_obj, &p, &r_obj, &f_obj, &spotx_obj, &spoty_obj, &spotradius_obj, &spotcontrast_obj, &planetangle_obj)) {
    PyErr_SetString(PyExc_ValueError, "Error parsing arguments.");
    return NULL;
  }

  /* Check argument dimensions and types. */
  if (PyArray_NDIM(planetx_obj) != 1 || PyArray_NDIM(planety_obj) != 1 || PyArray_NDIM(z_obj) != 1 || PyArray_NDIM(r_obj) != 1 || PyArray_NDIM(f_obj) != 1 || PyArray_NDIM(spotx_obj) != 1 || PyArray_NDIM(spoty_obj) != 1 || PyArray_NDIM(spotradius_obj) != 1 || PyArray_NDIM(spotcontrast_obj) != 1 || PyArray_NDIM(planetangle_obj) != 2 || PyArray_TYPE(planetx_obj) != PyArray_DOUBLE || PyArray_TYPE(planety_obj) != PyArray_DOUBLE || PyArray_TYPE(z_obj) != PyArray_DOUBLE || PyArray_TYPE(r_obj) != PyArray_DOUBLE || PyArray_TYPE(f_obj) != PyArray_DOUBLE || PyArray_TYPE(spotx_obj) != PyArray_DOUBLE || PyArray_TYPE(spoty_obj) != PyArray_DOUBLE || PyArray_TYPE(spotradius_obj) != PyArray_DOUBLE || PyArray_TYPE(spotcontrast_obj) != PyArray_DOUBLE || PyArray_TYPE(planetangle_obj) != PyArray_DOUBLE) {
    PyErr_SetString(PyExc_ValueError, "Argument dimensions or types not correct.");
    return NULL; 
  }

  /* Get dimensions. */
  int m = PyArray_DIM(planetx_obj, 0);
  int n = PyArray_DIM(r_obj, 0);
  int k = PyArray_DIM(spotx_obj, 0);

  /* Check argument shapes. */
  if (PyArray_DIM(planety_obj, 0) != m || PyArray_DIM(z_obj, 0) != m || PyArray_DIM(f_obj, 0) != n || PyArray_DIM(spoty_obj, 0) != k || PyArray_DIM(spotradius_obj, 0) != k || PyArray_DIM(spotcontrast_obj, 0) != k || PyArray_DIM(planetangle_obj, 0) != m || PyArray_DIM(planetangle_obj, 1) != n) {
    PyErr_SetString(PyExc_ValueError, "Argument shapes not correct.");
    return NULL; 
  }

  /* Interpret the input objects as numpy arrays. */
  PyObject *planetx_array = PyArray_FROM_OTF(planetx_obj, NPY_DOUBLE, NPY_IN_ARRAY);
  PyObject *planety_array = PyArray_FROM_OTF(planety_obj, NPY_DOUBLE, NPY_IN_ARRAY);
  PyObject *z_array = PyArray_FROM_OTF(z_obj, NPY_DOUBLE, NPY_IN_ARRAY);
  PyObject *r_array = PyArray_FROM_OTF(r_obj, NPY_DOUBLE, NPY_IN_ARRAY);
  PyObject *f_array = PyArray_FROM_OTF(f_obj, NPY_DOUBLE, NPY_IN_ARRAY);
  PyObject *spotx_array = PyArray_FROM_OTF(spotx_obj, NPY_DOUBLE, NPY_IN_ARRAY);
  PyObject *spoty_array = PyArray_FROM_OTF(spoty_obj, NPY_DOUBLE, NPY_IN_ARRAY);
  PyObject *spotradius_array = PyArray_FROM_OTF(spotradius_obj, NPY_DOUBLE, NPY_IN_ARRAY);
  PyObject *spotcontrast_array = PyArray_FROM_OTF(spotcontrast_obj, NPY_DOUBLE, NPY_IN_ARRAY);
  PyObject *planetangle_array = PyArray_FROM_OTF(planetangle_obj, NPY_DOUBLE, NPY_IN_ARRAY);

  /* If that didn't work, throw an exception. */
  if (planetx_array == NULL || planety_array == NULL || z_array == NULL || r_array == NULL || f_array == NULL || spotx_array == NULL || spoty_array == NULL || spotradius_array == NULL || spotcontrast_array == NULL || planetangle_array == NULL) {
    Py_XDECREF(planetx_array);
    Py_XDECREF(planety_array);
    Py_XDECREF(z_array);
    Py_XDECREF(r_array);
    Py_XDECREF(f_array);
    Py_XDECREF(spotx_array);
    Py_XDECREF(spoty_array);
    Py_XDECREF(spotradius_array);
    Py_XDECREF(spotcontrast_array);
    Py_XDECREF(planetangle_array);
    return NULL;
  }

  /* Get data pointers. */
  double *planetx = (double*)PyArray_DATA(planetx_array);
  double *planety = (double*)PyArray_DATA(planety_array);
  double *z = (double*)PyArray_DATA(z_array);
  double *r = (double*)PyArray_DATA(r_array);
  double *f = (double*)PyArray_DATA(f_array);
  double *spotx = (double*)PyArray_DATA(spotx_array);
  double *spoty = (double*)PyArray_DATA(spoty_array);
  double *spotradius = (double*)PyArray_DATA(spotradius_array);
  double *spotcontrast = (double*)PyArray_DATA(spotcontrast_array);
  double *planetangle = (double*)PyArray_DATA(planetangle_array);
  
  // Create answer numpy array, let Python allocate memory.
  PyArrayObject *answer = (PyArrayObject *)PyArray_FromDims(1, &m, NPY_DOUBLE);

  // Calculate answer.
  integratetransit(m, n, k, planetx, planety, z, p, r, f, spotx, spoty, spotradius, spotcontrast, planetangle, (double *)answer->data);

  /* Clean up. */
  Py_DECREF(planetx_array);
  Py_DECREF(planety_array);
  Py_DECREF(z_array);
  Py_DECREF(r_array);
  Py_DECREF(f_array);
  Py_DECREF(spotx_array);
  Py_DECREF(spoty_array);
  Py_DECREF(spotradius_array);
  Py_DECREF(spotcontrast_array);
  Py_DECREF(planetangle_array);

  // Return answer.
  return PyArray_Return(answer);
}

/* Wrapper function for elements. */
static PyObject *elements_wrapper(PyObject *self, PyObject *args, PyObject *kwds) {
  /* Input arguments. */
  PyObject *deltaT_obj;
  double period, a, k, h;

  // Keywords.
  static char *kwlist[] = {"deltaT", "period", "a", "k", "h", NULL};

  /* Parse the input tuple */
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "Odddd", kwlist, &deltaT_obj, &period, &a, &k, &h)) {
    PyErr_SetString(PyExc_ValueError, "Error parsing arguments.");
    return NULL;
  }

  /* Check argument dimensions and types. */
  if (PyArray_NDIM(deltaT_obj) != 1 || PyArray_TYPE(deltaT_obj) != PyArray_DOUBLE) {
    PyErr_SetString(PyExc_ValueError, "Argument dimensions or types not correct.");
    return NULL; 
  }

  /* Get dimensions. */
  int n = PyArray_DIM(deltaT_obj, 0);

  /* Interpret the input objects as numpy arrays. */
  PyObject *deltaT_array = PyArray_FROM_OTF(deltaT_obj, NPY_DOUBLE, NPY_IN_ARRAY);

  /* If that didn't work, throw an exception. */
  if (deltaT_array == NULL) {
    Py_XDECREF(deltaT_array);
    return NULL;
  }

  /* Get data pointers. */
  double *deltaT = (double*)PyArray_DATA(deltaT_array);
  
  // Create answer numpy arrays, let Python allocate memory.
  PyArrayObject *eta = (PyArrayObject *)PyArray_FromDims(1, &n, NPY_DOUBLE);
  PyArrayObject *xi = (PyArrayObject *)PyArray_FromDims(1, &n, NPY_DOUBLE);

  // Calculate answer.
  elements(deltaT, period, a, k, h, n, (double *)eta->data, (double *)xi->data);

  /* Clean up. */
  Py_DECREF(deltaT_array);

  /* Create answer tuple. */
  PyObject *answertuple = Py_BuildValue("(OO)", eta, xi);

  /* Now we have two references for the eta and xi numpy objects:
  one from when it was created, one when it was included in the tuple.
  As we only return the tuple, we need to decrement references. */
  Py_DECREF(eta);
  Py_DECREF(xi);

  // Return answer tuple.
  return answertuple;
}

/* Wrapper function for circleangle. */
static PyObject *circleangle_wrapper(PyObject *self, PyObject *args, PyObject *kwds) {
  /* Input arguments. */
  PyObject *r_obj;
  double p, z;

  // Keywords.
  static char *kwlist[] = {"r", "p", "z", NULL};

  /* Parse the input tuple */
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "Odd", kwlist, &r_obj, &p, &z)) {
    PyErr_SetString(PyExc_ValueError, "Error parsing arguments.");
    return NULL;
  }

  /* Check argument dimensions and types. */
  if (PyArray_NDIM(r_obj) != 1 || PyArray_TYPE(r_obj) != PyArray_DOUBLE) {
    PyErr_SetString(PyExc_ValueError, "Argument dimensions or types not correct.");
    return NULL; 
  }

  /* Get dimensions. */
  int n = PyArray_DIM(r_obj, 0);

  /* Interpret the input objects as numpy arrays. */
  PyObject *r_array = PyArray_FROM_OTF(r_obj, NPY_DOUBLE, NPY_IN_ARRAY);

  /* If that didn't work, throw an exception. */
  if (r_array == NULL) {
    Py_XDECREF(r_array);
    return NULL;
  }

  /* Get data pointers. */
  double *r = (double*)PyArray_DATA(r_array);
  
  // Create answer numpy arrays, let Python allocate memory.
  PyArrayObject *answer = (PyArrayObject *)PyArray_FromDims(1, &n, NPY_DOUBLE);

  // Calculate answer.
  circleangle(r, p, z, n, (double *)answer->data);

  /* Clean up. */
  Py_DECREF(r_array);

  // Return answer.
  return PyArray_Return(answer);
}

/* Wrapper function for ellipseangle. */
static PyObject *ellipseangle_wrapper(PyObject *self, PyObject *args, PyObject *kwds) {
  /* Input arguments. */
  PyObject *r_obj;
  double a, z;

  // Keywords.
  static char *kwlist[] = {"r", "a", "z", NULL};

  /* Parse the input tuple */
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "Odd", kwlist, &r_obj, &a, &z)) {
    PyErr_SetString(PyExc_ValueError, "Error parsing arguments.");
    return NULL;
  }

  /* Check argument dimensions and types. */
  if (PyArray_NDIM(r_obj) != 1 || PyArray_TYPE(r_obj) != PyArray_DOUBLE) {
    PyErr_SetString(PyExc_ValueError, "Argument dimensions or types not correct.");
    return NULL; 
  }

  /* Get dimensions. */
  int n = PyArray_DIM(r_obj, 0);

  /* Interpret the input objects as numpy arrays. */
  PyObject *r_array = PyArray_FROM_OTF(r_obj, NPY_DOUBLE, NPY_IN_ARRAY);

  /* If that didn't work, throw an exception. */
  if (r_array == NULL) {
    Py_XDECREF(r_array);
    return NULL;
  }

  /* Get data pointers. */
  double *r = (double*)PyArray_DATA(r_array);
  
  // Create answer numpy arrays, let Python allocate memory.
  PyArrayObject *answer = (PyArrayObject *)PyArray_FromDims(1, &n, NPY_DOUBLE);

  // Calculate answer.
  ellipseangle(r, a, z, n, (double *)answer->data);

  /* Clean up. */
  Py_DECREF(r_array);

  // Return answer.
  return PyArray_Return(answer);
}
