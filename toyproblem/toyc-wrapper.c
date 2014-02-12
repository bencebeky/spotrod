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
#include "toyc.h"

/* Docstrings */
static char module_docstring[] = 
"  This module is a fast C implementation of a toy problem.";
static char circleangleloop_docstring[] = 
"  circleangleloop(r, p, z)\n"
"  Calculate half central angle of the arc of circle of radius r\n"
"  that is inside a circle of radius p with separation of centers z.";
static char circleanglesorted_docstring[] = 
"  circleanglesorted(r, p, z)\n"
"  Calculate half central angle of the arc of circle of radius r\n"
"  that is inside a circle of radius p with separation of centers z.\n"
"  This version assumes r in increasing.";

/* Function wrappers for external use */
static PyObject *circleangleloop_wrapper(PyObject*, PyObject*, PyObject*);
static PyObject *circleanglesorted_wrapper(PyObject*, PyObject*, PyObject*);

/* Module specification */
static PyMethodDef module_methods[] = {
  {"circleangleloop", (PyCFunction)circleangleloop_wrapper, METH_VARARGS | METH_KEYWORDS, circleangleloop_docstring},
  {"circleanglesorted", (PyCFunction)circleanglesorted_wrapper, METH_VARARGS | METH_KEYWORDS, circleanglesorted_docstring},
  {NULL, NULL, 0, NULL}
};

/* Initialize the module */
PyMODINIT_FUNC inittoyc(void) {
  PyObject *m = Py_InitModule3("toyc", module_methods, module_docstring);
  if (m == NULL)
    return;
  /* Load numpy functionality. */
  import_array();
}

/* Wrapper function for circleangleloop. */
static PyObject *circleangleloop_wrapper(PyObject *self, PyObject *args,
                                         PyObject *kwds) {
  /* Input arguments. */
  double p, z;
  PyObject *r_obj;

  // Keywords.
  static char *kwlist[] = {"r", "p", "z", NULL};

  /* Parse the input tuple */
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "Odd", kwlist, 
                                   &r_obj, &p, &z))
    return NULL;

  /* Interpret the input object as a numpy array. */
  PyObject *r_array = PyArray_FROM_OTF(r_obj, NPY_DOUBLE, NPY_IN_ARRAY);

  /* If that didn't work, or the resulting array does not have the correct
   * number of dimensions or type, then abort. */
  if (r_array == NULL || PyArray_NDIM(r_array) != 1 ||
                           PyArray_TYPE(r_array) != PyArray_DOUBLE) {
    PyErr_SetString(PyExc_ValueError, 
                    "r cannot be converted to a suitable array.");
    return NULL; 
  }

  /* Read out dimensions and data pointers. */
  int n = (int)PyArray_DIM(r_array, 0);
  double *r_data = (double*)PyArray_DATA(r_array);

  /* Create answer numpy array, let Python allocate memory.
     Do not allocate memory manually and then use PyArray_FromDimsAndData! */
  PyArrayObject *answer_array = (PyArrayObject*)PyArray_FromDims(1, &n, NPY_DOUBLE);

  // Evaluate the model
  circleangleloop(r_data, p, z, n, (double*)PyArray_DATA(answer_array));

  /* Clean up. */
  Py_DECREF(r_array);

  // Return.
  return PyArray_Return(answer_array);
}

/* Wrapper function for circleanglesorted. */
static PyObject *circleanglesorted_wrapper(PyObject *self, PyObject *args, PyObject *kwds) {
  /* Input arguments. */
  double p, z;
  PyObject *r_obj;

  // Keywords.
  static char *kwlist[] = {"r", "p", "z", NULL};

  /* Parse the input tuple */
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "Odd", kwlist, &r_obj, &p, &z))
    return NULL;

  /* Check argument dimensions and types. */
  if (PyArray_NDIM(r_obj) != 1 || PyArray_TYPE(r_obj) != PyArray_DOUBLE) {
    PyErr_SetString(PyExc_ValueError, "Argument dimensions or types not correct.");
    return NULL; 
  }

  /* Interpret the input objects as numpy arrays. */
  PyObject *r_array = PyArray_FROM_OTF(r_obj, NPY_DOUBLE, NPY_IN_ARRAY);

  /* If that didn't work, throw an exception. */
  if (r_array == NULL) {
    Py_XDECREF(r_array);
    return NULL;
  }

  /* Read out dimensions and data pointers. */
  int n = (int)PyArray_DIM(r_array, 0);
  double *r_data = (double*)PyArray_DATA(r_array);

  /* Create answer numpy array, let Python allocate memory.
     Do not allocate memory manually and then use PyArray_FromDimsAndData! */
  PyArrayObject *answer = (PyArrayObject *)PyArray_FromDims(1, &n, NPY_DOUBLE);

  // Evaluate the model
  circleanglesorted(r_data, p, z, n, (double *)answer->data);

  /* Clean up. */
  Py_DECREF(r_array);

  // Return.
  return PyArray_Return(answer);
}

