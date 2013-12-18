#include <Python.h>
#include <numpy/arrayobject.h>
#include "toyc.h"

/* Docstrings */
static char module_docstring[] = "  This module is a fast C implementation of a toy problem.";
static char circleangle_docstring[] = 
"  circleangle(r, p, z)\n"
"  Calculate half central angle of the arc of circle of radius r\n"
"  that is inside a circle of radius p with separation of centers z.";

/* Function wrappers for external use */
static PyObject *circleangle_wrapper(PyObject*, PyObject*, PyObject*);

/* Module specification */
static PyMethodDef module_methods[] = {
  {"circleangle", (PyCFunction)circleangle_wrapper, METH_VARARGS | METH_KEYWORDS, circleangle_docstring},
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

/* Wrapper function for circleangle. */
static PyObject *circleangle_wrapper(PyObject *self, PyObject *args, PyObject *kwds) {
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
  circleangle(r_data, p, z, n, (double *)answer->data);

  /* Clean up. */
  Py_DECREF(r_array);

  // Return.
  return PyArray_Return(answer);
}

