from setuptools import setup, Extension
import numpy as np

name = "spotrod"
version = "1.0"
description = "A semi-analytic model for transits of spotted stars"
author = "Bence Béky"
author_email = "zsebkecske@gmail.com"
maintainer = author
maintainer_email = author_email
url = "https://github.com/bencebeky/spotrod"

module = Extension(
	name,
	sources=['spotrod.c', 'spotrod-python-wrapper.c'],
	include_dirs=[np.get_include()],
	extra_compile_args=['-Ofast']
)

setup(
	name=name,
	version=version,
	description=description,
	author=author,
	author_email=author_email,
	maintainer=maintainer,
	maintainer_email=maintainer_email,
	url=url,
	ext_modules=[module]
)