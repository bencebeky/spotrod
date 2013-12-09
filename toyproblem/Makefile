all: toyc.so toycython.so

toycython.so: toycython.pyx
	python toycython-setup.py build_ext --inplace

toyc.so: toyc.c toyc-wrapper.c toyc-setup.py
	python toyc-setup.py build_ext --inplace
