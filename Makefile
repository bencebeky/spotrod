all: spotrod.so

.PHONY: clean

clean:
	rm -f spotrod.so

spotrod.so: spotrod.c spotrod.h spotrod-python-wrapper.c setup.py
	python setup.py build_ext --inplace
