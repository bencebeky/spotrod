all: spotrod.so

.PHONY: clean

clean:
	rm -f spotrod.so

spotrod.so: spotrod.c spotrod.h spotrod-python-wrapper.c spotrod-setup.py
	python spotrod-setup.py build_ext --inplace
