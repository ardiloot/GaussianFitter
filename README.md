# README #

This library is developed for speedup of the SMolPhot software. It allows to speed-up the fitting of 2D Gaussian function significanly by eliminating the overhead to the Python fitting function. The code is written in Fortran because the fitting function lmdif from minpack is easily available in Fortran (the same function is also used by SciPy).

### Requirements ###

* fortran and c-compiler (e.g. [mingw-w64](http://sourceforge.net/projects/mingw-w64/files/))

### Setup ###

* (python setup.py install) or (pip install .)
