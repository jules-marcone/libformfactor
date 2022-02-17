# libformfactor

A C++ library for the efficient computation of scattering form factors
(Fourier shape transforms) of arbitrary polyhedra according to Wuttke,
[J Appl Cryst 54, 580-587 (2021)](https://doi.org/10.1107/S1600576721001710).

The library is in directory ff/. Tests are in directory test/.
To build the binaries and run the tests, do
```
mkdir build
cd build
cmake ..
make
ctest
make install
```
