# This file was generated on Fri Jan 31 00:04:09 2014 by awalther with the following command-line arguments: -hdf5root=/Users/awalther/lib/hdf5/ -with-ifort -hdflib=/Users/awalther/lib/hdf4//lib.
# System info: Darwin luna.local 13.0.0 Darwin Kernel Version 13.0.0: Thu Sep 19 22:22:27 PDT 2013; root:xnu-2422.1.72~6/RELEASE_X86_64 x86_64
fc = ifort
fflags = -O2 -assume byterecl
hdflibs = -L/Users/awalther/lib/hdf4//lib -lmfhdf -ldf -ljpeg -lz
hdf5libs = -I/Users/awalther/lib/hdf5/include/ -L/Users/awalther/lib/hdf5/lib/
hdf5links = -lhdf5_fortran -lhdf5 -lz
