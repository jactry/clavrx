# This file was generated on Thu Jun  5 12:17:00 2014 by awalther with the following command-line arguments: -hdf5root=/Users/awalther/lib/hdf5/ -with-ifort -hdflib=/Users/awalther/lib/hdf4//lib.
# System info: Darwin luna.ssec.wisc.edu 13.2.0 Darwin Kernel Version 13.2.0: Thu Apr 17 23:03:13 PDT 2014; root:xnu-2422.100.13~1/RELEASE_X86_64 x86_64
fc = ifort
fflags = -O2 -assume byterecl
hdflibs = -L/Users/awalther/lib/hdf4//lib -lmfhdf -ldf -ljpeg -lz
hdf5libs = -I/Users/awalther/lib/hdf5/include/ -L/Users/awalther/lib/hdf5/lib/
hdf5links = -lhdf5_fortran -lhdf5 -lz
