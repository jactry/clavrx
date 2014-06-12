# This file was generated on Thu Jun 12 09:55:22 2014 by awalther with the following command-line arguments: -hdf5root=/home/awalther/lib/hdf5/ -with-ifort -hdflib=/home/awalther/lib/hdf4//lib.
# System info: Linux odin.ssec.wisc.edu 2.6.18-53.1.4.el5 #1 SMP Wed Nov 14 10:37:27 EST 2007 x86_64 x86_64 x86_64 GNU/Linux
fc = gfortran
fflags = -O2 
hdflibs = -L/home/awalther/lib/hdf4//lib -lmfhdf -ldf -ljpeg -lz
hdf5libs = -I/home/awalther/lib/hdf5/include/ -L/home/awalther/lib/hdf5/lib/
hdf5links = -lhdf5_fortran -lhdf5 -lz
