# This file was generated on Thu May 23 16:23:29 2013 by awalther with the following command-line arguments: -with-ifort.
# System info: Darwin luna.ssec.wisc.edu 12.3.0 Darwin Kernel Version 12.3.0: Sun Jan  6 22:37:10 PST 2013; root:xnu-2050.22.13~1/RELEASE_X86_64 x86_64

fc = ifort -i-dynamic -mcmodel medium -shared-intel
fflags =  -C -O0 -g -m64 -warn unused 


hdflibs = -L/opt/hdf4/4.2.7-ifort/lib/ -lmfhdf -ldf -ljpeg -lz
hdf5libs = -I/opt/hdf5/1.8.9-ifort/include/ -L/opt/hdf5/1.8.9-ifort/lib/
hdf5links = -lhdf5_fortran -lhdf5hl_fortran
