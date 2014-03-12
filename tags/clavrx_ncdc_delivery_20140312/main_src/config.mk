# This file was generated on Fri Jan 31 00:04:09 2014 by awalther with the following command-line arguments: -hdf5root=/Users/awalther/lib/hdf5/ -with-ifort -hdflib=/Users/awalther/lib/hdf4//lib -hdfinc=/Users/awalther/lib/hdf4//include -nlcomp_dir=../cloud_team_nlcomp/ -dcomp_dir=../cloud_team_dcomp/ -acha_dir=../cloud_team_acha/.
# System info: Darwin luna.local 13.0.0 Darwin Kernel Version 13.0.0: Thu Sep 19 22:22:27 PDT 2013; root:xnu-2422.1.72~6/RELEASE_X86_64 x86_64
fc = ifort
fflags = -O2 -assume byterecl
fflags_pfast = -O2 -assume byterecl -fixed
fflags_sasrab_f77 = -O2 -assume byterecl -fixed -save
fflags_sasrab_f90 = -O2 -assume byterecl -save
ldflags = -O2 -assume byterecl
cpp = -cpp
cppflags = 
beconv = -convert big_endian
hdflibs = -L/Users/awalther/lib/hdf4//lib -lmfhdf -ldf -ljpeg -lz
hdfincs = -I/Users/awalther/lib/hdf4//include
hdf5libs = -I/Users/awalther/lib/hdf5/include/ -L/Users/awalther/lib/hdf5/lib/
hdf5links = -lhdf5_fortran -lhdf5 -lz
export CMASK=../jpssrr_nbvcm/
export CTYPE=../akh_cloud_type_repo/
export ACHA=../cloud_team_acha/
export DCOMP=../cloud_team_dcomp/
dcomplibs=-L$(DCOMP) -I$(DCOMP)
dcomplinks= -ldcomp
export NLCOMP=../cloud_team_nlcomp/
nlcomplibs=-L$(NLCOMP) -I$(NLCOMP)
nlcomplinks= -lnlcomp
