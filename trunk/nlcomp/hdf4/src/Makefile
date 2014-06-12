HDFDIR = $(HOME)/lib/hdf4
#COMP   = g95
#COMP   = lf95
COMP   = gfortran

HDFINC = $(HDFDIR)/include
ICAINC = ../include
ICALIB = ../lib
LIBPTH = $(ICALIB)/libicaf90hdf.a

CFLAGS = -w -I$(ICAINC) -I$(HDFINC)

LIBSRC = ica_f90_hdf_sds.f90 ica_f90_hdf_dim_sds.f90
LIBOBJ = ica_f90_hdf_sds.o   ica_f90_hdf_dim_sds.o
LIBMOD = ica_f90_hdf_sds.mod ica_f90_hdf_dim_sds.mod

all: lib

lib: $(LIBSRC)
	mkdir -p $(ICALIB)
	$(COMP) -c $(CFLAGS) $(LIBSRC) 
	ar -r $(LIBPTH) $(LIBOBJ)
	cp $(LIBMOD) $(ICAINC)

clean :
	rm -rf $(ICALIB) *.o *.mod
