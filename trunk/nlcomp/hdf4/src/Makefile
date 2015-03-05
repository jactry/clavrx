HDFDIR = $(HOME)/lib/hdf4
#COMP   = g95
#COMP   = lf95
COMP   = ifort

HDFINC = $(HDFDIR)/include
ICAINC = ../include
ICALIB = ../lib
LIBPTH = $(ICALIB)/libicaf90hdf.a

CFLAGS = -O2 -I$(ICAINC) -I$(HDFINC)

LIBSRC = ica_f90_hdf_sds.f90 ica_f90_hdf_dim_sds.f90
LIBOBJ = ica_f90_hdf_sds.o   ica_f90_hdf_dim_sds.o
LIBMOD = ica_f90_hdf_sds.mod ica_f90_hdf_dim_sds.mod

	ifdef COMP1
		COMP=$(COMP1) 
	endif  

all: lib

lib: $(LIBSRC)
	mkdir -p $(ICALIB)
	$(COMP) -c $(CFLAGS) $(LIBSRC) 
	ar -r $(LIBPTH) $(LIBOBJ)
	cp $(LIBMOD) $(ICAINC)

clean :
	rm -rf $(ICALIB) *.o *.mod
