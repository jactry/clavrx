HDFDIR = $(HOME)/lib/hdf4
#COMP   = g95
#COMP   = lf95
COMP   = ifort


   ifdef HDF4R
      HDFDIR=$(HDF4R)
   endif

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

.SUFFIXES:
.SUFFIXES: .f90 .F .o
   
#Compilation rule for f90 files
.f90.o:
	$(COMP) -c  $(CFLAGS)  $<
	

all: $(LIBPTH)

$(LIBPTH): $(LIBOBJ) 
	mkdir -p $(ICALIB)
	ar -r $(LIBPTH) $(LIBOBJ) 
	cp $(LIBMOD) $(ICAINC)  

clean :
	rm -rf $(ICALIB) *.o *.mod
