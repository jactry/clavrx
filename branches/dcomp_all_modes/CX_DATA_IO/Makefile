#
#

include config.mk

#fc=ifort -g 
fflags=-I./
DEPS=cx_sds_constants.F90
OBJ=cx_sds_test.o cx_sds_io_mod.o cx_hdf_read_mod.o cx_sds_type_definitions_mod.o cx_hdf_write_mod.o


HDFDIR = $(hdf4root)
HDF4_LIB = -L$(HDFDIR)/lib -I$(HDFDIR)/include  -lmfhdf -ldf -ljpeg -lz


.SUFFIXES: .F90 .F .o

.F90.o:
	$(fc) -c  $<  ${fflags} 

%.o: %.mod

.PHONY: all clean

all:cx_data_io libcx_data_io.a

cx_data_io: $(OBJ)
	@echo 'Building targets'
	$(fc) -o "cx_data_io" $(OBJ)  ${fflags}  ${HDF4_LIB}

libcx_data_io.a: ${OBJ}
	ar rvc libcx_data_io.a ${OBJ}
   
cx_sds_test.o:cx_sds_test.F90 cx_sds_io_mod.o cx_sds_type_definitions_mod.o
cx_sds_io_mod.o:cx_sds_io_mod.F90 cx_hdf_read_mod.o cx_sds_type_definitions_mod.o
cx_sds_type_definitions_mod.o:cx_sds_type_definitions_mod.F90
	$(fc) -c ${fflags}  cx_sds_type_definitions_mod.F90 
cx_hdf_read_mod.o:cx_hdf_read_mod.F90 cx_sds_type_definitions_mod.o
	$(fc) -c ${fflags} ${HDF4_LIB} cx_hdf_read_mod.F90 cx_sds_type_definitions_mod.o
cx_hdf_write_mod.o:cx_hdf_write_mod.f90 
	$(fc) -c ${fflags} ${HDF4_LIB} cx_hdf_write_mod.f90 cx_sds_type_definitions_mod.o  
clean:
	rm -f *.o *.mod
