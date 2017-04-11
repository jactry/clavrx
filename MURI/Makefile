# $Header:$

OBJDIR=./obj
RM := rm -rf

ARCH=gfortran
include arch/$(ARCH)

HDFDIR = $(hdf4root)
HDFINC = $(SSEC_HDF4_INC)
HDFLIB = $(SSEC_HDF4_LIB)
CX_READ_IO_DIR = ../CX_DATA_IO
CFLAGS = -O2 -w  -I$(HDFINC)
LFLAGS =  -L$(HDFLIB)  -lmfhdf -ldf -ljpeg -lz -lm -w

TEST_PATH = ./tests/
SUB_PATH = ./subtools/
PATH_FLAGS = -I$(TEST_PATH) -I$(SUB_PATH)

.SUFFIXES:
.SUFFIXES: .f90 .F .o

	ifdef COMP1
		$(fc)=COMP1 
	endif   

#Compilation rule for f90 files
.f90.o:
	$(fc) -c  $(CFLAGS) ${fflags} ${hdf5libs} ${cx_data_libs} ${cx_data_links}   $<
	
	
OBJS = muri_one_pixel_run.o muri_retrieval_mod.o muri_forward_mod.o muri_definitions_mod.o muri_array_loop_sub.o muri_interface_mod.o muri_lut_mod.o lib_array.o

.PHONY: all_muri clean

muri_clavrx: cx_read_io  muri_one_pixel_run libmuri.a
all_muri: cx_read_io  muri_one_pixel_run 

	     	
test_lut_mod: $(OBJS_TEST_LUT)
	@echo 'Building target for test lut'
	$(fc) ${fflags}  $(FFLAGS)    -o "test_lut" $(OBJS_TEST_LUT)  ${hdf5libs}  ${cx_data_libs} ${cx_data_links}  $(LFLAGS)
	
muri_one_pixel_run: $(OBJS)
	@echo 'Building target: $@'
	
	
	$(fc) ${fflags}  $(FFLAGS)    -o "muri_one_pixel_run" $(OBJS)  ${hdf5libs} ${cx_data_libs} ${cx_data_links} $(LFLAGS) $(PATH_FLAGS)
   
	@echo 'Finished building target: $@'
	@echo ' '	

libmuri.a: $(OBJS)
	ar rvc libmuri.a $(OBJS) 
   
cx_read_io:
	$(MAKE)  -C $(CX_READ_IO_DIR) ARCH=$(ARCH)

muri_array_loop_sub.o:muri_array_loop_sub.f90 muri_interface_mod.o muri_retrieval_mod.o muri_definitions_mod.o
muri_interface_mod.o:muri_interface_mod.f90
muri_one_pixel_run.o:muri_one_pixel_run.f90 muri_retrieval_mod.o muri_definitions_mod.o
muri_retrieval_mod.o:muri_retrieval_mod.f90 muri_definitions_mod.o muri_lut_mod.o lib_array.o
muri_definitions_mod.o:muri_definitions_mod.f90
muri_forward_mod.o:muri_forward_mod.f90 muri_lut_mod.o
muri_lut_mod.o:muri_lut_mod.f90
lib_array.o:lib_array.f90


clean:
	rm -f *.o *.mod
	rm -f muri_one_pixel_run
	rm -f nlcomp_one_pixel_run
	rm -f *.a
	rm -f test_lut
	rm -rf hdf4/lib
	
	-@echo ' '
