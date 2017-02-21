
# $Header$

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
	
	
OBJS =  \
dcomp_retrieval_mod.o \
dcomp_array_loop_sub.o \
dcomp_forward_mod.o \
dcomp_math_tools_mod.o \
dcomp_science_tools_mod.o \
M_kracken.o \
sensorname_from_wmoid.o \
dncomp_interface_def_mod.o \
dcomp_lut_mod.o \
view2d.o \
file_tools.o \
dcomp_lut_hdf_mod.o \
dncomp_trans_atmos_mod.o

OBJS_ARRAY_TEST := $(OBJS) dcomp_array_test.o
OBJS_ONE_TEST :=  $(OBJS) dcomp_one_pixel_run_program.o

OBJS_NLCOMP = dcomp_math_tools_mod.o \
dcomp_science_tools_mod.o \
M_kracken.o \
dcomp_lut_mod.o \
view2d.o \
file_tools.o \
dcomp_lut_hdf_mod.o \
nlcomp_array_loop_sub.o \
nlcomp_forward_mod.o \
dncomp_interface_def_mod.o \
nlcomp_one_pixel_run_program.o \
nlcomp_retrieval_mod.o


OBJS_TEST_LUT = dcomp_lut_mod.o \
test_dcomp_lut.o \
dcomp_math_tools_mod.o \
cx_hdf_read_mod.o \
file_tools.o

.PHONY: all_dncomp clean hdf4

all_dncomp:  cx_read_io dcomp_one_pixel_run nlcomp_one_pixel_run  libdncomp.a
test_dncomp:  cx_read_io  dcomp_array_test	
	     	
test_lut_mod: $(OBJS_TEST_LUT)
	@echo 'Building target for test lut'
	$(fc) ${fflags}  $(FFLAGS)    -o "test_lut" $(OBJS_TEST_LUT)  ${hdf5libs}  ${cx_data_libs} ${cx_data_links}  $(LFLAGS)
	
dcomp_one_pixel_run: $(OBJS_ONE_TEST)
	@echo 'Building target: $@'
	
	
	$(fc) ${fflags}  $(FFLAGS)    -o "dcomp_one_pixel_run" $(OBJS_ONE_TEST)  ${hdf5libs} ${cx_data_libs} ${cx_data_links} $(LFLAGS) $(PATH_FLAGS)
   
	@echo 'Finished building target: $@'
	@echo ' '	
   
dcomp_array_test:  $(OBJS_ARRAY_TEST)
	@echo 'Building target: $@'
	@echo $(OBJS_ARRAY_TEST)
	
	$(fc) ${fflags}  $(FFLAGS)    -o "dcomp_array_test" $(OBJS_ARRAY_TEST)  ${hdf5libs} ${cx_data_libs} ${cx_data_links} $(LFLAGS) $(PATH_FLAGS)
   
	@echo 'Finished building target: $@'
	@echo ' '	   

nlcomp_one_pixel_run: $(OBJS_NLCOMP)
	@echo 'Building target: $@'
	
	$(fc) ${fflags}  $(FFLAGS)    -o "nlcomp_one_pixel_run" $(OBJS_NLCOMP)   ${hdf5libs}   ${cx_data_libs} ${cx_data_links} $(LFLAGS)  $(PATH_FLAGS) 
	@echo 'Finished building target: $@'
	@echo ' '


libdncomp.a: $(OBJS)
	
	ar rvc libdncomp.a $(OBJS) $(OBJS_NLCOMP)


cx_read_io:
	$(MAKE)  -C $(CX_READ_IO_DIR) ARCH=$(ARCH)

dcomp_array_test.o: dcomp_array_test.f90 dcomp_science_tools_mod.o dcomp_array_loop_sub.o
dcomp_one_pixel_run_program.o: dcomp_one_pixel_run_program.f90  dcomp_retrieval_mod.o dcomp_science_tools_mod.o M_kracken.o
dcomp_retrieval_mod.o:dcomp_retrieval_mod.f90 dcomp_forward_mod.o dcomp_math_tools_mod.o dcomp_science_tools_mod.o
dcomp_array_loop_sub.o:dcomp_array_loop_sub.f90 dcomp_retrieval_mod.o dncomp_interface_def_mod.o view2d.o sensorname_from_wmoid.o dncomp_trans_atmos_mod.o
dcomp_forward_mod.o: dcomp_forward_mod.f90  dcomp_science_tools_mod.o   dcomp_lut_mod.o
dcomp_math_tools_mod.o: dcomp_math_tools_mod.f90
dcomp_lut_mod.o:dcomp_lut_mod.f90 dcomp_math_tools_mod.o file_tools.o dcomp_lut_hdf_mod.o
dcomp_science_tools_mod.o:dcomp_science_tools_mod.f90
dncomp_interface_def_mod.o:dncomp_interface_def_mod.f90
dncomp_trans_atmos_mod.o:dncomp_trans_atmos_mod.f90
dcomp_lut_hdf_mod.o:dcomp_lut_hdf_mod.f90
sensorname_from_wmoid.o:sensorname_from_wmoid.f90
nlcomp_one_pixel_run_program.o: nlcomp_one_pixel_run_program.f90  nlcomp_retrieval_mod.o dncomp_interface_def_mod.o
nlcomp_retrieval_mod.o:nlcomp_retrieval_mod.f90 nlcomp_forward_mod.o dcomp_math_tools_mod.o dcomp_science_tools_mod.o 
nlcomp_array_loop_sub.o:nlcomp_array_loop_sub.f90 nlcomp_retrieval_mod.o dncomp_interface_def_mod.o
nlcomp_forward_mod.o: nlcomp_forward_mod.f90  dcomp_math_tools_mod.o  M_kracken.o dcomp_lut_mod.o dcomp_science_tools_mod.o

M_kracken.o:$(SUB_PATH)/M_kracken.f90
	$(fc) -c ${fflags}  $(SUB_PATH)/M_kracken.f90

test_dcomp_lut.o:$(TEST_PATH)/test_dcomp_lut.f90 dcomp_lut_mod.o
	$(fc) -c ${fflags} ${hdflibs} $(COMOBJ) $(LFLAGS) $(TEST_PATH)/test_dcomp_lut.f90

view2d.o:$(SUB_PATH)/view2d.f90
	$(fc) -c ${fflags}  $(SUB_PATH)/view2d.f90
   
   
file_tools.o:file_tools.f90



clean:
	rm -f *.o *.mod
	rm -f dcomp_one_pixel_run
	rm -f nlcomp_one_pixel_run
	rm -f *.a
	rm -f test_lut
	rm -rf hdf4/lib
	$(MAKE) clean -C $(CX_READ_IO_DIR) 	
	-@echo ' '
