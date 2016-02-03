# $Header$


RM := rm -rf

include config.mk


HDFDIR = $(hdf4root)
HDFINC = $(HDFDIR)/include
HDFLIB = $(HDFDIR)/lib
ICAINC = hdf4/include
ICALIB = hdf4/lib
ICASRC = hdf4/src
CFLAGS = -O2 -w -I$(ICAINC) -I$(HDFINC)
LFLAGS = -L$(ICALIB) -L$(HDFLIB) -licaf90hdf -lmfhdf -ldf -ljpeg -lz -lm -w
COMSRC = hdf4/src/prntutil.f90
COMOBJ = prntutil.o
NCDF_LIB = /usr/local/lib


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
	$(fc) -c  $(CFLAGS) ${fflags} ${hdf5libs} ${hdflibs} $(COMSRC)  $<
	
	
OBJS = dcomp_one_pixel_run_program.o \
dcomp_retrieval_mod.o \
dcomp_array_loop_sub.o \
dcomp_forward_mod.o \
dcomp_math_tools_mod.o \
dcomp_science_tools_mod.o \
M_kracken.o \
dcomp_interface_types_mod.o \
dcomp_lut_mod.o \
view2d.o \
file_tools.o \
dcomp_lut_hdf_mod.o

OBJS_NLCOMP = dcomp_math_tools_mod.o \
dcomp_science_tools_mod.o \
M_kracken.o \
dcomp_interface_types_mod.o \
dcomp_lut_mod.o \
view2d.o \
file_tools.o \
dcomp_lut_hdf_mod.o \
nlcomp_array_loop_sub.o \
nlcomp_forward_mod.o \
nlcomp_interface_def_mod.o \
nlcomp_one_pixel_run_program.o \
nlcomp_retrieval_mod.o


OBJS_TEST_LUT = dcomp_lut_mod.o \
test_dcomp_lut.o \
dcomp_math_tools_mod.o \
file_tools.o

.PHONY: all_dncomp clean hdf4

all_dncomp: $(ICALIB)/libicaf90hdf.a dcomp_one_pixel_run nlcomp_one_pixel_run libdncomp.a
	
$(ICALIB)/libicaf90hdf.a:
	make COMP1="$(fc)" -C $(ICASRC)	
	     	
test_lut_mod: $(OBJS_TEST_LUT)
	@echo 'Building target for test lut'
	$(fc) ${fflags}  $(FFLAGS)    -o "test_lut" $(OBJS_TEST_LUT)  ${hdf5libs} ${hdflibs}  $(COMOBJ) $(LFLAGS)
	
dcomp_one_pixel_run: $(OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: Intel(R) Fortran Linker'
	
	$(fc) ${fflags}  $(FFLAGS)    -o "dcomp_one_pixel_run" $(OBJS)  ${hdf5libs} ${hdflibs} $(COMOBJ) $(LFLAGS) $(PATH_FLAGS)
   
	@echo 'Finished building target: $@'
	@echo ' '	

nlcomp_one_pixel_run: $(OBJS_NLCOMP)
	@echo 'Building target: $@'
	@echo 'Invoking: Intel(R) Fortran Linker'
	$(fc) ${fflags}  $(FFLAGS)    -o "nlcomp_one_pixel_run" $(OBJS_NLCOMP)  $(LFLAGS) ${hdf5libs} ${hdflibs} $(COMOBJ) $(PATH_FLAGS) 
	@echo 'Finished building target: $@'
	@echo ' '


libdncomp.a: $(OBJS)
	
	ar rvc libdncomp.a $(OBJS)
	
dcomp_one_pixel_run_program.o: dcomp_one_pixel_run_program.f90  dcomp_retrieval_mod.o dcomp_science_tools_mod.o M_kracken.o
dcomp_retrieval_mod.o:dcomp_retrieval_mod.f90 dcomp_forward_mod.o dcomp_math_tools_mod.o dcomp_science_tools_mod.o
dcomp_array_loop_sub.o:dcomp_array_loop_sub.f90 dcomp_retrieval_mod.o dcomp_interface_types_mod.o view2d.o
dcomp_forward_mod.o: dcomp_forward_mod.f90  dcomp_science_tools_mod.o   dcomp_lut_mod.o
dcomp_math_tools_mod.o: dcomp_math_tools_mod.f90
dcomp_lut_mod.o:dcomp_lut_mod.f90 dcomp_math_tools_mod.o file_tools.o dcomp_lut_hdf_mod.o
dcomp_science_tools_mod.o:dcomp_science_tools_mod.f90
dcomp_interface_types_mod.o:dcomp_interface_types_mod.f90
dcomp_lut_hdf_mod.o:dcomp_lut_hdf_mod.f90
nlcomp_one_pixel_run_program.o: nlcomp_one_pixel_run_program.f90  nlcomp_retrieval_mod.o nlcomp_interface_def_mod.o
nlcomp_retrieval_mod.o:nlcomp_retrieval_mod.f90 nlcomp_forward_mod.o dcomp_math_tools_mod.o dcomp_science_tools_mod.o 
nlcomp_array_loop_sub.o:nlcomp_array_loop_sub.f90 nlcomp_retrieval_mod.o nlcomp_interface_def_mod.o
nlcomp_forward_mod.o: nlcomp_forward_mod.f90  dcomp_math_tools_mod.o  M_kracken.o dcomp_lut_mod.o dcomp_science_tools_mod.o
nlcomp_interface_def_mod.o:nlcomp_interface_def_mod.f90

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
	rm -f *.a
	rm -f test_lut
	rm -rf hdf4/lib
	make clean  -C $(ICASRC)
	-@echo ' '
