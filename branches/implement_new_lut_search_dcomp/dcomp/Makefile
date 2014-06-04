# $Header:$


RM := rm -rf

include config.mk

# Where to find the RTM. Replace this with your own path to the object files and modules.
RTM_INC = /home/awalther/SOFTWARE/RTM/RTTOV/include
RTM_MOD = /home/awalther/SOFTWARE/RTM/RTTOV/mod 
RTM_LIB = /home/awalther/SOFTWARE/RTM/RTTOV/lib
RTMODEL = rttov10.2.0

HDFDIR = $(HOME)/lib/hdf4/
HDFINC = $(HDFDIR)/include
HDFLIB = $(HDFDIR)/lib
ICAINC = hdf4/include
ICALIB = hdf4/lib
ICASRC = hdf4/src
CFLAGS = -w -I$(ICAINC) -I$(HDFINC)
LFLAGS = -L$(ICALIB) -L$(HDFLIB) -licaf90hdf -lmfhdf -ldf -ljpeg -lz -lm
COMSRC = prntutil.f90
COMOBJ = prntutil.o
NCDF_LIB = /usr/local/lib

.SUFFIXES:
.SUFFIXES: .f90 .F .o

	


#Compilation rule for f90 files
.f90.o:
	$(fc) -c  $(CFLAGS) ${fflags} ${hdf5libs} ${hdflibs} $(COMSRC)  $<
	
	
OBJS = dcomp.o \
dcomp_forward_mod.o \
dcomp_tools.o \
M_kracken.o \
dcomp_mod.o \
get_planck_radiance_39um.o \
get_rad_refl_factor.o \
dcomp_array.o \
dcomp_interface_types_mod.o \
dcomp_data_pool_mod.o \
view2d.o \
file_tools.o 




OBJS_TEST_LUT = dcomp_data_pool_mod.o \
test_dcomp_lut.o \
dcomp_tools.o \
file_tools.o

.PHONY: clean hdf4

all: hdf4 dcomp 

hdf4:
	@echo 'Building target for hdf4: $@'
	make -C $(ICASRC)	
	     	
test_lut_mod: $(OBJS_TEST_LUT)
	@echo 'test lut'
	$(fc) ${fflags}  $(FFLAGS)    -o "test_lut" $(OBJS_TEST_LUT)  ${hdf5libs} ${hdflibs}  $(COMOBJ) $(LFLAGS)
	


dcomp: $(OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: Intel(R) Fortran Linker'
	
	$(fc) ${fflags}  $(FFLAGS)    -o "dcomp" $(OBJS)  ${hdf5libs} ${hdflibs} $(COMOBJ) $(LFLAGS)
   
	@echo 'Finished building target: $@'
	@echo ' '	

dcomp.o: dcomp.f90  dcomp_mod.o get_planck_radiance_39um.o
dcomp_forward_mod.o: dcomp_forward_mod.f90  dcomp_tools.o  M_kracken.o get_rad_refl_factor.o dcomp_data_pool_mod.o
dcomp_tools.o: dcomp_tools.f90
test_forward_computation.o:test_forward_computation.f90 dcomp_forward_mod.o
M_kracken.o:M_kracken.f90
dcomp_mod.o:dcomp_mod.f90 dcomp_forward_mod.o dcomp_tools.o get_rad_refl_factor.o 
get_planck_radiance_39um.o:get_planck_radiance_39um.f90
get_rad_refl_factor.o:get_rad_refl_factor.f90
dcomp_array.o:dcomp_array.f90 dcomp_mod.o dcomp_interface_types_mod.o view2d.o 
dcomp_interface_types_mod.o:dcomp_interface_types_mod.f90
dcomp_data_pool_mod.o:dcomp_data_pool_mod.f90 dcomp_tools.o file_tools.o
test_dcomp_lut.o:test_dcomp_lut.f90 dcomp_data_pool_mod.o
view2d.o:view2d.f90
test_array.o:test_array.f90 dcomp_interface_types_mod.o dcomp_array.o
file_tools.o:file_tools.f90



clean:
	rm -f *.o *.mod
	rm -f dcomp
	rm -f *.a
	rm -f test_lut
	rm -rf hdf4/lib
	-@echo ' '
