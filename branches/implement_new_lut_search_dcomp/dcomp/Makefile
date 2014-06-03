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
prntutil.o \
dcomp_data_pool_mod.o \
view2d.o \
file_tools.o 

OBJS_FULL = test_full.o \
dcomp_forward_mod.o \
dcomp_tools.o \
M_kracken.o \
dcomp_mod.o \
get_planck_radiance_39um.o \
get_rad_refl_factor.o \
dcomp_array.o \
dcomp_interface_types_mod.o \
prntutil.o \
dcomp_data_pool_mod.o \
view2d.o  \
file_tools.o

OBJS_TEST_ARRAY = test_array.o \
dcomp_forward_mod.o \
dcomp_tools.o \
M_kracken.o \
dcomp_mod.o \
get_planck_radiance_39um.o \
get_rad_refl_factor.o \
dcomp_array.o \
dcomp_interface_types_mod.o \
prntutil.o \
dcomp_data_pool_mod.o \
view2d.o  \
file_tools.o


OBJS_TEST = test_forward_computation.o \
dcomp_forward_mod.o \
dcomp_tools.o \
hdf_dcomp_awg.o \
dcomp_mod.o \
get_planck_radiance_39um.o \
get_rad_refl_factor.o \
file_tools.o




OBJS_TEST_LUT = dcomp_data_pool_mod.o \
test_dcomp_lut.o \
dcomp_tools.o \
file_tools.o

all: dcomp

test: $(OBJS_TEST)
	@echo 'test'
	$(fc) ${fflags}  $(FFLAGS)    -o "test_fwd" $(OBJS_TEST)  ${hdf5libs} ${hdflibs} 

test_full: $(OBJS_FULL)
	@echo 'test'
	$(fc) ${fflags}  $(FFLAGS)    -o "test_full" $(OBJS_FULL)  ${hdf5libs} ${hdflibs}   $(COMOBJ) $(LFLAGS)
   
test_array: $(OBJS_TEST_ARRAY)
	@echo 'test'
	$(fc) ${fflags}  $(FFLAGS)    -o "test_array" $(OBJS_TEST_ARRAY)  ${hdf5libs} ${hdflibs}   $(COMOBJ) $(LFLAGS)   
	
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
test_full.o: test_full.f90 dcomp_mod.o get_planck_radiance_39um.o
dcomp_forward_mod.o: dcomp_forward_mod.f90  dcomp_tools.o  M_kracken.o get_rad_refl_factor.o dcomp_data_pool_mod.o
dcomp_tools.o: dcomp_tools.f90
hdf_dcomp_awg.o: hdf_dcomp_awg.f90
test_forward_computation.o:test_forward_computation.f90 dcomp_forward_mod.o
M_kracken.o:M_kracken.f90
dcomp_mod.o:dcomp_mod.f90 dcomp_forward_mod.o dcomp_tools.o get_rad_refl_factor.o 
get_planck_radiance_39um.o:get_planck_radiance_39um.f90
get_rad_refl_factor.o:get_rad_refl_factor.f90
test_new_fwd.o:test_new_fwd.f90  dcomp_fwd_mod.o
dcomp_array.o:dcomp_array.f90 dcomp_mod.o dcomp_interface_types_mod.o view2d.o 
dcomp_interface_types_mod.o:dcomp_interface_types_mod.f90
dcomp_data_pool_mod.o:dcomp_data_pool_mod.f90 dcomp_tools.o file_tools.o
test_dcomp_lut.o:test_dcomp_lut.f90 dcomp_data_pool_mod.o
view2d.o:view2d.f90
test_array.o:test_array.f90 dcomp_interface_types_mod.o dcomp_array.o
file_tools.o:file_tools.f90



clean:
	rm -f *.o *.mod
	rm -f ../cr_start
	-@echo ' '
