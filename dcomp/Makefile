# $Id: Makefile 45 2014-01-30 17:55:23Z awalther $


RM := rm -rf

include config.mk

# Where to find the RTM. Replace this with your own path to the object files and modules.
RTM_INC = /home/awalther/SOFTWARE/RTM/RTTOV/include
RTM_MOD = /home/awalther/SOFTWARE/RTM/RTTOV/mod 
RTM_LIB = /home/awalther/SOFTWARE/RTM/RTTOV/lib
RTMODEL = rttov10.2.0

NCDF_LIB = /usr/local/lib

.SUFFIXES:
.SUFFIXES: .f90 .F .o

#Compilation rule for f90 files
.f90.o:
	$(fc) -c ${FFLAGS} ${fflags} ${hdf5libs} ${hdflibs}  $<
	
	
OBJS = dcomp.o \
dcomp_forward_mod.o \
dcomp_tools.o \
M_kracken.o \
dcomp_mod.o \
get_planck_radiance_39um.o \
get_rad_refl_factor.o \
dcomp_lut_def.o \
dcomp_lut_mod.o \
dcomp_array.o \
dcomp_interface_types_mod.o \
view2d.o


OBJS_TEST2 = test_dcomp.o \
dcomp_forward_mod.o \
dcomp_tools.o \
M_kracken.o \
dcomp_mod.o \
get_planck_radiance_39um.o \
get_rad_refl_factor.o \
dcomp_lut_mod.o \
dcomp_array.o \
dcomp_interface_types_mod.o

OBJS_TEST = test_forward_computation.o \
dcomp_forward_mod.o \
dcomp_tools.o \
hdf_dcomp_awg.o \
dcomp_mod.o \
get_planck_radiance_39um.o \
get_rad_refl_factor.o

OBJS_TEST_NEW = test_new_fwd.o \
dcomp_fwd_mod.o \
dcomp_lut_mod.o \
dcomp_tools.o  \
get_planck_radiance_39um.o \
get_rad_refl_factor.o 


OBJS_TEST_LUT = test_lut_mod.o \
dcomp_lut_mod.o

all: dcomp

test: $(OBJS_TEST)
	@echo 'test'
	$(fc) ${fflags}  $(FFLAGS)    -o "test_fwd" $(OBJS_TEST)  ${hdf5libs} ${hdflibs} 
   
test2: $(OBJS_TEST2)
	@echo 'test'
	$(fc) ${fflags}  $(FFLAGS)    -o "test2" $(OBJS_TEST2)  ${hdf5libs} ${hdflibs}   

test_new_fwd: $(OBJS_TEST_NEW)
	@echo 'fwd_test'
	$(fc) ${fflags}  $(FFLAGS)    -o "test_fwd_new" $(OBJS_TEST_NEW)  ${hdf5libs} ${hdflibs}
	
test_lut_mod: $(OBJS_TEST_LUT)
	@echo 'test lut'
	$(fc) ${fflags}  $(FFLAGS)    -o "test_lut" $(OBJS_TEST_LUT)  ${hdf5libs} ${hdflibs} 

dcomp: $(OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: Intel(R) Fortran Linker'
	$(fc) ${fflags}  $(FFLAGS)    -o "dcomp" $(OBJS)  ${hdf5libs} ${hdflibs} 
	@echo 'Finished building target: $@'
	@echo ' '	

dcomp.o: dcomp.f90 dcomp_mod.o get_planck_radiance_39um.o dcomp_lut_mod.o
test_dcomp.o: test_dcomp.f90 dcomp_mod.o get_planck_radiance_39um.o
dcomp_forward_mod.o: dcomp_forward_mod.f90  dcomp_tools.o  M_kracken.o get_rad_refl_factor.o dcomp_lut_def.o dcomp_lut_mod.o
dcomp_tools.o: dcomp_tools.f90
hdf_dcomp_awg.o: hdf_dcomp_awg.f90
test_forward_computation.o:test_forward_computation.f90 dcomp_forward_mod.o
M_kracken.o:M_kracken.f90
dcomp_mod.o:dcomp_mod.f90 dcomp_forward_mod.o dcomp_tools.o get_rad_refl_factor.o 
get_planck_radiance_39um.o:get_planck_radiance_39um.f90
get_rad_refl_factor.o:get_rad_refl_factor.f90
test_new_fwd.o:test_new_fwd.f90  dcomp_fwd_mod.o
dcomp_fwd_mod.o:dcomp_fwd_mod.f90 dcomp_tools.o get_planck_radiance_39um.o get_rad_refl_factor.o 
dcomp_lut_def.o:dcomp_lut_def.f90
dcomp_lut_mod.o:dcomp_lut_mod.f90 dcomp_lut_def.o
test_lut_mod.o:test_lut_mod.f90 dcomp_lut_mod.o
dcomp_array.o:dcomp_array.f90 dcomp_mod.o dcomp_interface_types_mod.o view2d.o
dcomp_interface_types_mod.o:dcomp_interface_types_mod.f90
view2d.o:view2d.f90

clean:
	rm -f *.o *.mod
	rm -f ../cr_start
	-@echo ' '
