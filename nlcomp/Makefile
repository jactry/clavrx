# $Header: https://svn.ssec.wisc.edu/repos/cloud_team_nlcomp/trunk/Makefile 8 2014-01-31 08:14:58Z awalther $


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
	
	
OBJS = nlcomp.o \
nlcomp_forward_mod.o \
dcomp_tools.o \
M_kracken.o \
nlcomp_mod.o \
get_planck_radiance_39um.o \
get_rad_refl_factor.o \
nlcomp_lut_mod.o \
nlcomp_array.o \
nlcomp_interface_def_mod.o

OBJS_TEST = test_forward_computation.o \
nlcomp_forward_mod.o \
dcomp_tools.o \
hdf_nlcomp_awg.o \
nlcomp_mod.o \
get_planck_radiance_39um.o \
get_rad_refl_factor.o

OBJS_TEST_NEW = test_new_fwd.o \
nlcomp_fwd_mod.o \
nlcomp_lut_mod.o \
nlcomp_tools.o  \
get_planck_radiance_39um.o \
get_rad_refl_factor.o 


OBJS_TEST_LUT = test_lut_mod.o \
nlcomp_lut_mod.o

all: nlcomp

test: $(OBJS_TEST)
	@echo 'test'
	$(fc) ${fflags}  $(FFLAGS)    -o "test_fwd" $(OBJS_TEST)  ${hdf5libs} ${hdflibs} 

test_new_fwd: $(OBJS_TEST_NEW)
	@echo 'fwd_test'
	$(fc) ${fflags}  $(FFLAGS)    -o "test_fwd_new" $(OBJS_TEST_NEW)  ${hdf5libs} ${hdflibs}
	
test_lut_mod: $(OBJS_TEST_LUT)
	@echo 'test lut'
	$(fc) ${fflags}  $(FFLAGS)    -o "test_lut" $(OBJS_TEST_LUT)  ${hdf5libs} ${hdflibs} 

nlcomp: $(OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: Intel(R) Fortran Linker'
	$(fc) ${fflags}  $(FFLAGS)    -o "nlcomp" $(OBJS)  ${hdf5libs} ${hdflibs} 
	@echo 'Finished building target: $@'
	@echo ' '	

nlcomp.o: nlcomp.f90  nlcomp_mod.o get_planck_radiance_39um.o nlcomp_interface_def_mod.o
nlcomp_forward_mod.o: nlcomp_forward_mod.f90  dcomp_tools.o  M_kracken.o get_rad_refl_factor.o nlcomp_lut_mod.o
dcomp_tools.o: dcomp_tools.f90
hdf_dcomp_awg.o: hdf_dcomp_awg.f90
test_forward_computation.o:test_forward_computation.f90 nlcomp_forward_mod.o
M_kracken.o:M_kracken.f90
nlcomp_mod.o:nlcomp_mod.f90 nlcomp_forward_mod.o dcomp_tools.o get_rad_refl_factor.o 
get_planck_radiance_39um.o:get_planck_radiance_39um.f90
get_rad_refl_factor.o:get_rad_refl_factor.f90
test_new_fwd.o:test_new_fwd.f90  nlcomp_fwd_mod.o
nlcomp_fwd_mod.o:nlcomp_fwd_mod.f90 nlcomp_lut_mod.o dcomp_tools.o get_planck_radiance_39um.o get_rad_refl_factor.o 
nlcomp_lut_mod.o:nlcomp_lut_mod.f90 
test_lut_mod.o:test_lut_mod.f90 nlcomp_lut_mod.o
nlcomp_array.o:nlcomp_array.f90 nlcomp_mod.o nlcomp_interface_def_mod.o
nlcomp_interface_def_mod.o:nlcomp_interface_def_mod.f90

clean:
	rm -f *.o *.mod
	rm -f ../cr_start
	-@echo ' '
