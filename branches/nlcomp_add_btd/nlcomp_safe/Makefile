# $Header: https://svn.ssec.wisc.edu/repos/cloud_team_nlcomp/trunk/Makefile 8 2014-01-31 08:14:58Z awalther $


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
	$(fc) -c ${CFLAGS} ${fflags} ${hdf5libs} ${hdflibs}  $(COMSRC)  $<
	
	
OBJS = nlcomp_one_pixel_run_program.o \
nlcomp_retrieval_mod.o \
nlcomp_array_loop_sub.o \
nlcomp_forward_mod.o \
dcomp_math_tools_mod.o \
dcomp_science_tools_mod.o \
M_kracken.o \
dcomp_lut_mod.o \
nlcomp_interface_def_mod.o \
file_tools.o \
clavrx_planck_mod.o  \
dcomp_lut_hdf_mod.o

.PHONY: clean hdf4

all: hdf4 nlcomp_one_pixel_run

hdf4:
	@echo  'kkf'$(LFLAGS)
	@echo 'Building target for hdf4: $@'
	make COMP1="$(fc)" -C $(ICASRC)	

nlcomp_one_pixel_run: $(OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: Intel(R) Fortran Linker'
	$(fc) ${fflags}  $(FFLAGS)    -o "nlcomp_one_pixel_run" $(OBJS)  $(LFLAGS) ${hdf5libs} ${hdflibs} $(COMOBJ) $(PATH_FLAGS) 
	@echo 'Finished building target: $@'
	@echo ' '	

nlcomp_one_pixel_run_program.o: nlcomp_one_pixel_run_program.f90  nlcomp_retrieval_mod.o nlcomp_interface_def_mod.o
nlcomp_retrieval_mod.o:nlcomp_retrieval_mod.f90 nlcomp_forward_mod.o dcomp_math_tools_mod.o dcomp_science_tools_mod.o 
nlcomp_array_loop_sub.o:nlcomp_array_loop_sub.f90 
nlcomp_forward_mod.o: nlcomp_forward_mod.f90  dcomp_math_tools_mod.o  M_kracken.o dcomp_lut_mod.o dcomp_science_tools_mod.o clavrx_planck_mod.o
dcomp_math_tools_mod.o: dcomp_math_tools_mod.f90
dcomp_lut_mod.o:dcomp_lut_mod.f90 dcomp_math_tools_mod.o file_tools.o dcomp_lut_hdf_mod.o
dcomp_science_tools_mod.o:dcomp_science_tools_mod.f90
nlcomp_interface_def_mod.o:nlcomp_interface_def_mod.f90
clavrx_planck_mod.o:clavrx_planck_mod.f90
dcomp_lut_hdf_mod.o:dcomp_lut_hdf_mod.f90

file_tools.o:file_tools.f90


M_kracken.o:$(SUB_PATH)/M_kracken.f90
	$(fc) -c ${fflags}  $(SUB_PATH)/M_kracken.f90

clean:
	rm -f *.o *.mod
	rm -f nlcomp_one_pixel_run
	rm -f *.a
	rm -rf hdf4/lib   
	make clean  -C $(ICASRC)
   
	-@echo ' '
