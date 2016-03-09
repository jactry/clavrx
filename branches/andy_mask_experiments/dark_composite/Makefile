optflags        = -O2 -assume byterecl 
f90     = ifort

#compile rules
#-------------
.SUFFIXES: .o .f90
.f90.o:
	$(f90) -c $(optflags) $<

#executables
#-----------
default: all

#all: batch_dark_ch1 dark_ch1
all: dark_ch1

batch_dark_ch1: batch_dark_ch1.o
	$(f90) -o batch_dark_ch1 batch_dark_ch1.o

#for this one, order of the .o's is significant because dependencies may not be complete
dark_ch1: constant.o goes_area_io_mod.o num_mod_gsip.o subs_new.o dark_ch1_gsip.o 
	$(f90) -o dark_ch1 constant.o dark_ch1_gsip.o num_mod_gsip.o subs_new.o goes_area_io_mod.o 

clean:
	-rm *.o *.mod

distclean: clean
	-rm batch_dark_ch1 dark_ch1

#dependencies for .o files
#-------------------------
batch_dark_ch1.o: batch_dark_ch1.f90

constant.o: constant.f90

dark_ch1_gsip.o: dark_ch1_gsip.f90 constant.o num_mod_gsip.o goes_area_io_mod.o

goes_area_io_mod.o: goes_area_io_mod.f90 constant.o

num_mod_gsip.o: num_mod_gsip.f90 constant.o

subs_new.o: subs_new.f90 constant.o goes_area_io_mod.o

