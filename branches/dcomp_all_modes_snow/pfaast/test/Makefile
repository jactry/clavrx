
fc = ifort -g -O0 -assume byterecl 

main_dir=../../main_src/
pfaast_dir = ../

.SUFFIXES: .f90 .F .o

 
#Compilation rule for f90 files
.f90.o:
	$(fc) -c  $<
   
OBJS = precision.o \
cx_pfaast_coef_mod.o \
cx_pfaast_test.o \
cx_string_tools_mod.o 

OBJS_RTM = precision.o \
cx_pfaast_coef_mod.o \
cx_pfaast_rtm_test.o \
cx_string_tools_mod.o \
cx_pfaast_mod.o \
cx_pfaast_tools_mod.o \
cx_pfaast_constants_mod.o



.PHONY: test clean

test: $(OBJS)
	@echo 'Building target: $@'
	$(fc) -o "test_pfaast" $(OBJS) -I../main_src/

test_rtm: $(OBJS_RTM)
	@echo 'Building target: $@'
	$(fc) -o "test_rtm" $(OBJS_RTM) -I../main_src/

cx_string_tools_mod.o: $(main_dir)/cx_string_tools_mod.f90  precision.o
	$(fc) -c $(main_dir)/cx_string_tools_mod.f90
precision.o: $(main_dir)/precision.f90
	$(fc) -c $(main_dir)/precision.f90    
cx_pfaast_coef_mod.o: $(pfaast_dir)/cx_pfaast_coef_mod.f90 cx_string_tools_mod.o
	$(fc) -c $(pfaast_dir)/cx_pfaast_coef_mod.f90 cx_string_tools_mod.o
cx_pfaast_test.o:cx_pfaast_test.f90 cx_pfaast_coef_mod.o
cx_pfaast_mod.o:$(pfaast_dir)/cx_pfaast_mod.f90 cx_pfaast_tools_mod.o cx_pfaast_constants_mod.o
	$(fc) -c $(pfaast_dir)/cx_pfaast_mod.f90 cx_pfaast_tools_mod.o cx_pfaast_constants_mod.o
cx_pfaast_tools_mod.o:$(pfaast_dir)/cx_pfaast_tools_mod.f90 cx_pfaast_constants_mod.o
	$(fc) -c $(pfaast_dir)/cx_pfaast_tools_mod.f90 cx_pfaast_constants_mod.o
cx_pfaast_constants_mod.o:$(pfaast_dir)/cx_pfaast_constants_mod.f90 
	$(fc) -c $(pfaast_dir)/cx_pfaast_constants_mod.f90      
cx_pfaast_rtm_test.o:cx_pfaast_rtm_test.f90 cx_pfaast_mod.o cx_pfaast_constants_mod.o


clean:
	rm -f *.o *.mod
	rm test_pfaast
