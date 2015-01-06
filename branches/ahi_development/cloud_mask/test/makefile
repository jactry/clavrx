FC=gfortran
FCFLAGS = -g

SRCS = ../main_src/file_tools.f90

PROGRAM = test_cloud_mask.o

all: $(PROGRAM)




test_cloud_mask.o: test_cloud_mask.f90 nb_cloud_mask_module.o	file_tools.o
	$(FC) test_cloud_mask.f90 nb_cloud_mask_module.o file_tools.o -o test_mask
file_tools.o: file_tools.f90
	$(FC) -c file_tools.f90 
nb_cloud_mask_module.o: ../nb_cloud_mask_module.f90 file_tools.o 
	$(FC) -c ../nb_cloud_mask_module.f90 file_tools.o 
	

clean:
	rm -f *.o *.mod test_mask
