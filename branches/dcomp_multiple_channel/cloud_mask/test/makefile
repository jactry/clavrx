FC=ifort
FCFLAGS = -g

SRCS = ../main_src/file_tools.f90

PROGRAM = test_cloud_mask.o

all: $(PROGRAM)




test_cloud_mask.o: test_cloud_mask.f90 naive_bayesian_cloud_mask_module.o	file_tools.o
	$(FC) test_cloud_mask.f90 naive_bayesian_cloud_mask_module.o ../../main_src/file_tools.o -o test_mask
file_tools.o: ../../main_src/file_tools.f90
	$(FC) -c ../../main_src/file_tools.f90 
naive_bayesian_cloud_mask_module.o: ../naive_bayesian_cloud_mask_module.f90 file_tools.o
	$(FC) -c ../naive_bayesian_cloud_mask_module.f90 file_tools.o

	

clean:
	rm -f *.o *.mod test_mask
