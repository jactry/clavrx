FC=gfortran
#FC=ifort
FCFLAGS = -g
NETCDFLIBS = -L/Users/awalther/lib/netcdf3/lib -I/Users/awalther/lib/netcdf3/include -lnetcdf

SRCS = ../main_src/file_tools.f90

PROGRAM = test_cloud_mask.o

all: $(PROGRAM)




test_cloud_mask.o: test_cloud_mask.f90 nb_cloud_mask_module.o nb_cloud_mask_addons_module.o nb_cloud_mask_clavrx_services_module.o file_tools.o
	$(FC)  test_cloud_mask.f90 -o test_mask file_tools.o nb_cloud_mask_clavrx_services_module.o nb_cloud_mask_addons_module.o nb_cloud_mask_module.o $(NETCDFLIBS)
file_tools.o: file_tools.f90
	$(FC) -c file_tools.f90 
nb_cloud_mask_clavrx_services_module.o: ../nb_cloud_mask_clavrx_services_module.f90
	$(FC) -c ../nb_cloud_mask_clavrx_services_module.f90
nb_cloud_mask_addons_module.o: ../nb_cloud_mask_addons_module.f90 nb_cloud_mask_clavrx_services_module.o
	$(FC) -c ../nb_cloud_mask_addons_module.f90
nb_cloud_mask_module.o: ../nb_cloud_mask_module.f90 nb_cloud_mask_clavrx_services_module.o file_tools.o 
	$(FC) -c $(NETCDFLIBS) ../nb_cloud_mask_module.f90 
	

clean:
	rm -f *.o *.mod test_mask
