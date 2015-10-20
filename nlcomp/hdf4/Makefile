@echo 'Building ....'

FC=ifort
   ifdef COMP1
		$(fc)=COMP1 
	endif  
   
   ifdef HDF4R
      $(HDF4)=HDF4R
   endif

all : lib 

lib :
	echo $(fc)
	cd src && make COMP1="$(fc)" HDF4R="$(HDF4)"

clean :
	cd src && make clean
	


