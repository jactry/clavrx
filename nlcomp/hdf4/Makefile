FC=gfortran
   ifdef COMP1
		$(fc)=COMP1 
	endif  

all : lib 

lib :
	echo $(fc)
	cd src && make COMP1="$(fc)"

clean :
	cd src && make clean
	


