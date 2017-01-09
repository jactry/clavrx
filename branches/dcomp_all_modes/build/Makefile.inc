


#EXT_BIN=.exe
#EXT_MOD=.mod
#EXT_OBJ=.o

#CLAVRX_VERSION=6.0.0
#EXT_LIB=.a
#MKINTF=$(TOPDIR)/build/mkintf.pl
#MVDMOD=$(TOPDIR)/build/mvdmod.pl
#MYPCPP=$(TOPDIR)/build/mypcpp.pl
#CPINCH=$(TOPDIR)/build/cpinch.pl

#FFLAG_MOD=-I
#FFLAG_D=-D
#CPP=cpp

#BLDDIR=$(TOPDIR)/$(INSTALLDIR)

#BINDIR=$(BLDDIR)/bin
#LIBDIR=$(BLDDIR)/lib
#OBJDIR=$(BLDDIR)/obj
#INCDIR=$(BLDDIR)/include
#MODDIR=$(BLDDIR)/mod

#TMPDIR=$(TOPDIR)/tmp-$(ARCH)/$(SUBDIR)
#SRCDIR=$(TOPDIR)/src/$(SUBDIR)



# default for compiler architecture macro is gfortran
# any value given on the command line gets higher priority
ARCH=gfortran
include $(TOPDIR)/build/arch/$(ARCH)

MAKEFILE_LOCAL=$(TOPDIR)/build/Makefile.local
# define local flags in MAKEFILE_LOCAL
include $(MAKEFILE_LOCAL)
