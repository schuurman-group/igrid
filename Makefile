#-----------------------------------------------------------------------
# Compiler flags
#-----------------------------------------------------------------------

#
# gfortran
#
F90	= gfortran
F77	= gfortran
CC	= gcc
F90OPTS = -cpp -g -ffixed-line-length-none -ffree-line-length-none -O3 -fbacktrace
CCOPTS  = -g -O0

#
# ifort
#
#F90	 = ifort
#F77	 = ifort
#CC	 = icc
#F90OPTS = -cpp -g -free -fopenmp -traceback -O3 -diag-disable 8290 -diag-disable 8291
#CCOPTS  = -g -O0

# External libraries
LIBS= -lblas -llapack
#LIBS = -mkl

#-----------------------------------------------------------------------
# Define object files
#-----------------------------------------------------------------------

INCLUDE = include/constants.o \
	include/channels.o \
	include/sysinfo.o

IOMODULES = iomodules/iomod.o \
	iomodules/parsemod.o \

IOQC = ioqc/ioqc.o

UTILITIES = utilities/timingmod.o \
	utilities/utils.o

MKCUT = mkcut/mkcutglobal.o \
	mkcut/mkcut.o

IGRID = igrid/igridglobal.o \
	igrid/eigen.o \
	igrid/spline.o \
	igrid/interpolation.o \
	igrid/igrid.o

OBJECTS_MKCUT = $(INCLUDE) \
	$(IOMODULES) \
	$(UTILITIES) \
	$(IOQC) \
	$(MKCUT)
OBJ_MKCUT = constants.o \
	channels.o \
	sysinfo.o \
	timingmod.o \
	iomod.o \
	parsemod.o \
	utils.o \
	ioqc.o \
	mkcutglobal.o \
	mkcut.o

OBJECTS_IGRID = $(INCLUDE) \
	$(IOMODULES) \
	$(UTILITIES) \
	$(IOQC) \
	$(IGRID)
OBJ_IGRID = constants.o \
	channels.o \
	sysinfo.o \
	timingmod.o \
	iomod.o \
	parsemod.o \
	utils.o \
	ioqc.o \
	eigen.o \
	igridglobal.o \
	spline.o \
	interpolation.o \
	igrid.o

#-----------------------------------------------------------------------
# Rules to create the program
#-----------------------------------------------------------------------
mkcut: $(OBJECTS_MKCUT)
	$(F90) $(F90OPTS) $(OBJ_MKCUT) $(LIBS) -o bin/mkcut
	rm -f *.o *~ *.mod 2>/dev/null

igrid: $(OBJECTS_IGRID)
	$(F90) $(F90OPTS) $(OBJ_IGRID) $(LIBS) -o bin/igrid
	rm -f *.o *~ *.mod 2>/dev/null

%.o: %.f90
	$(F90) -c $(F90OPTS) $<

%.o: %.f
	$(F77) -c $(F77OPTS) $<

%.o: %.c
	$(CC) $(CCOPTS)  -c $<

clean_all:
	rm -f *.o *~ *.mod
