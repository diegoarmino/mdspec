# Fortran compilation

HOME=/home/diegoa
BASEDIR=$(HOME)/dev/mdspec
SRCDIR=$(BASEDIR)

FC=ifort
FFLAGS= -O3 -ipo -xHost -FR 
FLIBS=

all: clean mdspec

mdspec: mdspec_lib.o mdspec_main.o
	$(FC) $(FFLAGS) -o $@ $^ $(FLIBS)

#----------------------------------------------------------------------------------------
# LIBRARIES COMPILATION
#----------------------------------------------------------------------------------------
mdspec_lib.o: mdspec_lib.f90
	$(FC) $(FFLAGS) -c $^ $(FLIBS)

#----------------------------------------------------------------------------------------
# MAIN PROGRAM COMPILATION
#----------------------------------------------------------------------------------------

mdspec_main.o: mdspec_lib.o mdspec_main.f90
	$(FC) $(FFLAGS) -c mdspec_main.f90 $(FLIBS) -o mdspec_main.o

#----------------------------------------------------------------------------------------

clean:
	rm -rf *.o *.mod $(SCRDIR)
