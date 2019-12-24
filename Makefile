# Fortran compilation

HOME=/home/diegoa
BASEDIR=$(HOME)/dev/mdspec
SRCDIR=$(BASEDIR)

FC=ifort
FFLAGS= -O3 -ipo -xHost -FR
FFLAGS=
INC= -I$(HOME)/progs/fftw/include -L$(HOME)/progs/fftw/lib -lfftw3
FLIBS=

all: clean mdspec

mdspec: mdspec_lib.o mdspec_main.o
	$(FC) $(FFLAGS) -o $@ $^ $(FLIBS) $(INC)

#----------------------------------------------------------------------------------------
# LIBRARIES COMPILATION
#----------------------------------------------------------------------------------------
mdspec_lib.o: mdspec_lib.f90
	$(FC) $(FFLAGS) $(INC) -c $^ $(FLIBS)

#----------------------------------------------------------------------------------------
# MAIN PROGRAM COMPILATION
#----------------------------------------------------------------------------------------

mdspec_main.o: mdspec_lib.o mdspec_main.f90
	$(FC) $(FFLAGS) $(INC) -c mdspec_main.f90 $(FLIBS) -o mdspec_main.o
#----------------------------------------------------------------------------------------

clean:
	rm -rf *.o *.mod $(SCRDIR)
