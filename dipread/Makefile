# Fortran compilation

FC=ifort
FFLAGS= -O3 -ipo -xHost -FR 
FLIBS=

all: dipread

dipread: dipread.f90
	$(FC) $(FFLAGS) -o $@ $^ $(FLIBS)
clean:
	rm -rf *.o *.mod $(SCRDIR)
