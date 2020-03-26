# Fortran compilation

HOME=/home/diegoa
BASEDIR=$(HOME)/dev/mdspec
SRCDIR=$(BASEDIR)

FC=gfortran
FFLAGS= -O3 -ipo -xHost -FR
FFLAGS= -cpp -DBINTRAJ -O3 -mtune=native 
INC= -I$(HOME)/progs/fftw/include -L$(HOME)/progs/fftw/lib -lfftw3
INC += -I$(AMBERHOME)/include -L$(AMBERHOME)/lib -lnetcdf -lnetcdff
FLIBS= 
##gfortran   -DBINTRAJ -DEMIL -DPUBFFT -DGNU_HACKS -O3 -mtune=native   -I/home/diegoa/progs/amber18/include -c AmberNetcdf.F90


all: clean mdspec

mdspec: mdspecNetcdf_mod.o pump_probe_mod.o linearSpectroscopy_mod.o mdspec_main.o
	$(FC) $(FFLAGS) -o $@ $^ $(FLIBS) $(INC)

#----------------------------------------------------------------------------------------
# NETCDF MODULE
#----------------------------------------------------------------------------------------
mdspecNetcdf_mod.o: mdspecNetcdf_mod.F90
	$(FC) $(FFLAGS) $(INC) -c $^ $(FLIBS)

#----------------------------------------------------------------------------------------
# PUMP-PROBE MODULE
#----------------------------------------------------------------------------------------
pump_probe_mod.o: mdspecNetcdf_mod.o pump_probe_mod.F90 
	$(FC) $(FFLAGS) $(INC) -c $^ $(FLIBS) 

#----------------------------------------------------------------------------------------
# LIBRARIES COMPILATION
#----------------------------------------------------------------------------------------
linearSpectroscopy_mod.o: linearSpectroscopy_mod.F90
	$(FC) $(FFLAGS) $(INC) -c $^ $(FLIBS)

#----------------------------------------------------------------------------------------
# MAIN PROGRAM COMPILATION
#----------------------------------------------------------------------------------------

mdspec_main.o: linearSpectroscopy_mod.o mdspec_main.F90
	$(FC) $(FFLAGS) $(INC) -c mdspec_main.F90 $(FLIBS) -o mdspec_main.o
#----------------------------------------------------------------------------------------

clean:
	rm -rf *.o *.mod $(SCRDIR)
