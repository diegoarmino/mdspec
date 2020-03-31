# Fortran compilation

HOME=/home/diegoa
BASEDIR=$(HOME)/dev/mdspec
SRCDIR=$(BASEDIR)

FC=gfortran
FFLAGS= -O3 -ipo -xHost -FR
FFLAGS= -cpp -DBINTRAJ -O3 -mtune=native -fno-backtrace
INC= -I$(HOME)/progs/fftw/include -L$(HOME)/progs/fftw/lib -lfftw3
INCNC= `nf-config --fflags --flibs --has-nc4`
INCNC= `nf-config --fflags --flibs`
## INC += -I$(AMBERHOME)/include -L$(AMBERHOME)/lib -lnetcdf -lnetcdff
FLIBS= 
##gfortran   -DBINTRAJ -DEMIL -DPUBFFT -DGNU_HACKS -O3 -mtune=native   -I/home/diegoa/progs/amber18/include -c AmberNetcdf.F90


all: clean mdspec

mdspec: nextprmtop_section_mod.o mdspecNetcdf_mod.o pump_probe_mod.o linearSpectroscopy_mod.o mdspec_main.o
	$(FC) $(FFLAGS) -o $@ $^ $(FLIBS) $(INC) $(INCNC)

#----------------------------------------------------------------------------------------
# NEXTPRMTOP SECTION MODULE
#----------------------------------------------------------------------------------------
nextprmtop_section_mod.o: nextprmtop_section_mod.F90
	$(FC) $(FFLAGS) $(INC) -c $^ $(FLIBS)

#----------------------------------------------------------------------------------------
# NETCDF MODULE
#----------------------------------------------------------------------------------------
mdspecNetcdf_mod.o: mdspecNetcdf_mod.F90
	$(FC) $(FFLAGS) $(INC) -c $^ $(FLIBS) `nf-config --fflags --flibs`

#----------------------------------------------------------------------------------------
# PUMP-PROBE MODULE
#----------------------------------------------------------------------------------------
pump_probe_mod.o: nextprmtop_section_mod.o mdspecNetcdf_mod.o pump_probe_mod.F90 
	$(FC) $(FFLAGS) $(INC) -c $^ $(FLIBS) `nf-config --fflags --flibs`


#----------------------------------------------------------------------------------------
# LIBRARIES COMPILATION
#----------------------------------------------------------------------------------------
linearSpectroscopy_mod.o: linearSpectroscopy_mod.F90
	$(FC) $(FFLAGS) $(INC) -c $^ $(FLIBS)

#----------------------------------------------------------------------------------------
# MAIN PROGRAM COMPILATION
#----------------------------------------------------------------------------------------

mdspec_main.o: pump_probe_mod.o linearSpectroscopy_mod.o mdspec_main.F90
	$(FC) $(FFLAGS) $(INC) -c mdspec_main.F90 $(FLIBS) -o mdspec_main.o `nf-config --fflags --flibs`


#----------------------------------------------------------------------------------------

clean:
	rm -rf *.o *.mod $(SCRDIR)
