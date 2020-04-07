#gfortran aka2.f90 -o aka2 `nf-config --fflags --flibs`
#gfortran aka3.f90 -o aka3 `/usr/bin/nf-config --fflags --flibs` -I${HOME}/progs/fftw/include -L${HOME}/progs/fftw/lib -lfftw3 nextprmtop_section_mod.o
#gfortran aka3.f90 -o aka3 -I/home/diegoa/progs/amber18/include -L/home/diegoa/progs/amber18/lib -lnetcdff -lnetcdf -Wl,-Bsymbolic-functions -Wl,-z,relro -I${HOME}/progs/fftw/include -L${HOME}/progs/fftw/lib -lfftw3 nextprmtop_section_mod.o
gfortran aka3.f90 -o aka3 `/usr/bin/nf-config --fflags --flibs` -O3 -I${HOME}/progs/fftw/include -L${HOME}/progs/fftw/lib -lfftw3 nextprmtop_section_mod.o

