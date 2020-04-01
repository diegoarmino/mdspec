#gfortran aka2.f90 -o aka2 `nf-config --fflags --flibs`
gfortran aka3.f90 -o aka3 `nf-config --fflags --flibs` -I${HOME}/progs/fftw/include -L${HOME}/progs/fftw/lib -lfftw3 nextprmtop_section_mod.o
