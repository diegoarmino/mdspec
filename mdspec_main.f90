   program mdspec
      use mdspec_lib

      implicit none
!     ------------------------------------------------------------------
      type(mdsin_type)       :: mdsin
      character(99)          :: inputfile
      character(99)          :: datafile
      character(99)          :: specfile
      integer                :: clcount
      integer                :: acdim
      integer                :: avgdim
      integer                :: ierr
      integer                :: i,j
      real*8                 :: intrmsd

      real*8,allocatable     :: acorr(:)
      real*8,allocatable     :: avgcorr(:)
      real*8,allocatable     :: rspec(:)
      real*8,allocatable     :: ispec(:)
      real*8,allocatable     :: newspec(:)
      real*8,allocatable     :: oldspec(:)
      real*8,allocatable     :: rmsdspec(:)
      real*8,allocatable     :: convplot(:)
!     ------------------------------------------------------------------

!     DEFAULT FILENAMES
      inputfile='input'
      datafile='data'
      specfile='spectra'

      clcount= command_argument_count()
      if (clcount /= 0 .AND. clcount /= 3) then
         write(*,'(A)') 'USAGE:'
         write(*,'(A)') 'mdspec requires either no arguments or three  '
         write(*,'(A)') 'as follows:'
         write(*,'(A)') ''
         write(*,'(A)') 'mdspec <input_file> <data_file> <spectra_file>'
         write(*,'(A)') ''
         write(*,'(A)') 'If no arguments are given, default filenames  '
         write(*,'(A)') 'are used: "input", "data" and "spectra".      '
      end if

!     READ FILE NAMES FROM COMMAND LINE 
      call get_command_argument(1,inputfile)
      call get_command_argument(2,datafile)
      call get_command_argument(3,specfile)


!     READ INPUT FILE
      call read_inputfile(inputfile,mdsin)

      acdim=mdsin%nsteps - 1
      avgdim=acdim-mdsin%navg+1
      allocate( acorr(acdim),avgcorr(avgdim),rspec(40000),ispec(40000) )

!     Initialize
      acorr   =0d0
      avgcorr =0d0
      rspec   =0d0
      ispec   =0d0

!     Calculate auto correlation function vacuum style.
      call vac_autocor(datafile,acorr,mdsin%nsteps,mdsin%ntraj)

!     Calculate running average of said autocorr.
      if (mdsin%navg > 1) then
         call run_avg(acorr,acdim,mdsin%navg,avgcorr)
      else
         avgcorr=acorr
      end if

!     Calculate real and imaginary Fourier transform (rspec and ispec).
      call fourier(avgcorr,mdsin%dodamp,mdsin%prefac,avgdim,mdsin%dt,& 
                   &mdsin%damp,rspec,ispec)

!     Print spectra
      call print_spectra(40000,rspec,ispec,specfile)

!     Convergence check
    write(*,'(A,2I)') 'convergence check option = ',mdsin%convergence_check
    IF (mdsin%convergence_check == 1) THEN
      allocate( newspec(40000),oldspec(40000),rmsdspec(40000),convplot(mdsin%ntraj) )

      do i=1,mdsin%ntraj
        !Initialize
        acorr   =0d0
        avgcorr =0d0
        rspec   =0d0
        ispec   =0d0
        newspec =0d0

!       Calculate auto correlation function vacuum style.
        call vac_autocor(datafile,acorr,mdsin%nsteps,i)

!       Calculate running average of said autocorr.
        if (mdsin%navg > 1) then
           call run_avg(acorr,acdim,mdsin%navg,avgcorr)
        else
           avgcorr=acorr
        end if

!       Calculate Fourier transform (rspec).
        if (i==1) then
           call fourier(avgcorr,mdsin%dodamp,mdsin%prefac,avgdim,mdsin%dt,&
                        &mdsin%damp,oldspec,ispec)
        else
           call fourier(avgcorr,mdsin%dodamp,mdsin%prefac,avgdim,mdsin%dt,&
                        &mdsin%damp,newspec,ispec)
        end if

!       Compute spectrum of RMSDs
        rmsdspec = Sqrt((newspec - oldspec)**2)

!       Integrate
        intrmsd = 0d0
        do j=1,40000
           intrmsd = intrmsd + rmsdspec(j)
        end do

        convplot(i) = intrmsd

        oldspec = newspec
        newspec = 0d0

        write(*,'(A,I5,A,I5)') 'STEP ', i,'/',mdsin%ntraj

     end do
    
     open(unit=4,file='convergence_plot.dat',iostat=ierr)
     if (ierr /= 0) then
        write(*,'(A,A)') 'ERROR OPENNING FILE ','convergence_plot.dat'
        STOP
     end if

     do i = 1,mdsin%ntraj
        write(4,'(I5,D17.8)') i, log10(convplot(i))
     end do

     close(unit=4,iostat=ierr)
     if (ierr /= 0) then
        write(*,'(A,A)') 'ERROR CLOSING FILE ','convergence_plot.dat'
        STOP
     end if

   END IF

     end program mdspec
