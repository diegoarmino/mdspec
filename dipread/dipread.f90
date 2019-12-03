      program dipread
!        CALCULATES THE AUTOCORRELATION FUNCTION FOR VACUUM SPECTRA.

         implicit none
        
!        --------------------------------------------------------------
         integer                  :: nsteps
         character(90)            :: datafile
         character(90)            :: outfile

         integer  :: i,j,k,ierr,clcount
         real*8   :: mu(3)
         real*8   :: amu
!        --------------------------------------------------------------

!        READ COMAND LINE ARGUMENTS
         clcount= command_argument_count()
         if (clcount /= 2) then
            write(*,'(A)') 'WHATS THE INPUT, DUDE???'
            write(*,'(A)') 'USAGE: dipread <datafile> <outfile>'
      
            STOP
         end if

         call get_command_argument(1,datafile)
         call get_command_argument(2,outfile)

!        ASK THE USER TO INTRODUCE THE LENGTH OF THE TRAJECTORY TO READ.
         write(*,'(A)') 'TOTAL NUMBER OF STEPS IN THE DIPOLE &
         &TRAJECTORY?'
         read(*,*) nsteps

!        OPENING INPUT FILE
         open(unit=1,file=datafile,iostat=ierr)
         if (ierr /= 0) then
            write(*,'(A,A)') 'ERROR OPENING INPUT FILE ',datafile
            STOP
         end if

         open(unit=2,file=outfile,iostat=ierr)
         if (ierr /= 0) then
            write(*,'(A,A)') 'ERROR OPENING OUTPUT FILE ',outfile
            STOP
         end if

!        COMPUTING AVERAGE OF DIPOLE MOMENT COMPONENTS.
         do i=1,nsteps
            read(1,*) (mu(j),j=1,3)
            amu = (mu(1) + mu(2) + mu(3))/3d0
            write(2,'(F15.10)') amu
         end do

      end program
