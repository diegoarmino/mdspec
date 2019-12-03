module mdspec_lib

   private
   public :: read_inputfile,vac_autocor,run_avg,fourier,print_spectra,mdsin_type
   
  
   type mdsin_type
      integer :: nsteps
      integer :: ntraj
      integer :: navg
      integer :: dodamp
      integer :: prefac
      integer :: convergence_check
      real*8  :: dt
      real*8  :: damp
   end type mdsin_type

   contains

      subroutine read_inputfile(inputfile,mdsin)
!        READS INPUT FILE CONTAINING NAMELIST WITH PARAMETERS 
         implicit none
       
!        --------------------------------------------------------------
         character(99),intent(in)      :: inputfile
         type(mdsin_type),intent(out)  :: mdsin

         integer :: nsteps
         integer :: ntraj
         integer :: navg
         integer :: dodamp
         integer :: prefac
         integer :: convergence_check
         real*8  :: dt
         real*8  :: damp

         namelist /mdspec/ nsteps,ntraj,navg,dt,damp,dodamp,prefac,convergence_check

         integer     :: ierr
!        --------------------------------------------------------------

         open(unit=2,file=inputfile,iostat=ierr)
         if (ierr /= 0) then
            write(*,'(A,A)') 'ERROR OPENING INPUT FILE ', inputfile
            STOP
         end if
        
         rewind 2
         read(2,nml=mdspec,iostat=ierr)
         if (ierr /= 0 ) then
            write(*,'(A)') 'ERROR READING NAMELIST'
            STOP
         end if

         mdsin%nsteps = nsteps
         mdsin%ntraj  = ntraj
         mdsin%navg   = navg
         mdsin%dodamp = dodamp
         mdsin%prefac = prefac
         mdsin%dt     = dt
         mdsin%damp   = damp
         mdsin%convergence_check   = convergence_check

         write(*,'(A)') 'INPUT PARAMETERS'
         write(*,'(A)') '----------------'
         write(*,'(A,I6)')'nsteps = ', mdsin%nsteps
         write(*,'(A,I6)')'ntraj = ', mdsin%ntraj
         write(*,'(A,I6)')'navg = ', mdsin%navg
         write(*,'(A,I6)')'dodamp = ', mdsin%dodamp
         write(*,'(A,I6)')'prefac = ', mdsin%prefac
         write(*,'(A,I6)')'convergence_check = ', mdsin%convergence_check
         write(*,'(A,F15.10)')'dt = ', mdsin%dt
         write(*,'(A,F15.10)')'damp = ', mdsin%damp
         write(*,'(A)') '----------------'
         
      end subroutine

      subroutine vac_autocor(datafile,acorr,nsteps,ntraj)
!        CALCULATES THE AUTOCORRELATION FUNCTION FOR VACUUM SPECTRA.

         implicit none
        
!        --------------------------------------------------------------
         integer,intent(in)       :: nsteps
         integer,intent(in)       :: ntraj
         character(90),intent(in) :: datafile
         real*8,intent(out)       :: acorr(nsteps-1)

         integer  :: i,j,k,ierr
         real*8   :: ref(ntraj)
         real*8   :: xdif(ntraj)
         real*8   :: xold(ntraj)
         real*8   :: xnew(ntraj)
         real*8   :: corr(nsteps-1)
!        --------------------------------------------------------------

!        OPENING INPUT FILE
         open(unit=1,file=datafile,iostat=ierr)
         if (ierr /= 0) then
            write(*,'(A,A)') 'ERROR OPENING INPUT FILE ',datafile
            STOP
         end if

!        FIRST VALUE OF AUTOCORRELATION FUNCTION
         corr=0d0
         write(*,'(A,I6)') 'READING DATAFILE FIRST TIME',1
         read(1,*) (xold(i),i=1,ntraj)
         read(1,*) (xnew(i),i=1,ntraj)
         ref = xnew - xold    ! difference
         do j=1,ntraj
            corr(1) = corr(1) + ref(j)**2  ! average over trajectories
         end do
         xold=xnew   ! updating x
         xnew=0

!        THE REST OF AUTOCORRELATION FUNCTION
         do i=2,nsteps-1
            write(*,'(A,I6)') 'READING DATAFILE ITERATIONS',i
!            read(1,'(F15.10)') (xnew(k),k=1,ntraj)
            read(1,*) (xnew(k),k=1,ntraj)
            xdif = xnew - xold   ! difference
            do j=1,ntraj
               corr(i) = corr(i) + xdif(j) * ref(j) ! average over trajectories
            end do
            xold = xnew
            xnew = 0d0
         end do
         acorr = corr/ntraj
         
!        CLOSING DATA FILE
         close(unit=1,iostat=ierr)
         if (ierr /= 0) then
            write(*,'(A,A)') 'ERROR CLOSING DATA FILE ',datafile
            STOP
         end if

         open(unit=10,file='autocorr.dat',iostat=ierr)
         if (ierr /= 0) then
            write(*,'(A,A)') 'ERROR OPENING INPUT FILE ',datafile
            STOP
         end if

         do i=1,nsteps-1
            write(10,'(I6,D20.10)') i,acorr(i)
         end do

      end subroutine

      subroutine run_avg(series,slen,navg,avgser)

!        CALCULATES A RUNNING AVERAGE OF A DATA SERIES.
         implicit none
         
!        --------------------------------------------------------------
         integer,intent(in)     :: slen         ! Length of the data series.
         integer,intent(in)     :: navg         ! # of steps to average over.
         real*8,intent(in)      :: series(slen) ! Data series.
         real*8,intent(out)     :: avgser(slen-navg+1)
 
         integer i,j
!        --------------------------------------------------------------

         do i=1,slen-navg+1
            do j=0,navg
               avgser(i) = avgser(i) + series(i+j)
            end do
         end do
         avgser = avgser / navg

      end subroutine

      subroutine fourier(series,dodamp,prefac,slen,dt,damp,rspec,ispec)
         implicit none

!        --------------------------------------------------------------
         integer,intent(in)        :: dodamp
         integer,intent(in)        :: prefac
         integer,intent(in)        :: slen
         real*8,intent(in)         :: dt
         real*8,intent(in)         :: damp
         real*8,intent(in)         :: series(slen)
         real*8,intent(out)        :: rspec(40000)
         real*8,intent(out)        :: ispec(40000)

         integer             :: i,j
         real*8              :: ddt
         real*8              :: dmp
         real*8              :: dp
         real*8              :: nu
         real*8              :: ftr
         real*8              :: fti
         real*8              :: t
         real*8,parameter    :: cc = 2.99792458d-5 ! speed of light cm/fs
         real*8,parameter    :: pi = 3.1415926535897932384626433832795d0
!        --------------------------------------------------------------

         ddt = dt 
         dmp = damp 

         do i = 1,10*4000
            nu = dble((0.1*i)) * cc
            ftr = 0d0
            fti = 0d0

            t=0d0
            do j=1,slen-1
               t = t + ddt

               dp = 1d0
               if (dodamp == 1) dp = dampfunc(damp,t)

               ftr = ftr + cos(2*pi*t*nu) * series(j) * dp
               fti = fti + sin(2*pi*t*nu) * series(j) * dp
            end do

            rspec(i) = ftr**2 + fti**2
!            rspec(i) = Sqrt(ftr**2 + fti**2)
!            rspec(i) = Abs(ftr) + Abs(fti)
            if (prefac == 1) rspec(i) = qmfac(nu,300d0) * rspec(i)
            ispec(i) = fti**2

         end do

      end subroutine

      real*8 function dampfunc(damp,t) 
         real*8  :: damp
         real*8  :: t
         dampfunc = exp(-t/damp)
         return
      end function

      real*8 function qmfac(nu,temp)
         real*8            :: nu
         real*8            :: temp
         real*8            :: beta
         real*8            :: num
         real*8            :: den
         real*8,parameter  :: kb = 0.69503476
         real*8,parameter  :: pi = 3.1415926535897932384626433832795d0

         beta = 1d0/(kb*temp)
         num = (beta*nu)**2
         den = pi*(1d0 - exp(-beta*nu))

         qmfac = num/den
         return
      end function
         
      subroutine print_spectra(slen,rspec,ispec,spectrafile)
         implicit none

!        --------------------------------------------------------------
         integer,intent(in)         :: slen
         real*8,intent(in)          :: rspec(slen)
         real*8,intent(in)          :: ispec(slen)
         character(99),intent(in)   :: spectrafile

         integer  :: ierr,i
         real*8   :: nu
!        --------------------------------------------------------------
 
         open(unit=3,file=spectrafile,iostat=ierr)
         if (ierr /= 0) then
            write(*,'(A,A)') 'ERROR OPENNING FILE ',spectrafile
            STOP
         end if

         do i = 1,10*4000
            nu = dble((0.1*i))
            write(3,'(F10.2,2D17.8)') nu, rspec(i), ispec(i)
         end do

         close(unit=3,iostat=ierr)
         if (ierr /= 0) then
            write(*,'(A,A)') 'ERROR CLOSING FILE ',spectrafile
            STOP
         end if

      end subroutine

end module mdspec_lib
