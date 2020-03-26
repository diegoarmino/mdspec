module linearSpectroscopy_mod

   private
   public :: read_inputfile,get_spectra_fast,vac_autocor,run_avg,fourier,print_spectra,mdsin_type


   type mdsin_type
      integer :: nsteps
      integer :: method
      integer :: ntraj
      integer :: navg
      integer :: dodamp
      integer :: prefac
      integer :: convergence_check
      double precision  :: dt
      double precision  :: damp
   end type mdsin_type

   contains

      subroutine read_inputfile(inputfile,mdsin)
!        READS INPUT FILE CONTAINING NAMELIST WITH PARAMETERS
         implicit none

!        --------------------------------------------------------------
         character(99),intent(in)      :: inputfile
         type(mdsin_type),intent(out)  :: mdsin

         integer :: nsteps
         integer :: method
         integer :: ntraj
         integer :: navg
         integer :: dodamp
         integer :: prefac
         integer :: convergence_check
         double precision  :: dt
         double precision  :: damp

         namelist /mdspec/ nsteps,ntraj,navg,dt,damp,dodamp,prefac,convergence_check,method

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
         mdsin%method = method
         mdsin%convergence_check   = convergence_check

         write(*,'(A)') 'INPUT PARAMETERS'
         write(*,'(A)') '----------------'
         write(*,'(A,I6)')'nsteps = ', mdsin%nsteps
         write(*,'(A,I6)')'ntraj = ', mdsin%ntraj
         write(*,'(A,I6)')'navg = ', mdsin%navg
         write(*,'(A,I6)')'dodamp = ', mdsin%dodamp
         write(*,'(A,I6)')'prefac = ', mdsin%prefac
         write(*,'(A,I6)')'convergence_check = ', mdsin%convergence_check
         write(*,'(A,I6)')'method = ', mdsin%method
         write(*,'(A,F15.10)')'dt = ', mdsin%dt
         write(*,'(A,F15.10)')'damp = ', mdsin%damp
         write(*,'(A)') '----------------'

      end subroutine

      subroutine get_spectra_fast(datafile,nsteps,ncoords,dt)
!        CALCULATES THE FOURIER TRANSFORM OF THE AUTOCORRELATION FUNCTION
!        AS THE SQUARE OF THE FOURIER TRANSFORM OF THE VELOCITY TRAJECTORY.

         use, intrinsic :: iso_c_binding
         implicit none

         include 'fftw3.f03'
!        --------------------------------------------------------------
         integer,intent(in)       :: nsteps
         integer,intent(in)       :: ncoords
         character(90),intent(in) :: datafile
         double precision         :: spectra(nsteps/2+1)
         double precision         :: dt

         character(20)      :: fmt
         integer  :: i,j,k,ierr,t,tau,icrd,maxlag
         complex*16         :: out(nsteps/2+1)
         double precision   :: cumul(nsteps/2+1)
         double precision   :: noise
         double precision   :: sample_rate
         double precision   :: nu_increment
         double precision   :: nu
         double precision   :: vel(nsteps)
         double precision   :: corr(nsteps-1)
         double precision   :: PI,nu1,nu2,nu3,tt

!        --------------------------------------------------------------
         integer*8 :: plan
         !type(C_PTR) :: plan
!        --------------------------------------------------------------

!        OPENING INPUT FILE
         open(unit=1,file=datafile,iostat=ierr)
         if (ierr /= 0) then
            write(*,'(A,A)') 'ERROR OPENING INPUT FILE ',datafile
            STOP
         end if

         open(unit=100,file='spectra.dat',iostat=ierr)
         if (ierr /= 0) then
            write(*,'(A,A)') 'ERROR OPENING INPUT FILE ',datafile
            STOP
         end if

         write(fmt,'("(",I6,"D24.15)")') nsteps-1
!        FOURIER TRANSFORM OF THE AUTOCORRELATION FUNCTION OF EACH COORDINATE
!        IN TURN...
         cumul=0d0
         do icrd=1,ncoords
            write(*,*) 'ANALYZING COORINATE No ', icrd
!           READ COORDINATE FROM FILE
            read(1,*) (vel(i),i=1,nsteps)
!           VELOCITIES WITH FOREWARD DIFFERENCES
            do i=1,nsteps-1
              vel(i)=vel(i)-vel(i+1)
            end do
            vel(nsteps)=0d0

!           @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!           DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG
!           @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!           CREATES FAKE DATA COMPOSED OF A SUM OF THREE SINE FUNCTIONS.
            !
            ! PI = 3.14159d0
            ! nu1= 0.5d0
            ! nu2= 1.0d0
            ! nu3= 2.22d0
            !
            ! do k=0,nsteps-1
            !    tt= real(k,8) * dt
            !    CALL RANDOM_NUMBER(noise)
            !    vel(k) = sin( 2d0 * PI * nu1 * tt ) + &
            !    &        sin( 2d0 * PI * nu2 * tt ) + &
            !    &        sin( 2d0 * PI * nu3 * tt ) !+ 10d0 * noise
            !    write(89,*) vel(k)
            ! end do
            !
!           @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!           DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG
!           @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

            if (icrd==1) then
               call dfftw_plan_dft_r2c_1d(plan,nsteps,vel,out,"FFTW_ESTIMATE")
            endif
            call dfftw_execute_dft_r2c(plan, vel, out)
            out = conjg(out)*out
            write(100,fmt) REAL(out)
            cumul = cumul + REAL(out)
         end do

!        Writing sum spectra to disk.
         open(unit=1000,file='sum_spectra.dat',iostat=ierr)
         if (ierr /= 0) then
            write(*,'(A,A)') 'ERROR OPENING INPUT FILE ',datafile
            STOP
         end if
         write(1000,fmt) cumul
         close(unit=1000,iostat=ierr)

!        Destroying fftw plan
         call dfftw_destroy_plan(plan)

!        CLOSING DATA AND COORD-SPECTRA FILE
         close(unit=1,iostat=ierr)
         close(unit=100,iostat=ierr)

!        CREATING A FREQUENCY AXIS AND WRITING IT TO FILE
         open(unit=1,file='freq_axis.dat',iostat=ierr)
         if (ierr /= 0) then
            write(*,'(A,A)') 'ERROR OPENING INPUT FILE ','freq_axis.dat'
            STOP
         end if

         sample_rate=1/dt
         sample_rate=sample_rate/1d12 ! dt in ps, freq in THz
         nu_increment = sample_rate/dble(nsteps)

         nu = 0d0
         do i = 1,nsteps/2+1
            write(1,'(D17.8)') nu
            nu = nu + nu_increment
         end do
         close(unit=1,iostat=ierr)

      end subroutine



      subroutine vac_autocor(datafile,acorr,nsteps,ncoords,dt,damp,dodamp)
!        CALCULATES THE AUTOCORRELATION FUNCTION FOR VACUUM SPECTRA.

         implicit none

!        --------------------------------------------------------------
         integer,intent(in)                 :: nsteps
         integer,intent(in)                 :: ncoords
         integer,intent(in)                 :: dodamp
         character(90),intent(in)           :: datafile
         double precision,intent(in)        :: dt
         double precision,intent(in)        :: damp
         double precision,intent(out)       :: acorr(nsteps-1)

         character(20)      :: fmt
         integer  :: i,j,k,ierr,t,tau,icrd,maxlag
         double precision   :: vel(nsteps)
         double precision   :: corr(nsteps-1)
         double precision   :: noise
         double precision   :: nu1,nu2,nu3,per1,per2,per3,tt,PI,dp,time
!        --------------------------------------------------------------

!        OPENING INPUT FILE
         open(unit=1,file=datafile,iostat=ierr)
         if (ierr /= 0) then
            write(*,'(A,A)') 'ERROR OPENING INPUT FILE ',datafile
            STOP
         end if

         open(unit=10,file='autocorr.dat',iostat=ierr)
         if (ierr /= 0) then
            write(*,'(A,A)') 'ERROR OPENING INPUT FILE ',datafile
            STOP
         end if

         write(fmt,'("(",I6,"D24.15)")') nsteps-1
!        FIRST VALUE OF AUTOCORRELATION FUNCTION
         maxlag=nsteps-2

         read(1,*) (vel(i),i=1,nsteps)
         !COMPUTE VELOCITIES WITH FOREWARD DIFFERENCES
         do i=1,nsteps-1
            vel(i)=vel(i)-vel(i+1)
         end do

!        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!        DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG
!        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         ! PI = 3.14159d0
         ! nu1= 0.5d0
         ! nu2= 1.0d0
         ! nu3= 2.22d0
         ! per1=1/nu1
         ! per2=1/nu2
         ! per3=1/nu3
         !
         ! do k=0,nsteps-1
         !    tt= real(k,8) * dt
         !    CALL RANDOM_NUMBER(noise)
         !    vel(k) = sin( 2d0 * PI * nu1 * tt ) + &
         !    &        sin( 2d0 * PI * nu2 * tt ) + &
         !    &        sin( 2d0 * PI * nu3 * tt ) !+ 10d0 * noise
         !    write(89,*) vel(k)
         ! end do
!        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!        DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG
!        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

         acorr=0d0
         corr=0d0
         do icrd=1,ncoords
            do tau=1,maxlag
               do t=1,nsteps-tau-1
                  corr(tau) = corr(tau) + vel(t) * vel(t+tau)
               end do
               dp = 1d0
               time = dt*tau
               if (dodamp == 1) dp = dampfunc(damp,time)
               corr(tau) = corr(tau) * dp
            end do
            acorr = corr/(nsteps-1)
            write(10,fmt) acorr
         end do

!        CLOSING DATA FILE
         close(unit=1,iostat=ierr)
         close(unit=10,iostat=ierr)

      end subroutine

      subroutine run_avg(series,slen,navg,avgser)

!        CALCULATES A RUNNING AVERAGE OF A DATA SERIES.
         implicit none

!        --------------------------------------------------------------
         integer,intent(in)     :: slen         ! Length of the data series.
         integer,intent(in)     :: navg         ! # of steps to average over.
         double precision,intent(in)      :: series(slen) ! Data series.
         double precision,intent(out)     :: avgser(slen-navg+1)

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
         double precision,intent(in)         :: dt
         double precision,intent(in)         :: damp
         double precision,intent(in)         :: series(slen)
         double precision,intent(out)        :: rspec(slen/2+1)
         double precision,intent(out)        :: ispec(slen/2+1)

         integer             :: i,j
         double precision              :: ddt
         double precision              :: dmp
         double precision              :: dp
         double precision              :: nu
         double precision              :: ftr
         double precision              :: fti
         double precision              :: t
         double precision              :: sample_rate
         double precision              :: nu_increment
         double precision,parameter    :: cc = 2.99792458d-5 ! speed of light cm/fs
         double precision,parameter    :: pi = 3.1415926535897932384626433832795d0
!        --------------------------------------------------------------

         ddt = dt
         dmp = damp
         sample_rate = 1/ddt
         nu_increment = sample_rate/dble(slen)

         nu = 0d0
         do i = 1,slen/2+1

            ftr = 0d0
            fti = 0d0

            t=0d0
            do j=1,slen-1
               t = t + ddt

               ftr = ftr + cos(2*pi*t*nu) * series(j)
               fti = fti + sin(2*pi*t*nu) * series(j)
            end do

            rspec(i) = ftr
            !rspec(i) = ftr**2 + fti**2
!            rspec(i) = Sqrt(ftr**2 + fti**2)
!            rspec(i) = Abs(ftr) + Abs(fti)
            if (prefac == 1) rspec(i) = qmfac(nu,300d0) * rspec(i)
            ispec(i) = fti
            nu = nu + nu_increment

         end do

      end subroutine

      double precision function dampfunc(damp,t)
         double precision  :: damp
         double precision  :: t
         dampfunc = exp(-t/damp)
         return
      end function

      double precision function qmfac(nu,temp)
         double precision            :: nu
         double precision            :: temp
         double precision            :: beta
         double precision            :: num
         double precision            :: den
         double precision,parameter  :: kb = 0.69503476
         double precision,parameter  :: pi = 3.1415926535897932384626433832795d0

         beta = 1d0/(kb*temp)
         num = (beta*nu)**2
         den = pi*(1d0 - exp(-beta*nu))

         qmfac = num/den
         return
      end function

      subroutine print_spectra(timestep,slen,rspec,ispec,spectrafile)
         implicit none

!        --------------------------------------------------------------
         integer,intent(in)         :: slen
         double precision,intent(in)          :: timestep
         double precision,intent(in)          :: rspec(slen)
         double precision,intent(in)          :: ispec(slen)
         character(99),intent(in)   :: spectrafile

         integer  :: ierr,i
         double precision   :: nu,sample_rate,nu_increment
!        --------------------------------------------------------------

         open(unit=3,file=spectrafile,iostat=ierr)
         if (ierr /= 0) then
            write(*,'(A,A)') 'ERROR OPENNING FILE ',spectrafile
            STOP
         end if

         sample_rate=1/timestep
         nu_increment = sample_rate/dble(slen)

         nu = 0d0
         do i = 1,slen/2+1
            write(3,'(3D17.8)') nu, rspec(i), ispec(i)
            nu = nu + nu_increment
         end do

         close(unit=3,iostat=ierr)
         if (ierr /= 0) then
            write(*,'(A,A)') 'ERROR CLOSING FILE ',spectrafile
            STOP
         end if

      end subroutine

end module linearSpectroscopy_mod
