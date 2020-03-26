module pump_probe_mod

   private
   public :: pump_probe


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


      subroutine pump_probe(datafile,nsteps,ncoords,dt)
!        CALCULATES THE TIME DEPENDENT POWER SPECTRA
!        AS THE SQUARE OF THE FOURIER TRANSFORM OF THE 
!        VELOCITY TRAJECTORY AT DIFFERENT TIME DELAYS FROM
!        THE FIRST FRAME.

         use netcdf
         use mdspecNetcdf_mod, only: NC_openRead,NC_setupCoordsVelo,NC_close, &
                                     NC_error,NC_readRestartBox,NC_setupRestart,&
                                     NC_setupMdcrd
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
         integer  :: i,j,k,ierr,t,tau,icrd,maxlag,iatm,xyz
         complex*16         :: out(nsteps/2+1)
         double precision   :: cumul(nsteps/2+1)
         double precision   :: noise
         double precision   :: sample_rate
         double precision   :: nu_increment
         double precision   :: nu
         double precision   :: vel(nsteps)
         double precision   :: corr(nsteps-1)
         double precision   :: PI,nu1,nu2,nu3,tt

         integer            :: coordVID, velocityVID, cellAngleVID, cellLengthVID, TempVID
         integer            :: timeVID, repidx_var_id, crdidx_var_id, remd_values_var_id
         integer                                       :: count(3)
         integer                                       :: start(3)
         integer                                       :: err
         integer                                       :: specframes
         integer                                       :: specdim
         character(len=80)                             :: title
         integer                                       :: natom
         integer                                       :: nframes
         double precision,allocatable, dimension(:,:)  :: Coords, Velo
         double precision,allocatable, dimension(:)    :: Veltraj
         double precision,allocatable, dimension(:,:)  :: tdspec
         double precision                              :: Time
         double precision, dimension(3)                :: box
         double precision                              :: alpha, beta, gamma
         logical                                       :: box_found, velocities_found

         integer :: ncid
!        --------------------------------------------------------------
         integer*8 :: plan
         !type(C_PTR) :: plan
!        --------------------------------------------------------------

!        OPENING INPUT FILE
         ! ---=== Open file
         write(*,*) ncid
         if (NC_openRead(datafile, ncid)) then
           write(*,'(a)') "read_nc_restart(): Could not open coordinate file."
           stop
         endif
         write(*,*) ncid
         ! Get number of atoms. err is used as a dummy arg for coordVID and velocityVID
         if (NC_setupCoordsVelo(ncid, natom, err, err)) stop
         ! Close file
         !call NC_close(ncid)
         write(*,*) natom


         ! Setup Trajectory
         if (NC_setupMdcrd(ncid, title, nframes, natom, &
                               coordVID, velocityVID, timeVID, &
                               cellLengthVID, cellAngleVID, TempVID)) stop

         allocate(Coords(3, natom), Velo(3, natom), veltraj(nframes))

!         ! Get coords
!         if (NC_error(nf90_get_var(ncid, coordVID, Coords(1:3,1:natom), &
!                            start = (/ 1, 1, 1 /), count = (/ 1, 3, natom /)),&
!               'reading restart coordinates')) stop

         ! Get velocities
         if (velocityVID.ne.-1) then
           velocities_found=.true.
           if (NC_error(nf90_get_var(ncid, velocityVID, Velo(1:3,1:natom), &
                                     start = (/ 1, 1, 1 /), count = (/ 1, 3, natom /)),&
                        "read_nc_restart(): Getting velocities")) stop
         endif

!         ! Get box information
!         if (NC_readRestartBox(ncid,box(1),box(2),box(3),alpha,beta,gamma)) then
!           box_found = .false.
!         else
!           box_found = .true.
!         endif

!        OPEN FILE HOLDING INDIVIUAL SPECTRA (ONE FOR EACH COORDINATE)
         open(unit=100,file='spectra.dat',iostat=ierr)
         if (ierr /= 0) then
            write(*,'(A,A)') 'ERROR OPENING INPUT FILE ',datafile
            STOP
         end if

         write(fmt,'("(",I6,"D24.15)")') nsteps-1
!        FOURIER TRANSFORM OF THE AUTOCORRELATION FUNCTION OF EACH COORDINATE
!        IN TURN...
         cumul=0d0

         specframes=nframes/2+1
         specdim=specframes/2+1
         allocate( tdspec(specdim,nframes/2+1) )
         tdspec=0d0

         do time0=1,nframes/2-2
         do iatm=1,natom
         do xyz=1,3
            start=(/ time0, xyz, natom /)
            count=(/ time0+specframes, xyz, natom /)
            write(*,*) 'ANALYZING ATOM No ', iatm
            ! Get velocities
            if (velocityVID.ne.-1) then
              velocities_found=.true.
              if (NC_error(nf90_get_var(ncid, velocityVID, Veltraj(1:nframes), &
                           start = start, count = count),&
                           "read_nc_restart(): Getting velocities")) stop
            endif

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

            if (iatm==1 .AND. xyz==1) then
               call dfftw_plan_dft_r2c_1d(plan,nsteps,veltraj,out,"FFTW_ESTIMATE")
            endif
            call dfftw_execute_dft_r2c(plan, veltraj, out)
            out = conjg(out)*out
            write(100,fmt) REAL(out)
            cumul = cumul + REAL(out)
         end do
         end do
         tdspec(:,time0) = cumul
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
         nu_increment = sample_rate/dble(specframes)

         nu = 0d0
         do i = 1,nsteps/2+1
            write(1,'(D17.8)') nu
            nu = nu + nu_increment
         end do
         close(unit=1,iostat=ierr)

      end subroutine






end module pump_probe_mod
