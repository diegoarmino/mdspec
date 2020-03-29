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
         use nextprmtop_section_mod
         use mdspecNetcdf_mod, only: NC_openRead,NC_setupCoordsVelo,NC_close, &
                                     NC_error,NC_readRestartBox,NC_setupRestart,&
                                     NC_setupMdcrd,checkNCerror
         use, intrinsic :: iso_c_binding
         implicit none

         include 'fftw3.f03'
!        --------------------------------------------------------------
         integer,intent(in)       :: nsteps
         integer,intent(in)       :: ncoords
         character(90),intent(in) :: datafile
         double precision         :: spectra(nsteps/2+1)
         double precision         :: dt

         character(20)      :: fmt,fmtin
         character(80)      :: type
         integer  :: i,j,k,ierr,t,tau,icrd,maxlag,iatm,xyz,time0
         complex*16         :: out(nsteps/2+1)
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
         double precision,allocatable, dimension(:)    :: veltraj
         double precision,allocatable, dimension(:)    :: veltraj_window
         double precision,allocatable, dimension(:)    :: atmass
         double precision,allocatable, dimension(:,:)  :: tdspec
         double precision,allocatable, dimension(:,:)  :: cumul
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

         write(*,*) nframes
         allocate(Coords(3, natom), Velo(3, natom), veltraj(nframes))

!         ! Get coords
!         if (NC_error(nf90_get_var(ncid, coordVID, Coords(1:3,1:natom), &
!                            start = (/ 1, 1, 499 /), count = (/ 3, natom, 1 /)),&
!               'reading trajectory coordinates')) stop
!         write(*,*) coords
!         stop

         ! Get velocities
!         if (velocityVID.ne.-1) then
!           velocities_found=.true.
!           if (NC_error(nf90_get_var(ncid, velocityVID, Velo(1:3,1:natom,1), &
!                                     start = (/ 1, 1, 500 /), count = (/ 3, natom, 1 /)),&
!                        "read_nc_restart(): Getting velocities")) stop
!         endif
!         write(*,*) velo
!         stop

!         ! Get box information
!         if (NC_readRestartBox(ncid,box(1),box(2),box(3),alpha,beta,gamma)) then
!           box_found = .false.
!         else
!           box_found = .true.
!         endif

!        OPEN FILE HOLDING INDIVIUAL SPECTRA (ONE FOR EACH COORDINATE)
!         open(unit=100,file='spectra.dat',iostat=ierr)
!         if (ierr /= 0) then
!            write(*,'(A,A)') 'ERROR OPENING INPUT FILE ',datafile
!            STOP
!         end if

         write(fmt,'("(",I6,"D24.15)")') nsteps-1
!        READ PRMTOP FILE AND READ ATOMIC MASSES
!        OPEN PRMTOP FILE
         allocate(atmass(natom))
         open(unit=101,file='top',iostat=ierr)
         if (ierr /= 0) then
            write(*,'(A,A)') 'ERROR OPENING INPUT FILE ',datafile
            STOP
         end if

         call nxtsec_reset()
         fmtin = '(5E16.8)'
         type = 'MASS'
         call nxtsec(101, 6, 0, fmtin, type, fmt, err)
       
         read(101, fmt) (atmass(i), i = 1, natom)
         close(unit=101)
         atmass=sqrt(atmass)


         err = nf90_close( ncid )
         write(*,*) err
!        FOURIER TRANSFORM OF THE AUTOCORRELATION FUNCTION OF EACH COORDINATE
!        IN TURN...

         specframes=nframes/2+1
         specdim=specframes/2+1
         allocate( tdspec(specdim,nframes/2+1), cumul(specdim,nframes/2+1)  )
         allocate( veltraj_window(specframes))
         cumul=0d0
         tdspec=0d0

         if (NC_error( nf90_open( datafile, NF90_NOWRITE, ncid ))) stop
         do iatm=1,natom
         do xyz=1,3
            write(*,*) ' xyz ',xyz,' iatm ', iatm
            start=(/ xyz, natom, 1 /)
            count=(/ 1, 1, nframes /)
            ! Get velocities
            if (NC_error(nf90_get_var(ncid, velocityVID, veltraj(1:nframes), &
                           start = start, count = count),&
                           "pump_probe(): Getting velocities")) stop

            veltraj=atmass(iatm)*veltraj

            do time0=1,nframes/2-2
               veltraj_window=veltraj(time0:time0+specframes)
               if (iatm==1 .AND. xyz==1) then
                  call dfftw_plan_dft_r2c_1d(plan,specframes,veltraj_window,out,"FFTW_ESTIMATE")
               endif
               call dfftw_execute_dft_r2c(plan, veltraj_window, out)
               out = conjg(out)*out
               cumul(:,time0) = REAL(out)
            end do
            err = nf90_close( ncid )
         end do
         tdspec = tdspec + cumul
         end do

         write(fmt,'("(",I6,"D24.15)")') specframes
         write(*,*) fmt
!        Writing sum spectra to disk.
         open(unit=1000,file='sum_spectra.dat',iostat=ierr)
         if (ierr /= 0) then
            write(*,'(A,A)') 'ERROR OPENING INPUT FILE ',datafile
            STOP
         end if
         do time0=1,nframes/2-2
            write(1000,fmt) tdspec(:,time0)
         end do
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
