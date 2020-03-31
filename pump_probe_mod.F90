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
         double precision,intent(in)         :: dt

         character(20)      :: fmt,fmtin
         character(80)      :: type
         integer  :: i,j,k,ierr,t,tau,icrd,maxlag,iatm,xyz,time0,i_frame
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
         integer                                       :: cnt(3)
         integer                                       :: strt(3)
         integer                                       :: err
         integer                                       :: specframes
         integer                                       :: specdim
         character(len=80)                             :: title
         integer                                       :: natom
         integer                                       :: nframes
!         double precision,allocatable, dimension(:,:)  :: Coords, Velo
         double precision,allocatable, dimension(:)    :: veltraj
         double precision,allocatable, dimension(:)    :: veltraj_window
         double precision,allocatable, dimension(:)    :: atmass
         !double precision,allocatable, dimension(:,:)  :: tdspec
         !double precision,allocatable, dimension(:,:)  :: cumul
         double precision,allocatable, dimension(:)  :: tdspec  ! DEBUG
         double precision,allocatable, dimension(:)  :: cumul  ! DEBUG
         double precision                              :: Time
         double precision, dimension(3)                :: box
         double precision                              :: alpha, beta, gamma
         logical                                       :: box_found, velocities_found

         integer :: ncid,spatial
!        --------------------------------------------------------------
         integer*8 :: plan
         !type(C_PTR) :: plan
!        --------------------------------------------------------------


         call getdims(datafile,nframes,spatial,natom,coordVID,velocityVID)
         write(*,*) nframes, spatial, natom, velocityVID, coordVID

         allocate(veltraj(nframes))

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
         !allocate( tdspec(specdim,specframes), cumul(specdim,specframes)  )
         allocate( tdspec(specdim), cumul(specdim)  )
         allocate( veltraj_window(specframes))
         cumul=0d0
         tdspec=0d0



         CALL check(nf90_open(datafile, nf90_nowrite, ncid)) 
         do iatm=1,natom
         do xyz=1,3
            write(*,*) ' xyz ',xyz,' iatm ', iatm
            strt=(/ xyz, iatm, 1 /)
            cnt=(/ 1, 1, nframes /)
            write(*,*)  'start = ', strt
            write(*,*)  'count = ',cnt
            write(*,*)  'ncid = ', ncid
            ! Get velocities
            call check(nf90_get_var(ncid,4,veltraj,start = strt,count = cnt)) 
               write(*,*) 'CACACACACC'
            veltraj=atmass(iatm)*veltraj

            !do time0=1,nframes/2-2
            do time0=1,1  ! DEBUG
               veltraj_window=veltraj(time0:time0+specframes)
               if (iatm==1 .AND. xyz==1) then
                  call dfftw_plan_dft_r2c_1d(plan,specframes,veltraj_window,out,"FFTW_ESTIMATE")
               endif
               call dfftw_execute_dft_r2c(plan, veltraj_window, out)
               out = conjg(out)*out
               !cumul(:,time0) = REAL(out)
               cumul(:) = REAL(out) ! DEBUG
            end do
            tdspec = tdspec + cumul
         end do
         end do
         CALL check(nf90_close( ncid ))


         write(fmt,'("(",I6,"D24.15)")') specframes
         write(*,*) fmt
!        Writing sum spectra to disk.
         open(unit=1000,file='sum_spectra.dat',iostat=ierr)
         if (ierr /= 0) then
            write(*,'(A,A)') 'ERROR OPENING INPUT FILE ',datafile
            STOP
         end if
         do time0=1,nframes/2-2
            !write(1000,fmt) tdspec(:,time0)
            write(1000,fmt) tdspec(:) ! DEBUG
         end do
         close(unit=1000,iostat=ierr)

!        Destroying fftw plan
         call dfftw_destroy_plan(plan)

!        CLOSING DATA AND COORD-SPECTRA FILE
         close(unit=1,iostat=ierr)

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

      SUBROUTINE getdims(infile,nsteps,spatial,natoms,coordVID,velVID)
      USE netcdf
      IMPLICIT NONE
      INTEGER(KIND=4), INTENT(OUT) :: nsteps,spatial,natoms,coordVID,velVID
      INTEGER(KIND=4) :: ncid
      CHARACTER(LEN=50), INTENT(IN) :: infile
      CHARACTER(LEN=50) :: xname, yname, zname,varnme
      
      CALL check(nf90_open(infile, nf90_nowrite, ncid))
      CALL check(nf90_inquire_dimension(ncid,1,xname,nsteps))
      CALL check(nf90_inquire_dimension(ncid,2,yname,spatial))
      CALL check(nf90_inquire_dimension(ncid,3,zname,natoms))
      call check(nf90_inq_varid(ncid,"coordinates",coordVID))
      call check(nf90_inq_varid(ncid,"velocities",velVID))
      CALL check(nf90_close(ncid))

      END SUBROUTINE getdims
      
      SUBROUTINE getandwritecoord(infile,nsteps,natoms,spatial,point,coordinate_evol)
      USE netcdf
      IMPLICIT NONE
      REAL (KIND=4) :: coordinate
      REAL (KIND=4), DIMENSION(NSTEPS) :: coordinate_evol
      INTEGER(KIND=4), INTENT(IN) :: natoms, nsteps, spatial
      !INTEGER(KIND=4), DIMENSION(2) :: dimids
      INTEGER(KIND=4) :: ncid, xtype, ndims, varid
      CHARACTER(LEN=50), INTENT(IN) :: infile
      CHARACTER(LEN=50) :: xname, vname
      integer, dimension(3) :: point
      integer :: i,j,k
      
      open(1, file = 'tseries.dat', status = 'unknown')  
      CALL check(nf90_open(infile, nf90_nowrite, ncid))
      do i = 1,natoms !natoms
        do j = 1,3 !3
          do k = 1,nsteps
            point = (/ j, i, k /)
            CALL check(nf90_get_var(ncid,3,coordinate,start = point))
            coordinate_evol(k)=coordinate
            enddo
          write(1,*) coordinate_evol(1:nsteps)
        enddo
          write(1,*) coordinate_evol(1:nsteps)
      enddo
      close(1)
      CALL check(nf90_close(ncid))
      END SUBROUTINE getandwritecoord
      
      SUBROUTINE check(istatus)
      USE netcdf
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: istatus
      IF (istatus /= nf90_noerr) THEN
      write(*,*) TRIM(ADJUSTL(nf90_strerror(istatus)))
      END IF
      END SUBROUTINE check
      
      
      
      end module pump_probe_mod
