program aka2
implicit none
character(len=50) :: infile
integer :: nsteps, spatial, natoms
double precision :: coordinate, dt
double precision, allocatable, dimension (:) :: coordinate_evol
integer, dimension (3) :: point

call get_command_argument(1,infile)
call getdims(infile,nsteps,spatial,natoms)
nsteps=nsteps

allocate(coordinate_evol(nsteps))
dt=1.0d-15
call pump_probe(infile,nsteps,natoms,spatial,point,coordinate_evol,dt)
   
contains

SUBROUTINE getdims(infile,nsteps,spatial,natoms)
USE netcdf
IMPLICIT NONE
INTEGER(KIND=4), INTENT(OUT) :: nsteps,spatial,natoms
INTEGER(KIND=4) :: ncid
CHARACTER(LEN=50), INTENT(IN) :: infile
CHARACTER(LEN=50) :: xname, yname, zname

CALL check(nf90_open(infile, nf90_nowrite, ncid))
CALL check(nf90_inquire_dimension(ncid,1,xname,nsteps))
CALL check(nf90_inquire_dimension(ncid,2,yname,spatial))
CALL check(nf90_inquire_dimension(ncid,3,zname,natoms))
CALL check(nf90_close(ncid))
END SUBROUTINE getdims

SUBROUTINE pump_probe(infile,nsteps,natoms,spatial,point,coordinate_evol,dt)
USE netcdf
IMPLICIT NONE
DOUBLE PRECISION :: coordinate,dt
DOUBLE PRECISION :: sample_rate,nu_increment,nu
DOUBLE PRECISION, DIMENSION(NSTEPS) :: coordinate_evol
!DOUBLE PRECISION, DIMENSION(NSTEPS/2+1) :: veltraj_window
DOUBLE PRECISION, DIMENSION(NSTEPS/2) :: veltraj_window
COMPLEX*16      , DIMENSION((NSTEPS/2+1)/2+1) :: out
DOUBLE PRECISION, DIMENSION((NSTEPS/2+1)/2+1) :: cumul
DOUBLE PRECISION, DIMENSION((NSTEPS/2+1)/2+1) :: tdspec

INTEGER(KIND=4), INTENT(IN) :: natoms, nsteps, spatial
INTEGER(KIND=4) :: ncid, xtype, ndims, varid
CHARACTER(LEN=50), INTENT(IN) :: infile
CHARACTER(LEN=50) :: xname, vname
CHARACTER(LEN=20) :: fmt
integer, dimension(3) :: point,endp
integer :: i,j,k,time0,ierr,specdim,specframes
integer*8 :: plan

specframes=int(nsteps/2)
specdim=int(specframes/2+1)
open(1, file = 'tseries.dat', status = 'unknown')  
CALL check(nf90_open(infile, nf90_nowrite, ncid))
write(*,*) natoms
tdspec=0d0
cumul=0d0
veltraj_window=0d0
do i = 1,natoms !natoms
  write(*,*) 'STARTING ATOM ',i,' OUT OF ',natoms
  do j = 1,3 !3
      point = (/ j, i, 1 /)
      endp = (/ 1,1,nsteps-1 /)
      CALL check(nf90_get_var(ncid,4,coordinate_evol,start = point,count = endp))
          write(*,*) 'after geting coords'
      do time0=1,1  ! DEBUG
       veltraj_window=coordinate_evol(time0:time0+nsteps/2-1)
       if (i==1 .AND. j==1) then
          call dfftw_plan_dft_r2c_1d(plan,specframes,veltraj_window,out,"FFTW_ESTIMATE")
          write(*,*) 'in plan creation'
       endif
       call dfftw_execute_dft_r2c(plan, veltraj_window, out)
          write(*,*) 'after fourier'
       out = conjg(out)*out
       cumul(:) = REAL(out) ! DEBUG
      end do
      tdspec=tdspec+cumul
  enddo
enddo
close(1)
CALL check(nf90_close(ncid))
write(*,*) tdspec


         write(fmt,'("(",I6,"D24.15)")') nsteps/2+1
         write(*,*) fmt
!        Writing sum spectra to disk.
         open(unit=1000,file='sum_spectra.dat',iostat=ierr)
         if (ierr /= 0) then
            write(*,'(A,A)') 'ERROR OPENING INPUT FILE ',infile
            STOP
         end if
        ! do time0=1,nsteps/2-2
            !write(1000,fmt) tdspec(:,time0)
            write(1000,fmt) tdspec(:) ! DEBUG
        ! end do
         close(unit=1000,iostat=ierr)

!        Destroying fftw plan
         call dfftw_destroy_plan(plan)

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
         close(1)



END SUBROUTINE pump_probe

SUBROUTINE check(istatus)
USE netcdf
IMPLICIT NONE
INTEGER, INTENT (IN) :: istatus
IF (istatus /= nf90_noerr) THEN
write(*,*) TRIM(ADJUSTL(nf90_strerror(istatus)))
END IF
END SUBROUTINE check


end program aka2

