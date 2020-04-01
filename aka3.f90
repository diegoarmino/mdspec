program aka2
implicit none
character(len=50) :: infile,topfile
integer :: nsteps, spatial, natoms
double precision :: coordinate, dt
double precision, allocatable, dimension (:,:) :: tdspec
double precision, allocatable, dimension (:)   :: atmass
integer                :: specframes
integer                :: specdim
INTEGER*8              :: plan

call get_command_argument(1,infile)
call get_command_argument(2,topfile)
call getdims(infile,nsteps,spatial,natoms)

specframes=int(nsteps/2)
specdim=int(specframes/2+1)

allocate( tdspec((NSTEPS/2+1)/2+1,NSTEPS/2), atmass(natoms) )

call get_at_mass_prmtop(natoms,topfile,atmass)

dt=1.0d-15
call make_fftw_plan(infile,nsteps,plan)
call td_spectrum(infile,atmass,plan,nsteps,natoms,spatial,dt,tdspec)
call print_results(nsteps,dt,tdspec)
   
!Destroying fftw plan
call dfftw_destroy_plan(plan)

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


SUBROUTINE get_at_mass_prmtop(natom,topfile,atmass)
!  READ ATOMIC MASSES FROM AMBER PRMTOP FILE.
   use nextprmtop_section_mod
   IMPLICIT NONE
   CHARACTER(LEN=50),                INTENT(IN)   :: topfile
   INTEGER,                          INTENT(IN)   :: natom
   DOUBLE PRECISION,DIMENSION(natom),INTENT(OUT)  :: atmass

   character(len=20)                                :: fmt, fmtin, type
   integer                                          :: err,i

!  OPEN PRMTOP FILE
   open(unit=101,file=topfile,iostat=err)
   if (err /= 0) then
      write(*,'(A,A)') 'ERROR OPENING INPUT FILE prmtop'
      STOP
   end if

   call nxtsec_reset()
   fmtin = '(5E16.8)'
   type = 'MASS'
   call nxtsec(101, 6, 0, fmtin, type, fmt, err)

   read(101, fmt) (atmass(i), i = 1, natom)
   close(unit=101)
   atmass=sqrt(atmass)

END SUBROUTINE get_at_mass_prmtop



SUBROUTINE make_fftw_plan(infile,nsteps,plan)
USE netcdf
USE,INTRINSIC :: iso_c_binding
IMPLICIT NONE
include 'fftw3.f03'
CHARACTER(LEN=50), INTENT(IN)  :: infile
INTEGER(KIND=4),   INTENT(IN)  :: nsteps
INTEGER*8,         INTENT(OUT) :: plan

INTEGER,          DIMENSION(3)                 :: point,endp
DOUBLE PRECISION, DIMENSION(NSTEPS)            :: coordinate_evol
DOUBLE PRECISION, DIMENSION(NSTEPS/2)          :: veltraj_window
COMPLEX*16      , DIMENSION((NSTEPS/2+1)/2+1)  :: out

INTEGER :: ncid,i,j,k,time0,specdim,specframes

write(*,*) 'MAKING FFTW PLAN...'

specframes=int(nsteps/2)
specdim=int(specframes/2+1)

! READ FIRST LINE OF DATASET AS A REPRESENTATIVE.
write(*,*) 'GETTING TRAJECTORY FROM NETCDF FILE'
CALL check(nf90_open(infile, nf90_nowrite, ncid))
point = (/ 1, 1, 1 /)
endp = (/ 1,1,nsteps-1 /)
CALL check(nf90_get_var(ncid,4,coordinate_evol,start = point, count = endp))
CALL check(nf90_close(ncid))
veltraj_window=coordinate_evol(1:int(nsteps/2))

! MEASURE THE BEST PLAN FOR FOURIER TRANSFORM. THIS IS MORE TIME-CONSUMING
! THAN THE ESTIMATE OPTION, BUT GIVEN THE AMOUNT OF CALLS TO FFTW IT IS 
! WORTH THE EFFORT.
write(*,*) 'FINDING BEST PLAN'
CALL dfftw_plan_dft_r2c_1d(plan,specframes,veltraj_window,out,"FFTW_MEASURE")

write(*,*) 'FINISHED FFTW PLAN!'

end subroutine make_fftw_plan

SUBROUTINE td_spectrum(infile,atmass,plan,nsteps,natoms,spatial,dt,tdspec)
USE netcdf
USE,INTRINSIC :: iso_c_binding
IMPLICIT NONE
include 'fftw3.f03'
CHARACTER(LEN=50), INTENT(IN) :: infile
INTEGER(KIND=4),   INTENT(IN) :: natoms, nsteps, spatial
INTEGER*8,         INTENT(IN) :: plan
DOUBLE PRECISION,  INTENT(IN) :: atmass(natoms)
DOUBLE PRECISION :: coordinate,dt
DOUBLE PRECISION, DIMENSION(NSTEPS) :: coordinate_evol
DOUBLE PRECISION, DIMENSION(NSTEPS/2) :: veltraj_window
COMPLEX*16      , DIMENSION((NSTEPS/2+1)/2+1) :: out
DOUBLE PRECISION, DIMENSION((NSTEPS/2+1)/2+1,NSTEPS/2),intent(out) :: tdspec

INTEGER(KIND=4) :: ncid, xtype, ndims, varid
CHARACTER(LEN=50) :: xname, vname
CHARACTER(LEN=20) :: fmt
integer, dimension(3) :: point,endp
integer :: i,j,k,time0,ierr,specdim,specframes

specframes=int(nsteps/2)
specdim=int(specframes/2+1)

CALL check(nf90_open(infile, nf90_nowrite, ncid))

tdspec=0d0
veltraj_window=0d0

do i = 1,natoms !natoms
  write(*,*) 'STARTING ATOM ',i,' OUT OF ',natoms
  do j = 1,3 !3
      point = (/ j, i, 1 /)
      endp = (/ 1,1,nsteps-1 /)
      CALL check(nf90_get_var(ncid,4,coordinate_evol,start = point,count = endp))
      do time0=1,specframes-1
         veltraj_window=coordinate_evol(time0:time0+nsteps/2-1)
         call dfftw_execute_dft_r2c(plan, veltraj_window, out)
         out = conjg(out)*out
         tdspec(:,time0) = tdspec(:,time0) + REAL(out) ! DEBUG
      end do
  end do
end do
close(1)
CALL check(nf90_close(ncid))

end subroutine td_spectrum

subroutine print_results(nsteps,dt,tdspec)
implicit none
INTEGER,                        INTENT(IN) :: nsteps
DOUBLE PRECISION,               INTENT(IN) :: dt
DOUBLE PRECISION,DIMENSION(:,:),INTENT(IN) :: tdspec

character(len=20)                          :: fmt
double precision                           :: nu
double precision                           :: nu_increment
double precision                           :: sample_rate
integer                                    :: ierr
integer                                    :: specframes
integer                                    :: specdim
integer                                    :: i,time0

specframes=int(nsteps/2)
specdim=int(specframes/2+1)

write(fmt,'("(",I6,"D24.15)")') nsteps/2+1
write(*,*) fmt
!Writing sum spectra to disk.
open(unit=1000,file='sum_spectra.dat',iostat=ierr)
if (ierr /= 0) then
   write(*,'(A,A)') 'ERROR OPENING INPUT FILE sum_spectra.dat'
   STOP
end if
do time0=1,specframes-1
   write(1000,fmt) tdspec(:,time0)
end do
close(unit=1000,iostat=ierr)


!CREATING A FREQUENCY AXIS AND WRITING IT TO FILE
open(unit=1,file='freq_axis.dat',iostat=ierr)
if (ierr /= 0) then
   write(*,'(A,A)') 'ERROR OPENING INPUT FILE ','freq_axis.dat'
   STOP
end if

sample_rate=1/dt
sample_rate=sample_rate/1d12 ! dt in ps, freq in THz
nu_increment = sample_rate/dble(specframes)

nu = 0d0
do i = 1,specdim
   write(1,'(D17.8)') nu
   nu = nu + nu_increment
end do
close(1)

end subroutine print_results


SUBROUTINE check(istatus)
USE netcdf
IMPLICIT NONE
INTEGER, INTENT (IN) :: istatus
IF (istatus /= nf90_noerr) THEN
write(*,*) TRIM(ADJUSTL(nf90_strerror(istatus)))
END IF
END SUBROUTINE check


end program aka2

