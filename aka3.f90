program aka2
implicit none
character(len=50) :: eqmd_file,neq1_file,neq2_file,topfile
integer :: nsteps, spatial, natoms
double precision :: coordinate, dt
double precision, allocatable, dimension (:,:) :: tdspec
double precision, allocatable, dimension (:)   :: atmass
integer                :: specframes
integer                :: specdim
INTEGER*8              :: plan

call get_command_argument(1,eqmd_file)
call get_command_argument(2,neq1_file)
call get_command_argument(3,neq2_file)
call get_command_argument(4,topfile)

!GET NETCDF INFORMATION
call getdims(eqmd_file,nsteps,spatial,natoms)
dt=1.0d-15
!specframes=int(nsteps/2)
specframes=int(40000)
specdim=int(specframes/2+1)

!allocate( tdspec((NSTEPS/2+1)/2+1,NSTEPS/2), atmass(natoms) )
allocate( tdspec(5000,10000), atmass(natoms) )

!GET ATOMIC MASSES FROM PRMTOP FILE AND SQRT THEM.
call get_at_mass_prmtop(natoms,topfile,atmass)

!MAKE PLAN FOR FOURIER TRANSFORM LIBRARY
call make_fftw_plan(eqmd_file,nsteps,plan)

!COMPUTE PUMP-PROBE SPECTRA
call td_spectrum(eqmd_file,neq1_file,neq2_file,atmass,plan,nsteps,natoms,spatial,dt,tdspec)
call print_results(nsteps,dt,tdspec)

! COMPUTE KINETIC ENERGY DIFFERENCES.
call kinetic_energy(eqmd_file,neq1_file,neq2_file,atmass,plan,nsteps,natoms,spatial,dt)

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

write(*,*) 'inside getdims. Attempting to open file ', infile
CALL check(nf90_open(infile, nf90_nowrite, ncid))
write(*,*) 'opened file ', infile, ncid
CALL check(nf90_inquire_dimension(ncid,1,xname,nsteps))
write(*,*) 'inquired nsteps ', xname, nsteps
CALL check(nf90_inquire_dimension(ncid,2,yname,spatial))
write(*,*) 'inquired spatial ', yname, spatial
CALL check(nf90_inquire_dimension(ncid,3,zname,natoms))
write(*,*) 'inquired natoms ', zname, natoms
CALL check(nf90_close(ncid))
write(*,*) 'closed file ', ncid
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
!DOUBLE PRECISION, DIMENSION(NSTEPS/2)          :: veltraj_window
DOUBLE PRECISION, DIMENSION(40000)          :: veltraj_window
!COMPLEX*16      , DIMENSION((NSTEPS/2+1)/2+1)  :: out
COMPLEX*16      , DIMENSION(40000/2+1)  :: out

INTEGER :: ncid,i,j,k,time0,specdim,specframes

write(*,*) 'MAKING FFTW PLAN...'

!specframes=int(nsteps/2)
specframes=int(40000)
specdim=int(specframes/2+1)

! READ FIRST LINE OF DATASET AS A REPRESENTATIVE.
write(*,*) 'GETTING TRAJECTORY FROM NETCDF FILE'
CALL check(nf90_open(infile, nf90_nowrite, ncid))
point = (/ 1, 1, 1 /)
endp = (/ 1,1,nsteps-1 /)
CALL check(nf90_get_var(ncid,4,coordinate_evol,start = point, count = endp))
CALL check(nf90_close(ncid))
!veltraj_window=coordinate_evol(1:int(nsteps/2))
veltraj_window=coordinate_evol(1:40000)

! MEASURE THE BEST PLAN FOR FOURIER TRANSFORM. THIS IS MORE TIME-CONSUMING
! THAN THE ESTIMATE OPTION, BUT GIVEN THE AMOUNT OF CALLS TO FFTW IT IS 
! WORTH THE EFFORT.
write(*,*) 'FINDING BEST PLAN'
CALL dfftw_plan_dft_r2c_1d(plan,specframes,veltraj_window,out,"FFTW_MEASURE")

write(*,*) 'FINISHED FFTW PLAN!'

end subroutine make_fftw_plan

SUBROUTINE td_spectrum(eqmd_file,neq1_file,neq2_file,atmass,plan,nsteps,natoms,spatial,dt,tdspec)
USE netcdf
USE,INTRINSIC :: iso_c_binding
IMPLICIT NONE
include 'fftw3.f03'
CHARACTER(LEN=50), INTENT(IN) :: eqmd_file,neq1_file,neq2_file
INTEGER(KIND=4),   INTENT(IN) :: natoms, nsteps, spatial
INTEGER*8,         INTENT(IN) :: plan
DOUBLE PRECISION,  INTENT(IN) :: atmass(natoms)
DOUBLE PRECISION :: coordinate,dt
DOUBLE PRECISION, DIMENSION(NSTEPS) :: eqmd_traj,neq1_traj,neq2_traj
!DOUBLE PRECISION, DIMENSION(NSTEPS/2) :: veltraj_window
DOUBLE PRECISION, DIMENSION(40000) :: veltraj_window
!COMPLEX*16      , DIMENSION((NSTEPS/2+1)/2+1) :: out
COMPLEX*16      , DIMENSION(40000/2+1) :: out
!DOUBLE PRECISION, DIMENSION((NSTEPS/2+1)/2+1,NSTEPS/2),intent(out) :: tdspec
DOUBLE PRECISION, DIMENSION(5000,10000),intent(out) :: tdspec

INTEGER(KIND=4) :: ncid1, ncid2, ncid3, xtype, ndims, varid
CHARACTER(LEN=50) :: xname, vname
CHARACTER(LEN=20) :: fmt
integer, dimension(3) :: point,endp
integer :: i,j,k,time0,ierr,specdim,specframes

!specframes=int(nsteps/2)
specframes=40000
specdim=int(specframes/2+1)

CALL check(nf90_open(eqmd_file, nf90_nowrite, ncid1))
CALL check(nf90_open(neq1_file, nf90_nowrite, ncid2))
CALL check(nf90_open(neq2_file, nf90_nowrite, ncid3))

tdspec=0d0
eqmd_traj=0d0
neq1_traj=0d0
neq2_traj=0d0

do i = 1,natoms !natoms
  write(*,*) 'STARTING ATOM ',i,' OUT OF ',natoms
  do j = 1,3 !3
      point = (/ j, i, 1 /)
      endp = (/ 1,1,nsteps-1 /)
      CALL check(nf90_get_var(ncid1,4,eqmd_traj,start = point,count = endp))
      CALL check(nf90_get_var(ncid2,4,neq1_traj,start = point,count = endp))
      CALL check(nf90_get_var(ncid3,4,neq2_traj,start = point,count = endp))
      eqmd_traj=eqmd_traj*atmass(i)
      neq1_traj=neq1_traj*atmass(i)
      neq2_traj=neq2_traj*atmass(i)
      do time0=1,10000
         veltraj_window=eqmd_traj(time0:time0+40000-1)
         call dfftw_execute_dft_r2c(plan, veltraj_window, out)
         out = conjg(out)*out
         tdspec(:,time0) = tdspec(:,time0) -       REAL(out(1:5000)) ! DEBUG

         veltraj_window=neq1_traj(time0:time0+40000-1)
         call dfftw_execute_dft_r2c(plan, veltraj_window, out)
         out = conjg(out)*out
         tdspec(:,time0) = tdspec(:,time0) + 0.5d0*REAL(out(1:5000)) ! DEBUG

         veltraj_window=neq2_traj(time0:time0+40000-1)
         call dfftw_execute_dft_r2c(plan, veltraj_window, out)
         out = conjg(out)*out
         tdspec(:,time0) = tdspec(:,time0) + 0.5d0*REAL(out(1:5000)) ! DEBUG
      end do
  end do
end do

! CLOSE NETCDF FILES
CALL check(nf90_close(ncid1))
CALL check(nf90_close(ncid2))
CALL check(nf90_close(ncid3))

end subroutine td_spectrum

SUBROUTINE kinetic_energy(eqmd_file,neq1_file,neq2_file,atmass,plan,nsteps,natoms,spatial,dt)
USE netcdf
USE,INTRINSIC :: iso_c_binding
IMPLICIT NONE
include 'fftw3.f03'
CHARACTER(LEN=50), INTENT(IN) :: eqmd_file,neq1_file,neq2_file
INTEGER(KIND=4),   INTENT(IN) :: natoms, nsteps, spatial
INTEGER*8,         INTENT(IN) :: plan
DOUBLE PRECISION,  INTENT(IN) :: atmass(natoms)
DOUBLE PRECISION :: coordinate,dt
DOUBLE PRECISION, DIMENSION(NSTEPS) :: eqmd_traj,neq1_traj,neq2_traj
DOUBLE PRECISION, DIMENSION(NSTEPS) :: eqmd_ke,neq1_ke,neq2_ke,dif_ke
DOUBLE PRECISION, DIMENSION(NSTEPS/2) :: veltraj_window
COMPLEX*16      , DIMENSION((NSTEPS/2+1)/2+1) :: out

INTEGER(KIND=4) :: ncid1, ncid2, ncid3
CHARACTER(LEN=50) :: xname, vname
CHARACTER(LEN=20) :: fmt
integer, dimension(3) :: point,endp
integer :: i,j,k,time0,err,specdim,specframes

specframes=int(nsteps/2)
specdim=int(specframes/2+1)

CALL check(nf90_open(eqmd_file, nf90_nowrite, ncid1))
CALL check(nf90_open(neq1_file, nf90_nowrite, ncid2))
CALL check(nf90_open(neq2_file, nf90_nowrite, ncid3))

tdspec    = 0d0
eqmd_traj = 0d0
neq1_traj = 0d0
neq2_traj = 0d0
eqmd_ke   = 0d0
neq1_ke   = 0d0
neq2_ke   = 0d0

do i = 1,natoms !natoms
  write(*,*) 'STARTING ATOM ',i,' OUT OF ',natoms
  do j = 1,3 !3
      point = (/ j, i, 1 /)
      endp = (/ 1,1,nsteps-1 /)
      CALL check(nf90_get_var(ncid1,4,eqmd_traj,start = point,count = endp))
      CALL check(nf90_get_var(ncid2,4,neq1_traj,start = point,count = endp))
      CALL check(nf90_get_var(ncid3,4,neq2_traj,start = point,count = endp))
      eqmd_traj=eqmd_traj*atmass(i)
      neq1_traj=neq1_traj*atmass(i)
      neq2_traj=neq2_traj*atmass(i)
      eqmd_ke = eqmd_ke + 0.5d0 * eqmd_traj * eqmd_traj
      neq1_ke = neq1_ke + 0.5d0 * neq1_traj * neq1_traj
      neq2_ke = neq2_ke + 0.5d0 * neq2_traj * neq2_traj
  end do
end do

! COMPUTE KINETIC ENERGY DIFFERENCE
dif_ke = 0.5d0 * (neq1_ke + neq2_ke) - eqmd_ke

! CLOSE NETCDF FILES
CALL check(nf90_close(ncid1))
CALL check(nf90_close(ncid2))
CALL check(nf90_close(ncid3))

!====================================================================
! PRINT RESULTS
!====================================================================
open(unit=123,file='diff_ke.dat',iostat=err)
do i=1,nsteps-1
   write(123,*) dif_ke(i)
end do
close(123)

open(unit=124,file='eqmd_ke.dat',iostat=err)
do i=1,nsteps-1
   write(124,*) eqmd_ke(i)
end do
close(124)

open(unit=124,file='neq1_ke.dat',iostat=err)
do i=1,nsteps-1
   write(124,*) neq1_ke(i)
end do
close(124)

open(unit=124,file='neq2_ke.dat',iostat=err)
do i=1,nsteps-1
   write(124,*) neq2_ke(i)
end do
close(124)

end subroutine kinetic_energy

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

!specframes=int(nsteps/2)
specframes=int(40000)
specdim=int(specframes/2+1)

write(fmt,'("(",I6,"D24.15)")') nsteps/2+1
write(*,*) fmt
!Writing sum spectra to disk.
open(unit=1000,file='sum_spectra.dat',iostat=ierr)
if (ierr /= 0) then
   write(*,'(A,A)') 'ERROR OPENING INPUT FILE sum_spectra.dat'
   STOP
end if
do time0=1,10000
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
do i = 1,10000
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

