program aka2
implicit none
character(len=50) :: infile
integer :: nsteps, spatial, natoms
real(4) :: coordinate 
real(4), allocatable, dimension (:) :: coordinate_evol
integer, dimension (3) :: point

call get_command_argument(1,infile)
call getdims(infile,nsteps,spatial,natoms)
nsteps=nsteps

allocate(coordinate_evol(nsteps))

call getandwritecoord(infile,nsteps,natoms,spatial,point,coordinate_evol)
   
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
integer, dimension(3) :: point,endp
integer :: i,j,k

open(1, file = 'tseries.dat', status = 'unknown')  
CALL check(nf90_open(infile, nf90_nowrite, ncid))
write(*,*) natoms
do i = 1,natoms !natoms
  do j = 1,3 !3
      point = (/ j, i, 1 /)
      endp = (/ 1,1,nsteps-1 /)
      CALL check(nf90_get_var(ncid,4,coordinate_evol,start = point,count = endp))
      write(1,*) coordinate_evol(1:nsteps)
  enddo
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


end program aka2

