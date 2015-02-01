MODULE nrutil
USE nrtype
IMPLICIT NONE
INTEGER (I4B), PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8
INTEGER (I4B), PARAMETER :: NPAR_GEOP=4,NPAR2_GEOP=2
INTEGER (I4B), PARAMETER :: NPAR_CUMSUM=16
INTEGER (I4B), PARAMETER :: NPAR_CUMPROD=8
INTEGER (I4B), PARAMETER :: NPAR_POLY=8
INTEGER (I4B), PARAMETER :: NPAR_POLYTERM=8

INTERFACE swap
MODULE PROCEDURE swap_i,swap_r,swap_rv,swap_c,swap_cv,swap_cm,swap_z,swap_zv,swap_zm,masked_swap_rs,masked_swap_rv,masked_swap_rm
END INTERFACE

INTERFACE assert_eq
MODULE PROCEDURE assert_eq2,assert_eq3,assert_eq4,assert_eqn
END INTERFACE

INTERFACE arth
MODULE PROCEDURE arth_r, arth_d, arth_i
END INTERFACE

INTERFACE outerprod
MODULE PROCEDURE outerprod_r,outerprod_d
END INTERFACE

INTERFACE outerdiff
MODULE PROCEDURE outerdiff_r,outerdiff_d,outerdiff_i
END INTERFACE

INTERFACE diagadd
MODULE PROCEDURE diagadd_rv,diagadd_r
END INTERFACE

CONTAINS

SUBROUTINE swap_i(a,b)
INTEGER (I4B), INTENT(INOUT) :: a,b
INTEGER (I4B) :: dum
dum=a
a=b
b=dum
END SUBROUTINE swap_i

SUBROUTINE swap_r(a,b)
REAL(SP), INTENT (INOUT) :: a,b
REAL(SP) :: dum
dum=a
a=b
b=dum
END SUBROUTINE swap_r

SUBROUTINE swap_rv(a,b)
REAL (SP), DIMENSION(:), INTENT (INOUT) :: a,b
REAL (SP), DIMENSION(SIZE(a)) :: dum
dum=a
a=b
b=dum
END SUBROUTINE swap_rv

SUBROUTINE swap_c(a,b)
COMPLEX (SPC), INTENT (INOUT) :: a,b
COMPLEX (SPC) :: dum
dum=a
a=b
b=dum
END SUBROUTINE swap_c

SUBROUTINE swap_cv(a,b)
COMPLEX (SPC), DIMENSION(:), INTENT (INOUT) :: a,b
COMPLEX (SPC), DIMENSION(SIZE(a)) :: dum
dum=a
a=b
b=dum

END SUBROUTINE swap_cv
SUBROUTINE swap_cm(a,b)
COMPLEX (SPC), DIMENSION(:,:), INTENT (INOUT) :: a,b
COMPLEX (SPC), DIMENSION(size(a,1),size(a,2)) :: dum
dum=a
a=b
b=dum
END SUBROUTINE swap_cm

SUBROUTINE swap_z(a,b)
COMPLEX (DPC), INTENT (INOUT) :: a,b
COMPLEX (DPC) :: dum
dum=a
a=b
b=dum
END SUBROUTINE swap_z

SUBROUTINE swap_zv(a,b)
COMPLEX (DPC), DIMENSION(:), INTENT (INOUT) :: a,b
COMPLEX (DPC), DIMENSION(SIZE(a)) :: dum
dum=a
a=b
b=dum
END SUBROUTINE swap_zv

SUBROUTINE swap_zm(a,b)
COMPLEX (DPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
COMPLEX(DPC), DIMENSION(size(a,1),size(a,2)) :: dum
dum=a
a=b
b=dum
END SUBROUTINE swap_zm

SUBROUTINE masked_swap_rs(a,b,mask)
REAL(SP), INTENT(INOUT) :: a,b
LOGICAL (LGT), INTENT(IN) :: mask
REAL(SP) :: swp
if (mask) then
	swp=a
	a=b
	b=swp
end if
END SUBROUTINE masked_swap_rs

SUBROUTINE masked_swap_rv(a,b,mask)
REAL (SP), DIMENSION(:), INTENT (INOUT) :: a,b
LOGICAL (LGT), DIMENSION(:), INTENT (IN) :: mask
REAL(SP), DIMENSION(size(a)) :: swp
where (mask)
	swp=a
	a=b
	b=swp
end where
END SUBROUTINE masked_swap_rv

SUBROUTINE masked_swap_rm(a,b,mask)
REAL (SP), DIMENSION(:,:), INTENT (INOUT) :: a,b
LOGICAL (LGT), DIMENSION(:,:), INTENT(IN) :: mask
REAL (SP), DIMENSION(size(a,1),size(a,2)) :: swp
where (mask)
	swp=a
	a=b
	b=swp
end where
END SUBROUTINE masked_swap_rm


!Routines for argument checking and wrror handling:
!SUBROUTINE assert1(n1,string)
!!Report and die if any logical is false (used for arg range checking).
!CHARACTER(LEN=*), INTENT(IN) :: string
!LOGICAL, INTENT(IN) :: n1
!if (.not. n1) then
!	write (*,*) 'nrerror: an assertion failed with this tag:', &
!	string
!	STOP 'program terminated by assert1'
!end if
!END SUBROUTINE assert1
!*******************************************************************
!SUBROUTINE assert2(n1,n2,string)
!CHARACTER (LEN=*), INTENT(IN) :: string
!LOGICAL, INTENT(IN) :: n1,n2
!if (.not. (n1 .and. n2)) then
!	write (*,*) 'nrerror: an assertion failed with this tag:', &
!	string
!	STOP 'program terminated by assert2'
!end if
!END SUBROUTINE assert2
!*******************************************************************
!SUBROUTINE assert3(n1,n2,n3,string)
!CHARACTER(LEN=*), INTENT(IN) :: string
!LOGICAL, INTENT(IN) :: n1,n2,n3
!if (.not. (n1 .and. n2 .and. n3)) then
!	write (*,*) 'nrerror: an assertion failed with this tag:', &
!	string
!	STOP 'program terminated by assert3'
!end if
!END SUBROUTINE assert3
!*******************************************************************
!SUBROUTINE assert4(n1,n2,n3,n4,string)
!CHARACTER(LEN=*), INTENT (IN) :: string
!LOGICAL, INTENT(IN) :: n1,n2,n3,n4
!if (.not. (n1 .and. n2 .and. n3 .and. n4)) then
!	write (*,*) 'nrerror: an assertion failed with this tag:', &
!	string
!	STOP 'program terminated by assert4'
!end if
!END SUBROUTINE assert4
!*******************************************************************
!SUBROUTINE assert_v(n,string)
!CHARACTER(LEN=*), INTENT(IN) :: string
!LOGICAL, DIMENSION(:), INTENT(IN) :: n
!if (.not. all(n)) then
!	write (*,*) 'nrerror: an assertion failed with this tag:', &
!	string
!	STOP 'program terminated by assert_v'
!end if
!END SUBROUTINE assert_v


!*******************************************************************
FUNCTION assert_eq2(n1,n2,string)
!Report and die if integers not all equal (used for size checking).
CHARACTER(LEN=*), INTENT (IN) :: string
INTEGER, INTENT(IN) :: n1,n2
INTEGER :: assert_eq2
if (n1 == n2) then
	assert_eq2=n1
else
	write (*,*) 'nrerror: an assert_eq failed with this tag:', &
	string
	STOP 'program terminated by assert_eq2'
end if
END FUNCTION assert_eq2
!*******************************************************************
FUNCTION assert_eq3(n1,n2,n3,string)
CHARACTER(LEN=*), INTENT(IN) :: string
INTEGER, INTENT(IN) :: n1,n2,n3
INTEGER :: assert_eq3
if (n1 == n2 .and. n2 == n3) then
	assert_eq3=n1
else
	write (*,*) 'nrerror: an assert_eq failed with this tag:', &
	string
STOP 'program terminated by assert_eq3'
end if
END FUNCTION assert_eq3
!*******************************************************************
FUNCTION assert_eq4(n1,n2,n3,n4,string)
CHARACTER(LEN=*), INTENT(IN) :: string
INTEGER, INTENT(IN) :: n1,n2,n3,n4
INTEGER :: assert_eq4
if (n1 == n2 .and. n2 == n3 .and. n3 == n4) then
	assert_eq4=n1
else
	write (*,*) 'nrerror: an assert_eq failed with this tag:', &
	string
	STOP 'program terminated by assert_eq4'
end if
END FUNCTION assert_eq4
!*******************************************************************
FUNCTION assert_eqn(nn,string)
CHARACTER(LEN=*), INTENT(IN) :: string
INTEGER, DIMENSION(:), INTENT(IN) :: nn
INTEGER :: assert_eqn
if (all(nn(2:) == nn(1))) then
	assert_eqn=nn(1)
else
	write (*,*) 'nrerror: an assert_eq failed with this tag:', &
	string
	STOP 'program terminated by assert_eqn'
end if
END FUNCTION assert_eqn
!*******************************************************************
SUBROUTINE nrerror(string)
!Report a message, then die.
CHARACTER(LEN=*), INTENT(IN) :: string
write (*,*) 'nrerror: ',string
STOP 'program terminated by nrerror'
END SUBROUTINE nrerror


!*******************************************************************
!Routines relating to polynomials and recurrences:
FUNCTION arth_r(first,increment,n)
!Array function returning an arithmetic progression.
REAL(SP), INTENT (IN) :: first,increment
INTEGER(I4B), INTENT (IN) :: n
REAL (SP), DIMENSION(n) :: arth_r
INTEGER (I4B) :: k,k2
REAL(SP) :: temp
if (n > 0) arth_r(1)=first
	if (n <= NPAR_ARTH) then
		do k=2,n
			arth_r(k)=arth_r(k-1)+increment
		end do
	else
		do k=2,NPAR2_ARTH
			arth_r(k)=arth_r(k-1)+increment
		end do
		temp=increment*NPAR2_ARTH
		k=NPAR2_ARTH
		do
			if (k >= n) exit
			k2=k+k
			arth_r(k+1:min(k2,n))=temp+arth_r(1:min(k,n-k))
			temp=temp+temp
			k=k2
		end do
	end if
END FUNCTION arth_r
!*******************************************************************
FUNCTION arth_d(first,increment,n)
REAL (DP), INTENT (IN) :: first,increment
INTEGER (I4B), INTENT(IN) :: n
REAL(DP), DIMENSION(n) :: arth_d
INTEGER (I4B) :: k,k2
REAL (DP) :: temp
if (n > 0) arth_d(1)=first
if (n <= NPAR_ARTH) then
	do k=2,n
		arth_d(k)=arth_d(k-1)+increment
	end do
else
	do k=2,NPAR2_ARTH
		arth_d(k)=arth_d(k-1)+increment
	end do
	temp=increment*NPAR2_ARTH
	k=NPAR2_ARTH
	do
		if (k >= n) exit
		k2=k+k
		arth_d(k+1:min(k2,n))=temp+arth_d(1:min(k,n-k))
		temp=temp+temp
		k=k2
	end do
end if
END FUNCTION arth_d
!*******************************************************************
FUNCTION arth_i(first,increment,n)
INTEGER (I4B), INTENT (IN) :: first,increment,n
INTEGER(I4B), DIMENSION(n) :: arth_i
INTEGER(I4B) :: k,k2,temp
if (n > 0) arth_i(1)=first
if (n <= NPAR_ARTH) then
	do k=2,n
		arth_i(k)=arth_i(k-1)+increment
	end do
else
	do k=2,NPAR2_ARTH
		arth_i(k)=arth_i(k-1)+increment
	end do
	temp=increment*NPAR2_ARTH
	k=NPAR2_ARTH
	do
		if (k >= n) exit
		k2=k+k
		arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
		temp=temp+temp
		k=k2
	end do
end if
END FUNCTION arth_i
!*******************************************************************
!Routines for skew operations on matrices:
SUBROUTINE diagadd_rv(dim,mat,diag)
!Adds vector or scalar diag to the diagonal of matrix mat.
integer(I4B)::dim
REAL(DP), DIMENSION(dim,dim), INTENT (INOUT) :: mat
REAL(DP), DIMENSION(dim), INTENT (IN) :: diag
INTEGER(I4B) :: j,n
n = assert_eq2(size(diag),min(size(mat,1),size(mat,2)),'diagadd_rv')
do j=1,n
	mat(j,j)=mat(j,j)+diag(j)
end do
END SUBROUTINE diagadd_rv
!*******************************************************************
SUBROUTINE diagadd_r(dim,mat,diag)
integer(I4B)::dim
REAL(DP), DIMENSION(dim,dim), INTENT (INOUT) :: mat
REAL(DP), INTENT (IN) :: diag
INTEGER (I4B) :: j,n
n = min(size(mat,1),size(mat,2))
do j=1,n
	mat(j,j)=mat(j,j)+diag
end do
END SUBROUTINE diagadd_r


!*******************************************************************
FUNCTION outerprod_r(a,b)
REAL (SP), DIMENSION(:), INTENT(IN) :: a,b
REAL(SP), DIMENSION(size(a),size(b)) :: outerprod_r
outerprod_r = spread(a,dim=2,ncopies=size(b)) * &
    spread(b,dim=1,ncopies=size(a))
END FUNCTION outerprod_r
!*******************************************************************
FUNCTION outerprod_d(a,b)
REAL (DP), DIMENSION(:), INTENT(IN) :: a,b
REAL (DP), DIMENSION(size(a),size(b)) :: outerprod_d
outerprod_d = spread(a,dim=2,ncopies=size(b)) * &
    spread(b,dim=1,ncopies=size(a))
END FUNCTION outerprod_d
!*******************************************************************
FUNCTION outerdiv(a,b)
REAL (SP), DIMENSION(:), INTENT(IN) :: a,b
REAL(SP), DIMENSION(size(a),size(b)) :: outerdiv
outerdiv = spread(a,dim=2,ncopies=size(b)) / &
    spread(b,dim=1,ncopies=size(a))
END FUNCTION outerdiv
!*******************************************************************
FUNCTION outersum(a,b)
REAL (SP), DIMENSION(:), INTENT(IN) :: a,b
REAL (SP), DIMENSION(size(a),size(b)) :: outersum
outersum = spread(a,dim=2,ncopies=size(b)) + &
spread(b,dim=1,ncopies=size(a))
END FUNCTION outersum
!*******************************************************************
FUNCTION outerdiff_r(a,b)
REAL(SP), DIMENSION(:), INTENT (IN) :: a,b
REAL(SP), DIMENSION(size(a),size(b)) :: outerdiff_r
outerdiff_r = spread(a,dim=2,ncopies=size(b)) - &
spread(b,dim=1,ncopies=size(a))
END FUNCTION outerdiff_r
!*******************************************************************
FUNCTION outerdiff_d(a,b)
REAL(DP), DIMENSION(:), INTENT(IN) :: a,b
REAL (DP), DIMENSION(size(a),size(b)) :: outerdiff_d
outerdiff_d = spread(a,dim=2,ncopies=size(b)) - &
spread(b,dim=1 ,ncopies=size (a))
END FUNCTION outerdiff_d
!*******************************************************************
FUNCTION outerdiff_i(a,b)
INTEGER(I4B), DIMENSION(:), INTENT(IN) :: a,b
INTEGER (I4B), DIMENSION(size(a),size(b)) :: outerdiff_i
outerdiff_i = spread(a,dim=2,ncopies=size(b)) - &
spread(b,dim=1,ncopies=size(a))
END FUNCTION outerdiff_i
!*******************************************************************
FUNCTION outerand(a,b)
LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: a,b
LOGICAL(LGT), DIMENSION(size(a),size(b)) :: outerand
outerand = spread(a,dim=2,ncopies=size(b)) .and. &
    spread(b,dim=1,ncopies=size(a))
END FUNCTION outerand

!*******************************************************************
FUNCTION upper_triangle(j,k,extra)
!Return an upper triangular logical mask.
INTEGER(I4B), INTENT (IN) :: j,k
INTEGER (I4B), OPTIONAL, INTENT(IN) :: extra
LOGICAL(LGT), DIMENSION(j,k) :: upper_triangle
INTEGER(I4B) :: n
n=0
if (present(extra)) n=extra
upper_triangle=(outerdiff(arth_i(1,1,j),arth_i(1,1,k)) < n)
ENDFUNCTION upper_triangle
!*******************************************************************
FUNCTION lower_triangle(j,k,extra)
!Return a lower triangular logical mask.
INTEGER(I4B), INTENT(IN) :: j,k
INTEGER(I4B), OPTIONAL, INTENT(IN) :: extra
LOGICAL(LGT), DIMENSION(j,k) :: lower_triangle
INTEGER (I4B) :: n
n=0
if (present(extra)) n=extra
lower_triangle=(outerdiff(arth_i(1,1,j),arth_i(1,1,k)) > -n)
END FUNCTION lower_triangle
!*******************************************************************
END MODULE nrutil
