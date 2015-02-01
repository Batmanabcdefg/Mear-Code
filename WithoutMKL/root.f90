subroutine zrhqr(order,a,rtr,rti)
use nrtype
use nrutil, only:assert_eq,nrerror
implicit none

integer::order
real(DP), dimension(order+1), intent(in)::a
real(DP), dimension(order), intent(out)::rtr,rti
integer(I4B)::k,m
integer(I4B), dimension(size(rtr))::indx
real(DP), dimension(size(a)-1,size(a)-1)::hess

m=assert_eq(size(rtr),size(rti),size(a)-1,'zrhqr')
if (a(m+1)==0.0) call &
	nrerror('zrhqr: Last value of array a must not 0')
hess(1,:)=-a(m:1:-1)/a(m+1)
hess(2:m,:)=0.0
do k=1,m-1
	hess(k+1,k)=1.0
end do
call balanc(order,hess)
call hqr(order,hess,rtr,rti)
call indexx_sp(order,rtr,indx)
rtr=rtr(indx)
rti=rti(indx)
end subroutine zrhqr
!*********************************************************
subroutine balanc(order,a)
use nrtype
use nrutil, only:assert_eq
implicit none

integer::order
real(DP), dimension(order,order), intent(inout)::a
real(DP), parameter::RADX=radix(a),SQRADX=RADX**2
integer(I4B)::i,last,ndum
real(DP)::c,f,g,r,s
ndum=assert_eq(size(a,1),size(a,2),'balanc')
do
	last=1
	do i=1,size(a,1)                         !calculate row and column norms.
		c=sum(abs(a(:,i)))-a(i,i)
		r=sum(abs(a(i,:)))-a(i,i)
		if (c/=0.0 .and. r/=0.0) then        !if both are nonzero
			g=r/RADX
			f=1.0
			s=c+r
			do                               !find the integer power of the
				if (c>=g) exit               !machine radix that comes closest
				f=f*RADX                     !to balanceing the matrix.
				c=c*SQRADX
			end do
			g=r*RADX
			do
				if (c<=g) exit
				f=f/RADX
				c=c/SQRADX
			end do
			if ((c+r)/f<0.95_sp*s) then
				last=0
				g=1.0_sp/f
				a(i,:)=a(i,:)*g              !apply similarity transformation.
				a(:,i)=a(:,i)*f
			end if
		end if
	end do
	if (last/=0) exit
end do
end subroutine balanc
!*********************************************************
SUBROUTINE hqr(order,a,wr,wi)
use nrtype
use nrutil, only:assert_eq,diagadd,nrerror,upper_triangle
IMPLICIT NONE

integer::order
REAL(DP), DIMENSION(order), INTENT (OUT) :: wr,wi
REAL(DP), DIMENSION(order,order), INTENT(INOUT) :: a

INTEGER (I4B) :: i,its,k,l,m,n,nn,mnnk

REAL(DP) :: anorm,p,q,r,s,t,u,v,w,x,y,z
REAL(DP) , DIMENSION(size(a,1)) :: pp
n=assert_eq(size(a,1),size(a,2),size(wr),size(wi),'hqr')
anorm=sum(abs(a),mask=upper_triangle(n,n,extra=2))
!Compute matrix norm for possibleuse in locatingsinglesmall subdiagonal element.
nn=n
t=0.0
do
	if (nn < 1) exit
	its=0
	iterate: do
		small: do l=nn,2,-1
			s=abs(a(l-1,l-1))+abs(a(l,l))
			if (s == 0.0) s=anorm
			if (abs(a(l,l-1))+s == s) then
				a(l,l-1)=0.0
				exit small
			end if
		end do small
		x=a(nn,nn)
		if (l == nn) then
			wr(nn)=x+t
			wi(nn)=0.0
			nn=nn-1
			exit iterate
		end if
		y=a(nn-1,nn-1)
		w=a(nn,nn-1)*a(nn-1,nn)
		if (l == nn-1) then
			p=0.5_sp*(y-x)
			q=p**2+w
			z=sqrt(abs(q))
			x=x+t
			if (q >= 0.0) then
				z=p+sign(z,p)
				wr(nn)=x+z
				wr(nn-1)=wr(nn)
				if (z /= 0.0) wr(nn)=x-w/z
				wi(nn)=0.0
				wi(nn-1)=0.0
			else
				wr(nn)=x+p
				wr(nn-1)=wr(nn)
				wi(nn)=z
				wi(nn-1) =-z
			end if
			nn=nn-2
			exit iterate
		end if
!No roots found. Continue iteration.
		if (its == 30) call nrerror('too many iterations in hqr')
		if (its == 10 .or. its == 20) then           !Form exceptional shift.
			t=t+x
			call diagadd(nn,a(1:nn,1:nn),-x)
			s=abs(a(nn,nn-1))+abs(a(nn-1,nn-2))
			x=0.75_sp*s
			y=x
			w=-0.4375_sp*s**2
		end if
		its=its+l
		do m=nn-2,l,-1
			z=a(m,m)
			r=x-z
			s=y-z
			p=(r*s-w)/a(m+1,m)+a(m,m+1)
			q=a(m+1,m+1)-z-r-s
			r=a(m+2,m+1)
			s=abs(p)+abs(q)+abs(r)
			p=p/s
			q=q/s
			r=r/s
			if (m == l) exit
			u=abs(a(m,m-l))*(abs(q)+abs(r))
			v=abs(p)*(abs(a(m-1,m-l))+abs(z)+abs(a(m+1,m+1)))
			if (u+v == v) exit
		end do
		do i=m+2,nn
			a(i,i-2)=0.0
			if (i /= m+2) a(i,i-3)=0.0
		end do
		do k=m,nn-l
			if (k /= m) then
				p=a(k,k-1)
				q=a(k+1,k-1)
				r=0.0
				if (k /= nn-1) r=a(k+2,k-1)
				x=abs(p)+abs(q)+abs(r)
				if (x /= 0.0) then
					p=p/x
					q=q/x
					r=r/x
				end if
			end if
			s=sign(sqrt(p**2+q**2+r**2),p)
			if (s /= 0.0) then
				if (k == m) then
					if (l /= m) a(k,k-1)=-a(k,k-1)
				else
					a(k,k-1)=-s*x
				end if
				p=p+s
				x=p/s
				y=q/s
				z=r/s
				q=q/p
				r=r/p
				pp(k:nn)=a(k,k:nn)+q*a(k+1,k:nn)
				if (k /= nn-1) then
					pp(k:nn)=pp(k:nn)+r*a(k+2,k:nn)
					a(k+2,k:nn)=a(k+2,k:nn)-pp(k:nn)*z
				end if
				a(k+1,k:nn)=a(k+1,k:nn)-pp(k:nn)*y
				a(k,k:nn)=a(k,k:nn)-pp(k:nn)*x
				mnnk=min(nn,k+3)                !Column modification.
				pp(l:mnnk)=x*a(l:mnnk,k)+y*a(l:mnnk,k+1)
				if (k /= nn-1) then
					pp(l:mnnk)=pp(l:mnnk)+z*a(l:mnnk,k+2)
					a(l:mnnk,k+2)=a(l:mnnk,k+2)-pp(l:mnnk)*r
				end if
				a(l:mnnk,k+1)=a(l:mnnk,k+1)-pp(l:mnnk)*q
				a(l:mnnk,k)=a(l:mnnk,k)-pp(l:mnnk)
			end if
		end do
	end do iterate
end do
END SUBROUTINE hqr
!*********************************************************
SUBROUTINE indexx_sp(order,arr,index)
use nrtype
use nrutil, only:arth,assert_eq,nrerror,swap
IMPLICIT NONE

integer::order
REAL(DP), DIMENSION(order), INTENT(IN) :: arr
INTEGER(I4B), DIMENSION(order), INTENT(OUT) :: index
INTEGER(I4B), PARAMETER :: NN=15, NSTACK=50
!Indexes an array arr, i.e., outputs the array index of length N such that arr(index(j))
!is in ascending order for j = 1,2,..., N. The input quantity arr is not changed.
REAL (DP) :: a
INTEGER(I4B) :: n,k,i,j,indext,jstack,l,r
INTEGER(I4B), DIMENSION(NSTACK) :: istack
n=assert_eq(size(index),size(arr),'indexx_sp')
index=arth(1,1,n)
jstack=0
l=1
r=n
do
	if (r-l < NN) then
		do j=l+1,r
			indext=index(j)
			a=arr(indext)
			do i=j-1,l,-1
				if (arr(index(i)) <= a) exit
				index(i+1)=index(i)
			end do
			index(i+1)=indext
		end do
		if (jstack == 0) RETURN
		r=istack(jstack)
		l=istack(jstack-1)
		jstack=jstack-2
	else
		k=(l+r)/2
		call swap(index(k),index(l+1))
		call icomp_xchg(index(l),index(r))
		call icomp_xchg(index(l+1),index(r))
		call icomp_xchg(index(l),index(l+1))
		i=l+1
		j=r
		indext=index(l+1)
		a=arr(indext)
		do
			do
				i=i+1
				if (arr(index(i)) >= a) exit
			end do
			do
				j=j-1
				if (arr(index(j)) <= a) exit
			end do
			if (j < i) exit
			call swap(index(i),index(j))
		end do
		index(l+1)=index(j)
		index(j)=indext
		jstack=jstack+2
		if (jstack > NSTACK) call nrerror('indexx: NSTACK too small')
		if (r-i+1 >= j-l) then
			istack(jstack)=r
			istack(jstack-1)=i
			r=j-1
		else
			istack(jstack)=j-1
			istack(jstack-1)=l
			l=i
		end if
	end if
end do
CONTAINS
SUBROUTINE icomp_xchg(i,j)
INTEGER(I4B), INTENT (INOUT) :: i,j
INTEGER(I4B) :: swp
if (arr(j) < arr(i)) then
	swp=i
	i=j
	j=swp
end if
END SUBROUTINE icomp_xchg
END SUBROUTINE indexx_sp
!*********************************************************

SUBROUTINE indexx_i4b(order,iarr,index)
use nrtype
use nrutil,only:arth,assert_eq,nrerror,swap
IMPLICIT NONE

integer::order
INTEGER(I4B), DIMENSION(order), INTENT (IN) :: iarr
INTEGER (I4B), DIMENSION(order), INTENT (OUT) :: index
INTEGER (I4B), PARAMETER :: NN=15, NSTACK=50
INTEGER (I4B) :: a
INTEGER(I4B) :: n,k,i,j,indext,jstack,l,r

INTEGER (I4B), DIMENSION (NSTACK) :: istack
n=assert_eq(size(index) ,size (iarr) ,'indexx_sp')
index=arth(1,1,n)
jstack=0
l=1
r=n
do
	if (r-l < NN) then
		do j=l+1,r
			indext=index(j)
			a=iarr(indext)
			do i=j-1,l,-1
				if (iarr(index(i)) <= a) exit
				index(i+1)=index(i)
			end do
			index (i+1)=indext
		end do
		if (jstack == 0) RETURN
		r=istack(jstack)
		l=istack(jstack-1)
		jstack=jstack-2
	else
		k=(l+r)/2
		call swap(index(k),index(l+1))
		call icomp_xchg(index(l),index(r))
		call icomp_xchg(index(l+1),index(r))
		call icomp_xchg(index(l),index(l+1))
		i=l+1
		j=r
		indext=index(l+1)
		a=iarr(indext)
		do
			do
				i=i+1
				if (iarr(index(i)) >= a) exit
			end do
			do
				j=j-1
				if (iarr(index(j)) <= a) exit
			end do
			if (j < i) exit
			call swap(index(i),index(j))
		end do
		index(l+1)=index(j)
		index(j)=indext
		jstack=jstack+2
		if (jstack > NSTACK) call nrerror('indexx: NSTACK too small')
		if (r-i+1 >= j-l) then
			istack(jstack)=r
			istack(jstack-1)=i
			r=j-1
		else
			istack(jstack)=j-1
			istack(jstack-1)=l
			l=i
		end if
	end if
end do
CONTAINS
SUBROUTINE icomp_xchg(i,j)
INTEGER (I4B), INTENT (INOUT) :: i,j
INTEGER (I4B) :: swp
if (iarr(j) < iarr(i)) then
	swp=i
	i=j
	j=swp
end if
END SUBROUTINE icomp_xchg
END SUBROUTINE indexx_i4b