!...Subroutine to compute the kernel CI2,GI2,UI2
!...see note for explanation
Subroutine Kernel2(dof,region,source,field,C2,G2,U2)
	USE DefinitionConstant
    USE GlobalData
    IMPLICIT NONE
    
	INTEGER, INTENT(IN)			:: region,dof
	REAL(KIND=DBL),INTENT(IN)	:: source(2),field(2)
    REAL(KIND=DBL),INTENT(OUT)	:: C2(dof,dof),G2(dof,dof),U2(dof,dof)

	!...local variables
	REAL(KIND=DBL)	:: r(2),r_mag,theta,temp,zeta,f(3)
    REAL(KIND=DBL)	:: ey1,pv1
	REAL(KIND=DBL)	:: CI2_table(dof,dof,181),GI2_table(dof,dof,181),UI2_table(dof,dof,181)
    INTEGER		:: i,j,itemp,itheta(3),quadrant
    INTEGER		:: matl_no,matl_type
	!...temporarily set value of n_theta and f_theta for general anisotropy
    INTEGER			:: n_theta
    REAL(KIND=DBL)	:: f_theta
    n_theta = 181
    f_theta = PI
    
    !...id # of material
    matl_no = region_mat_no(region)
    !...type of material
    matl_type = material_type(matl_no)
    
	!...position vector r = xi - y
    do i=1,2
      	r(i) = field(i) - source(i)
    enddo
    !...magnitude of r
    r_mag = dsqrt(r(1)*r(1) + r(2)*r(2))
    !...defensive programming
    if (r_mag.lt.SMALL_NUM) then
      	print*,'kernel2.f90: cannot compute theta when source point = field point'
        stop
    endif
    !...determine which quadrant the vector r is
    if (r(2).ge.SMALL_NUM) then
      	if (r(1).ge.SMALL_NUM) then
        	quadrant = 1
        else
          	quadrant = 2
        endif
    else
      	if (r(1).ge.SMALL_NUM) then
        	quadrant = 4
        else
          	quadrant = 3
        endif
    endif
    !...angle theta of vector r wrt 1-axis
    if ((quadrant.eq.1).or.(quadrant.eq.2)) then
      	theta = dacos(r(1)/r_mag)
    else
      	!...theta is pre-computed on (0-PI) since ln|r.e|=ln|r.(-e)|
      	theta = PI - dacos(r(1)/r_mag)
    endif
    !...debug
	!print'(a10,i1,a10,f6.3)','quadrant=',quadrant,'theta= ',theta
    !...compute CI2 kernel
	C2 = 0.0d0
    G2 = 0.0d0
    U2 = 0.0d0
    if (matl_type.eq.1) then
      	if (dof.eq.3) then
      		!...isotropy, get U2,G2 and C2 from closed form solution
        	ey1 = material(1,matl_no)
        	pv1 = material(2,matl_no)
      		!...compute CI2
        	temp = 0.25d0*ey1/(PI*(1.0d0-pv1*pv1))
        	C2(1,1) = -temp*dcos(theta)*dcos(theta)
        	C2(1,2) = -temp*dcos(theta)*dsin(theta)
        	C2(2,1) = C2(1,2)
        	C2(2,2) = -temp*dsin(theta)*dsin(theta)
        	!...compute GI2
        	temp = -0.25d0/(PI*(1.0d0-pv1))
        	G2(1,1) = temp*dcos(theta)*dsin(theta)
        	G2(1,2) = temp*dsin(theta)*dsin(theta)
        	G2(2,1) = -temp*dcos(theta)*dcos(theta)
        	G2(2,2) = -temp*dcos(theta)*dsin(theta)
        	!...compute UI2
        	temp = 0.25d0*(1.0d0+pv1)/(PI*ey1*(1.0d0-pv1))
        	U2(1,1) = temp*dcos(theta)*dcos(theta)
        	U2(1,2) = temp*dcos(theta)*dsin(theta)
        	U2(2,1) = U2(1,2)
        	U2(2,2) = temp*dsin(theta)*dsin(theta)
        else
          	print*,'not yet developed for isotropy of piezo/magneto case'
            stop
        endif
    else
      	!print*,'kernel CI2 is computed for anisotropy'
      	!...obtain 3 positions of theta in the table to be interpolated
        itemp = nint(dble(n_theta-1)*theta/f_theta)+1
        if (itemp.eq.1) then
          	itheta(1) = 1
            itheta(2) = 2
            itheta(3) = 3
        elseif (itemp.eq.n_theta) then
        	itheta(1) = n_theta - 2
            itheta(2) = n_theta - 1
            itheta(3) = n_theta
        else
          	itheta(1) = itemp - 1
            itheta(2) = itemp
            itheta(3) = itemp + 1
        endif
        !...obtain value of shape functions for interpolation
        zeta = dble(n_theta-1)*theta/f_theta - dble(itheta(1))
        f(1) = 0.5d0*zeta*(zeta-1.0d0)
        f(2) = (1.0d0-zeta)*(1.0d0+zeta)
        f(3) = 0.5d0*zeta*(zeta+1.0d0)
        !...interpolate value of C2,G2,U2
        do i=1,dof
          	do j=1,dof
            	C2(i,j) = CI2_table(i,j,itheta(1))*f(1) + CI2_table(i,j,itheta(2))*f(2) +&
                		CI2_table(i,j,itheta(3))*f(3)
                G2(i,j) = GI2_table(i,j,itheta(1))*f(1) + GI2_table(i,j,itheta(2))*f(2) +&
                		GI2_table(i,j,itheta(3))*f(3)
                U2(i,j) = UI2_table(i,j,itheta(1))*f(1) + UI2_table(i,j,itheta(2))*f(2) +&
                		UI2_table(i,j,itheta(3))*f(3)
            enddo
        enddo
    endif
End subroutine Kernel2