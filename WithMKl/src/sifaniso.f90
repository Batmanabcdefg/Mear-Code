SUBROUTINE CompSIFaniso(dof,elid,xele,uk,SIF,Emat)
	USE DefinitionConstant
  	USE GlobalData
    USE Remesh
	IMPLICIT NONE
    
	INTEGER, INTENT(IN) 		:: elid,dof
	REAL(KIND=DBL), INTENT(IN) 	:: xele(2,3),uk(dof),Emat(3,dof,dof,3)
    REAL(KIND=DBL), INTENT(OUT)	:: SIF(dof)
	! local variables
	INTEGER :: k,n,nve,ntt,ii,i1,i2,i3,i4,ierr
    REAL(KIND=DBL)	:: xitip
	REAL(KIND=DBL) 	:: dxds(2),ps(3),dps(3),eita(2),fa(2),fk(dof),&
                       dj, att(dof,dof),betaa, sinb, cosb, &
					   a(2),b(2),aa(dof,dof),bb(dof,dof),bbin(dof,dof),&
					   ab(dof,dof),BX(dof,dof),BTEMP(dof,dof),str(dof)
    REAL(KIND=DBL), ALLOCATABLE	:: zc(:),wt(:)

    !debug
!    REAL(KIND=DBL):: temp(dof,dof)

	!...Specify number of integration point for loop integration
    ntt = 100
    allocate (zc(ntt),wt(ntt),STAT= ierr)
    if (ierr.ne.0) then
      	print*,'zc,wt: allocation request denied'
        stop
    endif
    if (ntt.gt.max_nint) then
      	print*,'sifaniso.f90: ntt is greater than max_nint'
        stop
    endif
    do n = 1,ntt
      	zc(n) = xi(n,ntt)
        wt(n) = wi(n,ntt)
    enddo
	!CALL GaussIntPoint(ntt,zc,wt)
    
	!Obtain the local coordinates for crack front nodes 1
	if (elid.eq.CTIP1) then
		xitip = -1.0d0
    elseif (elid.eq.CTIP2) then
    	xitip = 1.0d0
    else
      	print*,'CompSIF: wrong element id for tip element'
        stop
    endif

	!...Initialize fk(3) and dxds(2)
    fk = 0.0d0
    dxds = 0.0d0

	!...Obtain the shape functions and their derivatives
    CALL reg_shape(elid, xitip, ps, dps)

	!...Compute dxds
    DO n = 1,NODE(elid)
    	DO  k = 1,2
        	dxds(k) = dxds(k) + xele(k,n)*dps(n)
        ENDDO
    ENDDO
    
	!...Compute dj
    dj = dsqrt(dxds(1)*dxds(1)+dxds(2)*dxds(2))

	!...Compute eita(2) and fa(2)
    DO nve = 1,2
		if (elid.eq.CTIP1) then
			!eita(nve) = dxds(nve,2)/dj2
			!a bug in above line(the sign of SIF), fixed in 2/20/98
			eita(nve) = -dxds(nve)/dj
        else
          	eita(nve) = dxds(nve)/dj
        endif
    END DO
    !...compute the UNIT normal vector fa(2)=-eita().cross_product.e3 (
    !fa() is the 2-axis of local coordinates)
    if (elid.eq.CTIP1) then
		fa(1) = eita(2)
		fa(2) = -eita(1)
    elseif (elid.eq.CTIP2) then
    	fa(1) = -eita(2)
		fa(2) = eita(1)
    else
      	print*,'error: id of tip element is out of range'
        stop
    endif
    
	!...Compute fk = uk of tip node in local coordinates
    !------------------------------------------------------------
    if (dof.eq.1) then
      	fk = uk
    endif
    !------------------------------------------------------------
    DO k = 1,2
    	fk(1) = fk(1) + eita(k)*uk(k) !Note: this is for mode II
    	fk(2) = fk(2) + fa(k)*uk(k)   !this is for mode I
    END DO
    !...for other media than elastic: u4/u5 is unchanged (as u3 of elastic media)
    if (dof.ge.3) then
    	do k = 3,dof
    		fk(k) = uk(k)
    	enddo
    endif
	!...Obtain transformation matrix from global to local: att(i,j) = e'_i.e_j, where e'_i are local coordinates
    !...for other media, axis 4,5...not changed -> same as e3 for the case of elastic media
    if (dof.eq.1) then
      	att = 1.0d0
    else
    	att(1,1) = eita(1)
    	att(1,2) = eita(2)
    	att(2,1) = fa(1)
    	att(2,2) = fa(2)
    endif
    if (dof.ge.3) then
      	!...complete the 1st and 2nd  row
      	do k = 3,dof
        	att(1,k) = 0.0d0
            att(2,k) = 0.0d0
        enddo
        !...for the 3rd/4th/5th row
        do k = 3,dof
          	do n = 1,dof
            	if (k.eq.n) then
                	att(k,n) = 1.0d0
                else
                  	att(k,n) = 0.0d0
                endif
            enddo
        enddo
    endif

	!...Initialize BX matrix
    BX = 0.0d0

	!...Loop for contour integral (reduce from 0->2pi to 0->pi because of symmetry)
	do ii = 1,ntt
		!......dummy angle for integration, beta, varies from 0 to pi
    	betaa=(zc(ii)+1.0d0)*PI/2.0d0
    	sinb=dsin(betaa)
		cosb=dcos(betaa)

		!......compute unit vector a and b based on global coordinate system
       	do i1=1,2
	    	a(i1)= cosb*eita(i1)+sinb*fa(i1)
	    	b(i1)= -sinb*eita(i1)+cosb*fa(i1)
	   	end do

		!......compute (a,a),(b,b) and using symmetry of (a,a) and (b,b)
	   	do i1 = 1,dof
			do i2 = i1,dof
		    	aa(i1,i2) = 0.0d0
		     	bb(i1,i2) = 0.0d0
             	!...indexes of a() and b() go from 1-2
		     	do i3 = 1,2
		        	do i4 = 1,2
		 	       		aa(i1,i2) = aa(i1,i2) + a(i3)*Emat(i3,i1,i2,i4)*a(i4)
				   		bb(i1,i2) = bb(i1,i2) + b(i3)*Emat(i3,i1,i2,i4)*b(i4)
					end do
			 	end do
		  	end do
	   	end do
        do i1 = 2,dof
          	do i2 = 1,i1-1
            	aa(i1,i2) = aa(i2,i1)
                bb(i1,i2) = bb(i2,i1)
            enddo
        enddo
		!...compute inverse of (b,b)
        if (dof.eq.1) then
          	bbin = 1.0d0/bb
        elseif (dof.eq.2) then
          	call matrix2_inverse(bb,bbin)
        elseif (dof.eq.3) then
        	call matrix3_inverse(bb,bbin)
        elseif (dof.eq.4) then
        	call matrix4_inverse(bb,bbin)
        elseif (dof.eq.5) then
        	call matrix5_inverse(bb,bbin)
        else
          	print*,'sifaniso.f90: cannot compute inverse of bb for case of dof =',dof
            stop
        endif
        !debug
!        if (ii.eq.1) then
!          	print*,'bb='
!        	do i1 = 1,dof
!          		print'(4(es20.10,2x))',(bb(i1,i2),i2=1,dof)
!        	enddo
!        endif
!	   	!debug
!        if (ii.eq.1) then
!          	print*,'bbin='
!        	do i1 = 1,dof
!          		print'(4(es20.10,2x))',(bbin(i1,i2),i2=1,dof)
!        	enddo
!            do i1 = 1,dof
!              	do i2 = 1,dof
!                	temp(i1,i2) = 0.0d0
!                	do i3 = 1,dof
!                    	temp(i1,i2) = temp(i1,i2) + bb(i1,i3)*bbin(i3,i2)
!                    enddo
!                enddo
!            enddo
!            print*,'bb*bbin='
!            do i1 = 1,dof
!          		print'(4(es20.10,2x))',(temp(i1,i2),i2=1,dof)
!        	enddo
!        endif
        
		!......compute (a,b) and use (a,b)=(b,a)^T
	   	do i1 = 1,dof
	    	do i2 = 1,dof
	        	ab(i1,i2) = 0.0d0
             	!...indexes of a() and b() go from 1-2
	   	     	do i3 = 1,2
		        	do i4 = 1,2
			       		ab(i1,i2) = ab(i1,i2) + a(i3)*Emat(i3,i1,i2,i4)*b(i4)
                    end do
		     	end do
		  	end do
	   	end do

		!......obtain loop integral for one integration point
       	do i1 = 1,dof
	      	do i2 = i1,dof
		     	BX(i1,i2) = BX(i1,i2) - aa(i1,i2)*(PI/2.0d0)*wt(ii)
             	do i3 = 1,dof
		        	do i4 = 1,dof
		           		BX(i1,i2) = BX(i1,i2) + ab(i1,i3)*bbin(i3,i4)*ab(i2,i4)*(PI/2.0d0)*wt(ii)
			    	end do
		     	end do
		  	end do
	   	end do
	end do !end loop of integration points

	!...adjust the constant (1/2pi) and the fact that contour 
	!integral reduced from 0->2pi to 0->pi
    do i1 = 1,dof
      	do i2 = i1,dof
        	BX(i1,i2) = -BX(i1,i2)/PI
        enddo
    enddo
    !...symmetry of BX
    do i1 = 2,dof
      	do i2 = 1,i1-1
        	BX(i1,i2) = BX(i2,i1)
        enddo
    enddo

	!...Transform BX to local coordinate system
    do i1 = 1,dof
		do i2 = i1,dof
	    	BTEMP(i1,i2) = 0.0d0
	      	do i3 = 1,dof
	   	    	do i4 = 1,dof
			    	BTEMP(i1,i2) = BTEMP(i1,i2) + att(i1,i3)*att(i2,i4)*BX(i3,i4)
		     	end do
		  	end do
	   	end do
	end do
    do i1 = 1,dof
      	do i2 = i1,dof
        	BX(i1,i2) = BTEMP(i1,i2)
        enddo
    enddo
    !...using symmetry
    do i1 = 2,dof
      	do i2 = 1,i1-1
        	BX(i1,i2) = BX(i2,i1)
        enddo
    enddo
	
	do i1 = 1,dof
		str(i1) = 0.0d0
	   	do i2 = 1,dof
	    	str(i1) = str(i1) + dsqrt(PI/2.0d0)*BX(i1,i2)*fk(i2)/dsqrt(dj)
	   	end do
	end do
	!...obtain the SIF in usual definition
    if (dof.eq.1) then
      	SIF = str
    else
    	SIF(1) = str(2)
    	SIF(2) = str(1)
    endif
    if (dof.ge.3) then
    	do k = 3,dof
    		SIF(k) = str(k)
    	enddo
    endif
END SUBROUTINE CompSIFaniso

              

