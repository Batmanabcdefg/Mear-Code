SUBROUTINE CompSIF(elid,xele,uk,SIF,ey,pv)
	USE DefinitionConstant
	USE GlobalData
	IMPLICIT NONE
    
	INTEGER ,INTENT(IN)				:: elid
	REAL(KIND=DBL),INTENT(IN) :: xele(2,3),uk(3),ey,pv
    REAL(KIND=DBL),INTENT(OUT):: SIF(3)

	! local variables
	INTEGER :: i,k,nve
	REAL(KIND=DBL) :: xitip,dxds(2),ps(3),dps(3),&
                          eita(2),fa(2),fk(3),alph,beta,dj

	! xitip: the local coordinates for crack front (node 1)
    if (elid.eq.CTIP1) then
		xitip = -1.0d0
    elseif (elid.eq.CTIP2) then
    	xitip = 1.0d0
    else
      	print*,'CompSIF: wrong element id for tip element'
        stop
    endif

    ! initialize fk(k,i)
    fk = 0.d0
    dxds = 0.d0
    !...compute shape functions at node 1 (for type 1 of tip element) or node 2 (for type 2)
    CALL reg_shape(elid,xitip,ps,dps)
	!...compute derivative of position vector r
    do i=1,NODE(elid)
      	do k=1,2
        	dxds(k) = dxds(k) + xele(k,i)*dps(i)
        enddo
    enddo
	dj = SQRT(dxds(1)*dxds(1)+dxds(2)*dxds(2))
	!...compute the UNIT tangential vector eita(2) (1-axis of local coordinates)
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
	!...compute ti(2) (for 2D ti(2) is the same as eita(2)

    !...compute displacement at node 1 or node 2 in LOCAL coordinate system
    do k = 1,2
    	fk(1) = fk(1) + fa(k)*uk(k)
    	fk(2) = fk(2) + eita(k)*uk(k)
    enddo
    fk(3) = uk(3)
    !...compute stress intensity factors
   	alph = 4.0d0*dsqrt(2.0d0/PI)*(1.0d0-pv*pv)/ey
   	beta = alph/(1.0d0-pv)
    SIF(1) = 2.0d0/alph*fk(1)/SQRT(dj)
    SIF(2) = 2.0d0/alph*fk(2)/SQRT(dj)
    SIF(3) = 2.0d0/beta*fk(3)/SQRT(dj)

END SUBROUTINE CompSIF
