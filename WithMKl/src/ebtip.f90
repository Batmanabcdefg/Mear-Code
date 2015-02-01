Subroutine ebtip(region,isrc,ifld,xy_src,xy_fld,eb)
!...Subroutine to calculate element vector eb3 at TIP node
	USE DefinitionConstant
	USE GlobalData
    IMPLICIT NONE
    INTEGER, INTENT(IN)			:: isrc,ifld,region
    REAL(KIND=DBL),INTENT(IN)	:: xy_src(:,:),xy_fld(:,:)
    REAL(KIND=DBL),INTENT(OUT)	:: eb(:,:)
    !...local variable

    if (isrc.eq.ifld) then
      	print*,'ebtip.f90: ebtip is never called for same element'
        stop
    else
      	!...separated elements: even with adjacent element, but source point is on tip, then it is separated elements
        !...obtain coordinates of source point y where sigma_u is evaluated
   		call ebtip_ij(region,elemid(isrc),elemid(ifld),eb,ng_in,xi(:,ng_in),&
            				wi(:,ng_in),xy_src,xy_fld)
    endif
	!--------------------------------------------------------------------------------
	CONTAINS
	Subroutine ebtip_ij(region,idtip,idi,eb,nt,xi,wi,xy_tip,xy_fld)
    !...Note: field element cannot be a tip element
    IMPLICIT NONE

    INTEGER, INTENT(IN)			:: idtip,idi,nt,region
    REAL(KIND=DBL), INTENT(OUT)	:: eb(:,:)
    REAL(KIND=DBL), INTENT(IN)	:: xi(:),wi(:)
    REAL(KIND=DBL), INTENT(IN)	:: xy_tip(:,:),xy_fld(:,:)

    !...local variables
    INTEGER			:: nei,li
    REAL(KIND=DBL)	:: ebbar(NDOF,3*NDOF),ebtwobar(NDOF,3*NDOF)
    REAL(KIND=DBL)	:: zeta
    REAL(KIND=DBL)	:: ps(3),dps(3)
    INTEGER			:: i,j,k,l,ierr
    REAL(KIND=DBL)	:: C2(NDOF,NDOF),G2(NDOF,NDOF),U2(NDOF,NDOF),UI1(NDOF,NDOF)
    REAL(KIND=DBL)	:: source(2),field(2),temp(2)
    REAL(KIND=DBL)	:: dxds(2),jacob,r

    !...get number of nodes of inner element
    nei = NODE(idi)

	!...compute ebbar and ebtwobar (regular integrals)
    ebbar = 0.0d0
    ebtwobar = 0.0d0
    !...get coordinates of source point (tip in this case)
    do k = 1,2
    	if (idtip.eq.CTIP1) then
        	!...tip node is (local) node 1
    		source(k) = xy_tip(k,1)
        elseif (idtip.eq.CTIP2) then
        	!...tip node is (local) node 2
        	source(k) = xy_tip(k,2)
        else
          	print*,'ebtip_ij(ebtip.f90): wrong id of (outer) tip elem:',idtip
            stop
        endif
    enddo
    !...loop over inner element
    do li=1,nt
       	!...Gauss point of regular Gauss integration
        zeta = xi(li)
        !...shape function & its derivatives at Gauss point
        call reg_shape(idi,zeta,ps,dps)
        !...compute derivative of position vector of (tip) element
       	dxds = 0.0d0
       	do i = 1,nei
           	do k = 1,2
               	dxds(k) = dxds(k) + dps(i)*xy_fld(k,i)
            enddo
       	enddo
       	!...compute jacobian of (tip) element
       	jacob = dsqrt(dxds(1)*dxds(1) + dxds(2)*dxds(2))
        !...position of field point
        field = 0.0d0
        do i=1,nei
        	do k=1,2
           		field(k) = field(k) + xy_fld(k,i)*ps(i)
           	enddo
        enddo
        !...compute distance r between source and field
        do k=1,2
           	temp(k) = field(k) - source(k)
        enddo
        r = dsqrt(temp(1)*temp(1) + temp(2)*temp(2))
		!...defensive programming
        if (r.le.SMALL_NUM) then
           	print*,'ebtip_ij: r is too small'
            stop
        endif
        !...compute kernels
        call Kernel2(NDOF,region,source,field,C2,G2,U2)
       	do j=1,nei
           	do k=1,NDOF
               	do l=1,NDOF
               		ebbar(k,NDOF*(j-1)+l) = ebbar(k,NDOF*(j-1)+l) + &
                        	UI1(l,k)*ps(j)*jacob*dlog(r)*wi(li)
                    ebtwobar(k,NDOF*(j-1)+l) = ebtwobar(k,NDOF*(j-1)+l) + &
                        U2(l,k)*ps(j)*jacob*wi(li)
                enddo
            enddo
        enddo
    enddo !over loop of inner integral
    eb = ebbar + ebtwobar
	End Subroutine ebtip_ij
End Subroutine ebtip
!--------------------------------------------------------------------------------






!--------------------------------------------------------------------------------
Subroutine eboth(region,isrc,ifld,xy_src,xy_fld,eb)
!...Subroutine to calculate element vector eb3 at SBL node
	USE DefinitionConstant
	USE GlobalData
    IMPLICIT NONE
    INTEGER, INTENT(IN)			:: isrc,ifld,region
    REAL(KIND=DBL),INTENT(IN)	:: xy_src(:,:),xy_fld(:,:)
    REAL(KIND=DBL),INTENT(OUT)	:: eb(:,:)
    !...local variable
    INTEGER						:: adjacent_type,SBLnode

	!...get the position of source point (SBL point)
    if (node_id(elnode(1,isrc)).eq.SBLP) then
      	SBLnode = 1
    elseif (node_id(elnode(2,isrc)).eq.SBLP) then
       	SBLnode = 2
    else
       	print*,'eboth is called for not a SBL crack element'
        stop
    endif
    if (isrc.eq.ifld) then
      	print*,'eboth (ebtip.f90) is never called for same element'
        stop
    elseif (adjacent(isrc,ifld).gt.0) then
    	!...defensive programming
        if ((node_id(elnode(1,isrc)).ne.SBLP).and.(node_id(elnode(2,isrc)).ne.SBLP)) then
          	print*,'ebtip.f90: eboth is never called for adjacent element where SOURCE element does not contain SBLP'
            stop
        endif
    	if (node_id(elnode(1,ifld)).eq.SBLP) then
        	!...source point is node 1 of field element
        	adjacent_type = 1
        elseif (node_id(elnode(2,ifld)).eq.SBLP) then
        	!...source point is node 2 of field element
            adjacent_type = 2
        else
          	!...this case only happen when field element is a crack element!
          	print*,'eboth (ebtip.f90) is never called for adjacent element where FIELD element does not contain SBLP'
            stop
        endif
        call eboth_adt(region,elemid(ifld),eb,ng_in,ng_ln,xi(:,ng_in),xi_log(:,ng_ln),&
            				wi(:,ng_in),wi_log(:,ng_ln),xy_fld,adjacent_type)
    else
      	!...separated elements
   		call eboth_ij(region,SBLnode,elemid(ifld),eb,ng_in,xi(:,ng_in),&
            				wi(:,ng_in),xy_src,xy_fld)
    endif
	!--------------------------------------------------------------------------------
	CONTAINS
    Subroutine eboth_adt(region,idi,eb,nt,ntln,xi,xiln,wi,wiln,xy_fld,adt_type)
    IMPLICIT NONE

    INTEGER, INTENT(IN)			:: idi,nt,ntln,region,adt_type
    REAL(KIND=DBL), INTENT(OUT)	:: eb(:,:)
    REAL(KIND=DBL), INTENT(IN)	:: xi(:),xiln(:),wi(:),wiln(:)
    REAL(KIND=DBL), INTENT(IN)	:: xy_fld(:,:)

    !...local variables
    INTEGER			:: nei,li
    REAL(KIND=DBL)	:: ebbar1(NDOF,3*NDOF),ebbar2(NDOF,3*NDOF),ebtwobar(NDOF,3*NDOF)
    REAL(KIND=DBL)	:: alpha,zeta
    REAL(KIND=DBL)	:: ps(3),dps(3)
    INTEGER			:: i,j,k,l,ierr
    REAL(KIND=DBL)	:: f(NDOF,NDOF,3)
    REAL(KIND=DBL)	:: C2(NDOF,NDOF),G2(NDOF,NDOF),U2(NDOF,NDOF),UI1(NDOF,NDOF)
    REAL(KIND=DBL)	:: source(2),field(2),temp(2)
    REAL(KIND=DBL)	:: dxds(2),jacob,r

    !...get number of nodes of inner element
    nei = NODE(idi)
	!...compute (weakly) singular integral ebbar2
	ebbar2 = 0.0d0
	!...begins outer integral loop (logarith integration)
	do li = 1,ntln
    	!...coordinate of integration point
		alpha = xiln(li)
        !...transform back to original variable
        if (adt_type.eq.1) then
          	!...singular at zeta=-1
          	zeta = 2.0d0*alpha - 1.0d0
        elseif (adt_type.eq.2) then
        	!...singular at zeta=1
        	zeta = 1.0d0 - 2.0d0*alpha
        else
          	print*,'eboth_adt(ebtip.f90): wrong adjacent type:',adt_type
            stop
        endif
		!...value of shape function derivatives at Gauss point
        call reg_shape(idi,zeta,ps,dps)
        !...compute derivative of position vector of (tip) element
       	dxds = 0.0d0
       	do i = 1,nei
           	do k = 1,2
               	dxds(k) = dxds(k) + dps(i)*xy_fld(k,i)
            enddo
       	enddo
       	!...compute jacobian of (tip) element
       	jacob = dsqrt(dxds(1)*dxds(1) + dxds(2)*dxds(2))
        !...compute value of f(alpha)
        do k = 1,NDOF
          	do i = 1,nei
            	do l = 1,NDOF
               		f(l,k,i) = UI1(l,k)*ps(i)*jacob
                enddo
            enddo
        enddo
        do k = 1,NDOF
          	do i = 1,nei
            	do l = 1,NDOF
                	ebbar2(k,NDOF*(i-1)+l) = ebbar2(k,NDOF*(i-1)+l) + &
                    	f(l,k,i)*wiln(li)
                enddo
            enddo
        enddo
    enddo !of loop over integral
    !...finally, multiply with 2.0d0 (Note: int(-ln(x)*f(x),x) = sum(w(i)*f(x(i)))
    ebbar2 = -2.0d0*ebbar2
    !----------------------------------------------------------------------
    !...compute ebbar1 and ebtwobar (regular integrals)
    ebbar1 = 0.0d0
    ebtwobar = 0.0d0
	!...get coordinates of source point y
    do k = 1,2
    	if (adt_type.eq.1) then
        	!...singular at zeta=-1 (source point is node 1)
    		source(k) = xy_fld(k,1)
        elseif (adt_type.eq.2) then
        	!...singular at zeta=1 (source point is node 2)
        	source(k) = xy_fld(k,2)
        endif
    enddo
    !...loop over inner element
    do li=1,nt
       	!...Gauss point of regular Gauss integration
        zeta = xi(li)
        !...shape function & its derivatives at Gauss point
        call reg_shape(idi,zeta,ps,dps)
        !...compute derivative of position vector of (tip) element
       	dxds = 0.0d0
       	do i = 1,nei
           	do k = 1,2
               	dxds(k) = dxds(k) + dps(i)*xy_fld(k,i)
            enddo
       	enddo
        !...compute jacobian of (tip) element
       	jacob = dsqrt(dxds(1)*dxds(1) + dxds(2)*dxds(2))
        !...position of field point
        field = 0.0d0
        do i=1,nei
        	do k=1,2
           		field(k) = field(k) + xy_fld(k,i)*ps(i)
           	enddo
        enddo
        !...compute distance r between source and field
        !...Note: since the Gauss point zeta never go to -1/1 then r never go to zero
        do k=1,2
           	temp(k) = field(k) - source(k)
        enddo
        r = dsqrt(temp(1)*temp(1) + temp(2)*temp(2))
		!...defensive programming
        if (r.le.SMALL_NUM) then
           	print*,'eboth_adt: r is too small'
            stop
        endif
        !...compute kernels
        call Kernel2(NDOF,region,source,field,C2,G2,U2)
       	do j=1,nei
           	do k=1,NDOF
               	do l=1,NDOF
                  	if (adt_type.eq.1) then
                   		ebbar1(k,NDOF*(j-1)+l) = ebbar1(k,NDOF*(j-1)+l) + &
                        	UI1(l,k)*ps(j)*jacob*dlog(2.0d0*r/(zeta+1.0d0))*wi(li)
                        
                    elseif (adt_type.eq.2) then
                   		ebbar1(k,NDOF*(j-1)+l) = ebbar1(k,NDOF*(j-1)+l) + &
                        	UI1(l,k)*ps(j)*jacob*dlog(2.0d0*r/(1.0d0-zeta))*wi(li)
                    endif
                    ebtwobar(k,NDOF*(j-1)+l) = ebtwobar(k,NDOF*(j-1)+l) + &
                        	U2(l,k)*ps(j)*jacob*wi(li)
                enddo
            enddo
        enddo
    enddo !over loop of inner integral
    eb = ebbar1 + ebbar2 + ebtwobar
	End Subroutine eboth_adt
    !----------------------------------------------------------------------
    Subroutine eboth_ij(region,sbnode,idi,eb,nt,xi,wi,xy_src,xy_fld)
    !...Note: there is no case of inner element is a tip element
    IMPLICIT NONE

    INTEGER, INTENT(IN)			:: sbnode,idi,nt,region
    REAL(KIND=DBL), INTENT(OUT)	:: eb(:,:)
    REAL(KIND=DBL), INTENT(IN)	:: xi(:),wi(:)
    REAL(KIND=DBL), INTENT(IN)	:: xy_src(:,:),xy_fld(:,:)

    !...local variables
    INTEGER			:: nei,li
    REAL(KIND=DBL)	:: ebbar(NDOF,3*NDOF),ebtwobar(NDOF,3*NDOF)
    REAL(KIND=DBL)	:: zeta
    REAL(KIND=DBL)	:: ps(3),dps(3)
    INTEGER			:: i,j,k,l,ierr
    REAL(KIND=DBL)	:: C2(NDOF,NDOF),G2(NDOF,NDOF),U2(NDOF,NDOF),UI1(NDOF,NDOF)
    REAL(KIND=DBL)	:: source(2),field(2),temp(2)
    REAL(KIND=DBL)	:: dxds(2),jacob,r

    !...get number of nodes of inner element
    nei = NODE(idi)

	!...compute ebbar and ebtwobar (regular integrals)
    ebbar = 0.0d0
    ebtwobar = 0.0d0
    !...get coordinates of source point (where y is evaluated)
    do k = 1,2
    	if (sbnode.eq.1) then
        	!...tip node is (local) node 1
    		source(k) = xy_src(k,1)
        elseif (sbnode.eq.2) then
        	!...tip node is (local) node 2
        	source(k) = xy_src(k,2)
        else
          	print*,'eboth_ij(ebtip.f90): wrong position of source point'
            stop
        endif
    enddo
    !...loop over inner element
    do li=1,nt
       	!...Gauss point of regular Gauss integration
        zeta = xi(li)
        !...shape function & its derivatives at Gauss point
        call reg_shape(idi,zeta,ps,dps)
        !...compute derivative of position vector of (tip) element
       	dxds = 0.0d0
       	do i = 1,nei
           	do k = 1,2
               	dxds(k) = dxds(k) + dps(i)*xy_fld(k,i)
            enddo
       	enddo
       	!...compute jacobian of (tip) element
       	jacob = dsqrt(dxds(1)*dxds(1) + dxds(2)*dxds(2))
        !...position of field point
        field = 0.0d0
        do i=1,nei
        	do k=1,2
           		field(k) = field(k) + xy_fld(k,i)*ps(i)
           	enddo
        enddo
        !...compute distance r between source and field
        do k=1,2
           	temp(k) = field(k) - source(k)
        enddo
        r = dsqrt(temp(1)*temp(1) + temp(2)*temp(2))
		!...defensive programming
        if (r.le.SMALL_NUM) then
           	print*,'eboth_ij: r is too small'
            stop
        endif
        !...compute kernels
        call Kernel2(NDOF,region,source,field,C2,G2,U2)
       	do j=1,nei
           	do k=1,NDOF
               	do l=1,NDOF
               		ebbar(k,NDOF*(j-1)+l) = ebbar(k,NDOF*(j-1)+l) + &
                        	UI1(l,k)*ps(j)*jacob*dlog(r)*wi(li)
                    ebtwobar(k,NDOF*(j-1)+l) = ebtwobar(k,NDOF*(j-1)+l) + &
                        U2(l,k)*ps(j)*jacob*wi(li)
                enddo
            enddo
        enddo
    enddo !over loop of inner integral
    eb = ebbar + ebtwobar
	End Subroutine eboth_ij
	!----------------------------------------------------------------------
    
End Subroutine eboth