Subroutine eb3tip(region,isrc,ifld,xy_src,xy_fld,eb)
!...Subroutine to calculate element vector eb(3,9) for evaluated sigma_u at TIP NODE
!...this subroutine is for tstress3.f90 (tstress with integration only on tip element)
	USE DefinitionConstant
	USE GlobalData
    IMPLICIT NONE
    INTEGER, INTENT(IN)			:: isrc,ifld,region
    REAL(KIND=DBL),INTENT(IN)	:: xy_src(:,:),xy_fld(:,:)
    REAL(KIND=DBL),INTENT(OUT)	:: eb(:,:)
    !...local variable
    INTEGER						:: i
    REAL(KIND=DBL)				:: src_pt(2)

    if (isrc.eq.ifld) then
      	print*,'eb3tip.f90: ebtip is never called for same crack element'
        stop
    else
      	!...separated elements: including adjacent element, but source point (y) is always far away from the common
        !boundary of the 2 adjacent elements
        !...obtain coordinates of source point y where sigma_u is evaluated
        do i = 1,2
        	if (elemid(isrc).eq.CTIP1) then
            	src_pt(i) = xy_src(i,1)
            elseif (elemid(isrc).eq.CTIP2) then
            	src_pt(i) = xy_src(i,2)
            else
            	print*,'eb3tip.f90: source element must be a crack tip element'
                stop
            endif
        enddo
   		call eb3tip_ij(region,elemid(ifld),eb,ng_in,xi(:,ng_in),&
            				wi(:,ng_in),src_pt,xy_fld)
    endif
	!--------------------------------------------------------------------------------
	CONTAINS
	Subroutine eb3tip_ij(region,idi,eb,nt,xi,wi,source,xy_fld)
    !...this is different from before: now coordinates of source (point where the value of sigma_u
    !is evaluated) become the argument of subroutine, and id of source element is not needed anymore
    IMPLICIT NONE

    INTEGER, INTENT(IN)			:: idi,nt,region
    REAL(KIND=DBL), INTENT(OUT)	:: eb(:,:)
    REAL(KIND=DBL), INTENT(IN)	:: xi(:),wi(:)
    REAL(KIND=DBL), INTENT(IN)	:: source(2),xy_fld(:,:)

    !...local variables
    INTEGER			:: nei,li
    REAL(KIND=DBL)	:: ebbar(NDOF,3*NDOF),ebtwobar(NDOF,3*NDOF)
    REAL(KIND=DBL)	:: zeta
    REAL(KIND=DBL)	:: ps(3),dps(3)
    INTEGER			:: i,j,k,l,ierr
    REAL(KIND=DBL)	:: C2(NDOF,NDOF),G2(NDOF,NDOF),U2(NDOF,NDOF),UI1(NDOF,NDOF)
    REAL(KIND=DBL)	:: field(2),temp(2)
    REAL(KIND=DBL)	:: dxds(2),jacob,r

    !...get number of nodes of inner element
    nei = NODE(idi)

	!...compute ebbar and ebtwobar (regular integrals)
    ebbar = 0.0d0
    ebtwobar = 0.0d0
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
	End Subroutine eb3tip_ij
End Subroutine eb3tip
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
Subroutine eb3oth(region,isrc,ifld,xy_src,xy_fld,eb)
!...Subroutine to calculate element vector eb(3,9) for evaluated sigma_u at OTHER boundary node (of tip element)
!...this subroutine is for tstress3.f90 (tstress with integration only on tip element)
	USE DefinitionConstant
	USE GlobalData
    IMPLICIT NONE
    INTEGER, INTENT(IN)			:: isrc,ifld,region
    REAL(KIND=DBL),INTENT(IN)	:: xy_src(:,:),xy_fld(:,:)
    REAL(KIND=DBL),INTENT(OUT)	:: eb(:,:)
    !...local variable
    INTEGER						:: adjacent_type,i
    REAL(KIND=DBL)				:: src_pt(2)

    if (isrc.eq.ifld) then
      	print*,'eb3oth.f90: eb3 is never called for same crack element'
        stop
    elseif (adjacent(isrc,ifld).gt.0) then
    	!...this case happens when tip element is also SBL element
    	!...defensive programming
        if ((node_id(elnode(1,isrc)).ne.SBLP).and.(node_id(elnode(2,isrc)).ne.SBLP)) then
          	print*,'ebtip.f90: ebtip is never called for adjacent element where SOURCE element does not contain SBLP'
            stop
        endif
    	if (node_id(elnode(1,ifld)).eq.SBLP) then
        	!...singular at node 1 of field element
        	adjacent_type = 1
        elseif (node_id(elnode(2,ifld)).eq.SBLP) then
        	!...singular at node 2 of field element
            adjacent_type = 2
        else
          	print*,'ebtip.f90: ebtip is never called for adjacent element where FIELD element does not contain SBLP'
            stop
        endif
        call eb3oth_adt(region,elemid(ifld),eb,ng_in,ng_ln,xi(:,ng_in),xi_log(:,ng_ln),&
            				wi(:,ng_in),wi_log(:,ng_ln),xy_fld,adjacent_type)
    else
      	!...separated elements
        !...obtain coordinates of source point y where sigma_u is evaluated
        do i = 1,2
        	if (elemid(isrc).eq.CTIP1) then
            	src_pt(i) = xy_src(i,1)
            elseif (elemid(isrc).eq.CTIP2) then
            	src_pt(i) = xy_src(i,2)
            else
            	print*,'eb3oth.f90: source element must be a crack tip element'
                stop
            endif
        enddo
   		call eb3oth_ij(region,elemid(ifld),eb,ng_in,xi(:,ng_in),&
            				wi(:,ng_in),src_pt,xy_fld)
    endif
	!--------------------------------------------------------------------------------
	CONTAINS
	Subroutine eb3oth_ij(region,idi,eb,nt,xi,wi,source,xy_fld)
    !...this subroutine is exactly the same as eb3tip_ij
    IMPLICIT NONE

    INTEGER, INTENT(IN)			:: idi,nt,region
    REAL(KIND=DBL), INTENT(OUT)	:: eb(:,:)
    REAL(KIND=DBL), INTENT(IN)	:: xi(:),wi(:)
    REAL(KIND=DBL), INTENT(IN)	:: source(2),xy_fld(:,:)

    !...local variables
    INTEGER			:: nei,li
    REAL(KIND=DBL)	:: ebbar(NDOF,3*NDOF),ebtwobar(NDOF,3*NDOF)
    REAL(KIND=DBL)	:: zeta
    REAL(KIND=DBL)	:: ps(3),dps(3)
    INTEGER			:: i,j,k,l,ierr
    REAL(KIND=DBL)	:: C2(NDOF,NDOF),G2(NDOF,NDOF),U2(NDOF,NDOF),UI1(NDOF,NDOF)
    REAL(KIND=DBL)	:: field(2),temp(2)
    REAL(KIND=DBL)	:: dxds(2),jacob,r

    !...get number of nodes of inner element
    nei = NODE(idi)

	!...compute ebbar and ebtwobar (regular integrals)
    ebbar = 0.0d0
    ebtwobar = 0.0d0
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
           	print*,'eb3oth_ij: r is too small'
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
	End Subroutine eb3oth_ij
	!----------------------------------------------------------------------
    Subroutine eb3oth_adt(region,idi,eb,nt,ntln,xi,xiln,wi,wiln,xy_fld,adt_type)
    IMPLICIT NONE

    INTEGER, INTENT(IN)			:: idi,nt,ntln,region,adt_type
    REAL(KIND=DBL), INTENT(OUT)	:: eb(:,:)
    REAL(KIND=DBL), INTENT(IN)	:: xi(:),xiln(:),wi(:),wiln(:)
    REAL(KIND=DBL), INTENT(IN)	:: xy_fld(:,:)

    !...local variables
    INTEGER			:: nei,li
    REAL(KIND=DBL)	:: ebbar1(NDOF,3*NDOF),ebbar2(NDOF,3*NDOF),ebtwobar(NDOF,3*NDOF)
    REAL(KIND=DBL)	:: alpha,zeta,zetapr
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
          	print*,'eb3oth_adt: wrong adjacent type:',adt_type
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
           	print*,'eb3oth_adt: r is too small'
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
                        	UI1(l,k)*ps(j)*jacob*dlog(2.0d0*r/(1.0d0-zetapr))*wi(li)
                    endif
                    ebtwobar(k,NDOF*(j-1)+l) = ebtwobar(k,NDOF*(j-1)+l) + &
                        	U2(l,k)*ps(j)*jacob*wi(li)
                enddo
            enddo
        enddo
    enddo !over loop of inner integral
    eb = ebbar1 + ebbar2 + ebtwobar
	End Subroutine eb3oth_adt
End Subroutine eb3oth