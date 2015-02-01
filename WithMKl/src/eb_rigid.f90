Subroutine eb_rigid(region,isrc,ifld,xy_src,xy_fld,local_src_node,eb)
!...Subroutine to calculate element vector ebb at fixed node on inner boundary
!...called by rigidinner1.f90
	USE DefinitionConstant
	USE GlobalData
    IMPLICIT NONE
    INTEGER, INTENT(IN)			:: isrc,ifld,region,local_src_node
    REAL(KIND=DBL),INTENT(IN)	:: xy_src(:,:),xy_fld(:,:)
    REAL(KIND=DBL),INTENT(OUT)	:: eb(:,:)
    !...local variable
    INTEGER						:: adjacent_type

    if (isrc.eq.ifld) then
      	!call ebb_ii(region,elemid(ifld),eb,ng_in,ng_ln,xi(:,ng_in),xi_log(:,ng_ln),&
        !			wi(:,ng_in),wi_log(:,ng_ln),xy_fld,local_src_node)
        !...this case is the same as adjacent elements
        call ebb_adt(region,elemid(ifld),eb,ng_in,ng_ln,xi(:,ng_in),xi_log(:,ng_ln),&
            				wi(:,ng_in),wi_log(:,ng_ln),xy_fld,local_src_node)
    elseif (elnode(1,ifld).eq.elnode(local_src_node,isrc)) then
       	!...source point is node 1 of field element
        adjacent_type = 1
        call ebb_adt(region,elemid(ifld),eb,ng_in,ng_ln,xi(:,ng_in),xi_log(:,ng_ln),&
            				wi(:,ng_in),wi_log(:,ng_ln),xy_fld,adjacent_type)
    elseif (elnode(2,ifld).eq.elnode(local_src_node,isrc)) then
       	!...source point is node 2 of field element
        adjacent_type = 2
        call ebb_adt(region,elemid(ifld),eb,ng_in,ng_ln,xi(:,ng_in),xi_log(:,ng_ln),&
            				wi(:,ng_in),wi_log(:,ng_ln),xy_fld,adjacent_type)
    else
      	!...separated elements
        !(this also includes adjacent element but source point is not a common node of 2 elements)
   		call ebb_ij(region,local_src_node,elemid(ifld),eb,ng_in,xi(:,ng_in),&
            				wi(:,ng_in),xy_src,xy_fld)
    endif
	!--------------------------------------------------------------------------------
	CONTAINS
    Subroutine ebb_adt(region,idi,eb,nt,ntln,xi,xiln,wi,wiln,xy_fld,adt_type)
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
        !...compute derivative of position vector
       	dxds = 0.0d0
       	do i = 1,nei
           	do k = 1,2
               	dxds(k) = dxds(k) + dps(i)*xy_fld(k,i)
            enddo
       	enddo
       	!...compute jacobian
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
        !...compute derivative of position vector
       	dxds = 0.0d0
       	do i = 1,nei
           	do k = 1,2
               	dxds(k) = dxds(k) + dps(i)*xy_fld(k,i)
            enddo
       	enddo
        !...compute jacobian
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
	End Subroutine ebb_adt
    
    !----------------------------------------------------------------------
    Subroutine ebb_ij(region,src_node,idi,eb,nt,xi,wi,xy_src,xy_fld)
    !...Note: there is no case of inner element is a tip element
    IMPLICIT NONE

    INTEGER, INTENT(IN)			:: src_node,idi,nt,region
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
    	if (src_node.eq.1) then
        	!...tip node is (local) node 1
    		source(k) = xy_src(k,1)
        elseif (src_node.eq.2) then
        	!...tip node is (local) node 2
        	source(k) = xy_src(k,2)
        else
          	print*,'ebb_ij(eb_rigid.f90): wrong position of source point'
            stop
        endif
    enddo
    !...loop over inner element
    do li=1,nt
       	!...Gauss point of regular Gauss integration
        if (idi.eq.CTIP1) then
          	print*,'ebb_ij (eb_rigid.f90): inner element is never tip1 element'
            stop
        elseif (idi.eq.CTIP2) then
        	print*,'ebb_ij (eb_rigid.f90): inner element is never tip2 element'
            stop
        else
        	zeta = xi(li)
        endif
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
           	print*,'ebb_ij: r is too small'
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
	End Subroutine ebb_ij
End Subroutine eb_rigid
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
Subroutine eb_rigid2(region,ifld,src_point,xy_fld,eb)
!...Subroutine to calculate element vector ebb at a node INSIDE inner boundary
!...called by rigidinner2.f90 and rigidinner4.f90
	USE DefinitionConstant
	USE GlobalData
    IMPLICIT NONE
    INTEGER, INTENT(IN)			:: ifld,region
    REAL(KIND=DBL),INTENT(IN)	:: src_point(:),xy_fld(:,:)
    REAL(KIND=DBL),INTENT(OUT)	:: eb(:,:)
  	!...separated elements
	call ebb_ij2(region,elemid(ifld),eb,ng_in,xi(:,ng_in),&
            				wi(:,ng_in),src_point,xy_fld)

	CONTAINS

    Subroutine ebb_ij2(region,idi,eb,nt,xi,wi,source,xy_fld)
    !...Note: there is no case of inner element is a tip element
    IMPLICIT NONE

    INTEGER, INTENT(IN)			:: idi,nt,region
    REAL(KIND=DBL), INTENT(OUT)	:: eb(:,:)
    REAL(KIND=DBL), INTENT(IN)	:: xi(:),wi(:)
    REAL(KIND=DBL), INTENT(IN)	:: source(:),xy_fld(:,:)

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
        if (idi.eq.CTIP1) then
          	print*,'ebb_ij2 (eb_rigid.f90): inner element is never tip1 element'
            stop
        elseif (idi.eq.CTIP2) then
        	print*,'ebb_ij2 (eb_rigid.f90): inner element is never tip2 element'
            stop
        else
        	zeta = xi(li)
        endif
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
           	print*,'ebb_ij2: r is too small'
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
	End Subroutine ebb_ij2
End Subroutine eb_rigid2