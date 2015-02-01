Subroutine eat_rigid(region,isrc,ifld,xy_src,xy_fld,local_src_node,eat)
!...Subroutine to calculate element vector eatt at fixed node on inner boundary
!...called by rigidinner1.f90
	USE DefinitionConstant
    USE GlobalData
    IMPLICIT NONE

    INTEGER, INTENT(IN)			:: isrc,ifld,region,local_src_node
    REAL(KIND=DBL), INTENT(IN)	:: xy_src(:,:),xy_fld(:,:)
    REAL(KIND=DBL), INTENT(OUT)	:: eat(:,:)

    !...local variables
    REAL(KIND=DBL)	:: eg(NDOF,3*NDOF),eh(NDOF,3*NDOF)
    INTEGER			:: adjacent_type

    if (isrc.eq.ifld) then
      	!...coincident elements
       	call egg_ii(region,elemid(isrc),eg,ng_in,ng_ln,xi(:,ng_in),xi_log(:,ng_ln),&
            				wi(:,ng_in),wi_log(:,ng_ln),xy_src,local_src_node)
        call ehh(local_src_node,elemid(ifld),eh,ng_in,xi(:,ng_in),wi(:,ng_in),xy_src,xy_fld)
   	elseif (elnode(1,ifld).eq.elnode(local_src_node,isrc)) then
       	!...source point is node 1 of field element
       	adjacent_type = 1
        call egg_adt(region,elemid(ifld),eg,ng_in,ng_ln,xi(:,ng_in),xi_log(:,ng_ln),&
        				wi(:,ng_in),wi_log(:,ng_ln),xy_fld,adjacent_type)
        call ehh(local_src_node,elemid(ifld),eh,ng_in,xi(:,ng_in),wi(:,ng_in),xy_src,xy_fld)
    elseif (elnode(2,ifld).eq.elnode(local_src_node,isrc)) then
        !...source point is node 2 of field element
        adjacent_type = 2
       	call egg_adt(region,elemid(ifld),eg,ng_in,ng_ln,xi(:,ng_in),xi_log(:,ng_ln),&
        				wi(:,ng_in),wi_log(:,ng_ln),xy_fld,adjacent_type)
        call ehh(local_src_node,elemid(ifld),eh,ng_in,xi(:,ng_in),wi(:,ng_in),xy_src,xy_fld)

    else
      	!...separated elements
        !(this also includes adjacent element but source point is not a common node of 2 elements)
       	call egg_ij(region,local_src_node,elemid(ifld),eg,ng_in,xi(:,ng_in),&
            				wi(:,ng_in),xy_src,xy_fld)
       	call ehh(local_src_node,elemid(ifld),eh,ng_in,xi(:,ng_in),wi(:,ng_in),xy_src,xy_fld)
    endif
    eat = eh - eg
    CONTAINS
	!--------------------------------------------------
    Subroutine egg_ii(region,idelem,eg,nt,ntln,xi,xiln,wi,wiln,xlelem,srcnode)
    IMPLICIT NONE

    INTEGER, INTENT(IN)			:: idelem,nt,ntln,region,srcnode
    REAL(KIND=DBL), INTENT(OUT)	:: eg(:,:)
    REAL(KIND=DBL), INTENT(IN)	:: xi(:),xiln(:),wi(:),wiln(:)
    REAL(KIND=DBL), INTENT(IN)	:: xlelem(:,:)

    !...local variables
    INTEGER			:: nei,li
    REAL(KIND=DBL)	:: egbar1(NDOF,3*NDOF),egbar2(NDOF,3*NDOF),egtwobar(NDOF,3*NDOF)
    REAL(KIND=DBL)	:: alpha,zeta
    REAL(KIND=DBL)	:: ps(3),dps(3)
    INTEGER			:: i,j,k,l
    REAL(KIND=DBL)	:: f(NDOF,NDOF,3)
    REAL(KIND=DBL)	:: C2(NDOF,NDOF),G2(NDOF,NDOF),U2(NDOF,NDOF),GI1(NDOF,NDOF)
    REAL(KIND=DBL)	:: source(2),field(2),temp(2)
    REAL(KIND=DBL)	:: dxds(2),jacob,r

    !...get number of nodes of inner element
    nei = NODE(idelem)

	!...compute (weakly) singular integral egbar2
	egbar2 = 0.0d0
	!...begins outer integral loop (logarith integration)
	do li = 1,ntln
    	!...coordinate of integration point
		alpha = xiln(li)
        !...transform back to original variable
        if (srcnode.eq.1) then
          	zeta = 2.0d0*alpha - 1.0d0
        elseif (srcnode.eq.2) then
        	zeta = 1.0d0 - 2.0d0*alpha
        else
          	print*,'egg_ii(eat_rigid.f90): wrong id source node:',srcnode
            stop
        endif
		!...value of (special) shape function derivatives at Gauss point
        call reg_shape(idelem,zeta,ps,dps)
        !...compute value of f(alpha)
        do k = 1,NDOF
          	do i = 1,nei
            	do l = 1,NDOF
               		f(l,k,i) = GI1(l,k)*dps(i)
                enddo
            enddo
        enddo
        do k = 1,NDOF
          	do i = 1,nei
            	do l = 1,NDOF
                	egbar2(k,NDOF*(i-1)+l) = egbar2(k,NDOF*(i-1)+l) + &
                    	f(l,k,i)*wiln(li)
                enddo
            enddo
        enddo
    enddo !of loop over integral
    !...finally, multiply with 2.0d0 (Note: int(-ln(x)*f(x),x) = sum(w(i)*f(x(i)))
    egbar2 = -2.0d0*egbar2
    
    !----------------------------------------------------------------------
    !...compute egbar1 and egtwobar (regular integrals)
    egbar1 = 0.0d0
    egtwobar = 0.0d0
	!...get coordinates of tip node
    do k = 1,2
    	if (srcnode.eq.1) then
        	!...source node is local node 1
    		source(k) = xlelem(k,1)
        elseif (srcnode.eq.2) then
        	!...source node is (local) node 2
        	source(k) = xlelem(k,2)
        endif
    enddo
    !...loop over inner element
    do li=1,nt
       	!...Gauss point of regular Gauss integration
        zeta = xi(li)
        !...shape function & its derivatives at Gauss point
        call reg_shape(idelem,zeta,ps,dps)
        !...compute derivative of position vector
       	dxds = 0.0d0
       	do i = 1,nei
           	do k = 1,2
               	dxds(k) = dxds(k) + dps(i)*xlelem(k,i)
            enddo
       	enddo
        !...position of field point
        field = 0.0d0
        do i=1,nei
        	do k=1,2
           		field(k) = field(k) + xlelem(k,i)*ps(i)
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
           	print*,'egg_ii: r is too small'
            stop
        endif
        !...compute kernels
        call Kernel2(NDOF,region,source,field,C2,G2,U2)
       	do j=1,nei
           	do k=1,NDOF
               	do l=1,NDOF
                  	if (srcnode.eq.1) then
                    	!singular when zeta = -1
                   		egbar1(k,NDOF*(j-1)+l) = egbar1(k,NDOF*(j-1)+l) + &
                        	GI1(l,k)*dps(j)*dlog(2.0d0*r/(zeta+1.0d0))*wi(li)
                        egtwobar(k,NDOF*(j-1)+l) = egtwobar(k,NDOF*(j-1)+l) + &
                        	G2(l,k)*dps(j)*wi(li)
                    elseif (srcnode.eq.2) then
                    	!singular when zeta = 1
                   		egbar1(k,NDOF*(j-1)+l) = egbar1(k,NDOF*(j-1)+l) + &
                        	GI1(l,k)*dps(j)*dlog(2.0d0*r/(1.0d0-zeta))*wi(li)
                        egtwobar(k,NDOF*(j-1)+l) = egtwobar(k,NDOF*(j-1)+l) + &
                        	G2(l,k)*dps(j)*wi(li)
                    endif
                enddo
            enddo
        enddo
    enddo !over loop of inner integral
    eg = egbar1 + egbar2 + egtwobar
    End Subroutine egg_ii

    !----------------------------------------------------------------------
    Subroutine egg_adt(region,idi,eg,nt,ntln,xi,xiln,wi,wiln,xy_fld,adt_type)
    IMPLICIT NONE

    INTEGER, INTENT(IN)			:: idi,nt,ntln,region,adt_type
    REAL(KIND=DBL), INTENT(OUT)	:: eg(:,:)
    REAL(KIND=DBL), INTENT(IN)	:: xi(:),xiln(:),wi(:),wiln(:)
    REAL(KIND=DBL), INTENT(IN)	:: xy_fld(:,:)

    !...local variables
    INTEGER			:: nei,li
    REAL(KIND=DBL)	:: egbar1(NDOF,3*NDOF),egbar2(NDOF,3*NDOF),egtwobar(NDOF,3*NDOF)
    REAL(KIND=DBL)	:: alpha,zeta
    REAL(KIND=DBL)	:: ps(3),dps(3)
    INTEGER			:: i,j,k,l
    REAL(KIND=DBL)	:: f(NDOF,NDOF,3)
    REAL(KIND=DBL)	:: C2(NDOF,NDOF),G2(NDOF,NDOF),U2(NDOF,NDOF),GI1(NDOF,NDOF)
    REAL(KIND=DBL)	:: source(2),field(2),temp(2)
    REAL(KIND=DBL)	:: r

    !...get number of nodes of inner element
    nei = NODE(idi)
	!...compute (weakly) singular integral ebbar2
	egbar2 = 0.0d0
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
          	print*,'egg_adt(eat_rigid.f90): wrong adjacent type:',adt_type
            stop
        endif
		!...value of shape function derivatives at Gauss point
        call reg_shape(idi,zeta,ps,dps)
        
        !...compute value of f(alpha)
        do k = 1,NDOF
          	do i = 1,nei
            	do l = 1,NDOF
               		f(l,k,i) = GI1(l,k)*dps(i)
                enddo
            enddo
        enddo
        do k = 1,NDOF
          	do i = 1,nei
            	do l = 1,NDOF
                	egbar2(k,NDOF*(i-1)+l) = egbar2(k,NDOF*(i-1)+l) + &
                    	f(l,k,i)*wiln(li)
                enddo
            enddo
        enddo
    enddo !of loop over integral
    !...finally, multiply with 2.0d0 (Note: int(-ln(x)*f(x),x) = sum(w(i)*f(x(i)))
    egbar2 = -2.0d0*egbar2
    !----------------------------------------------------------------------
    !...compute ebbar1 and ebtwobar (regular integrals)
    egbar1 = 0.0d0
    egtwobar = 0.0d0
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
           	print*,'egg_adt: r is too small'
            stop
        endif
        !...compute kernels
        call Kernel2(NDOF,region,source,field,C2,G2,U2)
       	do j=1,nei
           	do k=1,NDOF
               	do l=1,NDOF
                  	if (adt_type.eq.1) then
                   		egbar1(k,NDOF*(j-1)+l) = egbar1(k,NDOF*(j-1)+l) + &
                        	GI1(l,k)*dps(j)*dlog(2.0d0*r/(zeta+1.0d0))*wi(li)
                    elseif (adt_type.eq.2) then
                   		egbar1(k,NDOF*(j-1)+l) = egbar1(k,NDOF*(j-1)+l) + &
                        	GI1(l,k)*dps(j)*dlog(2.0d0*r/(1.0d0-zeta))*wi(li)
                    endif
                    egtwobar(k,NDOF*(j-1)+l) = egtwobar(k,NDOF*(j-1)+l) + &
                        	G2(l,k)*dps(j)*wi(li)
                enddo
            enddo
        enddo
    enddo !over loop of inner integral
    eg = egbar1 + egbar2 + egtwobar
	End Subroutine egg_adt
    !----------------------------------------------------------------------

    Subroutine egg_ij(region,srcnode,idi,eg,nt,xi,wi,xy_src,xy_fld)
    !...Note: this subroutine includes the case of inner element is a tip element
    !...In that case, special shape functions are used for inner integral
    IMPLICIT NONE

    INTEGER, INTENT(IN)			:: srcnode,idi,nt,region
    REAL(KIND=DBL), INTENT(OUT)	:: eg(:,:)
    REAL(KIND=DBL), INTENT(IN)	:: xi(:),wi(:)
    REAL(KIND=DBL), INTENT(IN)	:: xy_src(:,:),xy_fld(:,:)

    !...local variables
    INTEGER			:: nei,li
    REAL(KIND=DBL)	:: egbar(NDOF,3*NDOF),egtwobar(NDOF,3*NDOF)
    REAL(KIND=DBL)	:: zeta,zetapr
    REAL(KIND=DBL)	:: ps(3),dps(3),dpstip(3)
    INTEGER			:: i,j,k,l
    REAL(KIND=DBL)	:: C2(NDOF,NDOF),G2(NDOF,NDOF),U2(NDOF,NDOF),GI1(NDOF,NDOF)
    REAL(KIND=DBL)	:: source(2),field(2),temp(2)
    REAL(KIND=DBL)	:: dxds(2),r

    !...get number of nodes of inner element
    nei = NODE(idi)

	!...compute egbar and egtwobar (regular integrals)
    egbar = 0.0d0
    egtwobar = 0.0d0
	!...get coordinates of source point (where y is evaluated)
    do k = 1,2
    	if (srcnode.eq.1) then
        	!...tip node is (local) node 1
    		source(k) = xy_src(k,1)
        elseif (srcnode.eq.2) then
        	!...tip node is (local) node 2
        	source(k) = xy_src(k,2)
        else
          	print*,'egg_ij(eat_rigid.f90): wrong position of source point'
            stop
        endif
    enddo
    !...loop over inner element
    do li=1,nt
       	!...Gauss point of regular Gauss integration
        if (idi.eq.CTIP1) then
          	zetapr = xi(li)
            zeta = 0.5d0*(1.0d0+zetapr)*(1.0d0+zetapr) - 1.0d0
        elseif (idi.eq.CTIP2) then
        	zetapr = xi(li)
            zeta = 1.0d0 - 0.5d0*(1.0d0-zetapr)*(1.0d0-zetapr)
        else
        	zeta = xi(li)
        endif
        !...shape function & its derivatives at Gauss point
        call reg_shape(idi,zeta,ps,dps)
        if ((idi.eq.CTIP1).or.(idi.eq.CTIP2)) then
        	call tip_dshape(idi,zeta,dpstip)
        endif
        !...compute derivative of position vector of (tip) element
       	dxds = 0.0d0
       	do i = 1,nei
           	do k = 1,2
               	dxds(k) = dxds(k) + dps(i)*xy_fld(k,i)
            enddo
       	enddo
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
           	print*,'egg_ij(eat_rigid.f90): r is too small'
            stop
        endif
        !...compute kernels
        call Kernel2(NDOF,region,source,field,C2,G2,U2)
       	do j=1,nei
           	do k=1,NDOF
               	do l=1,NDOF
                  	if (idi.eq.CTIP1) then
               			egbar(k,NDOF*(j-1)+l) = egbar(k,NDOF*(j-1)+l) + &
                        	GI1(l,k)*dpstip(j)*(1.0d0+zetapr)*dlog(r)*wi(li)
                    	egtwobar(k,NDOF*(j-1)+l) = egtwobar(k,NDOF*(j-1)+l) + &
                        	G2(l,k)*dpstip(j)*(1.0d0+zetapr)*wi(li)
                    elseif (idi.eq.CTIP2) then
                    	egbar(k,NDOF*(j-1)+l) = egbar(k,NDOF*(j-1)+l) + &
                        	GI1(l,k)*dpstip(j)*(1.0d0-zetapr)*dlog(r)*wi(li)
                    	egtwobar(k,NDOF*(j-1)+l) = egtwobar(k,NDOF*(j-1)+l) + &
                        	G2(l,k)*dpstip(j)*(1.0d0-zetapr)*wi(li)
                    else
                      	egbar(k,NDOF*(j-1)+l) = egbar(k,NDOF*(j-1)+l) + &
                        	GI1(l,k)*dps(j)*dlog(r)*wi(li)
                    	egtwobar(k,NDOF*(j-1)+l) = egtwobar(k,NDOF*(j-1)+l) + &
                        	G2(l,k)*dps(j)*wi(li)
                    endif
                enddo
            enddo
        enddo
    enddo !over loop of inner integral
    eg = egbar + egtwobar
	End Subroutine egg_ij
    !----------------------------------------------------------------------
    Subroutine ehh(srcnode,idi,eh,nt,xi,wi,xy_src,xy_fld)
    !...Note: this subroutine is for all cases: coincident and separated elements
    IMPLICIT NONE

    INTEGER, INTENT(IN)			:: srcnode,idi,nt
    REAL(KIND=DBL), INTENT(OUT)	:: eh(:,:)
    REAL(KIND=DBL), INTENT(IN)	:: xi(:),wi(:)
    REAL(KIND=DBL), INTENT(IN)	:: xy_src(:,:),xy_fld(:,:)

    !...local variables
    INTEGER			:: nei,li
    REAL(KIND=DBL)	:: zeta,zetapr
    REAL(KIND=DBL)	:: HI(NDOF,NDOF)
    REAL(KIND=DBL)	:: ps(3),dps(3),pstip(3)
    INTEGER			:: i,j,k,l,ierr,alpha
    REAL(KIND=DBL)	:: source(2),field(2),temp(2)
    REAL(KIND=DBL)	:: dxds(2),jacob,r,normal(2)

    !...get number of nodes of inner element
    nei = NODE(idi)

	!...compute eh (regular integral)
    eh = 0.0d0

	!...get coordinates of tip node
    do k = 1,2
    	if (srcnode.eq.1) then
        	!...source point is node 1
    		source(k) = xy_src(k,1)
        elseif (srcnode.eq.2) then
        	!...source point is node 2
        	source(k) = xy_src(k,2)
        else
          	print*,'ehh(eat_rigid.f90): wrong id of source point'
            stop
        endif
    enddo
    !...loop over inner element
    do li = 1,nt
       	!...Gauss point of regular Gauss integration
        if (idi.eq.CTIP1) then
          	zetapr = xi(li)
            zeta = 0.5d0*(1.0d0+zetapr)*(1.0d0+zetapr) - 1.0d0
        elseif (idi.eq.CTIP2) then
        	zetapr = xi(li)
            zeta = 1.0d0 - 0.5d0*(1.0d0-zetapr)*(1.0d0-zetapr)
        else
        	zeta = xi(li)
        endif
        
        !...shape function & its derivatives at Gauss point
        call reg_shape(idi,zeta,ps,dps)
        if ((idi.eq.CTIP1).or.(idi.eq.CTIP2)) then
        	call tip_shape(idi,zeta,pstip)
        endif
        
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
           	print*,'ehh (eat_rigid.f90): r is too small'
            stop
        endif
        if (jacob.le.SMALL_NUM) then
        	print*,'ehh (eat_rigid.f90): jacobi is too small'
            stop
        endif
        !...compute normal vector n(xi) = s x e3
        temp(1) = dxds(2)
        temp(2) = -dxds(1)
        normal(1) = temp(1)/jacob
        normal(2) = temp(2)/jacob
        !...ready to calculate eh
       	do j=1,nei
           	do k=1,NDOF
               	do l=1,NDOF
                  	!...compute HI(l,k) = H(l,k,alpha).n(alpha)
                    HI = 0.0d0
                    if (k.eq.l) then
                       	do alpha=1,2
                       		HI(l,k) = HI(l,k) - (0.5d0/PI)*(field(alpha)-source(alpha))*normal(alpha)/(r*r)
                        enddo
                    endif
                    if (idi.eq.CTIP1) then
               			eh(k,NDOF*(j-1)+l) = eh(k,NDOF*(j-1)+l) + &
                        	HI(l,k)*pstip(j)*(1.0d0+zetapr)*jacob*wi(li)
                    elseif (idi.eq.CTIP2) then
                    	eh(k,NDOF*(j-1)+l) = eh(k,NDOF*(j-1)+l) + &
                        	HI(l,k)*pstip(j)*(1.0d0-zetapr)*jacob*wi(li)
                    else
                      	eh(k,NDOF*(j-1)+l) = eh(k,NDOF*(j-1)+l) + &
                        	HI(l,k)*ps(j)*jacob*wi(li)
                    endif
                enddo
            enddo
        enddo
    enddo !over loop of inner integral
	End Subroutine ehh
End Subroutine eat_rigid
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
Subroutine eat_rigid2(region,ifld,src_point,xy_fld,eat)
!...Subroutine to calculate element vector eatt at a node INSIDE inner boundary
!...called by rigidinner2.f90 and rigidinner4.f90
	USE DefinitionConstant
    USE GlobalData
    IMPLICIT NONE

    INTEGER, INTENT(IN)			:: ifld,region
    REAL(KIND=DBL), INTENT(IN)	:: xy_fld(:,:),src_point(:)
    REAL(KIND=DBL), INTENT(OUT)	:: eat(:,:)

    !...local variables
    REAL(KIND=DBL)	:: eg(NDOF,3*NDOF),eh(NDOF,3*NDOF)

   	!...separated elements
   	call egg_ij2(region,elemid(ifld),eg,ng_in,xi(:,ng_in),&
            				wi(:,ng_in),src_point,xy_fld)
   	call ehh2(elemid(ifld),eh,ng_in,xi(:,ng_in),wi(:,ng_in),src_point,xy_fld)
    eat = eh - eg
    
    CONTAINS
    Subroutine egg_ij2(region,idi,eg,nt,xi,wi,source,xy_fld)
    !...Note: this subroutine includes the case of inner element is a tip element
    !...In that case, special shape functions are used for inner integral
    IMPLICIT NONE
    INTEGER, INTENT(IN)			:: idi,nt,region
    REAL(KIND=DBL), INTENT(OUT)	:: eg(:,:)
    REAL(KIND=DBL), INTENT(IN)	:: xi(:),wi(:)
    REAL(KIND=DBL), INTENT(IN)	:: source(:),xy_fld(:,:)
    !...local variables
    INTEGER			:: nei,li
    REAL(KIND=DBL)	:: egbar(NDOF,3*NDOF),egtwobar(NDOF,3*NDOF)
    REAL(KIND=DBL)	:: zeta,zetapr
    REAL(KIND=DBL)	:: ps(3),dps(3),dpstip(3)
    INTEGER			:: i,j,k,l
    REAL(KIND=DBL)	:: C2(NDOF,NDOF),G2(NDOF,NDOF),U2(NDOF,NDOF),GI1(NDOF,NDOF)
    REAL(KIND=DBL)	:: field(2),temp(2)
    REAL(KIND=DBL)	:: dxds(2),r

    !...get number of nodes of inner element
    nei = NODE(idi)

	!...compute egbar and egtwobar (regular integrals)
    egbar = 0.0d0
    egtwobar = 0.0d0
    !...loop over inner element
    do li=1,nt
       	!...Gauss point of regular Gauss integration
        if (idi.eq.CTIP1) then
          	zetapr = xi(li)
            zeta = 0.5d0*(1.0d0+zetapr)*(1.0d0+zetapr) - 1.0d0
        elseif (idi.eq.CTIP2) then
        	zetapr = xi(li)
            zeta = 1.0d0 - 0.5d0*(1.0d0-zetapr)*(1.0d0-zetapr)
        else
        	zeta = xi(li)
        endif
        !...shape function & its derivatives at Gauss point
        call reg_shape(idi,zeta,ps,dps)
        if ((idi.eq.CTIP1).or.(idi.eq.CTIP2)) then
        	call tip_dshape(idi,zeta,dpstip)
        endif
        !...compute derivative of position vector of (tip) element
       	dxds = 0.0d0
       	do i = 1,nei
           	do k = 1,2
               	dxds(k) = dxds(k) + dps(i)*xy_fld(k,i)
            enddo
       	enddo
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
           	print*,'egg_ij2(eat_rigid.f90): r is too small'
            stop
        endif
        !...compute kernels
        call Kernel2(NDOF,region,source,field,C2,G2,U2)
       	do j=1,nei
           	do k=1,NDOF
               	do l=1,NDOF
                  	if (idi.eq.CTIP1) then
               			egbar(k,NDOF*(j-1)+l) = egbar(k,NDOF*(j-1)+l) + &
                        	GI1(l,k)*dpstip(j)*(1.0d0+zetapr)*dlog(r)*wi(li)
                    	egtwobar(k,NDOF*(j-1)+l) = egtwobar(k,NDOF*(j-1)+l) + &
                        	G2(l,k)*dpstip(j)*(1.0d0+zetapr)*wi(li)
                    elseif (idi.eq.CTIP2) then
                    	egbar(k,NDOF*(j-1)+l) = egbar(k,NDOF*(j-1)+l) + &
                        	GI1(l,k)*dpstip(j)*(1.0d0-zetapr)*dlog(r)*wi(li)
                    	egtwobar(k,NDOF*(j-1)+l) = egtwobar(k,NDOF*(j-1)+l) + &
                        	G2(l,k)*dpstip(j)*(1.0d0-zetapr)*wi(li)
                    else
                      	egbar(k,NDOF*(j-1)+l) = egbar(k,NDOF*(j-1)+l) + &
                        	GI1(l,k)*dps(j)*dlog(r)*wi(li)
                    	egtwobar(k,NDOF*(j-1)+l) = egtwobar(k,NDOF*(j-1)+l) + &
                        	G2(l,k)*dps(j)*wi(li)
                    endif
                enddo
            enddo
        enddo
    enddo !over loop of inner integral
    eg = egbar + egtwobar
	End Subroutine egg_ij2
    !----------------------------------------------------------------------
    Subroutine ehh2(idi,eh,nt,xi,wi,source,xy_fld)
    !...Note: this subroutine is for all cases: coincident and separated elements
    IMPLICIT NONE

    INTEGER, INTENT(IN)			:: idi,nt
    REAL(KIND=DBL), INTENT(OUT)	:: eh(:,:)
    REAL(KIND=DBL), INTENT(IN)	:: xi(:),wi(:)
    REAL(KIND=DBL), INTENT(IN)	:: source(:),xy_fld(:,:)

    !...local variables
    INTEGER			:: nei,li
    REAL(KIND=DBL)	:: zeta,zetapr
    REAL(KIND=DBL)	:: HI(NDOF,NDOF)
    REAL(KIND=DBL)	:: ps(3),dps(3),pstip(3)
    INTEGER			:: i,j,k,l,ierr,alpha
    REAL(KIND=DBL)	:: field(2),temp(2)
    REAL(KIND=DBL)	:: dxds(2),jacob,r,normal(2)

    !...get number of nodes of inner element
    nei = NODE(idi)

	!...compute eh (regular integral)
    eh = 0.0d0

	!...loop over inner element
    do li = 1,nt
       	!...Gauss point of regular Gauss integration
        if (idi.eq.CTIP1) then
          	zetapr = xi(li)
            zeta = 0.5d0*(1.0d0+zetapr)*(1.0d0+zetapr) - 1.0d0
        elseif (idi.eq.CTIP2) then
        	zetapr = xi(li)
            zeta = 1.0d0 - 0.5d0*(1.0d0-zetapr)*(1.0d0-zetapr)
        else
        	zeta = xi(li)
        endif
        
        !...shape function & its derivatives at Gauss point
        call reg_shape(idi,zeta,ps,dps)
        if ((idi.eq.CTIP1).or.(idi.eq.CTIP2)) then
        	call tip_shape(idi,zeta,pstip)
        endif
        
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
           	print*,'ehh2 (eat_rigid.f90): r is too small'
            stop
        endif
        if (jacob.le.SMALL_NUM) then
        	print*,'ehh2 (eat_rigid.f90): jacobi is too small'
            stop
        endif
        !...compute normal vector n(xi) = s x e3
        temp(1) = dxds(2)
        temp(2) = -dxds(1)
        normal(1) = temp(1)/jacob
        normal(2) = temp(2)/jacob
        !...ready to calculate eh
       	do j=1,nei
           	do k=1,NDOF
               	do l=1,NDOF
                  	!...compute HI(l,k) = H(l,k,alpha).n(alpha)
                    HI = 0.0d0
                    if (k.eq.l) then
                       	do alpha=1,2
                       		HI(l,k) = HI(l,k) - (0.5d0/PI)*(field(alpha)-source(alpha))*normal(alpha)/(r*r)
                        enddo
                    endif
                    if (idi.eq.CTIP1) then
               			eh(k,NDOF*(j-1)+l) = eh(k,NDOF*(j-1)+l) + &
                        	HI(l,k)*pstip(j)*(1.0d0+zetapr)*jacob*wi(li)
                    elseif (idi.eq.CTIP2) then
                    	eh(k,NDOF*(j-1)+l) = eh(k,NDOF*(j-1)+l) + &
                        	HI(l,k)*pstip(j)*(1.0d0-zetapr)*jacob*wi(li)
                    else
                      	eh(k,NDOF*(j-1)+l) = eh(k,NDOF*(j-1)+l) + &
                        	HI(l,k)*ps(j)*jacob*wi(li)
                    endif
                enddo
            enddo
        enddo
    enddo !over loop of inner integral
	End Subroutine ehh2
    !----------------------------------------------------------------------
End Subroutine eat_rigid2