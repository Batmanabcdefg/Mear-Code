Subroutine eat3tip(region,isrc,ifld,xy_src,xy_fld,eat)
!...Subroutine to calculate matrix eat(3,9) for evaluating sigma_u at boundary node which is TIP NODE
!...This subroutine is for tstress3.f90 (integrate over tip only)
!...but it has the same role as eattip.f90 for tstress2.f90 (integrate over the whole crack and only work for embeded crack)
	USE DefinitionConstant
    USE GlobalData
    IMPLICIT NONE

    INTEGER, INTENT(IN)			:: isrc,ifld,region
    REAL(KIND=DBL), INTENT(IN)	:: xy_src(:,:),xy_fld(:,:)
    REAL(KIND=DBL), INTENT(OUT)	:: eat(:,:)

    !...local variables
    REAL(KIND=DBL)	:: eg(NDOF,3*NDOF),eh(NDOF,3*NDOF)

    if (isrc.eq.ifld) then
      	!...coincident elements
       	call eg3tip_ii(region,elemid(isrc),eg,ng_in,ng_ln,xi(:,ng_in),xi_log(:,ng_ln),&
            				wi(:,ng_in),wi_log(:,ng_ln),xy_src)
        call eh3tip(elemid(isrc),elemid(ifld),eh,ng_in,xi(:,ng_in),wi(:,ng_in),xy_src,xy_fld)
    else
      	!...separated elements: can have case of field element is another tip element
       	call eg3tip_ij(region,elemid(isrc),elemid(ifld),eg,ng_in,xi(:,ng_in),&
            				wi(:,ng_in),xy_src,xy_fld)
       	call eh3tip(elemid(isrc),elemid(ifld),eh,ng_in,xi(:,ng_in),wi(:,ng_in),xy_src,xy_fld)
    endif
    eat = eh - eg
    CONTAINS
	!--------------------------------------------------
    Subroutine eg3tip_ii(region,idelem,eg,nt,ntln,xi,xiln,wi,wiln,xlelem)
    IMPLICIT NONE

    INTEGER, INTENT(IN)			:: idelem,nt,ntln,region
    REAL(KIND=DBL), INTENT(OUT)	:: eg(:,:)
    REAL(KIND=DBL), INTENT(IN)	:: xi(:),xiln(:),wi(:),wiln(:)
    REAL(KIND=DBL), INTENT(IN)	:: xlelem(:,:)

    !...local variables
    INTEGER			:: nei,li
    REAL(KIND=DBL)	:: egbar1(NDOF,3*NDOF),egbar2(NDOF,3*NDOF),egtwobar(NDOF,3*NDOF)
    REAL(KIND=DBL)	:: alpha,zeta,zetapr
    REAL(KIND=DBL)	:: ps(3),dps(3),dpstip(3)
    INTEGER			:: i,j,k,l,ierr
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
        if (idelem.eq.CTIP1) then
          	zetapr = 2.0d0*alpha - 1.0d0
            zeta = 0.5d0*(1.0d0+zetapr)*(1.0d0+zetapr) - 1.0d0
        elseif (idelem.eq.CTIP2) then
        	zetapr = 1.0d0 - 2.0d0*alpha
            zeta = 1.0d0 - 0.5d0*(1.0d0-zetapr)*(1.0d0-zetapr)
        else
          	print*,'egtip_ii(eattip.f90): wrong id of tip elem:',idelem
            stop
        endif
		!...value of (special) shape function derivatives at Gauss point
        call tip_dshape(idelem,zeta,dpstip)
        !...compute value of f(alpha)
        do k = 1,NDOF
          	do i = 1,nei
            	do l = 1,NDOF
                	if (idelem.eq.CTIP1) then
                		f(l,k,i) = GI1(l,k)*dpstip(i)*(1.0d0+zetapr)
                    elseif (idelem.eq.CTIP2) then
                    	f(l,k,i) = GI1(l,k)*dpstip(i)*(1.0d0-zetapr)
                    endif
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
    	if (idelem.eq.CTIP1) then
        	!...tip node is (local) node 1
    		source(k) = xlelem(k,1)
        elseif (idelem.eq.CTIP2) then
        	!...tip node is (local) node 2
        	source(k) = xlelem(k,2)
        else
          	print*,'egtip_ii: wrong id of tip elem'
            stop
        endif
    enddo
    !...loop over inner (tip) element
    do li=1,nt
       	!...Gauss point of regular Gauss integration
        zetapr = xi(li)
        !...transform back to original zeta
        if (idelem.eq.CTIP1) then
            zeta = 0.5d0*(1.0d0+zetapr)*(1.0d0+zetapr) - 1.0d0
        elseif (idelem.eq.CTIP2) then
            zeta = 1.0d0 - 0.5d0*(1.0d0-zetapr)*(1.0d0-zetapr)
        else
          	print*,'egtip_ii: wrong id of tip elem:',idelem
            stop
        endif
        !...shape function & its derivatives at Gauss point
        call reg_shape(idelem,zeta,ps,dps)
        call tip_dshape(idelem,zeta,dpstip)
        !...compute derivative of position vector of (tip) element
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
           	print*,'egtip_ii: r is too small'
            stop
        endif
        !...compute kernels
        call Kernel2(NDOF,region,source,field,C2,G2,U2)
       	do j=1,nei
           	do k=1,NDOF
               	do l=1,NDOF
                  	if (idelem.eq.CTIP1) then
                   		egbar1(k,NDOF*(j-1)+l) = egbar1(k,NDOF*(j-1)+l) + &
                        	GI1(l,k)*dpstip(j)*(1.0d0+zetapr)*dlog(4.0d0*r/((zetapr+1.0d0)*(zetapr+1.0d0)))*wi(li)
                        egtwobar(k,NDOF*(j-1)+l) = egtwobar(k,NDOF*(j-1)+l) + &
                        	G2(l,k)*dpstip(j)*(1.0d0+zetapr)*wi(li)
                    elseif (idelem.eq.CTIP2) then
                   		egbar1(k,NDOF*(j-1)+l) = egbar1(k,NDOF*(j-1)+l) + &
                        	GI1(l,k)*dpstip(j)*(1.0d0-zetapr)*dlog(4.0d0*r/((1.0d0-zetapr)*(1.0d0-zetapr)))*wi(li)
                        egtwobar(k,NDOF*(j-1)+l) = egtwobar(k,NDOF*(j-1)+l) + &
                        	G2(l,k)*dpstip(j)*(1.0d0-zetapr)*wi(li)
                    endif
                    
                enddo
            enddo
        enddo
    enddo !over loop of inner integral
    eg = egbar1 + 2.0d0*egbar2 + egtwobar
	End Subroutine eg3tip_ii
    !----------------------------------------------------------------------
    Subroutine eg3tip_ij(region,idtip,idi,eg,nt,xi,wi,xy_tip,xy_fld)
    !...Note: this subroutine includes the case of inner element is a tip element
    !...In that case, special shape functions are used for inner integral
    IMPLICIT NONE

    INTEGER, INTENT(IN)			:: idtip,idi,nt,region
    REAL(KIND=DBL), INTENT(OUT)	:: eg(:,:)
    REAL(KIND=DBL), INTENT(IN)	:: xi(:),wi(:)
    REAL(KIND=DBL), INTENT(IN)	:: xy_tip(:,:),xy_fld(:,:)

    !...local variables
    INTEGER			:: nei,li
    REAL(KIND=DBL)	:: egbar(NDOF,3*NDOF),egtwobar(NDOF,3*NDOF)
    REAL(KIND=DBL)	:: zeta,zetapr
    REAL(KIND=DBL)	:: ps(3),dps(3),dpstip(3)
    INTEGER			:: i,j,k,l,ierr
    REAL(KIND=DBL)	:: C2(NDOF,NDOF),G2(NDOF,NDOF),U2(NDOF,NDOF),GI1(NDOF,NDOF)
    REAL(KIND=DBL)	:: source(2),field(2),temp(2)
    REAL(KIND=DBL)	:: dxds(2),r

    !...get number of nodes of inner element
    nei = NODE(idi)

	!...compute egbar and egtwobar (regular integrals)
    egbar = 0.0d0
    egtwobar = 0.0d0
	!...get coordinates of tip node
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
           	print*,'egtip_ij: r is too small'
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
	End Subroutine eg3tip_ij
    !----------------------------------------------------------------------
    Subroutine eh3tip(idtip,idi,eh,nt,xi,wi,xy_tip,xy_fld)
    !...Note: this subroutine is for all cases: coincident and separated elements
    IMPLICIT NONE

    INTEGER, INTENT(IN)			:: idtip,idi,nt
    REAL(KIND=DBL), INTENT(OUT)	:: eh(:,:)
    REAL(KIND=DBL), INTENT(IN)	:: xi(:),wi(:)
    REAL(KIND=DBL), INTENT(IN)	:: xy_tip(:,:),xy_fld(:,:)

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
    	if (idtip.eq.CTIP1) then
        	!...tip node is (local) node 1
    		source(k) = xy_tip(k,1)
        elseif (idtip.eq.CTIP2) then
        	!...tip node is (local) node 2
        	source(k) = xy_tip(k,2)
        else
          	print*,'ehtip(ebtip.f90): wrong id of (outer) tip elem:',idtip
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
        	call tip_shape(idi,zeta,pstip)
        endif
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
           	print*,'ehtip: r is too small'
            stop
        endif
        if (jacob.le.SMALL_NUM) then
        	print*,'ehtip: jacobi is too small'
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
	End Subroutine eh3tip
End Subroutine eat3tip
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
Subroutine eat3oth(region,isrc,ifld,xy_src,xy_fld,eat)
!...Subroutine to calculate matrix eat(3,9) for evaluating sigma_u at BOUNDARY NODE of tip element
!...This subroutine is for tstress3.f90 (integrate over tip only)
!...but it has the same role as eattip.f90 for tstress2.f90 (integrate over the whole crack and only work for embeded crack)
	USE DefinitionConstant
    USE GlobalData
    IMPLICIT NONE

    INTEGER, INTENT(IN)			:: isrc,ifld,region
    REAL(KIND=DBL), INTENT(IN)	:: xy_src(:,:),xy_fld(:,:)
    REAL(KIND=DBL), INTENT(OUT)	:: eat(:,:)

    !...local variables
    REAL(KIND=DBL)	:: eg(NDOF,3*NDOF),eh(NDOF,3*NDOF)
    INTEGER		:: adjacent_type

    if (isrc.eq.ifld) then
      	!...coincident elements
       	call eg3oth_ii(region,elemid(isrc),eg,ng_in,ng_ln,xi(:,ng_in),xi_log(:,ng_ln),&
            				wi(:,ng_in),wi_log(:,ng_ln),xy_src)
        call eh3oth(elemid(isrc),elemid(ifld),eh,ng_in,xi(:,ng_in),wi(:,ng_in),xy_src,xy_fld)
    elseif (adjacent(isrc,ifld).gt.0) then
    	!...adjacent element
       	if (elnode(1,isrc).eq.elnode(2,ifld)) then
          	!...singular at node 2 of field element
            if (elemid(isrc).ne.CTIP2) then
              	print*,'eat3oth: error for connectivity of element (system#):',isrc
                stop
            endif
           	adjacent_type = 2
        elseif (elnode(2,isrc).eq.elnode(1,ifld)) then
           	!...singular at node 1 of field element
            if (elemid(isrc).ne.CTIP1) then
              	print*,'eat3oth: error for connectivity of element (system#):',isrc
                stop
            endif
           	adjacent_type = 1
        else
          	print*,'eat3oth: for adjacent elements, common node must be different in (local) number'
            stop
        endif
       	call eg3oth_adt(region,elemid(ifld),eg,ng_in,ng_ln,xi(:,ng_in),xi_log(:,ng_ln),&
        				wi(:,ng_in),wi_log(:,ng_ln),xy_fld,adjacent_type)
        call eh3oth(elemid(isrc),elemid(ifld),eh,ng_in,xi(:,ng_in),wi(:,ng_in),xy_src,xy_fld)
    else
      	!...separated elements: can have case of field element is another tip element
       	call eg3oth_ij(region,elemid(isrc),elemid(ifld),eg,ng_in,xi(:,ng_in),&
            				wi(:,ng_in),xy_src,xy_fld)
       	call eh3oth(elemid(isrc),elemid(ifld),eh,ng_in,xi(:,ng_in),wi(:,ng_in),xy_src,xy_fld)
    endif
    eat = eh - eg
    CONTAINS
	!--------------------------------------------------
    Subroutine eg3oth_ii(region,idelem,eg,nt,ntln,xi,xiln,wi,wiln,xlelem)
    IMPLICIT NONE

    INTEGER, INTENT(IN)			:: idelem,nt,ntln,region
    REAL(KIND=DBL), INTENT(OUT)	:: eg(:,:)
    REAL(KIND=DBL), INTENT(IN)	:: xi(:),xiln(:),wi(:),wiln(:)
    REAL(KIND=DBL), INTENT(IN)	:: xlelem(:,:)

    !...local variables
    INTEGER			:: nei,li
    REAL(KIND=DBL)	:: egbar1(NDOF,3*NDOF),egbar2(NDOF,3*NDOF),egtwobar(NDOF,3*NDOF)
    REAL(KIND=DBL)	:: alpha,zeta,zetapr
    REAL(KIND=DBL)	:: ps(3),dps(3),dpstip(3)
    INTEGER			:: i,j,k,l,ierr
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
        if (idelem.eq.CTIP1) then
          	!...singular at zeta = 1
            zetapr = 1.0d0 - 2.0d0*alpha
            zeta = 0.5d0*(1.0d0+zetapr)*(1.0d0+zetapr) - 1.0d0
        elseif (idelem.eq.CTIP2) then
        	zetapr = 2.0d0*alpha - 1.0d0
            zeta = 1.0d0 - 0.5d0*(1.0d0-zetapr)*(1.0d0-zetapr)
        else
          	print*,'eg3oth_ii: wrong id of tip elem:',idelem
            stop
        endif
		!...value of (special) shape function derivatives at Gauss point
        call tip_dshape(idelem,zeta,dpstip)
        !...compute value of f(alpha)
        do k = 1,NDOF
          	do i = 1,nei
            	do l = 1,NDOF
                	if (idelem.eq.CTIP1) then
                		f(l,k,i) = GI1(l,k)*dpstip(i)*(1.0d0+zetapr)
                    elseif (idelem.eq.CTIP2) then
                    	f(l,k,i) = GI1(l,k)*dpstip(i)*(1.0d0-zetapr)
                    endif
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
    egbar2 = -4.0d0*egbar2
    !----------------------------------------------------------------------
    !...compute egbar1 and egtwobar (regular integrals)
    egbar1 = 0.0d0
    egtwobar = 0.0d0
	!...get coordinates of boundary node
    do k = 1,2
    	if (idelem.eq.CTIP1) then
        	!...tip node is (local) node 1, thus boundary node is local node 2
    		source(k) = xlelem(k,2)
        elseif (idelem.eq.CTIP2) then
        	!...tip node is (local) node 2, thus boundary node is local node 1
        	source(k) = xlelem(k,1)
        else
          	print*,'eg3oth_ii: wrong id of tip elem'
            stop
        endif
    enddo
    !...loop over inner (tip) element
    do li=1,nt
       	!...Gauss point of regular Gauss integration
        zetapr = xi(li)
        !...transform back to original zeta
        if (idelem.eq.CTIP1) then
            zeta = 0.5d0*(1.0d0+zetapr)*(1.0d0+zetapr) - 1.0d0
        elseif (idelem.eq.CTIP2) then
            zeta = 1.0d0 - 0.5d0*(1.0d0-zetapr)*(1.0d0-zetapr)
        else
          	print*,'eg3oth_ii: wrong id of tip elem:',idelem
            stop
        endif
        !...shape function & its derivatives at Gauss point
        call reg_shape(idelem,zeta,ps,dps)
        call tip_dshape(idelem,zeta,dpstip)
        !...compute derivative of position vector of (tip) element
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
           	print*,'eg3oth_ii: r is too small'
            stop
        endif
        !...compute kernels
        call Kernel2(NDOF,region,source,field,C2,G2,U2)
       	do j=1,nei
           	do k=1,NDOF
               	do l=1,NDOF
                  	if (idelem.eq.CTIP1) then
                   		egbar1(k,NDOF*(j-1)+l) = egbar1(k,NDOF*(j-1)+l) + &
                        	GI1(l,k)*dpstip(j)*(1.0d0+zetapr)*dlog(4.0d0*r/((1.0d0-zetapr)*(1.0d0-zetapr)))*wi(li)
                        egtwobar(k,NDOF*(j-1)+l) = egtwobar(k,NDOF*(j-1)+l) + &
                        	G2(l,k)*dpstip(j)*(1.0d0+zetapr)*wi(li)
                    elseif (idelem.eq.CTIP2) then
                   		egbar1(k,NDOF*(j-1)+l) = egbar1(k,NDOF*(j-1)+l) + &
                        	GI1(l,k)*dpstip(j)*(1.0d0-zetapr)*dlog(4.0d0*r/((1.0d0+zetapr)*(1.0d0+zetapr)))*wi(li)
                        egtwobar(k,NDOF*(j-1)+l) = egtwobar(k,NDOF*(j-1)+l) + &
                        	G2(l,k)*dpstip(j)*(1.0d0-zetapr)*wi(li)
                    endif
                    
                enddo
            enddo
        enddo
    enddo !over loop of inner integral
    eg = egbar1 + egbar2 + egtwobar
	End Subroutine eg3oth_ii
    !----------------------------------------------------------------------
    Subroutine eg3oth_adt(region,idi,eg,nt,ntln,xi,xiln,wi,wiln,xy_fld,adt_type)
    IMPLICIT NONE

    INTEGER, INTENT(IN)			:: idi,nt,ntln,region,adt_type
    REAL(KIND=DBL), INTENT(OUT)	:: eg(:,:)
    REAL(KIND=DBL), INTENT(IN)	:: xi(:),xiln(:),wi(:),wiln(:)
    REAL(KIND=DBL), INTENT(IN)	:: xy_fld(:,:)

    !...local variables
    INTEGER			:: nei,li
    REAL(KIND=DBL)	:: egbar1(NDOF,3*NDOF),egbar2(NDOF,3*NDOF),egtwobar(NDOF,3*NDOF)
    REAL(KIND=DBL)	:: alpha,zeta,zetapr
    REAL(KIND=DBL)	:: ps(3),dps(3)
    INTEGER			:: i,j,k,l,ierr
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
          	print*,'eg3oth_adt: wrong adjacent type:',adt_type
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
        	!...singular at zeta=-1 (source point is node 1 of field element)
    		source(k) = xy_fld(k,1)
        elseif (adt_type.eq.2) then
        	!...singular at zeta=1 (source point is node 2 of field element)
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
           	print*,'eb3oth_adt: r is too small'
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
    enddo !over of inner integral
    eg = egbar1 + egbar2 + egtwobar
	End Subroutine eg3oth_adt
    !----------------------------------------------------------------------
    Subroutine eg3oth_ij(region,idtip,idi,eg,nt,xi,wi,xy_tip,xy_fld)
    !...Note: this subroutine includes the case of inner element is a tip element
    !(crack has 2 tip elements separately)
    !...In that case, special shape functions are used for inner integral
    IMPLICIT NONE

    INTEGER, INTENT(IN)			:: idtip,idi,nt,region
    REAL(KIND=DBL), INTENT(OUT)	:: eg(:,:)
    REAL(KIND=DBL), INTENT(IN)	:: xi(:),wi(:)
    REAL(KIND=DBL), INTENT(IN)	:: xy_tip(:,:),xy_fld(:,:)

    !...local variables
    INTEGER			:: nei,li
    REAL(KIND=DBL)	:: egbar(NDOF,3*NDOF),egtwobar(NDOF,3*NDOF)
    REAL(KIND=DBL)	:: zeta,zetapr
    REAL(KIND=DBL)	:: ps(3),dps(3),dpstip(3)
    INTEGER			:: i,j,k,l,ierr
    REAL(KIND=DBL)	:: C2(NDOF,NDOF),G2(NDOF,NDOF),U2(NDOF,NDOF),GI1(NDOF,NDOF)
    REAL(KIND=DBL)	:: source(2),field(2),temp(2)
    REAL(KIND=DBL)	:: dxds(2),r

    !...get number of nodes of inner element
    nei = NODE(idi)

	!...compute egbar and egtwobar (regular integrals)
    egbar = 0.0d0
    egtwobar = 0.0d0
	!...get coordinates of tip node
    do k = 1,2
    	if (idtip.eq.CTIP1) then
        	!...tip node is (local) node 1, thus the other boundary is local node 2
    		source(k) = xy_tip(k,2)
        elseif (idtip.eq.CTIP2) then
        	!...tip node is (local) node 2, thus the other boundary is local node 1
        	source(k) = xy_tip(k,1)
        else
          	print*,'eb3oth_ij: wrong id of (outer) tip elem:',idtip
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
           	print*,'eg3oth_ij: r is too small'
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
	End Subroutine eg3oth_ij
    !----------------------------------------------------------------------
    Subroutine eh3oth(idtip,idi,eh,nt,xi,wi,xy_tip,xy_fld)
    !...Note: this subroutine is for all cases: coincident and separated elements
    IMPLICIT NONE

    INTEGER, INTENT(IN)			:: idtip,idi,nt
    REAL(KIND=DBL), INTENT(OUT)	:: eh(:,:)
    REAL(KIND=DBL), INTENT(IN)	:: xi(:),wi(:)
    REAL(KIND=DBL), INTENT(IN)	:: xy_tip(:,:),xy_fld(:,:)

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
    	if (idtip.eq.CTIP1) then
        	!...tip node is (local) node 1, thus the other boundary is local node 2
    		source(k) = xy_tip(k,2)
        elseif (idtip.eq.CTIP2) then
        	!...tip node is (local) node 2, thus the other boundary is local node 1
        	source(k) = xy_tip(k,1)
        else
          	print*,'eh3oth: wrong id of (outer) tip elem:',idtip
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
        	call tip_shape(idi,zeta,pstip)
        endif
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
           	print*,'eh3oth: r is too small'
            stop
        endif
        if (jacob.le.SMALL_NUM) then
        	print*,'eh3oth: jacobi is too small'
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
	End Subroutine eh3oth
End Subroutine eat3oth