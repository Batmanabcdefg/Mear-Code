Subroutine ebts(region,isrc,ifld,xy_src,xy_fld,eb)
!...Subroutine to calculate element matrix eb for T-stress calculating
!...11/4/09: modify ekbts_adt_reg_reg with 2ND approach (improve integration for regular adjacent elements)

	USE DefinitionConstant
	USE GlobalData
    IMPLICIT NONE
    INTEGER, INTENT(IN)			:: isrc,ifld,region
    REAL(KIND=DBL),INTENT(IN)	:: xy_src(:,:),xy_fld(:,:)
    REAL(KIND=DBL),INTENT(OUT)	:: eb(:,:)
    
	!...local variables
    INTEGER						:: adjacent_type

    if (isrc.eq.ifld) then
      	!...coincident elements
        print*,'ebts.f90: ebts is never called for same element:',isrc
        stop
    elseif (adjacent(isrc,ifld).gt.0) then
    	!...adjacent elements, this happens when iscr is SBL
        !...continue to determine which type of adjacent
        if (elnode(1,ifld).eq.elnode(2,isrc)) then
          	adjacent_type = 1
        elseif (elnode(1,isrc).eq.elnode(2,ifld)) then
        	adjacent_type = 2
        !...add 2 cases for purpose of T-stress calculation
        elseif (elnode(1,isrc).eq.elnode(1,ifld)) then
        	!...singular at (eita = -1; zeta = -1)
            adjacent_type = 3
        elseif (elnode(2,isrc).eq.elnode(2,ifld)) then
        	!...sigular at (eita = 1; zeta = 1)
            adjacent_type = 4
        else
          	print*,'ebts.f90: error on connectivity of elements: ',elem_sys2user(isrc),elem_sys2user(ifld)
            stop
        endif
        if (ELTYPE(elemid(isrc)).eq.CTIP.or.ELTYPE(elemid(isrc)).eq.CREGULAR) then
          	!...no difference between tip/regular crack elements since test function
            !...tk uses regular shape function
          	if ((ELTYPE(elemid(ifld)).eq.BREGULAR).or.(ELTYPE(elemid(ifld)).eq.INTERFE)) then
            	!...this is the only case for crack elem adjacent to boundary elem at SBL
          		call ebts_adt_reg_reg(region,elemid(isrc),elemid(ifld),eb,ng_out,ng_ln,ng_out,xi(:,ng_out),xi_log(:,ng_ln),xi(:,ng_out),&
            				wi(:,ng_out),wi_log(:,ng_ln),wi(:,ng_out),xy_src,xy_fld,adjacent_type)
            else
              	print*,'ebts.f90: ebts never called for adjacent crack-crack elems:',isrc,ifld
                stop
            endif
        else
          	print*,'ebts.f90: ebts never called for source elem is boundary elem:',isrc
            stop
        endif
    else
      	!...separated elements
        if (ELTYPE(elemid(isrc)).eq.CTIP.or.ELTYPE(elemid(isrc)).eq.CREGULAR) then
          	if ((ELTYPE(elemid(ifld)).eq.BREGULAR).or.(ELTYPE(elemid(ifld)).eq.INTERFE)) then
          		call ebts_sed_reg_reg(region,elemid(isrc),elemid(ifld),eb,ng_out,ng_out,xi(:,ng_out),xi(:,ng_out),&
            				wi(:,ng_out),wi(:,ng_out),xy_src,xy_fld)
            else
              	print*,'ebts.f90: ebts never called for separated crack-crack elems:',isrc,ifld
                stop
            endif
        else
          	print*,'ebts.f90: ebts never called for source elem is boundary elem:',isrc
            stop
        endif
    endif

	!----------------------------------------------------------------------------------------------------
	CONTAINS
    
    Subroutine ebts_adt_reg_reg(region,ido,idi,ek,nto,ntoln,nti,xio,xioln,xii,wio,wioln,wii,xlo,xli,adt_type)
	!...11/4/09: use 2ND approach to improve integration of regular adjacent elements
    IMPLICIT NONE

    INTEGER, INTENT(IN)		:: ido,idi,nto,ntoln,nti,region
    INTEGER, INTENT(IN)		:: adt_type
    REAL(KIND=DBL), INTENT(OUT)	:: ek(:,:)
    REAL(KIND=DBL), INTENT(IN)	:: xio(:),xioln(:),xii(:),wio(:),wioln(:),wii(:)
    REAL(KIND=DBL), INTENT(IN)	:: xlo(:,:),xli(:,:)

    !...local variables
    INTEGER			:: neo,nei,lo,li
    REAL(KIND=DBL)	:: ekbar1(3*NDOF,3*NDOF),ek2bar(3*NDOF,3*NDOF)
    REAL(KIND=DBL)	:: ekbar21(3*NDOF,3*NDOF),ekbar22(3*NDOF,3*NDOF)
    REAL(KIND=DBL)	:: alpha,beta,eita,zeta,gamma,deltaa,rho,eitabar
    REAL(KIND=DBL)	:: pso(3),dpso(3),psi(3),dpsi(3),dxds(2)
    INTEGER			:: i,j,k,l
    REAL(KIND=DBL)	:: val_inner_integral(NDOF,NDOF,3,3),jacobo,jacobi
    REAL(KIND=DBL)	:: C2(NDOF,NDOF),G2(NDOF,NDOF),U2(NDOF,NDOF),UI1(NDOF,NDOF)
    REAL(KIND=DBL)	:: source(2),field(2)
    REAL(KIND=DBL)	:: temp(2),r

    !...get number of nodes of outer/inner elements
    neo = NODE(ido)
    nei = NODE(idi)

	!...compute matrix ekbar21 and ekbar22
	ekbar21 = 0.0d0
    ekbar22 = 0.0d0
	!...logarith integration: ekbar22 for adt_type = 1/3 and ekbar21 for adt_type = 2/4
	do lo = 1,ntoln
    	!...coordinate of integration point
        gamma = xioln(lo)
		!...initialize
        val_inner_integral = 0.0d0
		!...loop over inner integral (Gauss integration)
		do li =1,nti
			!...coordinate of Gauss point
        	deltaa = xii(li)
        	!...transformation back to original variables
            alpha = gamma*deltaa
            if ((adt_type.eq.1).or.(adt_type.eq.3)) then
            	beta = 1.0d0 - gamma
            elseif ((adt_type.eq.2).or.(adt_type.eq.4)) then
            	beta = gamma - 1.0d0
            else
              	print*,'ebts_adt_reg_reg: wrong adt_type:',adt_type
                stop
            endif
            !...distinguish adjacent types between 2/4 and 1/3
            if ((adt_type.eq.4).or.(adt_type.eq.3)) then
              	eitabar = alpha + beta
                eita = -1.0d0*eitabar
            else
            	eita = alpha + beta
            endif
            zeta = alpha - beta
            !...obtain value of shape function and its derivatives of 
            call reg_shape(ido,eita,pso,dpso)
            call reg_shape(idi,zeta,psi,dpsi)
        	!...jacobian of outer element
            dxds = 0.0d0
            do i = 1,neo
              	do k = 1,2
                	dxds(k) = dxds(k) + xlo(k,i)*dpso(i)
                enddo
            enddo
            jacobo = dsqrt(dxds(1)*dxds(1) + dxds(2)*dxds(2))
            !...jacobian of inner element
            dxds = 0.0d0
            do i = 1,nei
              	do k = 1,2
                	dxds(k) = dxds(k) + xli(k,i)*dpsi(i)
                enddo
            enddo
            jacobi = dsqrt(dxds(1)*dxds(1) + dxds(2)*dxds(2))
            !...calculate val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                	do k=1,NDOF
                    	do l=1,NDOF
                        	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + &
                            UI1(l,k)*dpso(i)*psi(j)*jacobi*wii(li)
                        enddo
                    enddo
                enddo
            enddo
        enddo !of loop over inner integral
        do i=1,neo
          	do j=1,nei
            	do k=1,NDOF
                	do l=1,NDOF
                    	if ((adt_type.eq.2).or.(adt_type.eq.4)) then
                        	ekbar21(NDOF*(i-1)+k,NDOF*(j-1)+l) = ekbar21(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		gamma*val_inner_integral(l,k,j,i)*wioln(lo)
                        elseif ((adt_type.eq.1).or.(adt_type.eq.3)) then
                    		ekbar22(NDOF*(i-1)+k,NDOF*(j-1)+l) = ekbar22(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		gamma*val_inner_integral(l,k,j,i)*wioln(lo)
                        else
                          	print*,'adt_type out of range'
                            stop
                        endif
                    enddo
                enddo
            enddo
        enddo
    enddo !of loop over outer integral
    !...finalize ekbar21/ekbar22
	if ((adt_type.eq.2).or.(adt_type.eq.4)) then
    	ekbar21 = -2.0d0*ekbar21
    elseif ((adt_type.eq.1).or.(adt_type.eq.3)) then
    	ekbar22 = -2.0d0*ekbar22
    endif
    
	!...regular Gauss integration: ekbar21 for adt_type = 1/3 and ekbar22 for adt_type = 2/4
	do lo = 1,nto
    	!...coordinate of integration point
        gamma = xio(lo)
		!...initialize
        val_inner_integral = 0.0d0
		!...loop over inner integral (Gauss integration)
		do li =1,nti
			!...coordinate of Gauss point
        	deltaa = xii(li)
        	!...transformation back to original variables
            if ((adt_type.eq.1).or.(adt_type.eq.3)) then
              	alpha = 0.5d0*(1.0d0 + gamma)*deltaa
            	beta = 0.5d0*(gamma - 1.0d0)
            elseif ((adt_type.eq.2).or.(adt_type.eq.4)) then
            	alpha = 0.5d0*(1.0d0-gamma)*deltaa
            	beta = 0.5d0*(gamma + 1.0d0)
            else
              	print*,'ebts_adt_reg_reg: wrong adt_type:',adt_type
                stop
            endif
            !...distinguish adjacent types between 2/4 and 1/3
            if ((adt_type.eq.4).or.(adt_type.eq.3)) then
              	eitabar = alpha + beta
                eita = -1.0d0*eitabar
            else
            	eita = alpha + beta
            endif
            zeta = alpha - beta
            !...obtain value of shape function and its derivatives of 
            call reg_shape(ido,eita,pso,dpso)
            call reg_shape(idi,zeta,psi,dpsi)
        	!...jacobian of outer element
            dxds = 0.0d0
            do i = 1,neo
              	do k = 1,2
                	dxds(k) = dxds(k) + xlo(k,i)*dpso(i)
                enddo
            enddo
            jacobo = dsqrt(dxds(1)*dxds(1) + dxds(2)*dxds(2))
            !...jacobian of inner element
            dxds = 0.0d0
            do i = 1,nei
              	do k = 1,2
                	dxds(k) = dxds(k) + xli(k,i)*dpsi(i)
                enddo
            enddo
            jacobi = dsqrt(dxds(1)*dxds(1) + dxds(2)*dxds(2))
            !...calculate val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                	do k=1,NDOF
                    	do l=1,NDOF
                        	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + &
                            UI1(l,k)*dpso(i)*psi(j)*jacobi*wii(li)
                        enddo
                    enddo
                enddo
            enddo
        enddo !of loop over inner integral
        do i=1,neo
          	do j=1,nei
            	do k=1,NDOF
                	do l=1,NDOF
                    	if ((adt_type.eq.2).or.(adt_type.eq.4)) then
                        	ekbar22(NDOF*(i-1)+k,NDOF*(j-1)+l) = ekbar22(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		(1.0d0-gamma)*(dlog(3.0d0+gamma)-dlog(2.0d0))*val_inner_integral(l,k,j,i)*wio(lo)
                        elseif ((adt_type.eq.1).or.(adt_type.eq.3)) then
                    		ekbar21(NDOF*(i-1)+k,NDOF*(j-1)+l) = ekbar21(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		(1.0d0+gamma)*(dlog(3.0d0-gamma)-dlog(2.0d0))*val_inner_integral(l,k,j,i)*wio(lo)
                        else
                          	print*,'adt_type out of range'
                            stop
                        endif
                    enddo
                enddo
            enddo
        enddo
    enddo !of loop over outer integral
    !finalize ekbar22/ekbar21
    if ((adt_type.eq.2).or.(adt_type.eq.4)) then
    	ekbar22 = 0.5d0*ekbar22
    elseif ((adt_type.eq.1).or.(adt_type.eq.3)) then
    	ekbar21 = 0.5d0*ekbar21
    endif

	!--------------------------------------------------
    !...compute matrix ekbar1
	ekbar1 = 0.0d0
	!...begins outer integral loop
	do lo = 1,nto
    	!...coordinate of Gauss point
    	if ((adt_type.eq.3).or.(adt_type.eq.4)) then
        	eitabar = xio(lo)
            eita = -1.0d0*eitabar
        elseif ((adt_type.eq.1).or.(adt_type.eq.2)) then
        	eita = xio(lo)
        else
          	print*,'ebts_adt_reg_reg: wrong adt_type'
            stop
        endif
        !...obtain value of shape function and its derivatives
        call reg_shape(ido,eita,pso,dpso)
        !...position of source point (y)
        source = 0.0d0
        do i=1,neo
        	do k=1,2
           		source(k) = source(k) + xlo(k,i)*pso(i)
           	enddo
        enddo
        !...jacobian of outer element
        dxds = 0.0d0
        do i = 1,neo
           	do k = 1,2
               	dxds(k) = dxds(k) + xlo(k,i)*dpso(i)
            enddo
        enddo
        jacobo = dsqrt(dxds(1)*dxds(1) + dxds(2)*dxds(2))
        
		!...initialize
        val_inner_integral = 0.0d0
		!...loop over inner integral (Gauss integration)
		do li =1,nti
			!...coordinate of Gauss point
        	zeta = xii(li)
            !...obtain value of shape function and its derivatives
            call reg_shape(idi,zeta,psi,dpsi)
            !...position of field point (xi)
        	field = 0.0d0
        	do i=1,nei
          		do k=1,2
            		field(k) = field(k) + xli(k,i)*psi(i)
            	enddo
        	enddo
            !...jacobian of inner element
            dxds = 0.0d0
            do i = 1,nei
              	do k = 1,2
                	dxds(k) = dxds(k) + xli(k,i)*dpsi(i)
                enddo
            enddo
            jacobi = dsqrt(dxds(1)*dxds(1) + dxds(2)*dxds(2))
            !...compute the vector r = xi - y
            do k=1,2
              	temp(k) = field(k) - source(k)
            enddo
            r = dsqrt(temp(1)*temp(1) + temp(2)*temp(2))
            !...defensive programming
            if (r.le.SMALL_NUM) then
              	print*,'ebts_adt_reg_reg: r is too small'
                stop
            endif
            !...calculate val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                	do k=1,NDOF
                    	do l=1,NDOF
                        	if (adt_type.eq.1) then
                        		val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + &
                            		UI1(l,k)*dpso(i)*psi(j)*jacobi*dlog(2.0d0*r/(2.0d0-eita+zeta))*wii(li)
                            elseif (adt_type.eq.3) then
                            	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + &
                            		UI1(l,k)*dpso(i)*psi(j)*jacobi*dlog(2.0d0*r/(2.0d0-eitabar+zeta))*wii(li)
                            elseif (adt_type.eq.2) then
                        		val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + &
                            		UI1(l,k)*dpso(i)*psi(j)*jacobi*dlog(2.0d0*r/(2.0d0-zeta+eita))*wii(li)
                            elseif (adt_type.eq.4) then
                            	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + &
                            	UI1(l,k)*dpso(i)*psi(j)*jacobi*dlog(2.0d0*r/(2.0d0-zeta+eitabar))*wii(li)
                            endif
                        enddo
                    enddo
                enddo
            enddo
        enddo !of loop over inner integral
        do i=1,neo
          	do j=1,nei
            	do k=1,NDOF
                	do l=1,NDOF
                    	ekbar1(NDOF*(i-1)+k,NDOF*(j-1)+l) = ekbar1(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		val_inner_integral(l,k,j,i)*wio(lo)
                    enddo
                enddo
            enddo
        enddo
    enddo !of loop over outer integral

    !--------------------------------------------------
    !...evaluate ek2bar (regular integral)
    ek2bar = 0.0d0
    do lo=1,nto
      	!...coordinate of Gauss point
      	eita = xio(lo)  
        !...shape function & its derivatives at Gauss point
        call reg_shape(ido,eita,pso,dpso)
        !...position of source point (y)
        source = 0.d0
        do i=1,neo
          	do k=1,2
            	source(k) = source(k) + xlo(k,i)*pso(i)
            enddo
        enddo
        !...jacobian of outer element
        dxds = 0.0d0
        do i = 1,neo
          	do k = 1,2
            	dxds(k) = dxds(k) + xlo(k,i)*dpso(i)
            enddo
        enddo
        jacobo = dsqrt(dxds(1)*dxds(1) + dxds(2)*dxds(2))
        
        !...loop over inner integral
        val_inner_integral = 0.0d0
        do li=1,nti
          	!...coordinate of Gauss point
            zeta = xii(li)
            !...shape function & its derivatives at Gauss point
        	call reg_shape(idi,zeta,psi,dpsi)
        	!...position of field point (xi)
        	field = 0.0d0
        	do i=1,nei
          		do k=1,2
            		field(k) = field(k) + xli(k,i)*psi(i)
            	enddo
        	enddo
            !...jacobian of inner element
            dxds = 0.0d0
            do i = 1,nei
              	do k = 1,2
                	dxds(k) = dxds(k) + xli(k,i)*dpsi(i)
                enddo
            enddo
            jacobi = dsqrt(dxds(1)*dxds(1) + dxds(2)*dxds(2))
            !...compute kernel
            call Kernel2(NDOF,region,source,field,C2,G2,U2)
			!...compute val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                	do k=1,NDOF
                    	do l=1,NDOF
                            val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + &
                            U2(l,k)*dpso(i)*psi(j)*jacobi*wii(li)
                        enddo
                    enddo
                enddo
            enddo
        enddo !over loop of inner integral
        do i=1,neo
          	do j=1,nei
            	do k=1,NDOF
                	do l=1,NDOF
                        ek2bar(NDOF*(i-1)+k,NDOF*(j-1)+l) = ek2bar(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        val_inner_integral(l,k,j,i)*wio(lo)
                    enddo
                enddo
            enddo
        enddo
    enddo !over loop of outer integral

    !...debugging
    !print*,'printing inside adt_reg_reg'
    !print*,'ekbar1:'
    !do k=1,3*neo
    !	print'(9(e12.6,1x))',(ekbar1(k,l),l=1,9)
    !enddo
    !print*,'ekbar2:'
    !do k=1,3*neo
    !	print'(9(e12.6,1x))',(ekbar2(k,l),l=1,9)
    !enddo
    !print*,'ek2bar:'
    !do k=1,3*neo
    !	print'(9(e12.6,1x))',(ek2bar(k,l),l=1,9)
    !enddo
    !...end of debugging
    
    !...element stiffness matrix
    ek = ekbar1 + ekbar21 + ekbar22 + ek2bar
	End Subroutine ebts_adt_reg_reg

    !----------------------------------------------------------------------------------------------------

	Subroutine ebts_sed_reg_reg(region,ido,idi,ek,nto,nti,xio,xii,wio,wii,xlo,xli)
    IMPLICIT NONE

    INTEGER			:: ido,idi,nto,nti,region
    REAL(KIND=DBL)	:: ek(:,:)
    REAL(KIND=DBL)	:: xio(:),xii(:),wio(:),wii(:)
    REAL(KIND=DBL)	:: xlo(:,:),xli(:,:)

    !...local variables
    INTEGER			:: neo,nei,lo,li
    REAL(KIND=DBL)	:: ekbar(3*NDOF,3*NDOF),ek2bar(3*NDOF,3*NDOF)
    REAL(KIND=DBL)	:: eita,zeta
    REAL(KIND=DBL)	:: pso(3),dpso(3),psi(3),dpsi(3)
    REAL(KIND=DBL)	:: jacobo,jacobi
    INTEGER					:: i,j,k,l,ierr
    REAL(KIND=DBL)	:: val_inner_integral1(NDOF,NDOF,3,3),val_inner_integral2(NDOF,NDOF,3,3)
    REAL(KIND=DBL)	:: C2(NDOF,NDOF),G2(NDOF,NDOF),U2(NDOF,NDOF),UI1(NDOF,NDOF)
    REAL(KIND=DBL)	:: source(2),field(2)
    REAL(KIND=DBL)	:: dxds(2)
    REAL(KIND=DBL)	:: temp(2),r

    !...get number of nodes of outer/inner elements
    neo = NODE(ido)
    nei = NODE(idi)

	!...compute ekbar (regular integral)
	ekbar = 0.0d0
    ek2bar = 0.0d0

	!...begins outer integral loop
	do lo = 1,nto
    	!...coordinate of integration point
		eita = xio(lo)
        !...obtain value of shape function and its derivatives of 
        call reg_shape(ido,eita,pso,dpso)
        !...jacobian of outer element
        dxds = 0.0d0
        do i = 1,neo
          	do k = 1,2
            	dxds(k) = dxds(k) + xlo(k,i)*dpso(i)
            enddo
        enddo
        jacobo = dsqrt(dxds(1)*dxds(1) + dxds(2)*dxds(2))

        !...position of source point (y)
        source = 0.0d0
        do i=1,neo
        	do k=1,2
           		source(k) = source(k) + xlo(k,i)*pso(i)
           	enddo
        enddo
		!...initialize
        val_inner_integral1 = 0.0d0
        val_inner_integral2 = 0.0d0
		!...loop over inner integral (Gauss integration)
		do li =1,nti
			!...coordinate of Gauss point
        	zeta = xii(li)
            !...obtain value of shape function and its derivatives of 
            call reg_shape(idi,zeta,psi,dpsi)
            !...jacobian of inner element
        	dxds = 0.0d0
        	do i = 1,nei
          		do k = 1,2
            		dxds(k) = dxds(k) + xli(k,i)*dpsi(i)
            	enddo
        	enddo
        	jacobi = dsqrt(dxds(1)*dxds(1) + dxds(2)*dxds(2))
            !...position of field point (xi)
        	field = 0.0d0
        	do i=1,nei
          		do k=1,2
            		field(k) = field(k) + xli(k,i)*psi(i)
            	enddo
        	enddo
            !...compute the vector r = xi - y
            do k=1,2
              	temp(k) = field(k) - source(k)
            enddo
            r = dsqrt(temp(1)*temp(1) + temp(2)*temp(2))
            !...defensive programming
            if (r.lt.SMALL_NUM) then
              	print*,'ebts.f90: r is too small for separated elements'
                stop
            endif
            !...compute kernel CI2
            call Kernel2(NDOF,region,source,field,C2,G2,U2)
            !...calculate val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                	do k=1,NDOF
                    	do l=1,NDOF
                        	val_inner_integral1(l,k,j,i) = val_inner_integral1(l,k,j,i) + &
                            	UI1(l,k)*dpso(i)*psi(j)*jacobi*dlog(r)*wii(li)
                            val_inner_integral2(l,k,j,i) = val_inner_integral2(l,k,j,i) + &
                            	U2(l,k)*dpso(i)*psi(j)*jacobi*wii(li)
                        enddo
                    enddo
                enddo
            enddo
        enddo !of loop over inner integral
        do i=1,neo
          	do j=1,nei
            	do k=1,NDOF
                	do l=1,NDOF
                    	ekbar(NDOF*(i-1)+k,NDOF*(j-1)+l) = ekbar(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        	val_inner_integral1(l,k,j,i)*wio(lo)
                        ek2bar(NDOF*(i-1)+k,NDOF*(j-1)+l) = ek2bar(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        	val_inner_integral2(l,k,j,i)*wio(lo)
                    enddo
                enddo
            enddo
        enddo
    enddo !of loop over outer integral
    
    !...element stiffness matrix
    ek = ekbar + ek2bar
	End Subroutine ebts_sed_reg_reg
End Subroutine ebts