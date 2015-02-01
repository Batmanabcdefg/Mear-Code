Subroutine eatts(region,isrc,ifld,xy_src,xy_fld,eat)
!...Subroutine to calculate element matrix eat for T-stress calculation
!11/4/09: the only difference with eatts.f90: change egts_adt_reg_reg using 2ND approach

	USE DefinitionConstant
    USE GlobalData
    IMPLICIT NONE

    INTEGER, INTENT(IN)			:: isrc,ifld,region
    REAL(KIND=DBL), INTENT(IN)	:: xy_src(:,:),xy_fld(:,:)
    REAL(KIND=DBL), INTENT(OUT)	:: eat(:,:)

    !...local variables
    INTEGER			:: adjacent_type,nno,nni
    REAL(KIND=DBL)	:: eg(3*NDOF,3*NDOF),eh(3*NDOF,3*NDOF)
    LOGICAL			:: debug
!    INTEGER			:: i,j
	
	!debug = .true.
    debug = .false.
    
    if (debug) then
      	print*,'eat called for elements: ',isrc,ifld
    endif
    nno = NODE(elemid(isrc))
    nni = NODE(elemid(ifld))
    
    if (isrc.eq.ifld) then
      	!...coincident elements
        if (ELTYPE(elemid(isrc)).eq.CTIP) then
          	!...even test function uses regular shape functions but in order to make
            !...ln(abs(r/eita-zeta)) be regular we still need to stretch out eita
           	call egts_ii_tip_tip(region,elemid(isrc),eg,ng_out,ng_ln,ng_in,xi(:,ng_out),xi_log(:,ng_ln),xi(:,ng_in),&
            				wi(:,ng_out),wi_log(:,ng_ln),wi(:,ng_in),xy_src)
            !...call eh_reg_tip since we use regular shape functions for outer tip element
            call ehts_reg_tip(elemid(isrc),elemid(isrc),eh,ng_out,ng_in,xi(:,ng_out),xi(:,ng_in),&
            				wi(:,ng_out),wi(:,ng_in),xy_src,xy_src)
            eat = eh - eg
            !debug
            !print*,'egts_ii_tip_tip called'
            !print*,'eg='
            !do i = 1,9
            !  	print'(9(f15.8,1x))',(eg(i,j),j=1,9)
            !enddo
            !print*,'ehts_reg_tip called'
            !print*,'eh='
            !do i = 1,9
            !  	print'(9(f15.8,1x))',(eh(i,j),j=1,9)
            !enddo
        elseif (ELTYPE(elemid(isrc)).eq.CREGULAR) then
        	call egts_ii_reg_reg(region,elemid(isrc),eg,ng_out,ng_ln,ng_in,xi(:,ng_out),xi_log(:,ng_ln),xi(:,ng_in),&
            				wi(:,ng_out),wi_log(:,ng_ln),wi(:,ng_in),xy_src)
            call ehts_reg_reg(elemid(isrc),elemid(isrc),eh,ng_out,ng_in,xi(:,ng_out),xi(:,ng_in),&
            				wi(:,ng_out),wi(:,ng_in),xy_src,xy_src)
            eat = eh - eg
            !debug
            !print*,'egts_ii_reg_reg called'
            !print*,'eg='
            !do i = 1,9
            !  	print'(9(f15.8,1x))',(eg(i,j),j=1,9)
            !enddo
            !print*,'ehts_reg_reg called'
            !print*,'eh='
            !do i = 1,9
            !  	print'(9(f15.8,1x))',(eh(i,j),j=1,9)
            !enddo
        !...4/28/09: add INTERFE taking account for interface elements
        elseif ((ELTYPE(elemid(isrc)).eq.BREGULAR).or.(ELTYPE(elemid(isrc)).eq.INTERFE)) then
        	!...reg-reg: T-stress calculation never goes into this case
            call egts_ii_reg_reg(region,elemid(isrc),eg,ng_out,ng_ln,ng_in,xi(:,ng_out),xi_log(:,ng_ln),xi(:,ng_in),&
            				wi(:,ng_out),wi_log(:,ng_ln),wi(:,ng_in),xy_src)
            !...use different order for eh
            call ehts_reg_reg(elemid(isrc),elemid(isrc),eh,ng_out,ng_in,xi(:,ng_out),xi(:,ng_in),&
            				wi(:,ng_out),wi(:,ng_in),xy_src,xy_src)
            eat = eh - eg
            !debug
            !print*,'egts_ii_reg_reg called'
            !print*,'eg='
            !do i = 1,9
            !  	print'(9(f15.8,1x))',(eg(i,j),j=1,9)
            !enddo
            !print*,'ehts_reg_reg called'
            !print*,'eh='
            !do i = 1,9
            !  	print'(9(f15.8,1x))',(eh(i,j),j=1,9)
            !enddo
        else
          	print*,'eat.f90: current version does not support (same element) with elem type: ',elemid(isrc)
            stop
        endif
    elseif (adjacent(isrc,ifld).gt.0) then
    	!...adjacent elements, continue to determine which type of adjacent
        if (elnode(1,ifld).eq.elnode(2,isrc)) then
          	!...singular at (eita = 1, zeta = -1)
          	adjacent_type = 1
        elseif (elnode(1,isrc).eq.elnode(2,ifld)) then
        	!...singular at (eita = -1, zeta = 1)
        	adjacent_type = 2
        !...add 2 more cases for the case of SBLP
        elseif (elnode(1,isrc).eq.elnode(1,ifld)) then
        	!...singular at (eita = -1; zeta = -1)
            adjacent_type = 3
        elseif (elnode(2,isrc).eq.elnode(2,ifld)) then
        	!...sigular at (eita = 1; zeta = 1)
            adjacent_type = 4
        else
          	print*,'eatts.f90: error on connectivity of elements: ',elem_sys2user(isrc),elem_sys2user(ifld)
            stop
        endif
        if (ELTYPE(elemid(isrc)).eq.CTIP.and.ELTYPE(elemid(ifld)).eq.CTIP) then
          	!...adjacent: tip-tip
           	call egts_adt_reg_tip(region,elemid(isrc),elemid(ifld),eg,ng_out,ng_out,xi(:,ng_out),xi(:,ng_out),&
            				wi(:,ng_out),wi(:,ng_out),xy_src,xy_fld,adjacent_type)
           	call ehts_reg_tip(elemid(isrc),elemid(ifld),eh,ng_out,ng_out,xi(:,ng_out),xi(:,ng_out),&
            				wi(:,ng_out),wi(:,ng_out),xy_src,xy_fld)
           	eat = eh - eg
            !debug
            !print*,'egts_adt_reg_tip called'
            !print*,'eg='
            !do i = 1,9
            !  	print'(9(f15.8,1x))',(eg(i,j),j=1,9)
            !enddo
            !print*,'ehts_reg_tip called'
            !print*,'eh='
            !do i = 1,9
            !  	print'(9(f15.8,1x))',(eh(i,j),j=1,9)
            !enddo
        elseif (ELTYPE(elemid(isrc)).eq.CTIP.and.ELTYPE(elemid(ifld)).ne.CTIP) then
        	!...adjacent tip-reg
           	call egts_adt_reg_reg(region,elemid(isrc),elemid(ifld),eg,ng_out,ng_ln,ng_out,xi(:,ng_out),xi_log(:,ng_ln),xi(:,ng_out),&
            				wi(:,ng_out),wi_log(:,ng_ln),wi(:,ng_out),xy_src,xy_fld,adjacent_type)
          	call ehts_reg_reg(elemid(isrc),elemid(ifld),eh,ng_out,ng_out,xi(:,ng_out),xi(:,ng_out),&
            				wi(:,ng_out),wi(:,ng_out),xy_src,xy_fld)
            eat = eh - eg
            !debug
            !print*,'egts_adt_reg_reg called'
            !print*,'eg='
            !do i = 1,9
            !  	print'(9(f15.8,1x))',(eg(i,j),j=1,9)
            !enddo
            !print*,'ehts_reg_reg called'
            !print*,'eh='
            !do i = 1,9
            !  	print'(9(f15.8,1x))',(eh(i,j),j=1,9)
            !enddo
		elseif (ELTYPE(elemid(isrc)).ne.CTIP.and.ELTYPE(elemid(ifld)).eq.CTIP) then
        	!...adjacent: reg-tip, use same integration order for both eg and eh
            call egts_adt_reg_tip(region,elemid(isrc),elemid(ifld),eg,ng_out,ng_out,xi(:,ng_out),xi(:,ng_out),&
            				wi(:,ng_out),wi(:,ng_out),xy_src,xy_fld,adjacent_type)
            call ehts_reg_tip(elemid(isrc),elemid(ifld),eh,ng_out,ng_out,xi(:,ng_out),xi(:,ng_out),&
            				wi(:,ng_out),wi(:,ng_out),xy_src,xy_fld)
            eat = eh - eg
            !debug
            !print*,'egts_adt_reg_tip called'
            !print*,'eg='
            !do i = 1,9
            !  	print'(9(f15.8,1x))',(eg(i,j),j=1,9)
            !enddo
            !print*,'ehts_reg_tip called'
            !print*,'eh='
            !do i = 1,9
            !  	print'(9(f15.8,1x))',(eh(i,j),j=1,9)
            !enddo
		elseif (ELTYPE(elemid(isrc)).ne.CTIP.and.ELTYPE(elemid(ifld)).ne.CTIP) then
        	!...adjacent: reg-reg, same integration order
            call egts_adt_reg_reg(region,elemid(isrc),elemid(ifld),eg,ng_out,ng_ln,ng_out,xi(:,ng_out),xi_log(:,ng_ln),xi(:,ng_out),&
            				wi(:,ng_out),wi_log(:,ng_ln),wi(:,ng_out),xy_src,xy_fld,adjacent_type)
            call ehts_reg_reg(elemid(isrc),elemid(ifld),eh,ng_out,ng_out,xi(:,ng_out),xi(:,ng_out),&
            				wi(:,ng_out),wi(:,ng_out),xy_src,xy_fld)
            eat = eh - eg
            !debug
            !print*,'egts_adt_reg_reg called'
            !print*,'eg='
            !do i = 1,9
            !  	print'(9(f15.8,1x))',(eg(i,j),j=1,9)
            !enddo
            !print*,'ehts_reg_reg called'
            !print*,'eh='
            !do i = 1,9
            !  	print'(9(f15.8,1x))',(eh(i,j),j=1,9)
            !enddo
        endif
    else
      	!...separated elements
        if (ELTYPE(elemid(isrc)).eq.CTIP.and.ELTYPE(elemid(ifld)).eq.CTIP) then
          	!...separated: tip-tip
           	call egts_sed_reg_tip(region,elemid(isrc),elemid(ifld),eg,ng_out,ng_out,xi(:,ng_out),xi(:,ng_out),&
            				wi(:,ng_out),wi(:,ng_out),xy_src,xy_fld)
          	call ehts_reg_tip(elemid(isrc),elemid(ifld),eh,ng_out,ng_out,xi(:,ng_out),xi(:,ng_out),&
            				wi(:,ng_out),wi(:,ng_out),xy_src,xy_fld)
            eat = eh - eg
            !debug
            !print*,'egts_sed_reg_tip called'
            !print*,'eg='
            !do i = 1,9
            !  	print'(9(f15.8,1x))',(eg(i,j),j=1,9)
            !enddo
            !print*,'ehts_reg_tip called'
            !print*,'eh='
            !do i = 1,9
            !  	print'(9(f15.8,1x))',(eh(i,j),j=1,9)
            !enddo
        elseif (ELTYPE(elemid(isrc)).eq.CTIP.and.ELTYPE(elemid(ifld)).ne.CTIP) then
        	!...separated tip-reg
           	call egts_sed_reg_reg(region,elemid(isrc),elemid(ifld),eg,ng_out,ng_out,xi(:,ng_out),xi(:,ng_out),&
            				wi(:,ng_out),wi(:,ng_out),xy_src,xy_fld)
           	call ehts_reg_reg(elemid(isrc),elemid(ifld),eh,ng_out,ng_out,xi(:,ng_out),xi(:,ng_out),&
            				wi(:,ng_out),wi(:,ng_out),xy_src,xy_fld)
           	eat = eh - eg
            !debug
            !print*,'egts_sed_reg_reg called'
            !print*,'eg='
            !do i = 1,9
            !  	print'(9(f15.8,1x))',(eg(i,j),j=1,9)
            !enddo
            !print*,'ehts_reg_reg called'
            !print*,'eh='
            !do i = 1,9
            !  	print'(9(f15.8,1x))',(eh(i,j),j=1,9)
            !enddo
		elseif (ELTYPE(elemid(isrc)).ne.CTIP.and.ELTYPE(elemid(ifld)).eq.CTIP) then
        	!...separate: reg-tip, same integration order
            call egts_sed_reg_tip(region,elemid(isrc),elemid(ifld),eg,ng_out,ng_out,xi(:,ng_out),xi(:,ng_out),&
            				wi(:,ng_out),wi(:,ng_out),xy_src,xy_fld)
            call ehts_reg_tip(elemid(isrc),elemid(ifld),eh,ng_out,ng_out,xi(:,ng_out),xi(:,ng_out),&
            				wi(:,ng_out),wi(:,ng_out),xy_src,xy_fld)
            eat = eh - eg
            !debug
            !print*,'egts_sed_reg_tip called'
            !print*,'eg='
            !do i = 1,9
            !  	print'(9(f15.8,1x))',(eg(i,j),j=1,9)
            !enddo
            !print*,'ehts_reg_tip called'
            !print*,'eh='
            !do i = 1,9
            !  	print'(9(f15.8,1x))',(eh(i,j),j=1,9)
            !enddo
		elseif (ELTYPE(elemid(isrc)).ne.CTIP.and.ELTYPE(elemid(ifld)).ne.CTIP) then
        	!...separate: reg-reg, same integration order
            call egts_sed_reg_reg(region,elemid(isrc),elemid(ifld),eg,ng_out,ng_out,xi(:,ng_out),xi(:,ng_out),&
            				wi(:,ng_out),wi(:,ng_out),xy_src,xy_fld)
            call ehts_reg_reg(elemid(isrc),elemid(ifld),eh,ng_out,ng_out,xi(:,ng_out),xi(:,ng_out),&
            				wi(:,ng_out),wi(:,ng_out),xy_src,xy_fld)
            eat = eh - eg
            !debug
            !print*,'egts_sed_reg_reg called'
            !print*,'eg='
            !do i = 1,9
            !  	print'(9(f15.8,1x))',(eg(i,j),j=1,9)
            !enddo
            !print*,'ehts_reg_reg called'
            !print*,'eh='
            !do i = 1,9
            !  	print'(9(f15.8,1x))',(eh(i,j),j=1,9)
            !enddo
        endif
    endif
    CONTAINS
	!--------------------------------------------------
    Subroutine egts_ii_reg_reg(region,idelem,eg,nto,ntoln,nti,xio,xioln,xii,wio,wioln,wii,xlelem)
    IMPLICIT NONE

    INTEGER, INTENT(IN)			:: idelem,nto,ntoln,nti,region
    REAL(KIND=DBL), INTENT(OUT)	:: eg(:,:)
    REAL(KIND=DBL), INTENT(IN)	:: xio(:),xioln(:),xii(:),wio(:),wioln(:),wii(:)
    REAL(KIND=DBL), INTENT(IN)	:: xlelem(:,:)

    !...local variables
    INTEGER			:: neo,nei,lo,li
    REAL(KIND=DBL)	:: egbar1(3*NDOF,3*NDOF),egbar2(3*NDOF,3*NDOF),egtwobar(3*NDOF,3*NDOF)
    REAL(KIND=DBL)	:: alpha,alphapr,beta,eita,zeta
    REAL(KIND=DBL)	:: pso(3),dpso(3),psi(3),dpsi(3)
    INTEGER			:: i,j,k,l,ierr
    REAL(KIND=DBL)	:: val_inner_integral(NDOF,NDOF,3,3)
    REAL(KIND=DBL)	:: val_inner_integral1(NDOF,NDOF,3,3),val_inner_integral2(NDOF,NDOF,3,3)
    REAL(KIND=DBL)	:: f1(NDOF,NDOF,3,3),f2(NDOF,NDOF,3,3)
    REAL(KIND=DBL)	:: C2(NDOF,NDOF),G2(NDOF,NDOF),U2(NDOF,NDOF),GI1(NDOF,NDOF)
    REAL(KIND=DBL)	:: source(2),field(2)
    REAL(KIND=DBL)	:: dxdso(2),jacobo
    REAL(KIND=DBL)	:: temp(2),r

    !...get number of nodes of outer/inner elements
    neo = NODE(idelem)
    nei = NODE(idelem)

	!...compute (weakly) singular integral egbar2
	egbar2 = 0.0d0
	!...begins outer integral loop (logarith integration)
	do lo = 1,ntoln
    	!...coordinate of integration point
		beta = xioln(lo)
		!...initialize
        val_inner_integral = 0.0d0
		!...loop over inner integral (Gauss integration)
		do li =1,nti
			!...coordinate of Gauss point
        	alphapr = xii(li)
        	!...transformation back to original variables
            alpha = (1.0d0 - beta)*alphapr
            eita = alpha + beta
            zeta = alpha - beta
            !...obtain values and derivatives of REGULAR shape functions of
            !(original) outer integration for f(alphapr,betapr)
            call reg_shape(idelem,eita,pso,dpso)
            !...compute derivative of position vector of outer element
            dxdso = 0.0d0
            do i = 1,neo
              	do k = 1,2
                	dxdso(k) = dxdso(k) + dpso(i)*xlelem(k,i)
                enddo
            enddo
            !...compute jacobian of outer element
            jacobo = dsqrt(dxdso(1)*dxdso(1) + dxdso(2)*dxdso(2))
            !...obtain derivatives of shape functions of
            !(original) inner integration for f(alphapr,betapr)
            call reg_shape(idelem,zeta,psi,dpsi)
            !...compute value of f(alphapr,betapr)
            do i=1,neo
              	do j=1,nei
                	do k=1,NDOF
                    	do l=1,NDOF
                       	    f1(l,k,j,i) = dpso(i)*dpsi(j)*GI1(l,k)
                        enddo
                    enddo
                enddo
            enddo
            !...repeat with beta = -beta
            eita = alpha - beta
            zeta = alpha + beta
            
            !...obtain values and derivatives of REGULAR shape functions of
            !(original) outer integration for f_klij(alpha,-beta)
            call reg_shape(idelem,eita,pso,dpso)
            !...compute derivatives of position vector of outer element
            dxdso = 0.0d0
            do i = 1,neo
              	do k = 1,2
                	dxdso(k) = dxdso(k) + dpso(i)*xlelem(k,i)
                enddo
            enddo
            !...compute jacobian of outer element
            jacobo = dsqrt(dxdso(1)*dxdso(1) + dxdso(2)*dxdso(2))
            !...obtain values of shape functions and its derivatives of
            !(original) inner integration for f_klij(alpha,-beta)
            call reg_shape(idelem,zeta,psi,dpsi)
            !...compute value of f(alphapr,-beta)
            do i=1,neo
              	do j=1,nei
                	do k=1,NDOF
                    	do l=1,NDOF
                       	    f2(l,k,j,i) = dpso(i)*dpsi(j)*GI1(l,k)
                        enddo
                    enddo
                enddo
            enddo
            !...calculate val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                	do k=1,NDOF
                    	do l=1,NDOF
                        	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + &
                            (f1(l,k,j,i)+f2(l,k,j,i))*wii(li)
                        enddo
                    enddo
                enddo
            enddo
        enddo !of loop over inner integral
        do i=1,neo
          	do j=1,nei
            	do k=1,NDOF
                	do l=1,NDOF
                    	egbar2(NDOF*(i-1)+k,NDOF*(j-1)+l) = egbar2(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                            (1.0d0-beta)*val_inner_integral(l,k,j,i)*wioln(lo)
                    enddo
                enddo
            enddo
        enddo
    enddo !of loop over outer integral
    !...finally, multiply with 2.0d0 (Note: int(-ln(x)*f(x),x) = sum(w(i)*f(x(i)))
    egbar2 = -2.0d0*egbar2
    !...debugging
    !print*,'egbar2 = '
    !do k=1,3*neo
    !	print'(9(e12.6,1x))',(egbar2(k,l),l=1,9)
    !enddo
    !--------------------------------------------------
    !...evaluate egbar1 and egtwobar (regular integrals)
    egbar1 = 0.0d0
    egtwobar = 0.0d0
    do lo=1,nto
      	!...coordinate of Gauss point
      	eita = xio(lo)
        !...shape function & its derivatives at Gauss point
        call reg_shape(idelem,eita,pso,dpso)
        !...position of source point (y)
        source = 0.0d0
        do i=1,neo
          	do k=1,2
            	source(k) = source(k) + xlelem(k,i)*pso(i)
            enddo
        enddo
        !...compute derivative of position vector of outer element
        dxdso = 0.0d0
        do i = 1,neo
           	do k = 1,2
               	dxdso(k) = dxdso(k) + dpso(i)*xlelem(k,i)
            enddo
        enddo
        !...compute jacobian of outer element
        jacobo = dsqrt(dxdso(1)*dxdso(1) + dxdso(2)*dxdso(2))
        !...loop over inner integral
        val_inner_integral1 = 0.0d0   !...for egbar1
        val_inner_integral2 = 0.0d0   !...for egtwobar
        do li=1,nti
          	!...coordinate of Gauss point
            zeta = xii(li)

            !...shape function & its derivatives at Gauss point
        	call reg_shape(idelem,zeta,psi,dpsi)

        	!...position of field point (xi)
        	field = 0.0d0
        	do i=1,nei
          		do k=1,2
            		field(k) = field(k) + xlelem(k,i)*psi(i)
            	enddo
        	enddo
            !...compute r
            do k=1,2
              	temp(k) = field(k) - source(k)
            enddo
            r = dsqrt(temp(1)*temp(1) + temp(2)*temp(2))
			!...defensive programming
            if (r.le.SMALL_NUM) then
              	print*,'eatts.f90: r is too small'
                stop
            endif
            !...compute kernels
            call Kernel2(NDOF,region,source,field,C2,G2,U2)
			!...compute val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                	do k=1,NDOF
                    	do l=1,NDOF
                           	!...not include jacobian of outer element yet
                       		val_inner_integral1(l,k,j,i) = val_inner_integral1(l,k,j,i) + &
                            		GI1(l,k)*dpso(i)*dpsi(j)*dlog(2.0d0*r/dabs(zeta-eita))*wii(li)
                           	val_inner_integral2(l,k,j,i) = val_inner_integral2(l,k,j,i) + &
                            		G2(l,k)*dpso(i)*dpsi(j)*wii(li)
                        enddo
                    enddo
                enddo
            enddo
        enddo !over loop of inner integral
        do i=1,neo
          	do j=1,nei
            	do k=1,NDOF
                	do l=1,NDOF
                    	egbar1(NDOF*(i-1)+k,NDOF*(j-1)+l) = egbar1(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        	val_inner_integral1(l,k,j,i)*wio(lo)
                        egtwobar(NDOF*(i-1)+k,NDOF*(j-1)+l) = egtwobar(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        	val_inner_integral2(l,k,j,i)*wio(lo)
                    enddo
                enddo
            enddo
        enddo
    enddo !over loop of outer integral
    !...debugging
    !print*,'egbar1 = '
    !do k=1,3*neo
    !	print'(9(e12.6,1x))',(egbar1(k,l),l=1,9)
    !enddo
    !print*,'egtwobar = '
    !do k=1,3*neo
    !	print'(9(e12.6,1x))',(egtwobar(k,l),l=1,9)
    !enddo
    !...element stiffness matrix
    eg = egbar1 + egbar2 + egtwobar
	End Subroutine egts_ii_reg_reg
    !--------------------------------------------------
	Subroutine egts_ii_tip_tip(region,idelem,eg,nto,ntoln,nti,xio,xioln,xii,wio,wioln,wii,xlelem)
    !...Note: even test function uses regular shape function on tip element, but in order to
    !make ln(dabs(2r/zetapr-eitapr)) to be regular, we still need to stretch out eita to eitapr
    !Note: just stretch out eita, NOT use tip shape function (Not call tip_shape(...)...)
    !...10/5/09: fix the bug in ln(2r/(zetapr-eitapr))
    IMPLICIT NONE

    INTEGER, INTENT(IN)			:: idelem,nto,ntoln,nti,region
    REAL(KIND=DBL), INTENT(OUT)	:: eg(:,:)
    REAL(KIND=DBL), INTENT(IN)	:: xio(:),xioln(:),xii(:),wio(:),wioln(:),wii(:)
    REAL(KIND=DBL), INTENT(IN)	:: xlelem(:,:)

    !...local variables
    INTEGER			:: neo,nei,lo,li
    REAL(KIND=DBL)	:: egbar1(3*NDOF,3*NDOF),egbar2(3*NDOF,3*NDOF),egtwobar(3*NDOF,3*NDOF)
    REAL(KIND=DBL)	:: egbar31(3*NDOF,3*NDOF),egbar32(3*NDOF,3*NDOF)
    REAL(KIND=DBL)	:: alpha,alphapr,beta,eita,zeta,eitapr,zetapr
    REAL(KIND=DBL)	:: rho,gamma,deltaa
    REAL(KIND=DBL)	:: pso(3),dpso(3),psi(3),dpsi(3)
    REAL(KIND=DBL)	:: dtipi(3)
    INTEGER			:: i,j,k,l,ierr
    REAL(KIND=DBL)	:: val_inner_integral(NDOF,NDOF,3,3)
    REAL(KIND=DBL)	:: val_inner_integral1(NDOF,NDOF,3,3),val_inner_integral2(NDOF,NDOF,3,3)
    REAL(KIND=DBL)	:: f1(NDOF,NDOF,3,3),f2(NDOF,NDOF,3,3)
    REAL(KIND=DBL)	:: C2(NDOF,NDOF),G2(NDOF,NDOF),U2(NDOF,NDOF),GI1(NDOF,NDOF)
    REAL(KIND=DBL)	:: source(2),field(2)
    REAL(KIND=DBL)	:: dxdso(2),jacobo
    REAL(KIND=DBL)	:: temp(2),r

    !...get number of nodes of outer/inner elements
    neo = NODE(idelem)
    nei = NODE(idelem)

	!...compute (weakly) singular integral egbar2
	egbar2 = 0.0d0

	!...begins outer integral loop (logarith integration)
	do lo = 1,ntoln
    	!...coordinate of integration point
		beta = xioln(lo)
		!...initialize
        val_inner_integral = 0.0d0
		!...loop over inner integral (Gauss integration)
		do li =1,nti
			!...coordinate of Gauss point
        	alphapr = xii(li)
        	!...transformation back to original variables
            alpha = (1.0d0 - beta)*alphapr
            eitapr = alpha + beta
            zetapr = alpha - beta
            if (idelem.eq.CTIP1) then
              	eita = 0.5d0*(eitapr+1.0d0)*(eitapr+1.0d0)-1.0d0
                zeta = 0.5d0*(zetapr+1.0d0)*(zetapr+1.0d0)-1.0d0
            elseif (idelem.eq.CTIP2) then
            	eita = 1.0d0-0.5d0*(1.0d0-eitapr)*(1.0d0-eitapr)
                zeta = 1.0d0-0.5d0*(1.0d0-zetapr)*(1.0d0-zetapr)
            else
                print*,'egts_ii_tip_tip: id of element is out of range'
                stop
            endif
            !...obtain values and derivatives of REGULAR shape functions of
            !(original) outer integration for f(alphapr,betapr)
            call reg_shape(idelem,eita,pso,dpso)
            !...compute derivative of position vector of outer element
            !dxdso = 0.0d0
            !do l = 1,neo
            !  	do k = 1,2
            !    	dxdso(k) = dxdso(k) + dpso(l)*xlelem(k,l)
            !    enddo
            !enddo
            !...compute jacobian of outer element
            !jacobo = dsqrt(dxdso(1)*dxdso(1) + dxdso(2)*dxdso(2))
            !...obtain derivatives of SPECIAL shape functions of
            !(original) inner integration for f(alphapr,betapr)
            call tip_dshape(idelem,zeta,dpsi)
            !...compute value of f(alphapr,betapr)
            do i=1,neo
              	do j=1,nei
                	do k=1,NDOF
                    	do l=1,NDOF
                    	    if (idelem.eq.CTIP1) then
                        	    f1(l,k,j,i) = (eitapr+1.0d0)*(zetapr+1.0d0)*dpso(i)*dpsi(j)*GI1(l,k)
                        	else
                        	    f1(l,k,j,i) = (1.0d0-eitapr)*(1.0d0-zetapr)*dpso(i)*dpsi(j)*GI1(l,k)
                        	endif
                        enddo
                    enddo
                enddo
            enddo
            !...repeat with beta = -beta
            eitapr = alpha - beta
            zetapr = alpha + beta
            if (idelem.eq.CTIP1) then
              	eita = 0.5d0*(eitapr+1.0d0)*(eitapr+1.0d0)-1.0d0
                zeta = 0.5d0*(zetapr+1.0d0)*(zetapr+1.0d0)-1.0d0
            elseif (idelem.eq.CTIP2) then
            	eita = 1.0d0-0.5d0*(1.0d0-eitapr)*(1.0d0-eitapr)
                zeta = 1.0d0-0.5d0*(1.0d0-zetapr)*(1.0d0-zetapr)
            else
                print*,'egts_ii_tip_tip: id of element is out of range'
                stop
            endif
            !...obtain values and derivatives of REGULAR shape functions of
            !(original) outer integration for f(alpha,-beta)
            call reg_shape(idelem,eita,pso,dpso)
            !...compute derivatives of position vector of outer element
            !dxdso = 0.0d0
            !do l = 1,neo
            !  	do k = 1,2
            !    	dxdso(k) = dxdso(k) + dpso(l)*xlelem(k,l)
            !    enddo
            !enddo
            !...compute jacobian of outer element
            !jacobo = dsqrt(dxdso(1)*dxdso(1) + dxdso(2)*dxdso(2))
            !...obtain values of SPECIAL shape functions and its derivatives of
            !(original) inner integration for f_klij(alpha,-beta)
            call tip_dshape(idelem,zeta,dpsi)
            !...compute value of f(alphapr,-beta)
            do i=1,neo
              	do j=1,nei
                	do k=1,NDOF
                    	do l=1,NDOF
                    	    if (idelem.eq.CTIP1) then
                        	    f2(l,k,j,i) = (eitapr+1.0d0)*(zetapr+1.0d0)*dpso(i)*dpsi(j)*GI1(l,k)
                        	else
                        	    f2(l,k,j,i) = (1.0d0-eitapr)*(1.0d0-zetapr)*dpso(i)*dpsi(j)*GI1(l,k)
                        	endif
                        enddo
                    enddo
                enddo
            enddo
            !...calculate val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                	do k=1,NDOF
                    	do l=1,NDOF
                        	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + &
                            (f1(l,k,j,i)+f2(l,k,j,i))*wii(li)
                        enddo
                    enddo
                enddo
            enddo
        enddo !of loop over inner integral
        do i=1,neo
          	do j=1,nei
            	do k=1,NDOF
                	do l=1,NDOF
                    	egbar2(NDOF*(i-1)+k,NDOF*(j-1)+l) = egbar2(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                            (1.0d0-beta)*val_inner_integral(l,k,j,i)*wioln(lo)
                    enddo
                enddo
            enddo
        enddo
    enddo !of loop over outer integral
    !...finally, multiply with 2.0d0 (Note: int(-ln(x)*f(x),x) = sum(w(i)*f(x(i)))
    egbar2 = -2.0d0*egbar2
    !...debugging
    !print*,'egbar2 = '
    !do k=1,3*neo
    !	print'(9(e12.6,1x))',(egbar2(k,l),l=1,9)
    !enddo
    !--------------------------------------------------
    !...evaluate egbar1 and egtwobar (regular integrals)
    egbar1 = 0.0d0
    egtwobar = 0.0d0
    do lo=1,nto
      	!...coordinate of Gauss point
      	eitapr = xio(lo)
        !...transform back to original eita
        if (idelem.eq.CTIP1) then
        	eita = 0.5d0*(eitapr+1.0d0)*(eitapr+1.0d0)-1.0d0
        else
            eita = 1.0d0-0.5d0*(1.0d0-eitapr)*(1.0d0-eitapr)
        endif
        !...shape function & its derivatives at Gauss point
        call reg_shape(idelem,eita,pso,dpso)
        !...position of source point (y)
        source = 0.0d0
        do i=1,neo
          	do k=1,2
            	source(k) = source(k) + xlelem(k,i)*pso(i)
            enddo
        enddo
        !...compute derivative of position vector of outer element
        dxdso = 0.0d0
        do l = 1,neo
           	do k = 1,2
               	dxdso(k) = dxdso(k) + dpso(l)*xlelem(k,l)
            enddo
        enddo
        !...compute jacobian of outer element
        jacobo = dsqrt(dxdso(1)*dxdso(1) + dxdso(2)*dxdso(2))
        !...loop over inner integral
        val_inner_integral1 = 0.0d0   !...for egbar1
        val_inner_integral2 = 0.0d0   !...for egtwobar
        do li=1,nti
          	!...coordinate of Gauss point
            zetapr = xii(li)
            !...transform back to original zeta
            if (idelem.eq.CTIP1) then
                zeta = 0.5d0*(zetapr+1.0d0)*(zetapr+1.0d0)-1.0d0
            else
                zeta = 1.0d0-0.5d0*(1.0d0-zetapr)*(1.0d0-zetapr)
            endif
            !...shape function & its derivatives at Gauss point
        	call reg_shape(idelem,zeta,psi,dpsi)
            call tip_dshape(idelem,zeta,dtipi)
        	!...position of field point (xi)
        	field = 0.0d0
        	do i=1,nei
          		do k=1,2
            		field(k) = field(k) + xlelem(k,i)*psi(i)
            	enddo
        	enddo
            !...compute r
            do k=1,2
              	temp(k) = field(k) - source(k)
            enddo
            r = dsqrt(temp(1)*temp(1) + temp(2)*temp(2))
			!...defensive programming
            if (r.le.SMALL_NUM) then
              	print*,'egts_ii_tip_tip: r is too small'
                stop
            endif
            !...compute kernels
            call Kernel2(NDOF,region,source,field,C2,G2,U2)
			!...compute val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                	do k=1,NDOF
                    	do l=1,NDOF
                        	if (idelem.eq.CTIP1) then
                            	!...jacobian of outer element yet: cancel out with D-operator
                        		val_inner_integral1(l,k,j,i) = val_inner_integral1(l,k,j,i) + &
                            		GI1(l,k)*dpso(i)*dtipi(j)*(eitapr+1.0d0)*(zetapr+1.0d0)*dlog(4.0d0*r/dabs(zetapr-eitapr)/dabs(2.0d0+zetapr+eitapr))*wii(li)
                            	val_inner_integral2(l,k,j,i) = val_inner_integral2(l,k,j,i) + &
                            		G2(l,k)*dpso(i)*dtipi(j)*(eitapr+1.0d0)*(zetapr+1.0d0)*wii(li)
                            else
                              	!...jacobian of outer element yet: cancle out with D-operator
                              	val_inner_integral1(l,k,j,i) = val_inner_integral1(l,k,j,i) + &
                            		GI1(l,k)*dpso(i)*dtipi(j)*(1.0d0-eitapr)*(1.0d0-zetapr)*dlog(4.0d0*r/dabs(zetapr-eitapr)/dabs(2.0d0-zetapr-eitapr))*wii(li)
                            	val_inner_integral2(l,k,j,i) = val_inner_integral2(l,k,j,i) + &
                            		G2(l,k)*dpso(i)*dtipi(j)*(1.0d0-eitapr)*(1.0d0-zetapr)*wii(li)
                            endif
                        enddo
                    enddo
                enddo
            enddo
        enddo !over loop of inner integral
        do i=1,neo
          	do j=1,nei
            	do k=1,NDOF
                	do l=1,NDOF
                    	egbar1(NDOF*(i-1)+k,NDOF*(j-1)+l) = egbar1(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        	val_inner_integral1(l,k,j,i)*wio(lo)
                        egtwobar(NDOF*(i-1)+k,NDOF*(j-1)+l) = egtwobar(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        	val_inner_integral2(l,k,j,i)*wio(lo)
                    enddo
                enddo
            enddo
        enddo
    enddo !over loop of outer integral
    !...debugging
    !print*,'egbar1 = '
    !do k=1,3*neo
    !	print'(9(e12.6,1x))',(egbar1(k,l),l=1,9)
    !enddo
    !print*,'egtwobar = '
    !do k=1,3*neo
    !	print'(9(e12.6,1x))',(egtwobar(k,l),l=1,9)
    !enddo

    !...calculate egbar31
    egbar31 = 0.0d0
	!...begins outer integral loop
	do lo = 1,nto
    	!...coordinate of integration point
        if (idelem.eq.CTIP1) then
			rho = xio(lo)
        	gamma = 1.0d0 - 0.5d0*(1.0d0-rho)*(1.0d0-rho)
        elseif (idelem.eq.CTIP2) then
        	gamma = xio(lo)
        else
          	print*,'egts_ii_tip_tip: wrong id of tip element'
            stop
        endif
        val_inner_integral = 0.0d0
		!...loop over inner integral (Gauss integration)
		do li =1,nti
			!...coordinate of Gauss point
        	deltaa = xii(li)
        	!...transformation back to original variables
            alpha = 0.5d0*(1.0d0-gamma)*deltaa
            beta = 0.5d0*(1.0d0+gamma)
            eitapr = -1.0d0*(alpha + beta)
            zetapr = alpha - beta
            !...transform to get rid of singularity of tip shape function
            if (idelem.eq.CTIP1) then
            	eita = 0.5d0*(eitapr+1.0d0)*(eitapr+1.0d0)-1.0d0
                zeta = 0.5d0*(zetapr+1.0d0)*(zetapr+1.0d0)-1.0d0
            elseif (idelem.eq.CTIP2) then
            	eita = 1.0d0-0.5d0*(1.0d0-eitapr)*(1.0d0-eitapr)
                zeta = 1.0d0-0.5d0*(1.0d0-zetapr)*(1.0d0-zetapr)
            else
              	print*,'eg_ii_tip_tip: id of tip element is out of range'
                stop
            endif
            !...obtain value of shape function and its derivatives of
            call reg_shape(idelem,eita,pso,dpso)            
            call tip_dshape(idelem,zeta,dtipi)
        	
            !...calculate val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                	do k=1,NDOF
                    	do l=1,NDOF
                        	if (idelem.eq.CTIP1) then
                        		val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + &
                            		GI1(l,k)*dpso(i)*dtipi(j)*(1.0d0+eitapr)*(1.0d0+zetapr)*wii(li)
                            elseif (idelem.eq.CTIP2) then
                              	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + &
                            		GI1(l,k)*dpso(i)*dtipi(j)*(1.0d0-eitapr)*(1.0d0-zetapr)*wii(li)
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
                    	if (idelem.eq.CTIP1) then
                       		egbar31(NDOF*(i-1)+k,NDOF*(j-1)+l) = egbar31(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		((1.0d0-rho)**3.0d0)*dlog(0.5d0*(1.0d0-rho))*val_inner_integral(l,k,j,i)*wio(lo)
                        elseif (idelem.eq.CTIP2) then
                        	egbar31(NDOF*(i-1)+k,NDOF*(j-1)+l) = egbar31(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		(1.0d0-gamma)*dlog(0.5d0*(3.0d0+gamma))*val_inner_integral(l,k,j,i)*wio(lo)
                        endif
                    enddo
                enddo
            enddo
        enddo
    enddo !of loop over outer integral
    egbar31 = 0.5d0*egbar31
    !...for ekbar32
    egbar32 = 0.0d0
    do lo = 1,nto
    	!...coordinate of integration point
        if (idelem.eq.CTIP1) then
			gamma = xio(lo)
        elseif (idelem.eq.CTIP2) then
        	rho = xio(lo)
            gamma = 0.5d0*(1.0d0+rho)*(1.0d0+rho)-1.0d0
        endif
        val_inner_integral = 0.0d0
		!...loop over inner integral (Gauss integration)
		do li =1,nti
			!...coordinate of Gauss point
        	deltaa = xii(li)
        	!...transformation back to original variables
            alpha = 0.5d0*(1.0d0+gamma)*deltaa
            beta = 0.5d0*(gamma-1.0d0)
            eitapr = -1.0d0*(alpha + beta)
            zetapr = alpha - beta
            !...transform to get rid of singularity of tip shape function
            if (idelem.eq.CTIP1) then
            	eita = 0.5d0*(eitapr+1.0d0)*(eitapr+1.0d0)-1.0d0
                zeta = 0.5d0*(zetapr+1.0d0)*(zetapr+1.0d0)-1.0d0
            elseif (idelem.eq.CTIP2) then
            	eita = 1.0d0-0.5d0*(1.0d0-eitapr)*(1.0d0-eitapr)
                zeta = 1.0d0-0.5d0*(1.0d0-zetapr)*(1.0d0-zetapr)
            else
              	print*,'egts_ii_tip_tip: id of tip element is out of range'
                stop
            endif
            !...obtain value of shape function and its derivatives of 
            call reg_shape(idelem,eita,pso,dpso)
            call tip_dshape(idelem,zeta,dtipi)
        	
            !...calculate val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                	do k=1,NDOF
                    	do l=1,NDOF
                        	if (idelem.eq.CTIP1) then
                        		val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + &
                            		GI1(l,k)*dpso(i)*dtipi(j)*(1.0d0+eitapr)*(1.0d0+zetapr)*wii(li)
                            elseif (idelem.eq.CTIP2) then
                              	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + &
                            		GI1(l,k)*dpso(i)*dtipi(j)*(1.0d0-eitapr)*(1.0d0-zetapr)*wii(li)
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
                    	if (idelem.eq.CTIP1) then
                       		egbar32(NDOF*(i-1)+k,NDOF*(j-1)+l) = egbar32(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		(1.0d0+gamma)*dlog(0.5d0*(3.0d0-gamma))*val_inner_integral(l,k,j,i)*wio(lo)
                        elseif (idelem.eq.CTIP2) then
                        	egbar32(NDOF*(i-1)+k,NDOF*(j-1)+l) = egbar32(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		((1.0d0+rho)**3.0d0)*dlog(0.5d0*(1.0d0+rho))*val_inner_integral(l,k,j,i)*wio(lo)
                        endif
                    enddo
                enddo
            enddo
        enddo
    enddo !of loop over outer integral
    egbar32 = 0.5d0*egbar32
    
    !...element stiffness matrix
    eg = egbar1 + egbar2 + egtwobar + egbar31 + egbar32
    
	End Subroutine egts_ii_tip_tip
    !-----------------------------------------------
	Subroutine egts_adt_reg_tip(region,ido,idi,eg,nto,nti,xio,xii,wio,wii,xlo,xli,adt_type)
    IMPLICIT NONE

    INTEGER, INTENT(IN)		:: ido,idi,nto,nti,region
    INTEGER, INTENT(IN)		:: adt_type
    REAL(KIND=DBL), INTENT(OUT)	:: eg(:,:)
    REAL(KIND=DBL), INTENT(IN)	:: xio(:),xii(:),wio(:),wii(:)
    REAL(KIND=DBL), INTENT(IN)	:: xlo(:,:),xli(:,:)

    !...local variables
    INTEGER			:: neo,nei,lo,li
    REAL(KIND=DBL)	:: egbar1(3*NDOF,3*NDOF),egbar2(3*NDOF,3*NDOF),egtwobar(3*NDOF,3*NDOF)
    REAL(KIND=DBL)	:: alpha,beta,eita,zeta,gamma,deltaa,zetapr,rho,eitabar
    REAL(KIND=DBL)	:: pso(3),dpso(3),psi(3),dpsi(3),dtipi(3)
    INTEGER			:: i,j,k,l,ierr
    REAL(KIND=DBL)	:: val_inner_integral(NDOF,NDOF,3,3)
    REAL(KIND=DBL)	:: C2(NDOF,NDOF),G2(NDOF,NDOF),U2(NDOF,NDOF),GI1(NDOF,NDOF)
    REAL(KIND=DBL)	:: source(2),field(2),dxdso(2),jacobo
    REAL(KIND=DBL)	:: temp(2),r

    !...get number of nodes of outer/inner elements
    neo = NODE(ido)
    nei = NODE(idi)

	!...compute matrix egbar2: singular if adt_type=2 or 4
	egbar2=0.0d0

	!...begins outer integral loop
	do lo = 1,nto
    	!...coordinate of integration point, add the case of adt_type = 4 (3/23/09)
		if ((adt_type.eq.2).or.(adt_type.eq.4)) then
          	!...need stretching of gamma to rho to improve ln integration
            rho = xio(lo)
            gamma = 0.5d0*(rho+1.0d0)*(rho+1.0d0)-1.0d0
        elseif ((adt_type.eq.1).or.(adt_type.eq.3)) then
        	!...no need to stretch gamma
			gamma = xio(lo)
        else
          	print*,'error in egts_adt_reg_tip: adt_type out of range'
            stop
        endif
		!...initialize
        val_inner_integral = 0.0d0
		!...loop over inner integral (Gauss integration)
		do li =1,nti
			!...coordinate of Gauss point
        	deltaa = xii(li)
        	!...transformation back to original variables
            alpha = 0.5d0*(gamma+1.0d0)*deltaa
            beta = 0.5d0*(gamma-1.0d0)
            !...now, distinguish adjacent types between 2/4 and 1/3
            if ((adt_type.eq.4).or.(adt_type.eq.3)) then
              	eitabar = alpha + beta
                eita = -1.0d0*eitabar
            else
            	eita = alpha + beta
            endif
            zetapr = alpha - beta
            !...transform to get rid of singularity of tip shape function
            if (idi.eq.CTIP1) then
                zeta = 0.5d0*(zetapr+1.0d0)*(zetapr+1.0d0)-1.0d0
            elseif (idi.eq.CTIP2) then
                zeta = 1.0d0-0.5d0*(1.0d0-zetapr)*(1.0d0-zetapr)
            else
              	print*,'eg_adt_reg_tip: id of inner element is out of range'
                stop
            endif
            !...obtain value of shape function and its derivatives
            call reg_shape(ido,eita,pso,dpso)
            call reg_shape(idi,zeta,psi,dpsi)
            call tip_dshape(idi,zeta,dtipi)
            !...compute derivative of position vector of outer element
        	dxdso = 0.0d0
        	do i = 1,neo
           		do k = 1,2
               		dxdso(k) = dxdso(k) + dpso(i)*xlo(k,i)
            	enddo
        	enddo
        	!...compute jacobian of outer element
        	jacobo = dsqrt(dxdso(1)*dxdso(1) + dxdso(2)*dxdso(2))
        	!...position of source point (y)
        	source = 0.0d0
        	do i=1,neo
          		do k=1,2
            		source(k) = source(k) + xlo(k,i)*pso(i)
            	enddo
        	enddo
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
            if (r.le.SMALL_NUM) then
              	print*,'egts_adt_reg_tip: r is too small'
                stop
            endif
            !...calculate val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                	do k=1,NDOF
                    	do l=1,NDOF
                        	if (idi.eq.CTIP2) then
                        		val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + &
                            		GI1(l,k)*dpso(i)*dtipi(j)*(1.0d0-zetapr)*dlog(r)*wii(li)
                            elseif (idi.eq.CTIP1) then
                              	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + &
                            		GI1(l,k)*dpso(i)*dtipi(j)*(1.0d0+zetapr)*dlog(r)*wii(li)
                            else
                              	print*,'eg_adt_reg_tip: inner element must be a tip element'
                                stop
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
                    	!...add one more case of adt_type = 4
                    	if ((adt_type.eq.2).or.(adt_type.eq.4)) then
                        	egbar2(NDOF*(i-1)+k,NDOF*(j-1)+l) = egbar2(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		((1.0d0+rho)**3.0d0)*val_inner_integral(l,k,j,i)*wio(lo)
                        elseif ((adt_type.eq.1).or.(adt_type.eq.3)) then
                    		egbar2(NDOF*(i-1)+k,NDOF*(j-1)+l) = egbar2(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		(1.0d0+gamma)*val_inner_integral(l,k,j,i)*wio(lo)
                        else
                          	print*,'adt_type out of range'
                            stop
                        endif
                    enddo
                enddo
            enddo
        enddo
    enddo !of loop over outer integral
    !...finalize ekbar2
    if ((adt_type.eq.2).or.(adt_type.eq.4)) then
      	egbar2 = 0.25d0*egbar2
    elseif ((adt_type.eq.1).or.(adt_type.eq.3)) then
    	egbar2 = 0.5d0*egbar2
    else
      	print*,'adt_type out of range'
        stop
    endif
	!--------------------------------------------------
    !...compute matrix ekbar1: singular if adt_type=1 or 3
	egbar1=0.0d0
	!...begins outer integral loop
	do lo = 1,nto
    	!...coordinate of integration point
		if ((adt_type.eq.1).or.(adt_type.eq.3)) then
          	!...need to stretch gamma to rho to improve ln integration
            rho = xio(lo)
            gamma = 1.0d0 - 0.5d0*(1.0d0-rho)*(1.0d0-rho)
        elseif ((adt_type.eq.2).or.(adt_type.eq.4)) then
        	!...no need to stretch gamma
			gamma = xio(lo)
        else
          	print*,'error in egts_adt_reg_tip: adt_type out of range'
            stop
        endif
		!...initialize
        val_inner_integral = 0.0d0
		!...loop over inner integral (Gauss integration)
		do li =1,nti
			!...coordinate of Gauss point
        	deltaa = xii(li)
        	!...transformation back to original variables
            alpha = 0.5d0*(1.0d0-gamma)*deltaa
            beta = 0.5d0*(gamma+1.0d0)
            !...distinguish adjacent types between 2/4 and 1/3
            if ((adt_type.eq.3).or.(adt_type.eq.4)) then
              	eitabar = alpha + beta
                eita = -1.0d0*eitabar
            else
            	eita = alpha + beta
            endif
            zetapr = alpha - beta
            !...transform to get rid of singularity of tip shape function
            if (idi.eq.CTIP1) then
                zeta = 0.5d0*(zetapr+1.0d0)*(zetapr+1.0d0)-1.0d0
            elseif (idi.eq.CTIP2) then
                zeta = 1.0d0-0.5d0*(1.0d0-zetapr)*(1.0d0-zetapr)
            else
              	print*,'eg_adt_reg_tip: id of inner element must be a tip element'
                stop
            endif
            !...obtain value of shape function and its derivatives of 
            call reg_shape(ido,eita,pso,dpso)
            call reg_shape(idi,zeta,psi,dpsi)
            call tip_dshape(idi,zeta,dtipi)
            !...compute derivative of position vector of outer element
        	dxdso = 0.0d0
        	do i = 1,neo
           		do k = 1,2
               		dxdso(k) = dxdso(k) + dpso(i)*xlo(k,i)
            	enddo
        	enddo
        	!...compute jacobian of outer element
        	jacobo = dsqrt(dxdso(1)*dxdso(1) + dxdso(2)*dxdso(2))
        	!...position of source point (y)
        	source = 0.0d0
        	do i=1,neo
          		do k=1,2
            		source(k) = source(k) + xlo(k,i)*pso(i)
            	enddo
        	enddo
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
            if (r.le.SMALL_NUM) then
              	print*,'r is too small'
                stop
            endif
            !...calculate val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                	do k=1,NDOF
                    	do l=1,NDOF
                        	if (idi.eq.CTIP2) then
                        		val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + &
                            		GI1(l,k)*dpso(i)*dtipi(j)*(1.0d0-zetapr)*dlog(r)*wii(li)
                            elseif (idi.eq.CTIP1) then
                              	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + &
                            		GI1(l,k)*dpso(i)*dtipi(j)*(1.0d0+zetapr)*dlog(r)*wii(li)
                            else
                              	print*,'egts_adt_reg_tip: inner element must be a tip element'
                                stop
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
                    	if ((adt_type.eq.1).or.(adt_type.eq.3)) then
                        	egbar1(NDOF*(i-1)+k,NDOF*(j-1)+l) = egbar1(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		((1.0d0-rho)**3.0d0)*val_inner_integral(l,k,j,i)*wio(lo)
                        elseif ((adt_type.eq.2).or.(adt_type.eq.4)) then
                    		egbar1(NDOF*(i-1)+k,NDOF*(j-1)+l) = egbar1(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		(1.0d0-gamma)*val_inner_integral(l,k,j,i)*wio(lo)
                        else
                          	print*,'adt_type out of range'
                            stop
                        endif
                    enddo
                enddo
            enddo
        enddo
    enddo !of loop over outer integral
    !...finalize ekbar1
    if ((adt_type.eq.1).or.(adt_type.eq.3)) then
      	egbar1 = 0.25d0*egbar1
    elseif ((adt_type.eq.2).or.(adt_type.eq.4)) then
    	egbar1 = 0.5d0*egbar1
    else
      	print*,'adt_type out of range'
        stop
    endif
    !--------------------------------------------------
    !...evaluate egtwobar (regular integrals)
    egtwobar = 0.0d0
    do lo=1,nto
      	!...coordinate of Gauss point
        eita = xio(lo)
        !...shape function & its derivatives at Gauss point
        call reg_shape(ido,eita,pso,dpso)
        !...compute derivative of position vector of outer element
        dxdso = 0.0d0
        do l = 1,neo
        	do k = 1,2
           		dxdso(k) = dxdso(k) + dpso(l)*xlo(k,l)
           	enddo
        enddo
        !...compute jacobian of outer element
        jacobo = dsqrt(dxdso(1)*dxdso(1) + dxdso(2)*dxdso(2))
        !...position of source point (y)
        source = 0.0d0
        do i=1,neo
          	do k=1,2
            	source(k) = source(k) + xlo(k,i)*pso(i)
            enddo
        enddo
        
        !...loop over inner integral
        val_inner_integral = 0.0d0
        do li=1,nti
          	!...coordinate of Gauss point
            zetapr = xii(li)
            !...transform to get rid of singularity at tip shape functions
        	if (idi.eq.CTIP1) then
        		zeta = 0.5d0*(zetapr+1.0d0)*(zetapr+1.0d0)-1.0d0
        	elseif (idi.eq.CTIP2) then
           		zeta = 1.0d0-0.5d0*(1.0d0-zetapr)*(1.0d0-zetapr)
        	else
           		print*,'eg_adt_reg_tip: id of inner element is out of range'
            	stop
        	endif
            !...shape function & its derivatives at Gauss point
        	call reg_shape(idi,zeta,psi,dpsi)
            call tip_dshape(idi,zeta,dtipi)
        	!...position of field point (xi)
        	field = 0.0d0
        	do i=1,nei
          		do k=1,2
            		field(k) = field(k) + xli(k,i)*psi(i)
            	enddo
        	enddo
        	
            !...compute kernels
            call Kernel2(NDOF,region,source,field,C2,G2,U2)
			!...compute val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                	do k=1,NDOF
                    	do l=1,NDOF
                        	if (idi.eq.CTIP2) then
                            	!...jacobian of outer element: cancel out with D-operator
                            	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + &
                            		G2(l,k)*dpso(i)*dtipi(j)*(1.0d0-zetapr)*wii(li)
                            elseif (idi.eq.CTIP1) then
                              	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + &
                            		G2(l,k)*dpso(i)*dtipi(j)*(1.0d0+zetapr)*wii(li)
                            else
                              	print*,'eg_adt_tip_tip: inner element must be a tip element'
                                stop
                            endif
                        enddo
                    enddo
                enddo
            enddo
        enddo !over loop of inner integral
        do i=1,neo
          	do j=1,nei
            	do k=1,NDOF
                	do l=1,NDOF
                        egtwobar(NDOF*(i-1)+k,NDOF*(j-1)+l) = egtwobar(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        val_inner_integral(l,k,j,i)*wio(lo)
                    enddo
                enddo
            enddo
        enddo
    enddo !over loop of outer integral

    !...element stiffness matrix
    eg = egbar1 + egbar2 + egtwobar
	End Subroutine egts_adt_reg_tip

	!-----------------------------------------------

    Subroutine egts_adt_reg_reg(region,ido,idi,ek,nto,ntoln,nti,xio,xioln,xii,wio,wioln,wii,xlo,xli,adt_type)
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
    REAL(KIND=DBL)	:: val_inner_integral(NDOF,NDOF,3,3),jacobo
    REAL(KIND=DBL)	:: C2(NDOF,NDOF),G2(NDOF,NDOF),U2(NDOF,NDOF),GI1(NDOF,NDOF)
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
              	print*,'egts_adt_reg_reg: wrong adt_type:',adt_type
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
            
            !...calculate val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                	do k=1,NDOF
                    	do l=1,NDOF
                        	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + &
                            GI1(l,k)*dpso(i)*dpsi(j)*wii(li)
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
              	print*,'egts_adt_reg_reg: wrong adt_type:',adt_type
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
            
            !...calculate val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                	do k=1,NDOF
                    	do l=1,NDOF
                        	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + &
                            GI1(l,k)*dpso(i)*dpsi(j)*wii(li)
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
          	print*,'eg_adt_reg_reg: wrong adt_type'
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
            
            !...compute the vector r = xi - y
            do k=1,2
              	temp(k) = field(k) - source(k)
            enddo
            r = dsqrt(temp(1)*temp(1) + temp(2)*temp(2))
            !...defensive programming
            if (r.le.SMALL_NUM) then
              	print*,'egts_adt_reg_reg: r is too small'
                stop
            endif
            !...calculate val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                	do k=1,NDOF
                    	do l=1,NDOF
                        	if (adt_type.eq.1) then
                        		val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + &
                            		GI1(l,k)*dpso(i)*dpsi(j)*dlog(2.0d0*r/(2.0d0-eita+zeta))*wii(li)
                            elseif (adt_type.eq.3) then
                            	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + &
                            		GI1(l,k)*dpso(i)*dpsi(j)*dlog(2.0d0*r/(2.0d0-eitabar+zeta))*wii(li)
                            elseif (adt_type.eq.2) then
                        		val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + &
                            		GI1(l,k)*dpso(i)*dpsi(j)*dlog(2.0d0*r/(2.0d0-zeta+eita))*wii(li)
                            elseif (adt_type.eq.4) then
                            	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + &
                            	GI1(l,k)*dpso(i)*dpsi(j)*dlog(2.0d0*r/(2.0d0-zeta+eitabar))*wii(li)
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
            
            !...compute kernel
            call Kernel2(NDOF,region,source,field,C2,G2,U2)
			!...compute val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                	do k=1,NDOF
                    	do l=1,NDOF
                            val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + &
                            G2(l,k)*dpso(i)*dpsi(j)*wii(li)
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
	End Subroutine egts_adt_reg_reg
    !-----------------------------------------------
    
    Subroutine egts_sed_reg_tip(region,ido,idi,eg,nto,nti,xio,xii,wio,wii,xlo,xli)
    IMPLICIT NONE

    INTEGER, INTENT(IN)		:: ido,idi,nto,nti,region
    REAL(KIND=DBL), INTENT(OUT)	:: eg(:,:)
    REAL(KIND=DBL), INTENT(IN)	:: xio(:),xii(:),wio(:),wii(:)
    REAL(KIND=DBL), INTENT(IN)	:: xlo(:,:),xli(:,:)

    !...local variables
    INTEGER			:: neo,nei,lo,li
    REAL(KIND=DBL)	:: egbar(3*NDOF,3*NDOF),egtwobar(3*NDOF,3*NDOF)
    REAL(KIND=DBL)	:: eita,zetapr,zeta
    REAL(KIND=DBL)	:: pso(3),dpso(3),psi(3),dpsi(3),dtipi(3)
    INTEGER			:: i,j,k,l,ierr
    REAL(KIND=DBL)	:: val_inner_integral1(NDOF,NDOF,3,3),val_inner_integral2(NDOF,NDOF,3,3)
    REAL(KIND=DBL)	:: C2(NDOF,NDOF),G2(NDOF,NDOF),U2(NDOF,NDOF),GI1(NDOF,NDOF)
    REAL(KIND=DBL)	:: source(2),field(2)
    REAL(KIND=DBL)	:: dxdso(2),jacobo
    REAL(KIND=DBL)	:: temp(2),r

    !...get number of nodes of outer/inner elements
    neo = NODE(ido)
    nei = NODE(idi)

	!...compute egbar and egtwobar (regular integrals)
	egbar = 0.0d0
    egtwobar = 0.0d0

	!...begins outer integral loop
	do lo = 1,nto
    	!...coordinate of integration point
		eita = xio(lo)
        !...obtain value of shape function and its derivatives of 
        call reg_shape(ido,eita,pso,dpso)
        !...compute derivative of position vector of outer element
        dxdso = 0.0d0
        do i = 1,neo
        	do k = 1,2
           		dxdso(k) = dxdso(k) + dpso(i)*xlo(k,i)
           	enddo
        enddo
        !...compute jacobian of outer element
        jacobo = dsqrt(dxdso(1)*dxdso(1) + dxdso(2)*dxdso(2))
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
        	zetapr = xii(li)
            !...transform to remove singularity at tip shape functions
            if (idi.eq.CTIP1) then
            	zeta = 0.5d0*(zetapr+1.0d0)*(zetapr+1.0d0)-1.0d0
            elseif (idi.eq.CTIP2) then
            	zeta = 1.0d0-0.5d0*(1.0d0-zetapr)*(1.0d0-zetapr)
            else
              	print*,'eg_sed_tip_tip: id of inner element is out of range'
                stop
            endif
            !...obtain value of shape function and its derivatives of 
            call reg_shape(idi,zeta,psi,dpsi)
            call tip_dshape(idi,zeta,dtipi)
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
              	print*,'egts_sed_reg_tip: r is too small for separated elements'
                stop
            endif
            !...compute kernels
            call Kernel2(NDOF,region,source,field,C2,G2,U2)
            !...calculate val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                	do k=1,NDOF
                    	do l=1,NDOF
                        	if (idi.eq.CTIP1) then
                            	!...jacobian of outer element: cancel out with D-operator
                        		val_inner_integral1(l,k,j,i) = val_inner_integral1(l,k,j,i) + &
                            		GI1(l,k)*dpso(i)*dtipi(j)*(1.0d0+zetapr)*dlog(r)*wii(li)
                            	val_inner_integral2(l,k,j,i) = val_inner_integral2(l,k,j,i) + &
                            		G2(l,k)*dpso(i)*dtipi(j)*(1.0d0+zetapr)*wii(li)
                            elseif (idi.eq.CTIP2) then
                            	val_inner_integral1(l,k,j,i) = val_inner_integral1(l,k,j,i) + &
                            		GI1(l,k)*dpso(i)*dtipi(j)*(1.0d0-zetapr)*dlog(r)*wii(li)
                            	val_inner_integral2(l,k,j,i) = val_inner_integral2(l,k,j,i) + &
                            		G2(l,k)*dpso(i)*dtipi(j)*(1.0d0-zetapr)*wii(li)
                            else
                              	print*,'egts_sed_reg_tip: id of inner element is out of range'
                                stop
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
                   		egbar(NDOF*(i-1)+k,NDOF*(j-1)+l) = egbar(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		val_inner_integral1(l,k,j,i)*wio(lo)
                       	egtwobar(NDOF*(i-1)+k,NDOF*(j-1)+l) = egtwobar(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		val_inner_integral2(l,k,j,i)*wio(lo)
                    enddo
                enddo
            enddo
        enddo
    enddo !of loop over outer integral
    
    !...element stiffness matrix
    eg = egbar + egtwobar
	End Subroutine egts_sed_reg_tip
    !-----------------------------------------------
    Subroutine egts_sed_reg_reg(region,ido,idi,eg,nto,nti,xio,xii,wio,wii,xlo,xli)
    
    IMPLICIT NONE

    INTEGER, INTENT(IN)		:: ido,idi,nto,nti,region
    REAL(KIND=DBL), INTENT(OUT)	:: eg(:,:)
    REAL(KIND=DBL), INTENT(IN)	:: xio(:),xii(:),wio(:),wii(:)
    REAL(KIND=DBL), INTENT(IN)	:: xlo(:,:),xli(:,:)

    !...local variables
    INTEGER			:: neo,nei,lo,li
    REAL(KIND=DBL)	:: egbar(3*NDOF,3*NDOF),egtwobar(3*NDOF,3*NDOF)
    REAL(KIND=DBL)	:: eita,zeta
    REAL(KIND=DBL)	:: pso(3),dpso(3),psi(3),dpsi(3)
    INTEGER			:: i,j,k,l,ierr
    REAL(KIND=DBL)	:: val_inner_integral1(NDOF,NDOF,3,3),val_inner_integral2(NDOF,NDOF,3,3)
    REAL(KIND=DBL)	:: C2(NDOF,NDOF),G2(NDOF,NDOF),U2(NDOF,NDOF),GI1(NDOF,NDOF)
    REAL(KIND=DBL)	:: source(2),field(2)
    REAL(KIND=DBL)	:: dxdso(2),jacobo
    REAL(KIND=DBL)	:: temp(2),r

    !...get number of nodes of outer/inner elements
    neo = NODE(ido)
    nei = NODE(idi)

	!...compute egbar and egtwobar (regular integrals)
	egbar = 0.0d0
    egtwobar = 0.0d0

	!...begins outer integral loop
	do lo = 1,nto
    	!...coordinate of integration point
		eita = xio(lo)
        !...obtain value of shape function and its derivatives of 
        call reg_shape(ido,eita,pso,dpso)
        !...compute derivative of position vector of outer element
        dxdso = 0.0d0
        do i = 1,neo
        	do k = 1,2
           		dxdso(k) = dxdso(k) + dpso(i)*xlo(k,i)
           	enddo
        enddo
        !...compute jacobian of outer element
        jacobo = dsqrt(dxdso(1)*dxdso(1) + dxdso(2)*dxdso(2))
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
              	print*,'egts_sed_reg_reg: r is too small for separated elements'
                stop
            endif
            !...compute kernels
            call Kernel2(NDOF,region,source,field,C2,G2,U2)
            !...calculate val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                	do k=1,NDOF
                    	do l=1,NDOF
                           	!...jacobian of outer element: cancel out with D-operator
                       		val_inner_integral1(l,k,j,i) = val_inner_integral1(l,k,j,i) + &
                            		GI1(l,k)*dpso(i)*dpsi(j)*dlog(r)*wii(li)
                           	val_inner_integral2(l,k,j,i) = val_inner_integral2(l,k,j,i) + &
                            		G2(l,k)*dpso(i)*dpsi(j)*wii(li)
                        enddo
                    enddo
                enddo
            enddo
        enddo !of loop over inner integral
        do i=1,neo
          	do j=1,nei
            	do k=1,NDOF
                	do l=1,NDOF
                   		egbar(NDOF*(i-1)+k,NDOF*(j-1)+l) = egbar(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		val_inner_integral1(l,k,j,i)*wio(lo)
                       	egtwobar(NDOF*(i-1)+k,NDOF*(j-1)+l) = egtwobar(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		val_inner_integral2(l,k,j,i)*wio(lo)
                    enddo
                enddo
            enddo
        enddo
    enddo !of loop over outer integral
    
    !...element stiffness matrix
    eg = egbar + egtwobar
	End Subroutine egts_sed_reg_reg
!-----------------------------------------------

	Subroutine ehts_reg_tip(ido,idi,eh,nto,nti,xio,xii,wio,wii,xlo,xli)
    !...This subroutine is used for all cases: coincident/adjacent/separated elements
    !...For the case of coincident element => call subroutine with different Gauss orders for outer/inner integrals
    IMPLICIT NONE

    INTEGER, INTENT(IN)			:: ido,idi,nto,nti
    REAL(KIND=DBL), INTENT(OUT)	:: eh(:,:)
    REAL(KIND=DBL), INTENT(IN)	:: xio(:),xii(:),wio(:),wii(:)
    REAL(KIND=DBL), INTENT(IN)	:: xlo(:,:),xli(:,:)

    !...local variables
    INTEGER			:: neo,nei,lo,li
    REAL(KIND=DBL)	:: eita,zeta,zetapr
    REAL(KIND=DBL)	:: pso(3),dpso(3),psi(3),dpsi(3),psitip(3)
    INTEGER			:: i,j,k,l,alpha
    REAL(KIND=DBL)	:: val_inner_integral(NDOF,NDOF,3,3)
    REAL(KIND=DBL)	:: HI(NDOF,NDOF)
    REAL(KIND=DBL)	:: source(2),field(2)
    REAL(KIND=DBL)	:: dxdso(2),dxdsi(2),jacobo,jacobi
    REAL(KIND=DBL)	:: temp(2),r,normal(2)
    LOGICAL			:: debug

	!debug = .true.
    debug = .false.
    
    !...get number of nodes of outer/inner elements
    neo = NODE(ido)
    nei = NODE(idi)
    
    !...evaluate eh (regular integrals)
    eh = 0.0d0
    do lo=1,nto
      	!...coordinate of Gauss point
      	eita = xio(lo)
        !...shape function & its derivatives at Gauss point
        call reg_shape(ido,eita,pso,dpso)
        !...compute derivative of position vector of outer element
        dxdso = 0.0d0
        do i = 1,neo
           	do k = 1,2
               	dxdso(k) = dxdso(k) + dpso(i)*xlo(k,i)
            enddo
        enddo
        !...compute jacobian of outer element
        jacobo = dsqrt(dxdso(1)*dxdso(1) + dxdso(2)*dxdso(2))
        !...position of source point (y)
        source = 0.0d0
        do i=1,neo
          	do k=1,2
            	source(k) = source(k) + xlo(k,i)*pso(i)
            enddo
        enddo
        val_inner_integral = 0.0d0
        !...loop over inner integral
        do li=1,nti
          	!...coordinate of Gauss point
            zetapr = xii(li)
            !...transform back to original zeta
            if (idi.eq.CTIP1) then
                zeta = 0.5d0*(zetapr+1.0d0)*(zetapr+1.0d0)-1.0d0
            elseif (idi.eq.CTIP2) then
                zeta = 1.0d0-0.5d0*(1.0d0-zetapr)*(1.0d0-zetapr)
            else
              	print*,'ehts_reg_tip: inner element must be a tip element'
                stop
            endif
            !...shape function & its derivatives at Gauss point
        	call reg_shape(idi,zeta,psi,dpsi)
            !...also need tip shape functions
            call tip_shape(idi,zeta,psitip)
            
            !...compute derivative of position vector of inner element
        	dxdsi = 0.0d0
        	do i = 1,nei
           		do k = 1,2
               		dxdsi(k) = dxdsi(k) + dpsi(i)*xli(k,i)
            	enddo
        	enddo
        	!...compute jacobian of inner element
        	jacobi = dsqrt(dxdsi(1)*dxdsi(1) + dxdsi(2)*dxdsi(2))
        	!...position of field point (xi)
        	field = 0.0d0
        	do i=1,nei
          		do k=1,2
            		field(k) = field(k) + xli(k,i)*psi(i)
            	enddo
        	enddo
            !...compute r
            do k=1,2
              	temp(k) = field(k) - source(k)
            enddo
            r = dsqrt(temp(1)*temp(1) + temp(2)*temp(2))
			!...defensive programming
            if (r.le.SMALL_NUM) then
              	print*,'ehts_reg_tip: r is too small'
                stop
            endif
            if (jacobi.le.SMALL_NUM) then
              	print*,'ehts_reg_tip: jacobi is too small'
                stop
            endif
            !...compute normal vector n(xi) = s x e3
            temp(1) = dxdsi(2)
            temp(2) = -dxdsi(1)
            normal(1) = temp(1)/jacobi
            normal(2) = temp(2)/jacobi
            if (debug) then
              	print'(a5,f12.8,a1,f12.8,a1)','n = (',normal(1),';',normal(2),')'
            endif
			!...compute val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                	do k=1,NDOF
                    	do l=1,NDOF
                        	!...compute HI(l,k)=H(l,k,alpha).n(alpha)
                            HI = 0.0d0
                            if (k.eq.l) then
                              	do alpha=1,2
                              		HI(l,k) = HI(l,k) - (0.5d0/PI)*(field(alpha)-source(alpha))*normal(alpha)/(r*r)
                                enddo
                        	endif
                        	if (idi.eq.CTIP1) then
                            	!...jacobian of outer element: cancel out with D-operator
                        		val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + &
                            		HI(l,k)*dpso(i)*psitip(j)*(zetapr+1.0d0)*jacobi*wii(li)
                            elseif (idi.eq.CTIP2) then
                              	!...jacobian of outer element: cancel out with D-operator
                              	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + &
                            		HI(l,k)*dpso(i)*psitip(j)*(1.0d0-zetapr)*jacobi*wii(li)
                            else
                              	print*,'ehts_tip_tip: inner element must be a tip element'
                                stop
                            endif
                        enddo
                    enddo
                enddo
            enddo
        enddo !over loop of inner integral
        do i=1,neo
          	do j=1,nei
            	do k=1,NDOF
                	do l=1,NDOF
                    	eh(NDOF*(i-1)+k,NDOF*(j-1)+l) = eh(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        	val_inner_integral(l,k,j,i)*wio(lo)
                    enddo
                enddo
            enddo
        enddo
    enddo !over loop of outer integral
	End Subroutine ehts_reg_tip
!-----------------------------------------------

	Subroutine ehts_reg_reg(ido,idi,eh,nto,nti,xio,xii,wio,wii,xlo,xli)
    !...This subroutine is used for both adjacent/separated elements
    IMPLICIT NONE

    INTEGER, INTENT(IN)		:: ido,idi,nto,nti
    REAL(KIND=DBL), INTENT(OUT)	:: eh(:,:)
    REAL(KIND=DBL), INTENT(IN)	:: xio(:),xii(:),wio(:),wii(:)
    REAL(KIND=DBL), INTENT(IN)	:: xlo(:,:),xli(:,:)

    !...local variables
    INTEGER			:: neo,nei,lo,li
    REAL(KIND=DBL)	:: eita,zeta
    REAL(KIND=DBL)	:: pso(3),dpso(3),psi(3),dpsi(3)
    INTEGER			:: i,j,k,l,alpha
    REAL(KIND=DBL)	:: val_inner_integral(NDOF,NDOF,3,3)
    REAL(KIND=DBL)	:: HI(NDOF,NDOF)
    REAL(KIND=DBL)	:: source(2),field(2)
    REAL(KIND=DBL)	:: dxdso(2),dxdsi(2),jacobo,jacobi
    REAL(KIND=DBL)	:: temp(2),r,normal(2)
    LOGICAL			:: debug

    !debug = .true.
    debug = .false.

    !...get number of nodes of outer/inner elements
    neo = NODE(ido)
    nei = NODE(idi)
    
    !...evaluate eh (regular integrals)
    eh = 0.0d0
    do lo=1,nto
      	!...coordinate of Gauss point
      	eita = xio(lo)
        !...shape function & its derivatives at Gauss point
        call reg_shape(ido,eita,pso,dpso)
        !...compute derivative of position vector of outer element
        dxdso = 0.0d0
        do i = 1,neo
           	do k = 1,2
               	dxdso(k) = dxdso(k) + dpso(i)*xlo(k,i)
            enddo
        enddo
        !...compute jacobian of outer element
        jacobo = dsqrt(dxdso(1)*dxdso(1) + dxdso(2)*dxdso(2))
        !...position of source point (y)
        source = 0.0d0
        do i=1,neo
          	do k=1,2
            	source(k) = source(k) + xlo(k,i)*pso(i)
            enddo
        enddo
        val_inner_integral = 0.0d0
        !...loop over inner integral
        do li=1,nti
          	!...coordinate of Gauss point
            zeta = xii(li)
            !...shape function & its derivatives at Gauss point
        	call reg_shape(idi,zeta,psi,dpsi)
            !...compute derivative of position vector of inner element
        	dxdsi = 0.0d0
        	do i = 1,nei
           		do k = 1,2
               		dxdsi(k) = dxdsi(k) + dpsi(i)*xli(k,i)
            	enddo
        	enddo
        	!...compute jacobian of inner element
        	jacobi = dsqrt(dxdsi(1)*dxdsi(1) + dxdsi(2)*dxdsi(2))
        	!...position of field point (xi)
        	field = 0.0d0
        	do i=1,nei
          		do k=1,2
            		field(k) = field(k) + xli(k,i)*psi(i)
            	enddo
        	enddo
            !...compute r
            do k=1,2
              	temp(k) = field(k) - source(k)
            enddo
            r = dsqrt(temp(1)*temp(1) + temp(2)*temp(2))
			!...defensive programming
            if (r.le.SMALL_NUM) then
              	print*,'ehts_reg_reg: r is too small'
                stop
            endif
            if (jacobi.le.SMALL_NUM) then
              	print*,'ehts_reg_reg: jacobi is too small'
                stop
            endif
            !...compute normal vector n(xi) = s x e3
            temp(1) = dxdsi(2)
            temp(2) = -dxdsi(1)
            normal(1) = temp(1)/jacobi
            normal(2) = temp(2)/jacobi
            if (debug) then
              	print'(a5,f12.8,a1,f12.8,a1)','n = (',normal(1),';',normal(2),')'
            endif
			!...compute val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                	do k=1,NDOF
                    	do l=1,NDOF
                        	!...compute H.n
                            HI = 0.0d0
                            if (k.eq.l) then
                              	do alpha=1,2
                              		HI(l,k) = HI(l,k) - (0.5d0/PI)*(field(alpha)-source(alpha))*normal(alpha)/(r*r)
                                enddo
                        	endif
                           	!...jacobian of outer element: cancel out with D-operator
                       		val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + &
                            		HI(l,k)*dpso(i)*psi(j)*jacobi*wii(li)
                        enddo
                    enddo
                enddo
            enddo
        enddo !over loop of inner integral
        do i=1,neo
          	do j=1,nei
            	do k=1,NDOF
                	do l=1,NDOF
                    	eh(NDOF*(i-1)+k,NDOF*(j-1)+l) = eh(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        	val_inner_integral(l,k,j,i)*wio(lo)
                    enddo
                enddo
            enddo
        enddo
    enddo !over loop of outer integral
	End Subroutine ehts_reg_reg
End Subroutine eatts
