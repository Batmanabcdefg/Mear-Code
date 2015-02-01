Subroutine ek_d(region,isrc,ifld,xy_src,xy_fld,ekd)

!...This subroutine is changed for the formulation, Saumik

!...In the original Han code, the kernal is divided into 2 parts 
!...1) the term that multiplies with ln(r), CI1
!...2) the terms that does not multiply with ln(r), C2
!...3) both these terms multiply with derivatives of shape functions

!...(CI1*ln(r)+C2)*D \delta u *D \delta u

!...In this formulation, Saumik 
!...1) CI1 gets replaced by B1 and C2 gets replaced with B2. Both these terms multiply with derivatives of shape functions.
!...2) There is E1 which multiplies with ln(r). There is E2 which does not. So the exact same terms that multiply with B1 get multiplied with E1. 
!...The exact same terms that multiply with B2 get multiplied with E2. And finally both these terms multiply with shape functions.

!...(B1*ln(r)+B2)*D \delta u *D \delta u + (E1*ln(r)+E2)*\delta u *\delta u

!...In this formulation, Saumik
!...1) 3 dimensional kronecker delta terms are used, E1 and E2
!...2) The node normals n_1,n_2 are used, E1 and E2
!...3) w_1,w_2,x_1,x_2,y,z_1,z_2 are functions of r and ht. They are used in B1,B2,E1,E2

!...The only difference with ek_d.f90: fix the bug in subroutine ekd_ii_tip_tip for ln(2r/(zeta-eita))
!...10/26/09: difference with ek_E1.f90: use 2ND approach (see note...) for ekd_adt_reg_reg
!...Subroutine to calculate element matrix ekd
!...All sub-subroutines are from unbounded code, just add the word "ekd_" at the begining
!...9/30/10: different with ek_E2.f90: improve adj_reg_tip and adj_tip_reg to avoid regular Gauss 
!...integration for ln r	
	
	USE DefinitionConstant
	USE GlobalData
    IMPLICIT NONE
    
	INTEGER, INTENT(IN)			:: isrc,ifld,region
    REAL(KIND=DBL),INTENT(IN)	:: xy_src(:,:),xy_fld(:,:)
	
	REAL(KIND=DBL),INTENT(OUT)	:: ekd(:,:)
    !...local variables
    INTEGER		                :: adjacent_type,matl_no,matl_type
	REAL(KIND=DBL)              :: ey,nu,phi
	LOGICAL		                :: debug
    
	!debug = .true.
    debug = .false.
    
	if (debug) then
    	print*,'ek_d called for elements: ',isrc,ifld
    endif

	!...we need to feed this everytime ek_d is invoked, since region number might be different, and we don't know when that changes
	!...id # of material
    matl_no = region_mat_no(region)
    !...type of material
    matl_type = material_type(matl_no)
	!...set material properties and necessary information

	!...we haven't developed ek_d.f90 for anything other than isotropy,so it is assumed that material is isotropic, two independent material 
	!...constants
	if ((NDOF.eq.3).or.(NDOF.eq.2)) then
        !...for isotropy, plane strain and antiplane are uncoupled, and the closed form of
        !...kernels are independent
        	ey = material(1,matl_no)   !...Young's modulus
        	nu = material(2,matl_no)   !...Poisson ratio
    else
          	print*,'not develop yet isotropy for other media than elasticity'
            stop
    endif

	phi = ey/(8.0d0*pi*(1.0d0-nu*nu))

!----------------------------------------------------    
	if (isrc.eq.ifld) then
      	!...coincident elements
		!...can happen D_{cc} in (3)
        if (ELTYPE(elemid(isrc)).eq.CTIP) then
          	!...tip-tip
            call ekd_ii_tip_tip(region,elemid(isrc),ekd,ng_out,ng_ln,ng_in,xi(:,ng_out),xi_log(:,ng_ln),xi(:,ng_in),&
            				wi(:,ng_out),wi_log(:,ng_ln),wi(:,ng_in),xy_src,phi,ht,nu)
							!...added 2 arguments
            !ekd = -1.0d0*ekd
        !...4.0d0/28.0d0/09: add INTERFE
        elseif ((ELTYPE(elemid(isrc)).eq.CREGULAR).or.(ELTYPE(elemid(isrc)).eq.BREGULAR).or.&
        (ELTYPE(elemid(isrc)).eq.INTERFE)) then
        	!...reg-reg
            call ekd_ii_reg_reg(region,elemid(isrc),ekd,ng_out,ng_ln,ng_in,xi(:,ng_out),xi_log(:,ng_ln),xi(:,ng_in),&
            				wi(:,ng_out),wi_log(:,ng_ln),wi(:,ng_in),xy_src,phi,ht,nu)
            !ekd = -1.0d0*ekd
        else
          	print*,'ek_d.f90: current version does not support elem type: ',elemid(isrc)
            stop
        endif
!----------------------------------------------------------------------------------------------------
    elseif (adjacent(isrc,ifld).gt.0) then
    	!...adjacent elements, continue to determine which type of adjacent
        if (elnode(1,ifld).eq.elnode(2,isrc)) then
          	!...singular at (eita = 1; zeta = -1)
          	adjacent_type = 1
        elseif (elnode(1,isrc).eq.elnode(2,ifld)) then
        	!...singular at (eita = -1; zeta = 1)
        	adjacent_type = 2
        !...add 2 more cases for the case of SBLP
        elseif (elnode(1,isrc).eq.elnode(1,ifld)) then
        	!...singular at (eita = -1; zeta = -1)
            adjacent_type = 3
        elseif (elnode(2,isrc).eq.elnode(2,ifld)) then
        	!...singular at (eita = 1; zeta = 1)
            adjacent_type = 4
        else
          	print*,'ek_d.f90: error on connectivity of elements: ',elem_sys2user(isrc),elem_sys2user(ifld)
            stop
        endif
		!...adj_type = 1 or 4 - singular at eita = 1
		!...adj_type = 2 or 3 - singular at zeta = 1
        if (ELTYPE(elemid(isrc)).eq.CTIP.and.ELTYPE(elemid(ifld)).eq.CTIP) then
          	!...adjacent: tip-tip, same integration order
			!...adjacent elements, tip elements
			!...this still implies we are working with D_{cc}, now E_{cc}
            call ekd_adt_tip_tip(region,elemid(isrc),elemid(ifld),ekd,ng_out,ng_out,xi(:,ng_out),xi(:,ng_out),&
            				wi(:,ng_out),wi(:,ng_out),xy_src,xy_fld,adjacent_type,phi,ht,nu)
            !ekd = -1.0d0*ekd
        elseif (ELTYPE(elemid(isrc)).eq.CTIP.and.ELTYPE(elemid(ifld)).ne.CTIP) then
        	!...adjacent: tip-reg, same integration order
            call ekd_adt_tip_reg(region,elemid(isrc),elemid(ifld),ekd,ng_in,ng_ln,ng_out,xi(:,ng_in),xi_log(:,ng_ln),xi(:,ng_out),&
            				wi(:,ng_in),wi_log(:,ng_ln),wi(:,ng_out),xy_src,xy_fld,adjacent_type,phi,ht,nu)
			!...Module 3.3.3.5 in Han's thesis,crack tip source element,regular field element
            !ekd = -1.0d0*ekd
		elseif (ELTYPE(elemid(isrc)).ne.CTIP.and.ELTYPE(elemid(ifld)).eq.CTIP) then
        	!...adjacent: reg-tip, same integration order
            call ekd_adt_reg_tip(region,elemid(isrc),elemid(ifld),ekd,ng_out,ng_ln,ng_in,xi(:,ng_out),xi_log(:,ng_ln),xi(:,ng_in),&
            				wi(:,ng_out),wi_log(:,ng_ln),wi(:,ng_in),xy_src,xy_fld,adjacent_type,phi,ht,nu)
            !ekd = -1.0d0*ekd
		elseif (ELTYPE(elemid(isrc)).ne.CTIP.and.ELTYPE(elemid(ifld)).ne.CTIP) then
        	!...adjacent: reg-reg, same integration order
            call ekd_adt_reg_reg(region,elemid(isrc),elemid(ifld),ekd,ng_out,ng_ln,ng_out,xi(:,ng_out),xi_log(:,ng_ln),xi(:,ng_out),&
            				wi(:,ng_out),wi_log(:,ng_ln),wi(:,ng_out),xy_src,xy_fld,adjacent_type,phi,ht,nu)
            !ekd = -1.0d0*ekd
        endif
!----------------------------------------------------------------------------------------------------
    else
      	!...separated elements
        if (ELTYPE(elemid(isrc)).eq.CTIP.and.ELTYPE(elemid(ifld)).eq.CTIP) then
          	!...separate: tip-tip, same integration order
            call ekd_sed_tip_tip(region,elemid(isrc),elemid(ifld),ekd,ng_out,ng_out,xi(:,ng_out),xi(:,ng_out),&
            				wi(:,ng_out),wi(:,ng_out),xy_src,xy_fld,phi,ht,nu)
            !ekd = -1.0d0*ekd
        elseif (ELTYPE(elemid(isrc)).eq.CTIP.and.ELTYPE(elemid(ifld)).ne.CTIP) then
        	!...separate: tip-reg, same integration order
            call ekd_sed_tip_reg(region,elemid(isrc),elemid(ifld),ekd,ng_out,ng_out,xi(:,ng_out),xi(:,ng_out),&
            				wi(:,ng_out),wi(:,ng_out),xy_src,xy_fld,phi,ht,nu)
            !ekd = -1.0d0*ekd
		elseif (ELTYPE(elemid(isrc)).ne.CTIP.and.ELTYPE(elemid(ifld)).eq.CTIP) then
        	!...separate: reg-tip, same integration order
            call ekd_sed_reg_tip(region,elemid(isrc),elemid(ifld),ekd,ng_out,ng_out,xi(:,ng_out),xi(:,ng_out),&
            				wi(:,ng_out),wi(:,ng_out),xy_src,xy_fld,phi,ht,nu)
            !ekd = -1.0d0*ekd
		elseif (ELTYPE(elemid(isrc)).ne.CTIP.and.ELTYPE(elemid(ifld)).ne.CTIP) then
        	!...separate: reg-reg, same integration order
            call ekd_sed_reg_reg(region,elemid(isrc),elemid(ifld),ekd,ng_out,ng_out,xi(:,ng_out),xi(:,ng_out),&
            				wi(:,ng_out),wi(:,ng_out),xy_src,xy_fld,phi,ht,nu)
            !ekd = -1.0d0*ekd
        endif
    endif
    !----------------------------------------------------------------------------------------------------
	CONTAINS
	Subroutine ekd_ii_tip_tip(region,idelem,ek,nto,ntoln,nti,xio,xioln,xii,wio,wioln,wii,xlelem,phi,ht,nu)
    IMPLICIT NONE

	INTEGER, INTENT(IN)		    :: idelem,nto,ntoln,nti,region
    REAL(KIND=DBL), INTENT(OUT)	:: ek(:,:)
    REAL(KIND=DBL), INTENT(IN)	:: xio(:),xioln(:),xii(:),wio(:),wioln(:),wii(:)
	REAL(KIND=DBL), INTENT(IN)	:: xlelem(:,:)
	REAL(KIND=DBL), INTENT(IN)	:: phi,ht,nu
	!...local variables
    INTEGER		    :: neo,nei,lo,li
    REAL(KIND=DBL)	:: ekbar1(3*NDOF,3*NDOF),ekbar2(3*NDOF,3*NDOF),ek2bar(3*NDOF,3*NDOF)
    REAL(KIND=DBL)	:: ekbar31(3*NDOF,3*NDOF),ekbar32(3*NDOF,3*NDOF)
    REAL(KIND=DBL)	:: alpha,alphapr,beta,eita,zeta,eitapr,zetapr
    REAL(KIND=DBL)  :: rho,gamma,deltaa
	REAL(KIND=DBL)	:: tipo(3),dtipo(3),tipi(3),dtipi(3),pso(3),dpso(3),psi(3),dpsi(3)
    REAL(KIND=DBL)	:: tipo1(3),tipi1(3),tipo2(3),tipi2(3),pso1(3),dpso1(3),psi1(3),dpsi1(3)
    REAL(KIND=DBL)	:: dtipo1(3),dtipi1(3),dtipo2(3),dtipi2(3),pso2(3),dpso2(3),psi2(3),dpsi2(3)
	INTEGER			:: i,j,k,l
    REAL(KIND=DBL)	:: B1(NDOF,NDOF),B2(NDOF,NDOF)
	REAL(KIND=DBL)	:: val_inner_integral(NDOF,NDOF,3,3)
    REAL(KIND=DBL)	:: val_inner_integral1(NDOF,NDOF,3,3),val_inner_integral2(NDOF,NDOF,3,3)
    REAL(KIND=DBL)	:: f1(NDOF,NDOF,3,3),f2(NDOF,NDOF,3,3)
	REAL(KIND=DBL)	:: E1(NDOF,NDOF),E2(NDOF,NDOF)
	REAL(KIND=DBL)	:: field(2),source(2),r,temp(2)
    
	
	 
	!------------------------------------------------------------------------------------------------
	neo = NODE(idelem)
	!...no of nodes of outer element
    nei = NODE(idelem)
	!...no of nodes of inner element
    !------------------------------------------------------------------------------------------------
	!...compute (weakly) singular integral ekbar2
	!...Evaluating term 2 in Eqn (3.76)
	ekbar2 = 0.0d0
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
				!...eita on the source element
                zeta = 1.0d0-0.5d0*(1.0d0-zetapr)*(1.0d0-zetapr)
				!...zeta on the field element
            else
                print*,'ekd_ii_tip_tip: id of element is out of range'
                stop
            endif
            
			call tip_dshape(idelem,eita,dtipo1)
		        call tip_shape(idelem,eita,tipo1)
			call reg_shape(idelem,eita,pso1,dpso1)
			call tip_dshape(idelem,zeta,dtipi1)
		        call tip_shape(idelem,zeta,tipi1)
			call reg_shape(idelem,zeta,psi1,dpsi1)
			


			call findB1(phi,B1,xlelem,xlelem,pso1,dpso1,psi1,dpsi1)
			call findE1(phi,nu,E1,xlelem,xlelem,pso1,dpso1,psi1,dpsi1)
			
                        f1 = 0.0d0
			!...this initialization wasn't there before
			do i=1,neo
              	do j=1,nei
                	do k=1,NDOF
                    	do l=1,NDOF
                            if (idelem.eq.CTIP1) then
                                 f1(l,k,j,i) = f1(l,k,j,i) + (eitapr+1.0d0)*(zetapr+1.0d0)*(dtipo1(i)*dtipi1(j)*B1(l,k) + tipo1(i)*tipi1(j)*E1(l,k))

			    else
                                 f1(l,k,j,i) = f1(l,k,j,i) + (-eitapr+1.0d0)*(-zetapr+1.0d0)*(dtipo1(i)*dtipi1(j)*B1(l,k) + tipo1(i)*tipi1(j)*E1(l,k))
								 
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
                print*,'ekd_ii_tip_tip: id of element is out of range'
                stop
            endif

			call tip_dshape(idelem,eita,dtipo2)
			call tip_shape(idelem,eita,tipo2)
			call reg_shape(idelem,eita,pso2,dpso2)
			call tip_dshape(idelem,zeta,dtipi2)
			call tip_shape(idelem,zeta,tipi2)
			call reg_shape(idelem,zeta,psi2,dpsi2)
			

			
			call findB1(phi,B1,xlelem,xlelem,pso2,dpso2,psi2,dpsi2)
			call findE1(phi,nu,E1,xlelem,xlelem,pso2,dpso2,psi2,dpsi2)
			
                        f2 = 0.0d0
			!...this initialization wasn't there before
			!...compute value of f(alphapr,-beta)
            do i=1,neo
              	do j=1,nei
                	do k=1,NDOF
                    	do l=1,NDOF
                            if (idelem.eq.CTIP1) then
                        	    f2(l,k,j,i) = f2(l,k,j,i) + (eitapr+1.0d0)*(zetapr+1.0d0)*(dtipo2(i)*dtipi2(j)*B1(l,k) + tipo2(i)*tipi2(j)*E1(l,k))
			    else
                        	    f2(l,k,j,i) = f2(l,k,j,i) + (-eitapr+1.0d0)*(-zetapr+1.0d0)*(dtipo2(i)*dtipi2(j)*B1(l,k) + tipo2(i)*tipi2(j)*E1(l,k))
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
                    	ekbar2(NDOF*(i-1)+k,NDOF*(j-1)+l) = ekbar2(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                            (1.0d0-beta)*val_inner_integral(l,k,j,i)*wioln(lo)
                    enddo
                enddo
            enddo
        enddo
    enddo !of loop over outer integral
    !...finally, multiply with 2.0 (Note: int(-ln(x)*f(x),x) = sum(w(i)*f(x(i)))
    ekbar2 = -2.0d0*ekbar2
	!...Evaluated term 2 in Eqn (3.76)
    !--------------------------------------------------
    !...evaluate ekbar1 and ek2bar (regular integrals)
	!...evaluating term 1 in Eqn (3.76)
    ekbar1 = 0.0d0
    ek2bar = 0.0d0
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
        call tip_dshape(idelem,eita,dtipo)
	call tip_shape(idelem,eita,tipo)
        call reg_shape(idelem,eita,pso,dpso)
        

				
	!...loop over inner integral
        val_inner_integral1 = 0.0d0
        val_inner_integral2 = 0.0d0
        
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
                        call tip_dshape(idelem,zeta,dtipi)
			call tip_shape(idelem,zeta,tipi)
        	        call reg_shape(idelem,zeta,psi,dpsi)
        	

			call findB1(phi,B1,xlelem,xlelem,pso,dpso,psi,dpsi)
			call findB2(phi,B2,xlelem,xlelem,pso,dpso,psi,dpsi)
			
                        call findE1(phi,nu,E1,xlelem,xlelem,pso,dpso,psi,dpsi)
                        call findE2(phi,nu,E2,xlelem,xlelem,pso,dpso,psi,dpsi)

                        !...position of source point (y)
                        source = 0.0d0
                        do i=1,neo
          	            do k=1,2
            	                source(k) = source(k) + xlelem(k,i)*pso(i)
                            enddo
                         enddo

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
			
			
			!...compute val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                	do k=1,NDOF
                    	do l=1,NDOF
                            if (idelem.eq.CTIP1) then
                        	val_inner_integral1(l,k,j,i) = val_inner_integral1(l,k,j,i) + (B1(l,k)*dtipo(i)*dtipi(j) + E1(l,k)*tipo(i)*tipi(j))*(eitapr+1.0d0)*(zetapr+1.0d0)*dlog(4.0d0*r/dabs(zetapr-eitapr)/dabs(2.0d0+zetapr+eitapr))*wii(li)

								
				val_inner_integral2(l,k,j,i) = val_inner_integral2(l,k,j,i) + (B2(l,k)*dtipo(i)*dtipi(j) + E2(l,k)*tipo(i)*tipi(j))*(eitapr+1.0d0)*(zetapr+1.0d0)*wii(li)

										
			    else
                              	val_inner_integral1(l,k,j,i) = val_inner_integral1(l,k,j,i) + (B1(l,k)*dtipo(i)*dtipi(j) + E1(l,k)*tipo(i)*tipi(j))*(1.0d0-eitapr)*(1.0d0-zetapr)*dlog(4.0d0*r/dabs(zetapr-eitapr)/dabs(2.0d0-zetapr-eitapr))*wii(li)

								
				val_inner_integral2(l,k,j,i) = val_inner_integral2(l,k,j,i) + (B2(l,k)*dtipo(i)*dtipi(j) + E2(l,k)*tipo(i)*tipi(j))*(1.0d0-eitapr)*(1.0d0-zetapr)*wii(li)

								
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
                    	ekbar1(NDOF*(i-1)+k,NDOF*(j-1)+l) = ekbar1(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        	val_inner_integral1(l,k,j,i)*wio(lo)
                        ek2bar(NDOF*(i-1)+k,NDOF*(j-1)+l) = ek2bar(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        	val_inner_integral2(l,k,j,i)*wio(lo)
                    enddo
                enddo
            enddo
        enddo
    enddo !over loop of outer integral
	!...evaluated term 1 in Eqn (3.76)
	!...calculate ekbar31 and ekbar32
	!...evaluating term 1 in Eqn (3.78.0d0)
	ekbar31 = 0.0d0
	!...begins outer integral loop
	do lo = 1,nto
    	!...coordinate of integration point
        if (idelem.eq.CTIP1) then
			rho = xio(lo)
        	gamma = 1.0d0 - 0.5d0*(1.0d0-rho)*(1.0d0-rho)
        elseif (idelem.eq.CTIP2) then
        	gamma = xio(lo)
        else
          	print*,'ek_E1.f90,ii_tip_tip: wrong id of tip element'
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
              	print*,'ek_E1.f90,ii_tip_tip: id of tip element is out of range'
                stop
            endif
            
			call tip_shape(idelem,eita,tipo)
			call tip_shape(idelem,zeta,tipi)
			call tip_dshape(idelem,eita,dtipo)
                        call tip_dshape(idelem,zeta,dtipi)
        	        call reg_shape(idelem,eita,pso,dpso)
			call reg_shape(idelem,zeta,psi,dpsi)

			
			
			call findB1(phi,B1,xlelem,xlelem,pso,dpso,psi,dpsi)
			call findE1(phi,nu,E1,xlelem,xlelem,pso,dpso,psi,dpsi)
			
			!...calculate val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                	do k=1,NDOF
                    	do l=1,NDOF
                            if (idelem.eq.CTIP1) then
                        	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + (B1(l,k)*dtipo(i)*dtipi(j) + E1(l,k)*tipo(i)*tipi(j))*(1.0d0+eitapr)*(1.0d0+zetapr)*wii(li)

								
			    elseif (idelem.eq.CTIP2) then
                              	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + (B1(l,k)*dtipo(i)*dtipi(j) + E1(l,k)*tipo(i)*tipi(j))*(1.0d0-eitapr)*(1.0d0-zetapr)*wii(li)

								
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
                       		ekbar31(NDOF*(i-1)+k,NDOF*(j-1)+l) = ekbar31(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		((1.0d0-rho)**3.0d0)*dlog(0.5d0*(1.0d0-rho))*val_inner_integral(l,k,j,i)*wio(lo)
                        elseif (idelem.eq.CTIP2) then
                        	ekbar31(NDOF*(i-1)+k,NDOF*(j-1)+l) = ekbar31(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		(1.0d0-gamma)*dlog(0.5d0*(3.0d0+gamma))*val_inner_integral(l,k,j,i)*wio(lo)
                        endif
                    enddo
                enddo
            enddo
        enddo
    enddo !of loop over outer integral
    ekbar31 = 0.5d0*ekbar31
	!...evaluated term 1 in Eqn (3.78.0d0)
    !----------for ekbar32
	!...evaluating term 2 in Eqn (3.78.0d0)
    ekbar32 = 0.0d0
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
              	print*,'ek_E1.f90,ii_tip_tip: id of tip element is out of range'
                stop
            endif
            
			call tip_shape(idelem,eita,tipo)
			call tip_shape(idelem,zeta,tipi)
			call reg_shape(idelem,eita,pso,dpso)
			call tip_dshape(idelem,eita,dtipo)
                        call tip_dshape(idelem,zeta,dtipi)
			call reg_shape(idelem,zeta,psi,dpsi)
			

			
			call findB1(phi,B1,xlelem,xlelem,pso,dpso,psi,dpsi)
			call findE1(phi,nu,E1,xlelem,xlelem,pso1,dpso1,psi1,dpsi1)
		    
			!...calculate val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                	do k=1,NDOF
                    	do l=1,NDOF
                            if (idelem.eq.CTIP1) then
                        	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + (B1(l,k)*dtipo(i)*dtipi(j) + E1(l,k)*tipo(i)*tipi(j))*(1.0d0+eitapr)*(1.0d0+zetapr)*wii(li)

								
			    elseif (idelem.eq.CTIP2) then
                              	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + (B1(l,k)*dtipo(i)*dtipi(j) + E1(l,k)*tipo(i)*tipi(j))*(1.0d0-eitapr)*(1.0d0-zetapr)*wii(li)

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
                       		ekbar32(NDOF*(i-1)+k,NDOF*(j-1)+l) = ekbar32(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		(1.0d0+gamma)*dlog(0.5d0*(3.0d0-gamma))*val_inner_integral(l,k,j,i)*wio(lo)
                        elseif (idelem.eq.CTIP2) then
                        	ekbar32(NDOF*(i-1)+k,NDOF*(j-1)+l) = ekbar32(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		((1.0d0+rho)**3.0d0)*dlog(0.5d0*(1.0d0+rho))*val_inner_integral(l,k,j,i)*wio(lo)
                        endif
                    enddo
                enddo
            enddo
        enddo
    enddo !of loop over outer integral
    ekbar32 = 0.5d0*ekbar32
    !...evaluated term 2 in Eqn (3.78.0d0)
    !...element stiffness matrix
    ek = ekbar1 + ekbar2 + ek2bar + ekbar31 + ekbar32

	End Subroutine ekd_ii_tip_tip
	!----------------------------------------------------------------------------------------------------
	Subroutine ekd_ii_reg_reg(region,idelem,ek,nto,ntoln,nti,xio,xioln,xii,wio,wioln,wii,xlelem,phi,ht,nu)
    IMPLICIT NONE

    INTEGER, INTENT(IN)		    :: idelem,nto,ntoln,nti,region
    REAL(KIND=DBL), INTENT(OUT)	:: ek(:,:)
    REAL(KIND=DBL), INTENT(IN)	:: xio(:),xioln(:),xii(:),wio(:),wioln(:),wii(:)
	REAL(KIND=DBL), INTENT(IN)	:: xlelem(:,:)
	REAL(KIND=DBL), INTENT(IN)	:: phi,ht,nu
	!...local variables
    INTEGER		:: neo,nei,lo,li
	REAL(KIND=DBL)	:: ekbar1(3*NDOF,3*NDOF),ekbar2(3*NDOF,3*NDOF),ek2bar(3*NDOF,3*NDOF)
    REAL(KIND=DBL)	:: alpha,alphapr,beta,eita,zeta
    REAL(KIND=DBL)	:: pso(3),dpso(3),psi(3),dpsi(3)
    REAL(KIND=DBL)	:: pso1(3),dpso1(3),psi1(3),dpsi1(3)
    REAL(KIND=DBL)	:: pso2(3),dpso2(3),psi2(3),dpsi2(3)
    INTEGER			:: i,j,k,l
    REAL(KIND=DBL)	:: B1(NDOF,NDOF),B2(NDOF,NDOF),B11(NDOF,NDOF),B12(NDOF,NDOF)
    REAL(KIND=DBL)	:: val_inner_integral(NDOF,NDOF,3,3)
    REAL(KIND=DBL)	:: val_inner_integral1(NDOF,NDOF,3,3),val_inner_integral2(NDOF,NDOF,3,3)
   
	REAL(KIND=DBL)	:: E1(NDOF,NDOF),E2(NDOF,NDOF),E11(NDOF,NDOF),E12(NDOF,NDOF)
    REAL(KIND=DBL)	:: field(2),source(2),r,temp(2)
	
	neo = NODE(idelem)
    nei = NODE(idelem)
	!...compute (weakly) singular integral ekbar2
	ekbar2 = 0.0d0

	!...begins outer integral loop (logarith integration)
	do lo = 1,ntoln
    	!...coordinate of integration point
	      beta = xioln(lo)
		!...initialize
              val_inner_integral = 0.0d0
		!...loop over inner integral (Gauss integration)
	      do li = 1,nti
			!...coordinate of Gauss point
        	    alphapr = xii(li)
        	!...transformation back to original variables
                    alpha = (1.0d0 - beta)*alphapr
            
		    eita = alpha + beta
                    zeta = alpha - beta
            
			!...obtain value of shape function and its derivatives of (original) outer integration for f_klij(alpha,beta)
                    call reg_shape(idelem,eita,pso1,dpso1)
            	
            !...obtain value of shape function and its derivatives of (original) inner integration for f_klij(alpha,beta)
                    call reg_shape(idelem,zeta,psi1,dpsi1)
		

			
	    !...this is for beta=beta
                    call findB1(phi,B11,xlelem,xlelem,pso1,dpso1,psi1,dpsi1)
	            call findE1(phi,nu,E11,xlelem,xlelem,pso1,dpso1,psi1,dpsi1)
	    !-------------------------------------
                    eita = alpha - beta
                    zeta = alpha + beta
            
			!...obtain value of shape function and its derivatives of (original) outer integration for f_klij(alpha,-beta)
                    call reg_shape(idelem,eita,pso2,dpso2)
            	
            !...obtain value of shape function and its derivatives of (original) inner integration for f_klij(alpha,-beta)
                    call reg_shape(idelem,zeta,psi2,dpsi2)
            	
            !...calculate val_inner_integral(l,k,j,i)
            

			
	            call findB1(phi,B12,xlelem,xlelem,pso2,dpso2,psi2,dpsi2)
	            call findE1(phi,nu,E12,xlelem,xlelem,pso2,dpso2,psi2,dpsi2)
			
	            do i=1,neo
              	        do j=1,nei
                            do k=1,NDOF
                    	        do l=1,NDOF
                                    val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + (B11(l,k)*dpso1(i)*dpsi1(j) + B12(l,k)*dpso2(i)*dpsi2(j) + E11(l,k)*pso1(i)*psi1(j) + & 
							          E12(l,k)*pso2(i)*psi2(j))*wii(li)

			        enddo
                            enddo
                        enddo
                    enddo
               enddo !of loop over inner integral
               do i=1,neo
          	   do j=1,nei
            	       do k=1,NDOF
                	   do l=1,NDOF
                    	       ekbar2(NDOF*(i-1)+k,NDOF*(j-1)+l) = ekbar2(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                                                                   (1.0d0-beta)*val_inner_integral(l,k,j,i)*wioln(lo)
                           enddo
                       enddo
                   enddo
               enddo
         enddo !of loop over outer integral
    !...finally, multiply with 2.0 (Note: int(-ln(x)*f(x),x) = sum(w(i)*f(x(i)))
    ekbar2 = -2.0d0*ekbar2
    !--------------------------------------------------
    !...evaluate ekbar1 and ek2bar (regular integrals)
    ekbar1 = 0.0d0
    ek2bar = 0.0d0
    do lo=1,nto
      	!...coordinate of Gauss point
      	eita = xio(lo)
        
        !...shape function & its derivatives at Gauss point
        call reg_shape(idelem,eita,pso,dpso)
        

				
		!...loop over inner integral
        val_inner_integral1 = 0.0d0
        val_inner_integral2 = 0.0d0
        do li=1,nti
          	!...coordinate of Gauss point
            zeta = xii(li)

            !...shape function & its derivatives at Gauss point
             call reg_shape(idelem,zeta,psi,dpsi)
                

			
	     call findB1(phi,B1,xlelem,xlelem,pso,dpso,psi,dpsi)
			
	     call findB2(phi,B2,xlelem,xlelem,pso,dpso,psi,dpsi)
	     call findE1(phi,nu,E1,xlelem,xlelem,pso,dpso,psi,dpsi)
             call findE2(phi,nu,E2,xlelem,xlelem,pso,dpso,psi,dpsi)

             !...position of source point (y)
                        source = 0.0d0
                        do i=1,neo
          	            do k=1,2
            	                source(k) = source(k) + xlelem(k,i)*pso(i)
                            enddo
                         enddo

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
			
			!...compute val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                    do k=1,NDOF
                    	do l=1,NDOF
                            val_inner_integral1(l,k,j,i) = val_inner_integral1(l,k,j,i) + (B1(l,k)*dpso(i)*dpsi(j) + E1(l,k)*pso(i)*psi(j))*dlog(2.0d0*r/dabs(zeta-eita))*wii(li)
                            val_inner_integral2(l,k,j,i) = val_inner_integral2(l,k,j,i) + (B2(l,k)*dpso(i)*dpsi(j) + E2(l,k)*pso(i)*psi(j))*wii(li)

			enddo
                    enddo
                enddo
            enddo
        enddo !over loop of inner integral
        do i=1,neo
          	do j=1,nei
            	do k=1,NDOF
                	do l=1,NDOF
                    	ekbar1(NDOF*(i-1)+k,NDOF*(j-1)+l) = ekbar1(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        	val_inner_integral1(l,k,j,i)*wio(lo)
                        ek2bar(NDOF*(i-1)+k,NDOF*(j-1)+l) = ek2bar(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        	val_inner_integral2(l,k,j,i)*wio(lo)
                    enddo
                enddo
            enddo
        enddo
    enddo !over loop of outer integral

    !...element stiffness matrix
    ek = ekbar1 + ekbar2 + ek2bar	  
	End Subroutine ekd_ii_reg_reg
	!----------------------------------------------------------------------------------------------------
    Subroutine ekd_adt_reg_reg(region,ido,idi,ek,nto,ntoln,nti,xio,xioln,xii,wio,wioln,wii,xlo,xli,adt_type,phi,ht,nu)
	!...10/27/09: use 2ND approach to improve integration of regular adjacent elements
    IMPLICIT NONE

    INTEGER, INTENT(IN)		:: ido,idi,nto,ntoln,nti,region
    INTEGER, INTENT(IN)		:: adt_type
    REAL(KIND=DBL), INTENT(OUT)	:: ek(:,:)
    REAL(KIND=DBL), INTENT(IN)	:: xio(:),xioln(:),xii(:),wio(:),wioln(:),wii(:)
    REAL(KIND=DBL), INTENT(IN)	:: xlo(:,:),xli(:,:)
	REAL(KIND=DBL), INTENT(IN)	:: phi,ht,nu
	!...local variables
    INTEGER			:: neo,nei,lo,li
    REAL(KIND=DBL)	:: ekbar1(3*NDOF,3*NDOF),ek2bar(3*NDOF,3*NDOF)
    REAL(KIND=DBL)	:: ekbar21(3*NDOF,3*NDOF),ekbar22(3*NDOF,3*NDOF)
    REAL(KIND=DBL)	:: alpha,beta,eita,zeta,gamma,deltaa,rho,eitabar
    REAL(KIND=DBL)	:: pso(3),dpso(3),psi(3),dpsi(3)
    INTEGER			:: i,j,k,l
    REAL(KIND=DBL)	:: val_inner_integral(NDOF,NDOF,3,3)
    
    REAL(KIND=DBL)	:: B1(NDOF,NDOF),B2(NDOF,NDOF)
	REAL(KIND=DBL)	:: E1(NDOF,NDOF),E2(NDOF,NDOF)
	REAL(KIND=DBL)	:: field(2),source(2),r,temp(2)
	
	neo = NODE(ido)
    nei = NODE(idi)
	!...compute matrix ekbar21 and ekbar22
	ekbar21 = 0.0d0
    ekbar22 = 0.0d0
	!...logarith integration: ekbar22 for adt_type = 1/3 and ekbar21 for adt_type = 2/4.0d0
	!...here we are evaluating I_4.0d0 in Eqn (3.91), we use log Gauss integration for outer element
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
			!...Eqn (3.90)(2)
                if ((adt_type.eq.1).or.(adt_type.eq.3)) then
            	beta = 1.0d0 - gamma
                elseif ((adt_type.eq.2).or.(adt_type.eq.4)) then
            	beta = gamma - 1.0d0
				!...Eqn (3.90)(1)
                else
              	print*,'ekd_adt_reg_reg: wrong adt_type:',adt_type
                stop
                endif
            !...distinguish adjacent types between 2/4.0d0 and 1/3
                if ((adt_type.eq.4).or.(adt_type.eq.3)) then
              	eitabar = alpha + beta
                eita = -1.0d0*eitabar
                else
            	eita = alpha + beta
				!...Eqn (3.8.0d06)(1)
                endif
                zeta = alpha - beta
			!...Eqn (3.8.0d06)(2)
            !...obtain value of shape function and its derivatives of 
                call reg_shape(ido,eita,pso,dpso)
            

						
                call reg_shape(idi,zeta,psi,dpsi)

			!...field element
		

		call findB1(phi,B1,xlo,xli,pso,dpso,psi,dpsi)
		call findE1(phi,nu,E1,xlo,xli,pso,dpso,psi,dpsi)
            
			!...calculate val_inner_integral(l,k,j,i)
                do i=1,neo
              	    do j=1,nei
                        do k=1,NDOF
                    	    do l=1,NDOF
                                val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + (B1(l,k)*dpso(i)*dpsi(j) + E1(l,k)*pso(i)*psi(j))*wii(li)

                            
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
	if ((adt_type.eq.2).or.(adt_type.eq.4)) then
    	ekbar21 = -2.0d0*ekbar21
    elseif ((adt_type.eq.1).or.(adt_type.eq.3)) then
    	ekbar22 = -2.0d0*ekbar22
    endif
	!...In this way, I_4.0d0 in Eqn (3.91) is evaluated
    
	!...regular Gauss integration: ekbar21 for adt_type = 1/3 and ekbar22 for adt_type = 2/4.0d0
	!...Now we are evaluating I_3 in Eqn (3.8.0d09)
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
				!...Eqn (3.8.0d08.0d0)(1)
            	beta = 0.5d0*(gamma - 1.0d0)
				!...Eqn (3.8.0d08.0d0)(2)
            elseif ((adt_type.eq.2).or.(adt_type.eq.4)) then
            	alpha = 0.5d0*(1.0d0-gamma)*deltaa
            	beta = 0.5d0*(gamma + 1.0d0)
            else
              	print*,'ekd_adt_reg_reg: wrong adt_type:',adt_type
                stop
            endif
            !...distinguish adjacent types between 2/4.0d0 and 1/3
            if ((adt_type.eq.4).or.(adt_type.eq.3)) then
              	eitabar = alpha + beta
                eita = -1.0d0*eitabar
            else
            	eita = alpha + beta
				!...again, Eqn (3.8.0d06)(1)
            endif
            zeta = alpha - beta
			!...Eqn (3.8.0d06)(2)
            !...obtain value of shape function and its derivatives of 
            call reg_shape(ido,eita,pso,dpso)

						
			call reg_shape(idi,zeta,psi,dpsi)


			call findB1(phi,B1,xlo,xli,pso,dpso,psi,dpsi)
			call findE1(phi,nu,E1,xlo,xli,pso,dpso,psi,dpsi)
            
			!...calculate val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                    do k=1,NDOF
                    	do l=1,NDOF
                            val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + (B1(l,k)*dpso(i)*dpsi(j) + E1(l,k)*pso(i)*psi(j))*wii(li)

                           
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
    if ((adt_type.eq.2).or.(adt_type.eq.4)) then
    	ekbar22 = 0.5d0*ekbar22
    elseif ((adt_type.eq.1).or.(adt_type.eq.3)) then
    	ekbar21 = 0.5d0*ekbar21
    endif
	!...In this way, I_3 in Eqn (3.8.0d09) is evaluated
    !--------------------------------------------------
    !...compute matrix ekbar1
	!...Now we are evaluating I_1 in Eqn (3.8.0d05)
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
          	print*,'ekd_adt_reg_reg: wrong adt_type'
            stop
        endif
        !...obtain value of shape function and its derivatives
        call reg_shape(ido,eita,pso,dpso)
        

				
		!...initialize
        val_inner_integral = 0.0d0
		!...loop over inner integral (Gauss integration)
		do li =1,nti
			!...coordinate of Gauss point
        	zeta = xii(li)
            !...obtain value of shape function and its derivatives
            call reg_shape(idi,zeta,psi,dpsi)
        


			call findB1(phi,B1,xlo,xli,pso,dpso,psi,dpsi)
			call findE1(phi,nu,E1,xlo,xli,pso,dpso,psi,dpsi)

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
        
                         !...compute r
                         do k=1,2
              	             temp(k) = field(k) - source(k)
                         enddo
                         r = dsqrt(temp(1)*temp(1) + temp(2)*temp(2))
            
			!...calculate val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                    do k=1,NDOF
                    	do l=1,NDOF
                            if (adt_type.eq.1) then
                        	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + (B1(l,k)*dpso(i)*dpsi(j) + E1(l,k)*pso(i)*psi(j))*dlog(2.0d0*r/(2.0d0-eita+zeta))*wii(li)

								
			    elseif (adt_type.eq.3) then
                            	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + (B1(l,k)*dpso(i)*dpsi(j) + E1(l,k)*pso(i)*psi(j))*dlog(2.0d0*r/(2.0d0-eitabar+zeta))*wii(li)

								
			    elseif (adt_type.eq.2) then
                        	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + (B1(l,k)*dpso(i)*dpsi(j) + E1(l,k)*pso(i)*psi(j))*dlog(2.0d0*r/(2.0d0-zeta+eita))*wii(li)

								
			    elseif (adt_type.eq.4) then
                            	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + (B1(l,k)*dpso(i)*dpsi(j) + E1(l,k)*pso(i)*psi(j))*dlog(2.0d0*r/(2.0d0-zeta+eitabar))*wii(li)

			    endif
                        enddo
                    enddo
                enddo
            enddo
        enddo  !of loop over inner integral
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
	!...In this way, I_1 in Eqn (3.8.0d05) is evaluated
    
	!--------------------------------------------------
    !...evaluate ek2bar (regular integral)
	!...We evaluate \hat{K}_{J\beta}^{I\alpha} part of C_l^k in Eqn (3.61)
    ek2bar = 0.0d0
    do lo=1,nto
      	!...coordinate of Gauss point
      	eita = xio(lo)  
        !...shape function & its derivatives at Gauss point
        call reg_shape(ido,eita,pso,dpso)
        

				
        !...loop over inner integral
        val_inner_integral = 0.0d0
        do li=1,nti
          	!...coordinate of Gauss point
            zeta = xii(li)
            !...shape function & its derivatives at Gauss point
        	call reg_shape(idi,zeta,psi,dpsi)
        	


			call findB2(phi,B2,xlo,xli,pso,dpso,psi,dpsi)
			call findE2(phi,nu,E2,xlo,xli,pso,dpso,psi,dpsi)
			
			!...compute val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                    do k=1,NDOF
                    	do l=1,NDOF
                            val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + (B2(l,k)*dpso(i)*dpsi(j) + E2(l,k)*pso(i)*psi(j))*wii(li)

							
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
	!...Here we evaluate the part of stiffness which is a function of direction of r
    
    !...element stiffness matrix
    ek = ekbar1 + ekbar21 + ekbar22 + ek2bar        
	End Subroutine ekd_adt_reg_reg
	!-----------------------------------------------
	Subroutine ekd_adt_reg_tip(region,ido,idi,ek,nto,ntoln,nti,xio,xioln,xii,wio,wioln,wii,xlo,xli,adt_type,phi,ht,nu)
	!...this subroutine is new comparing with ek_E2.f90
    IMPLICIT NONE
    INTEGER, INTENT(IN)			:: ido,idi,nto,ntoln,nti,region
    INTEGER, INTENT(IN)			:: adt_type
    REAL(KIND=DBL), INTENT(OUT)	:: ek(:,:)
    REAL(KIND=DBL), INTENT(IN)	:: xio(:),xioln(:),xii(:),wio(:),wioln(:),wii(:)
    REAL(KIND=DBL), INTENT(IN)	:: xlo(:,:),xli(:,:)
	REAL(KIND=DBL), INTENT(IN)	:: phi,ht,nu
	!...local variables
    INTEGER		:: neo,nei,lo,li
    REAL(KIND=DBL)	:: ekbar1a(3*NDOF,3*NDOF),ekbar1b1(3*NDOF,3*NDOF),ekbar1b2(3*NDOF,3*NDOF)
    REAL(KIND=DBL)	:: ekbar2(3*NDOF,3*NDOF),ekbar3(3*NDOF,3*NDOF),ek2bar(3*NDOF,3*NDOF)
    REAL(KIND=DBL)	:: alpha,beta,eita,zeta,gamma,deltaa,zetapr,zetaprpr,rho,eitabar
    REAL(KIND=DBL)	:: pso(3),dpso(3),dtipi(3),tipi(3),psi(3),dpsi(3)
    INTEGER			:: i,j,k,l
    REAL(KIND=DBL)	:: val_inner_integral(NDOF,NDOF,3,3)
	REAL(KIND=DBL)	:: B1(NDOF,NDOF),B2(NDOF,NDOF)
    REAL(KIND=DBL)	:: E1(NDOF,NDOF),E2(NDOF,NDOF)
	REAL(KIND=DBL)	:: field(2),source(2),r,temp(2)
    
	REAL(KIND=DBL),parameter	:: avalue = -0.4d0, bvalue = 0.4d0 !...for subinterval splitting
	
	neo = NODE(ido)
    nei = NODE(idi)
	!...computing the I_1 equivalent of ekd_adt_reg_reg
	!...compute matrix ekbar1a
	ekbar1a=0.0d0
	!...begins outer integration (Gauss int.)
	do lo = 1,nto
        !...coordinate of integration point
        if ((adt_type.eq.4).or.(adt_type.eq.3)) then
        	eitabar = xio(lo)
            eita = -1.0d0*eitabar
        elseif ((adt_type.eq.2).or.(adt_type.eq.1)) then
          	eita = xio(lo)
        else
          	print*,'error in ekd_adt_reg_tip (ekbar1a): wrong adjacent type'
            stop
        endif
        !...obtain value of shape function and its derivatives
    
	    call reg_shape(ido,eita,pso,dpso)


				
		!...initialize
        val_inner_integral = 0.0d0
		!...loop over inner integral (Gauss integration)
		do li =1,nti
			!...coordinate of Gauss point
        	zetapr = xii(li)
        	!...transformation back to zeta
            if ((adt_type.eq.2).or.(adt_type.eq.4)) then
              	zeta = 0.5d0*(1.0d0 - bvalue)*zetapr + 0.5d0*(1.0d0 + bvalue)
				!...This is used instead of Eqn (3.73)
            elseif ((adt_type.eq.1).or.(adt_type.eq.3)) then
            	zeta = 0.5d0*(avalue + 1.0d0)*zetapr + 0.5d0*(avalue - 1.0d0)
            endif
            !...obtain value of shape function and its derivatives
            call tip_shape(idi,zeta,tipi)
			call tip_dshape(idi,zeta,dtipi)
			call reg_shape(idi,zeta,psi,dpsi)


			
			call findB1(phi,B1,xlo,xli,pso,dpso,psi,dpsi)
			call findE1(phi,nu,E1,xlo,xli,pso,dpso,psi,dpsi)

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
        
                         !...compute r
                         do k=1,2
              	             temp(k) = field(k) - source(k)
                         enddo
                         r = dsqrt(temp(1)*temp(1) + temp(2)*temp(2))
            
			!...calculate val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                    do k=1,NDOF
                    	do l=1,NDOF
                            if (adt_type.eq.1) then
                        	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + (B1(l,k)*dpso(i)*dtipi(j) + E1(l,k)*pso(i)*tipi(j))*dlog(2.0d0*r/(2.0d0-eita+zetapr))*wii(li)

								
			    elseif (adt_type.eq.3) then
                            	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + (B1(l,k)*dpso(i)*dtipi(j) + E1(l,k)*pso(i)*tipi(j))*dlog(2.0d0*r/(2.0d0-eitabar+zetapr))*wii(li)

								 
			    elseif (adt_type.eq.2) then
                        	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + (B1(l,k)*dpso(i)*dtipi(j) + E1(l,k)*pso(i)*tipi(j))*dlog(2.0d0*r/(2.0d0-zetapr+eita))*wii(li)

								
			    elseif (adt_type.eq.4) then
                            	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + (B1(l,k)*dpso(i)*dtipi(j) + E1(l,k)*pso(i)*tipi(j))*dlog(2.0d0*r/(2.0d0-zetapr+eitabar))*wii(li)

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
                    	ekbar1a(NDOF*(i-1)+k,NDOF*(j-1)+l) = ekbar1a(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		val_inner_integral(l,k,j,i)*wio(lo)
                    enddo
                enddo
            enddo
        enddo
    enddo !of loop over outer integral
    !...finalize ekbar1a
    if ((adt_type.eq.2).or.(adt_type.eq.4)) then
      	ekbar1a = ekbar1a*0.5d0*(1.0d0-bvalue)
    elseif ((adt_type.eq.1).or.(adt_type.eq.3)) then
    	ekbar1a = ekbar1a*0.5d0*(avalue+1.0d0)
    endif
	!...I_1 equivalent of ekd_adt_reg_reg
	!--------------------------------------------------
    !...I_3 equivalent of ekd_adt_reg_reg
	!...compute matrix ekbar1b1 (regular)
	ekbar1b1 = 0.0d0
	!...begins outer integral loop
	do lo = 1,nto
    	!...coordinate of integration point
        gamma = xio(lo)
		!...initialize
        val_inner_integral = 0.0d0
		!...loop over inner integral (Gauss integration)
		do li =1,nti
			!...coordinate of Gauss point
        	deltaa = xii(li)
        	!...transformation back to alpha and beta
            if ((adt_type.eq.1).or.(adt_type.eq.3)) then
            	alpha = 0.5d0*(1.0d0+gamma)*deltaa
            	beta = 0.5d0*(gamma-1.0d0)
            elseif ((adt_type.eq.2).or.(adt_type.eq.4)) then
            	alpha = 0.5d0*(1.0d0-gamma)*deltaa
            	beta = 0.5d0*(gamma+1.0d0)
            else
              	print*,'ekd_adt_reg_tip (ekbar1b1): wrong adt_type'
                stop
            endif
            !...transform to original zeta and eita
            if ((adt_type.eq.3).or.(adt_type.eq.4)) then
              	eitabar = alpha + beta
                eita = -1.0d0*eitabar
            else
            	eita = alpha + beta
            endif
            zetapr = alpha - beta
            !...transform to original zeta
            if ((adt_type.eq.1).or.(adt_type.eq.3)) then
                zeta = 0.5d0*(avalue+1.0d0)*zetapr + 0.5d0*(avalue-1.0d0)
            elseif ((adt_type.eq.2).or.(adt_type.eq.4)) then
                zeta = 0.5d0*(1.0d0-bvalue)*zetapr + 0.5d0*(1.0d0+bvalue)
            endif
            !...obtain value of shape function and its derivatives of 
            call reg_shape(ido,eita,pso,dpso)

						
            call tip_shape(idi,zeta,tipi)
			!...field element
            call tip_dshape(idi,zeta,dtipi)
			!...field element
			call reg_shape(idi,zeta,psi,dpsi)



			call findB1(phi,B1,xlo,xli,pso,dpso,psi,dpsi)
			call findE1(phi,nu,E1,xlo,xli,pso,dpso,psi,dpsi)

			!...calculate val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                    do k=1,NDOF
                    	do l=1,NDOF
                            val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + (B1(l,k)*dpso(i)*dtipi(j) + E1(l,k)*pso(i)*tipi(j))*wii(li)

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
                        	ekbar1b1(NDOF*(i-1)+k,NDOF*(j-1)+l) = ekbar1b1(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		(1.0d0+gamma)*dlog(0.5d0*(3.0d0-gamma))*val_inner_integral(l,k,j,i)*wio(lo)
                        elseif ((adt_type.eq.2).or.(adt_type.eq.4)) then
                    		ekbar1b1(NDOF*(i-1)+k,NDOF*(j-1)+l) = ekbar1b1(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		(1.0d0-gamma)*dlog(0.5d0*(3.0d0+gamma))*val_inner_integral(l,k,j,i)*wio(lo)
                        endif
                    enddo
                enddo
            enddo
        enddo
    enddo !of loop over outer integral
    !...finalize ekbar1b1
    if ((adt_type.eq.1).or.(adt_type.eq.3)) then
      	ekbar1b1 = 0.25d0*(avalue+1.0d0)*ekbar1b1
    elseif ((adt_type.eq.2).or.(adt_type.eq.4)) then
    	ekbar1b1 = 0.25d0*(1.0d0-bvalue)*ekbar1b1
    endif
    !...I_3 equivalent of ekd_adt_reg_reg
	!--------------------------------------------------
    !...I_4.0d0 equivalent of ekd_adt_reg_reg
	!...compute matrix ekbar1b2 (singular)
	ekbar1b2 = 0.0d0
	!...begins outer integral loop (logarith integration)
	do lo = 1,ntoln
    	!...coordinate of integration point
        gamma = xioln(lo)
		!...initialize
        val_inner_integral = 0.0d0
		!...loop over inner integral (Gauss integration)
		do li =1,nti
			!...coordinate of Gauss point
        	deltaa = xii(li)
        	!...transformation back to alpha and beta
            if ((adt_type.eq.1).or.(adt_type.eq.3)) then
            	alpha = deltaa*gamma
            	beta = 1.0d0-gamma
            elseif ((adt_type.eq.2).or.(adt_type.eq.4)) then
            	alpha = deltaa*gamma
            	beta = gamma-1.0d0
            else
              	print*,'ekd_adt_reg_tip (ekbar1b2): wrong adt_type'
                stop
            endif
            !...transform to original zeta and eita
            if ((adt_type.eq.3).or.(adt_type.eq.4)) then
              	eitabar = alpha + beta
                eita = -1.0d0*eitabar
            else
            	eita = alpha + beta
            endif
            zetapr = alpha - beta
            !...transform back to original zeta
            if ((adt_type.eq.1).or.(adt_type.eq.3)) then
                zeta = 0.5d0*(avalue+1.0d0)*zetapr + 0.5d0*(avalue-1.0d0)
            elseif ((adt_type.eq.2).or.(adt_type.eq.4)) then
                zeta = 0.5d0*(1.0d0-bvalue)*zetapr + 0.5d0*(1.0d0+bvalue)
            endif
            !...obtain value of shape function and its derivatives of 
            call reg_shape(ido,eita,pso,dpso)


						
            call tip_shape(idi,zeta,tipi)
			!...field element
            call tip_dshape(idi,zeta,dtipi)
			!...field element
			call reg_shape(idi,zeta,psi,dpsi)


			call findB1(phi,B1,xlo,xli,pso,dpso,psi,dpsi)
			call findE1(phi,nu,E1,xlo,xli,pso,dpso,psi,dpsi)
            
			!...calculate val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                    do k=1,NDOF
                    	do l=1,NDOF
                            val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + (B1(l,k)*dpso(i)*dtipi(j) + E1(l,k)*pso(i)*tipi(j))*wii(li)

							
			enddo
                    enddo
                enddo
            enddo
        enddo !of loop over inner integral
        do i=1,neo
          	do j=1,nei
            	do k=1,NDOF
                	do l=1,NDOF
                       	ekbar1b2(NDOF*(i-1)+k,NDOF*(j-1)+l) = ekbar1b2(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		gamma*val_inner_integral(l,k,j,i)*wioln(lo)
                    enddo
                enddo
            enddo
        enddo
    enddo !of loop over outer integral
    !...finalize ekbar1b2
    if ((adt_type.eq.1).or.(adt_type.eq.3)) then
      	ekbar1b2 = -1.0d0*(avalue+1.0d0)*ekbar1b2
    elseif ((adt_type.eq.2).or.(adt_type.eq.4)) then
    	ekbar1b2 = -1.0d0*(1.0d0-bvalue)*ekbar1b2
    endif
    !...I_4.0d0 equivalent of ekd_adt_reg_reg
	!--------------------------------------------------
    !...evaluate ekbar2 (regular integral)
    ekbar2 = 0.0d0
    do lo=1,nto
      	!...coordinate of Gauss point
        if ((adt_type.eq.3).or.(adt_type.eq.4)) then
          	eitabar = xio(lo)
            eita = -1.0d0*eitabar
        elseif ((adt_type.eq.1).or.(adt_type.eq.2)) then
      		eita = xio(lo)
        else
          	print*,'ekd_adt_reg_tip (ekbar2): wrong adjacent type'
            stop
        endif
        !...shape function & its derivatives at Gauss point
        call reg_shape(ido,eita,pso,dpso)



				
        !...loop over inner integral
        val_inner_integral = 0.0d0
        do li=1,nti
          	!...coordinate of Gauss point
            zetapr = xii(li)
            !...transform to original zeta
            zeta = 0.5d0*(bvalue-avalue)*zetapr + 0.5d0*(bvalue+avalue)
            !...shape function & its derivatives at Gauss point
        	call tip_shape(idi,zeta,tipi)
            call tip_dshape(idi,zeta,dtipi)
			call reg_shape(idi,zeta,psi,dpsi)


			call findB1(phi,B1,xlo,xli,pso,dpso,psi,dpsi)
			call findE1(phi,nu,E1,xlo,xli,pso,dpso,psi,dpsi)

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
        
                         !...compute r
                         do k=1,2
              	             temp(k) = field(k) - source(k)
                         enddo
                         r = dsqrt(temp(1)*temp(1) + temp(2)*temp(2))
			
			!...calculate val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                    do k=1,NDOF
                    	do l=1,NDOF
                            val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + (B1(l,k)*dpso(i)*dtipi(j) + E1(l,k)*pso(i)*tipi(j))*dlog(r)*wii(li)

			enddo
                    enddo
                enddo
            enddo
        enddo !of loop over inner integral
        do i=1,neo
          	do j=1,nei
            	do k=1,NDOF
                	do l=1,NDOF
                       	ekbar2(NDOF*(i-1)+k,NDOF*(j-1)+l) = ekbar2(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		val_inner_integral(l,k,j,i)*wio(lo)
                    enddo
                enddo
            enddo
        enddo
    enddo !of loop over outer integral
    !...finalize ekbar2
   	ekbar2 = 0.5d0*(bvalue-avalue)*ekbar2
    !--------------------------------------------------
    !...evaluate ekbar3 (singular at tip node)
    ekbar3 = 0.0d0
    do lo=1,nto
      	!...coordinate of Gauss point
        if ((adt_type.eq.3).or.(adt_type.eq.4)) then
          	eitabar = xio(lo)
            eita = -1.0d0*eitabar
        elseif ((adt_type.eq.1).or.(adt_type.eq.2)) then
      		eita = xio(lo)
        else
          	print*,'ekd_adt_reg_tip (ekbar3): wrong adjacent type'
            stop
        endif
        !...shape function & its derivatives at Gauss point
        call reg_shape(ido,eita,pso,dpso)


				
        !...loop over inner integral
        val_inner_integral = 0.0d0
        do li=1,nti
          	!...coordinate of Gauss point
            zetaprpr = xii(li)
            !...transform back to zetapr
            if ((adt_type.eq.1).or.(adt_type.eq.3)) then
              	zetapr = 0.5d0*dsqrt(2.0d0*(1.0d0-bvalue))*zetaprpr + 1.0d0-0.5d0*dsqrt(2.0d0*(1.0d0-bvalue))
            elseif ((adt_type.eq.2).or.(adt_type.eq.4)) then
            	zetapr = 0.5d0*dsqrt(2.0d0*(avalue+1.0d0))*zetaprpr + 0.5d0*dsqrt(2.0d0*(avalue+1.0d0))-1.0d0
            endif
            !...transform back to zeta to get rid of singularity of tip shape functions
            if ((adt_type.eq.1).or.(adt_type.eq.3)) then
              	zeta = 1.0d0 - 0.5d0*(1.0d0-zetapr)*(1.0d0-zetapr)
            elseif ((adt_type.eq.2).or.(adt_type.eq.4)) then
            	zeta = 0.5d0*(1.0d0+zetapr)*(1.0d0+zetapr) - 1.0d0
            endif
            !...shape function & its derivatives at Gauss point
        	call tip_shape(idi,zeta,tipi)
            call tip_dshape(idi,zeta,dtipi)
			call reg_shape(idi,zeta,psi,dpsi)


			
			call findB1(phi,B1,xlo,xli,pso,dpso,psi,dpsi)
			call findE1(phi,nu,E1,xlo,xli,pso,dpso,psi,dpsi)

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
        
                         !...compute r
                         do k=1,2
              	             temp(k) = field(k) - source(k)
                         enddo
                         r = dsqrt(temp(1)*temp(1) + temp(2)*temp(2))
            
			!...calculate val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                    do k=1,NDOF
                    	do l=1,NDOF
                            if ((adt_type.eq.1).or.(adt_type.eq.3)) then
                       		val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + (B1(l,k)*dpso(i)*dtipi(j) + E1(l,k)*pso(i)*tipi(j))*dlog(r)*(1.0d0-zetaprpr)*wii(li)

							
			    elseif ((adt_type.eq.2).or.(adt_type.eq.4)) then
                       		val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + (B1(l,k)*dpso(i)*dtipi(j) + E1(l,k)*pso(i)*tipi(j))*dlog(r)*(1.0d0+zetaprpr)*wii(li)

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
                       	ekbar3(NDOF*(i-1)+k,NDOF*(j-1)+l) = ekbar3(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		val_inner_integral(l,k,j,i)*wio(lo)
                    enddo
                enddo
            enddo
        enddo
    enddo !of loop over outer integral
    !...finalize ekbar3
    if ((adt_type.eq.1).or.(adt_type.eq.3)) then
   		ekbar3 = 0.5d0*(1.0d0-bvalue)*ekbar3
    elseif ((adt_type.eq.2).or.(adt_type.eq.4)) then
   		ekbar3 = 0.5d0*(avalue+1.0d0)*ekbar3
    endif
    !--------------------------------------------------
    !...evaluate ek2bar (regular integral, associated with the sencond part of the kernel decomposition)
    ek2bar = 0.0d0
    do lo=1,nto
      	!...coordinate of Gauss point
      	eita = xio(lo)  
        !...shape function & its derivatives at Gauss point
        call reg_shape(ido,eita,pso,dpso)


        !...loop over inner integral
        val_inner_integral = 0.0d0
        do li=1,nti
          	!...coordinate of Gauss point
            zetapr = xii(li)
            !...transform to get rid of singularity of tip shape function
            if (idi.eq.CTIP1) then
              	zeta = 0.5d0*(zetapr+1.0d0)*(zetapr+1.0d0)-1.0d0
            elseif (idi.eq.CTIP2) then
            	zeta = 1.0d0-0.5d0*(1.0d0-zetapr)*(1.0d0-zetapr)
            else
              	print*,'ek_d3.f90,adt_reg_tip (ek2bar): wrong id of tip element'
                stop
            endif
            !...shape function & its derivatives at Gauss point
        	call tip_shape(idi,zeta,tipi)
            call tip_dshape(idi,zeta,dtipi)
			call reg_shape(idi,zeta,psi,dpsi)


			call findB2(phi,B2,xlo,xli,pso,dpso,psi,dpsi)
			call findE2(phi,nu,E2,xlo,xli,pso,dpso,psi,dpsi)

			!...compute val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                    do k=1,NDOF
                    	do l=1,NDOF
                            if (idi.eq.CTIP1) then
                            	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + (B2(l,k)*dpso(i)*dtipi(j) + E2(l,k)*pso(i)*tipi(j))*(1.0d0+zetapr)*wii(li)

			    elseif (idi.eq.CTIP2) then
                            	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + (B2(l,k)*dpso(i)*dtipi(j) + E2(l,k)*pso(i)*tipi(j))*(1.0d0-zetapr)*wii(li)

								
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
                        ek2bar(NDOF*(i-1)+k,NDOF*(j-1)+l) = ek2bar(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        val_inner_integral(l,k,j,i)*wio(lo)
                    enddo
                enddo
            enddo
        enddo
    enddo !over loop of outer integral
    !...combine all together to have element stiffness matrix
    ek = ekbar1a + ekbar1b1 + ekbar1b2 + ekbar2 + ekbar3 + ek2bar
	End Subroutine ekd_adt_reg_tip
    
!----------------------------------------------------------------------------------------------------
	Subroutine ekd_adt_tip_reg(region,ido,idi,ek,nto,ntoln,nti,xio,xioln,xii,wio,wioln,wii,xlo,xli,adt_type,phi,ht,nu)
	!...this subroutine is new comparing with ek_E2.f90
    IMPLICIT NONE
    INTEGER, INTENT(IN)			:: ido,idi,nto,ntoln,nti,region
    INTEGER, INTENT(IN)			:: adt_type
    REAL(KIND=DBL), INTENT(OUT)	:: ek(:,:)
    REAL(KIND=DBL), INTENT(IN)	:: xio(:),xioln(:),xii(:),wio(:),wioln(:),wii(:)
    REAL(KIND=DBL), INTENT(IN)	:: xlo(:,:),xli(:,:)
	REAL(KIND=DBL), INTENT(IN)	:: phi,ht,nu
    !...local variables
    INTEGER			:: neo,nei,lo,li
    REAL(KIND=DBL)	:: ekbar1a(3*NDOF,3*NDOF),ekbar1b1(3*NDOF,3*NDOF),ekbar1b2(3*NDOF,3*NDOF)
    REAL(KIND=DBL)	:: ekbar2(3*NDOF,3*NDOF),ekbar3(3*NDOF,3*NDOF),ek2bar(3*NDOF,3*NDOF)
    REAL(KIND=DBL)	:: alpha,beta,eita,eitapr,eitaprpr,zeta,gamma,deltaa,zetabar,rho
    REAL(KIND=DBL)	:: tipo(3),dtipo(3),psi(3),dpsi(3),pso(3),dpso(3)
    INTEGER			:: i,j,k,l
    REAL(KIND=DBL)	:: val_inner_integral(NDOF,NDOF,3,3)
    REAL(KIND=DBL)	:: B1(NDOF,NDOF),B2(NDOF,NDOF)
	REAL(KIND=DBL)	:: E1(NDOF,NDOF),E2(NDOF,NDOF)
	REAL(KIND=DBL)	:: field(2),source(2),r,temp(2)
   
	REAL(KIND=DBL),parameter	:: avalue = -0.4d0,bvalue = 0.4d0 !...for subinterval splitting
	
    !...adj_type = 1 or 4.0d0 - singular at eita = 1
	!...adj_type = 2 or 3 - singular at zeta = 1
	
	!...get number of nodes of outer/inner elements
    neo = NODE(ido)
	!...source element
    nei = NODE(idi)
	!...field element
	!...compute matrix ekbar1a
	ekbar1a=0.0d0
	!...begins outer integration (Gauss int.),source element
	do lo = 1,nto
        !...coordinate of integration point
        eitapr = xio(lo)
        !...transform back to original eita
        if ((adt_type.eq.1).or.(adt_type.eq.4)) then
          	eita = 0.5d0*(1.0d0-bvalue)*eitapr + 0.5d0*(1.0d0+bvalue)
        elseif ((adt_type.eq.2).or.(adt_type.eq.3)) then
        	eita = 0.5d0*(avalue+1.0d0)*eitapr + 0.5d0*(avalue-1.0d0)
        else
          	print*,'ekd_adt_tip_reg (ekbar1a): wrong adjacent type'
            stop
        endif
        !...obtain value of shape function and its derivatives
        call tip_shape(ido,eita,tipo)
        call tip_dshape(ido,eita,dtipo)
		call reg_shape(ido,eita,pso,dpso)

				
		!...initialize
        val_inner_integral = 0.0d0
		!...loop over inner integral (Gauss integration)
		do li =1,nti
			!...looping over the Gauss points for inner integration
            if ((adt_type.eq.3).or.(adt_type.eq.4)) then
        		zetabar = xii(li)
                zeta = -1.0d0*zetabar
            elseif ((adt_type.eq.1).or.(adt_type.eq.2)) then
            	zeta = xii(li)
            endif
            !...obtain value of shape function and its derivatives
            call reg_shape(idi,zeta,psi,dpsi)

			
		    call findB1(phi,B1,xlo,xli,pso,dpso,psi,dpsi)
			call findE1(phi,nu,E1,xlo,xli,pso,dpso,psi,dpsi)


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
        
                         !...compute r
                         do k=1,2
              	             temp(k) = field(k) - source(k)
                         enddo
                         r = dsqrt(temp(1)*temp(1) + temp(2)*temp(2))

			!...calculate val_inner_integral(l,k,j,i)
            do i=1,neo
			!...number of nodes in source element
              	do j=1,nei
				!...number of nodes in field element
                	do k=1,NDOF
					!...DX_k(y)
                    	do l=1,NDOF
						!...DY_l(\xi)
                            if (adt_type.eq.1) then
                        	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + (B1(l,k)*dtipo(i)*dpsi(j) + E1(l,k)*tipo(i)*psi(j))*dlog(2.0d0*r/(2.0d0+zeta-eitapr))*wii(li)

								
			    elseif (adt_type.eq.4) then
                            	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + (B1(l,k)*dtipo(i)*dpsi(j) + E1(l,k)*tipo(i)*psi(j))*dlog(2.0d0*r/(2.0d0+zetabar-eitapr))*wii(li)

			    elseif (adt_type.eq.2) then
                        	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + (B1(l,k)*dtipo(i)*dpsi(j) + E1(l,k)*tipo(i)*psi(j))*dlog(2.0d0*r/(2.0d0-zeta+eitapr))*wii(li)

			    elseif (adt_type.eq.3) then
                            	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + (B1(l,k)*dtipo(i)*dpsi(j) + E1(l,k)*tipo(i)*psi(j))*dlog(2.0d0*r/(2.0d0-zetabar+eitapr))*wii(li)

								
			    endif
							!...evaluating term 1 in (3.68.0d0)
                        enddo
                    enddo
                enddo
            enddo
        enddo !of loop over inner integral
        do i=1,neo
          	do j=1,nei
            	do k=1,NDOF
                	do l=1,NDOF
                    	ekbar1a(NDOF*(i-1)+k,NDOF*(j-1)+l) = ekbar1a(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		val_inner_integral(l,k,j,i)*wio(lo)
                    enddo
                enddo
            enddo
        enddo
    enddo !of loop over outer integral
    !...finalize ekbar1a
    if ((adt_type.eq.1).or.(adt_type.eq.4)) then
      	ekbar1a = ekbar1a*0.5d0*(1.0d0-bvalue)
    elseif ((adt_type.eq.2).or.(adt_type.eq.3)) then
    	ekbar1a = ekbar1a*0.5d0*(avalue+1.0d0)
    endif

	!--------------------------------------------------
    !...compute matrix ekbar1b1 (regular)
	ekbar1b1 = 0.0d0
	!...begins outer integral loop
	do lo = 1,nto
    	!...coordinate of integration point
        gamma = xio(lo)
		!...initialize
        val_inner_integral = 0.0d0
		!...loop over inner integral (Gauss integration)
		do li =1,nti
			!...coordinate of Gauss point
        	deltaa = xii(li)
        	!...transformation back to alpha and beta
            if ((adt_type.eq.1).or.(adt_type.eq.4)) then
            	alpha = 0.5d0*(1.0d0+gamma)*deltaa
            	beta = 0.5d0*(gamma-1.0d0)
            elseif ((adt_type.eq.2).or.(adt_type.eq.3)) then
            	alpha = 0.5d0*(1.0d0-gamma)*deltaa
            	beta = 0.5d0*(gamma+1.0d0)
            else
              	print*,'ekd_adt_tip_reg (ekbar1b1): wrong adt_type'
                stop
            endif
            !...transform to original zeta and eita
            if ((adt_type.eq.3).or.(adt_type.eq.4)) then
              	zetabar = alpha - beta
                zeta = -1.0d0*zetabar
            else
            	zeta = alpha - beta
            endif
            eitapr = alpha + beta
            !...transform back to original eita
            if ((adt_type.eq.1).or.(adt_type.eq.4)) then
                eita = 0.5d0*(1.0d0-bvalue)*eitapr + 0.5d0*(1.0d0+bvalue)
            elseif ((adt_type.eq.2).or.(adt_type.eq.3)) then
                eita = 0.5d0*(avalue+1.0d0)*eitapr + 0.5d0*(avalue-1.0d0)
            endif
            
            !...obtain value of shape function and its derivatives
            call tip_shape(ido,eita,tipo)
			call tip_dshape(ido,eita,dtipo)
			call reg_shape(ido,eita,pso,dpso)

						
		
			call reg_shape(idi,zeta,psi,dpsi)


			call findB1(phi,B1,xlo,xli,pso,dpso,psi,dpsi)
			call findE1(phi,nu,E1,xlo,xli,pso,dpso,psi,dpsi)
			!...calculate val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                    do k=1,NDOF
                    	do l=1,NDOF
                            val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + (B1(l,k)*dtipo(i)*dpsi(j) + E1(l,k)*tipo(i)*psi(j))*wii(li)

							
			enddo
                    enddo
                enddo
            enddo
        enddo !of loop over inner integral
        do i=1,neo
          	do j=1,nei
            	do k=1,NDOF
                	do l=1,NDOF
                    	if ((adt_type.eq.1).or.(adt_type.eq.4)) then
                        	ekbar1b1(NDOF*(i-1)+k,NDOF*(j-1)+l) = ekbar1b1(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		(1.0d0+gamma)*dlog(0.5d0*(3.0d0-gamma))*val_inner_integral(l,k,j,i)*wio(lo)
                        elseif ((adt_type.eq.2).or.(adt_type.eq.3)) then
                    		ekbar1b1(NDOF*(i-1)+k,NDOF*(j-1)+l) = ekbar1b1(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		(1.0d0-gamma)*dlog(0.5d0*(3.0d0+gamma))*val_inner_integral(l,k,j,i)*wio(lo)
                        endif
                    enddo
                enddo
            enddo
        enddo
    enddo !of loop over outer integral
    !...finalize ekbar1b1
    if ((adt_type.eq.1).or.(adt_type.eq.4)) then
      	ekbar1b1 = 0.25d0*(1.0d0-bvalue)*ekbar1b1
    elseif ((adt_type.eq.2).or.(adt_type.eq.3)) then
    	ekbar1b1 = 0.25d0*(avalue+1.0d0)*ekbar1b1
    endif
    
    !--------------------------------------------------
    !...compute matrix ekbar1b2 (singular)
	ekbar1b2 = 0.0d0
	!...begins outer integral loop (logarith integration)
	do lo = 1,ntoln
    	!...coordinate of integration point
        gamma = xioln(lo)
		!...initialize
        val_inner_integral = 0.0d0
		!...loop over inner integral (Gauss integration)
		do li =1,nti
			!...coordinate of Gauss point
        	deltaa = xii(li)
        	!...transformation back to alpha and beta
            if ((adt_type.eq.1).or.(adt_type.eq.4)) then
            	alpha = deltaa*gamma
            	beta = 1.0d0-gamma
            elseif ((adt_type.eq.2).or.(adt_type.eq.3)) then
            	alpha = deltaa*gamma
            	beta = gamma-1.0d0
            else
              	print*,'ekd_adt_tip_reg (ekbar1b2): wrong adt_type'
                stop
            endif
            !...transform to original zeta and eita
            if ((adt_type.eq.3).or.(adt_type.eq.4)) then
              	zetabar = alpha - beta
                zeta = -1.0d0*zetabar
            else
            	zeta = alpha - beta
            endif
            eitapr = alpha + beta
            !...transform to original eita
            if ((adt_type.eq.1).or.(adt_type.eq.4)) then
                eita = 0.5d0*(1.0d0-bvalue)*eitapr + 0.5d0*(1.0d0+bvalue)
            elseif ((adt_type.eq.2).or.(adt_type.eq.3)) then
                eita = 0.5d0*(avalue+1.0d0)*eitapr + 0.5d0*(avalue-1.0d0)
            endif
            !...obtain value of shape function and its derivatives of 
            call tip_shape(ido,eita,tipo)
			call tip_dshape(ido,eita,dtipo)
			call reg_shape(ido,eita,pso,dpso)

						
        
			call reg_shape(idi,zeta,psi,dpsi)



			call findB1(phi,B1,xlo,xli,pso,dpso,psi,dpsi)
			call findE1(phi,nu,E1,xlo,xli,pso,dpso,psi,dpsi)
            !...calculate val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                    do k=1,NDOF
                    	do l=1,NDOF
                            val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + (B1(l,k)*dtipo(i)*dpsi(j) + E1(l,k)*tipo(i)*psi(j))*wii(li)

							
			enddo
                    enddo
                enddo
            enddo
        enddo !of loop over inner integral
        do i=1,neo
          	do j=1,nei
            	do k=1,NDOF
                	do l=1,NDOF
                       	ekbar1b2(NDOF*(i-1)+k,NDOF*(j-1)+l) = ekbar1b2(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		gamma*val_inner_integral(l,k,j,i)*wioln(lo)
                    enddo
                enddo
            enddo
        enddo
    enddo !of loop over outer integral
    !...finalize ekbar1b2
    if ((adt_type.eq.1).or.(adt_type.eq.4)) then
      	ekbar1b2 = -1.0d0*(1.0d0-bvalue)*ekbar1b2
    elseif ((adt_type.eq.2).or.(adt_type.eq.3)) then
    	ekbar1b2 = -1.0d0*(avalue+1.0d0)*ekbar1b2
    endif
    !--------------------------------------------------
    !...evaluate ekbar2 (regular integral)
    ekbar2 = 0.0d0
    do lo=1,nto
      	!...coordinate of Gauss point
        eitapr = xio(lo)
        !...transform back to original eita
        eita = 0.5d0*(bvalue-avalue)*eitapr + 0.5d0*(bvalue+avalue)
        
        !...shape function & its derivatives at Gauss point
        call tip_shape(ido,eita,tipo)
        call tip_dshape(ido,eita,dtipo)
        call reg_shape(ido,eita,pso,dpso)


				
        !...loop over inner integral
        val_inner_integral = 0.0d0
        do li=1,nti
          	!...coordinate of Gauss point
            if ((adt_type.eq.3).or.(adt_type.eq.4)) then
            	zetabar = xii(li)
                zeta = -1.0d0*zetabar
            elseif ((adt_type.eq.1).or.(adt_type.eq.2)) then
            	zeta = xii(li)
            else
              	print*,'ek_d3.f90,adt_tip_reg(ekbar2): wrong adjacent type'
                stop
            endif
            
            !...shape function & its derivatives at Gauss point
        	call reg_shape(idi,zeta,psi,dpsi)


			
			call findB1(phi,B1,xlo,xli,pso,dpso,psi,dpsi)
			call findE1(phi,nu,E1,xlo,xli,pso,dpso,psi,dpsi)

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
        
                         !...compute r
                         do k=1,2
              	             temp(k) = field(k) - source(k)
                         enddo
                         r = dsqrt(temp(1)*temp(1) + temp(2)*temp(2))

            !...calculate val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                    do k=1,NDOF
                    	do l=1,NDOF
                            val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + (B1(l,k)*dtipo(i)*dpsi(j) + E1(l,k)*tipo(i)*psi(j))*dlog(r)*wii(li)

							
			enddo
                    enddo
                enddo
            enddo
        enddo !of loop over inner integral
        do i=1,neo
          	do j=1,nei
            	do k=1,NDOF
                	do l=1,NDOF
                       	ekbar2(NDOF*(i-1)+k,NDOF*(j-1)+l) = ekbar2(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		val_inner_integral(l,k,j,i)*wio(lo)
                    enddo
                enddo
            enddo
        enddo
    enddo !of loop over outer integral
    !...finalize ekbar2
   	ekbar2 = 0.5d0*(bvalue-avalue)*ekbar2

    !--------------------------------------------------
    !...evaluate ekbar3 (singular at tip node)
    ekbar3 = 0.0d0
    do lo=1,nto
      	!...coordinate of Gauss point
        eitaprpr = xio(lo)
        !...transform back to eitapr
        if ((adt_type.eq.1).or.(adt_type.eq.4)) then
          	eitapr = dsqrt(0.5d0*(avalue+1.0d0))*eitaprpr + dsqrt(0.5d0*(avalue+1.0d0)) - 1.0d0
        elseif ((adt_type.eq.2).or.(adt_type.eq.3)) then
      		eitapr = dsqrt(0.5d0*(1.0d0-bvalue))*eitaprpr + 1.0d0 - dsqrt(0.5d0*(1.0d0-bvalue))
        else
          	print*,'ekd_adt_tip_reg (ekbar3): wrong adjacent type'
            stop
        endif
        !...transform back to eita (to get rid of singularity of derivative of crack-tip shape function
        if ((adt_type.eq.1).or.(adt_type.eq.4)) then
          	eita = 0.5d0*(1.0d0+eitapr)*(1.0d0+eitapr) - 1.0d0
        elseif ((adt_type.eq.2).or.(adt_type.eq.3)) then
        	eita = 1.0d0 - 0.5d0*(1.0d0-eitapr)*(1.0d0-eitapr)
        endif
        !...shape function & its derivatives at Gauss point
        call tip_shape(ido,eita,tipo)
        call tip_dshape(ido,eita,dtipo)
        call reg_shape(ido,eita,pso,dpso)

        !...loop over inner integral
        val_inner_integral = 0.0d0
        do li=1,nti
          	!...coordinate of Gauss point
            if ((adt_type.eq.3).or.(adt_type.eq.4)) then
            	zetabar = xii(li)
                zeta = -1.0d0*zetabar
            elseif ((adt_type.eq.1).or.(adt_type.eq.2)) then
            	zeta = xii(li)
            endif
            !...shape function & its derivatives at Gauss point
        	call reg_shape(idi,zeta,psi,dpsi)



			call findB1(phi,B1,xlo,xli,pso,dpso,psi,dpsi)
			call findE1(phi,nu,E1,xlo,xli,pso,dpso,psi,dpsi)

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
        
                         !...compute r
                         do k=1,2
              	             temp(k) = field(k) - source(k)
                         enddo
                         r = dsqrt(temp(1)*temp(1) + temp(2)*temp(2))

            !...calculate val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                    do k=1,NDOF
                    	do l=1,NDOF
                            if ((adt_type.eq.1).or.(adt_type.eq.4)) then
                       		val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + (B1(l,k)*dtipo(i)*dpsi(j) + E1(l,k)*tipo(i)*psi(j))*dlog(r)*wii(li)

								
			    elseif ((adt_type.eq.2).or.(adt_type.eq.3)) then
                       		val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + (B1(l,k)*dtipo(i)*dpsi(j) + E1(l,k)*tipo(i)*psi(j))*dlog(r)*wii(li)

								
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
                    	if ((adt_type.eq.1).or.(adt_type.eq.4)) then
                       		ekbar3(NDOF*(i-1)+k,NDOF*(j-1)+l) = ekbar3(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		val_inner_integral(l,k,j,i)*(eitaprpr+1.0d0)*wio(lo)
                        elseif ((adt_type.eq.2).or.(adt_type.eq.3)) then
                       		ekbar3(NDOF*(i-1)+k,NDOF*(j-1)+l) = ekbar3(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		val_inner_integral(l,k,j,i)*(1.0d0-eitaprpr)*wio(lo)
                        endif
                    enddo
                enddo
            enddo
        enddo
    enddo !of loop over outer integral
    !...finalize ekbar3
    if ((adt_type.eq.1).or.(adt_type.eq.4)) then
   		ekbar3 = 0.5d0*(avalue+1.0d0)*ekbar3
    elseif ((adt_type.eq.2).or.(adt_type.eq.3)) then
   		ekbar3 = 0.5d0*(1.0d0-bvalue)*ekbar3
    endif
    !--------------------------------------------------
    !...evaluate ek2bar (regular integral, associated with the sencond part of the kernel decomposition)
    ek2bar = 0.0d0
    do lo=1,nto
      	!...coordinate of Gauss point
        eitapr = xio(lo)
        !...transform to get rid of singularity at tip shape functions
        if (ido.eq.CTIP1) then
        	eita = 0.5d0*(eitapr+1.0d0)*(eitapr+1.0d0)-1.0d0
        elseif (ido.eq.CTIP2) then
           	eita = 1.0d0-0.5d0*(1.0d0-eitapr)*(1.0d0-eitapr)
        else
           	print*,'ek_d3.f90,adt_tip_reg (ek2bar): id of outer element is out of range'
            stop
        endif
        !...shape function & its derivatives at Gauss point
        call tip_shape(ido,eita,tipo)
        call tip_dshape(ido,eita,dtipo)
		call reg_shape(ido,eita,pso,dpso)


				
        !...loop over inner integral
        val_inner_integral = 0.0d0
        do li=1,nti
          	!...coordinate of Gauss point
            zeta = xii(li)
            !...shape function & its derivatives at Gauss point
        	call reg_shape(idi,zeta,psi,dpsi)


			call findB2(phi,B2,xlo,xli,pso,dpso,psi,dpsi)
			call findE2(phi,nu,E2,xlo,xli,pso,dpso,psi,dpsi)
			!...compute val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                    do k=1,NDOF
                    	do l=1,NDOF
                            if (ido.eq.CTIP1) then
                            	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + (B2(l,k)*dtipo(i)*dpsi(j) + E2(l,k)*tipo(i)*psi(j))*(1.0d0+eitapr)*wii(li)
								

								
			    else
                              	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + (B2(l,k)*dtipo(i)*dpsi(j) + E2(l,k)*tipo(i)*psi(j))*(1.0d0-eitapr)*wii(li)

								
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
                        ek2bar(NDOF*(i-1)+k,NDOF*(j-1)+l) = ek2bar(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        	val_inner_integral(l,k,j,i)*wio(lo)
                    enddo
                enddo
            enddo
        enddo
    enddo !over loop of outer integral
	!-----------------------------------------------
    !...combine all together to have element stiffness matrix
    ek = ekbar1a + ekbar1b1 + ekbar1b2 + ekbar2 + ekbar3 + ek2bar
	End Subroutine ekd_adt_tip_reg

	!-----------------------------------------------
	Subroutine ekd_adt_tip_tip(region,ido,idi,ek,nto,nti,xio,xii,wio,wii,xlo,xli,adt_type,phi,ht,nu)
    !...10/10/10: this subroutine is temporarily kept here, not improved yet, so result of 2 crack elements will not be good
	
    IMPLICIT NONE

    INTEGER, INTENT(IN)		:: ido,idi,nto,nti,region
    INTEGER, INTENT(IN)		:: adt_type
    REAL(KIND=DBL), INTENT(OUT)	:: ek(:,:)
    REAL(KIND=DBL), INTENT(IN)	:: xio(:),xii(:),wio(:),wii(:)
    REAL(KIND=DBL), INTENT(IN)	:: xlo(:,:),xli(:,:)
	REAL(KIND=DBL), INTENT(IN)	:: phi,ht,nu
	!...local variables
    INTEGER		:: neo,nei,lo,li
    REAL(KIND=DBL)	:: ekbar1(3*NDOF,3*NDOF),ekbar2(3*NDOF,3*NDOF),ek2bar(3*NDOF,3*NDOF)
    REAL(KIND=DBL)	:: alpha,beta,eita,zeta,gamma,deltaa,eitapr,zetapr,eitabar
    REAL(KIND=DBL)  :: rho
    REAL(KIND=DBL)	:: tipo(3),tipi(3),dtipo(3),dtipi(3),pso(3),dpso(3),psi(3),dpsi(3)
    INTEGER			:: i,j,k,l
    REAL(KIND=DBL)	:: val_inner_integral(NDOF,NDOF,3,3)
    REAL(KIND=DBL)	:: B1(NDOF,NDOF),B2(NDOF,NDOF)
	REAL(KIND=DBL)	:: E1(NDOF,NDOF),E2(NDOF,NDOF)
	REAL(KIND=DBL)	:: field(2),source(2),r,temp(2)
   
	
	!...get number of nodes of outer/inner elements
    neo = NODE(ido)
    nei = NODE(idi)
	!...compute matrix ekbar2: singular if adt_type=2
	ekbar2=0.0d0

	!...begins outer integral loop
	do lo = 1,nto
    	!...coordinate of integration point
		if ((adt_type.eq.2).or.(adt_type.eq.4)) then
          	!...need stretching of gamma to rho to improve ln integration
            rho = xio(lo)
            gamma = 0.5d0*(rho+1.0d0)*(rho+1.0d0)-1.0d0
        elseif ((adt_type.eq.1).or.(adt_type.eq.3)) then
        	!...no need to stretch gamma
			gamma = xio(lo)
        else
          	print*,'error in ekd_adt_tip_tip: adt_type out of range'
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
            !...distinguish adjacent types between 2/4.0d0 and 1/3
            if ((adt_type.eq.3).or.(adt_type.eq.4)) then
              	eitabar = alpha + beta
                eitapr = -1.0d0*eitabar
            else
            	eitapr = alpha + beta
            endif
            zetapr = alpha - beta
            !...transform to get rid of singularity of tip shape function
            if (ido.eq.CTIP1) then
            	eita = 0.5d0*(eitapr+1.0d0)*(eitapr+1.0d0)-1.0d0
            elseif (ido.eq.CTIP2) then
            	eita = 1.0d0-0.5d0*(1.0d0-eitapr)*(1.0d0-eitapr)
            else
              	print*,'ekd_adt_tip_tip: id of outer element is out of range'
                stop
            endif
            if (idi.eq.CTIP1) then
                zeta = 0.5d0*(zetapr+1.0d0)*(zetapr+1.0d0)-1.0d0
            elseif (idi.eq.CTIP2) then
                zeta = 1.0d0-0.5d0*(1.0d0-zetapr)*(1.0d0-zetapr)
            else
              	print*,'ekd_adt_tip_tip: id of inner element is out of range'
                stop
            endif
            !...obtain value of shape function and its derivatives of 
            call tip_shape(ido,eita,tipo)
            call tip_shape(idi,zeta,tipi)
            call reg_shape(ido,eita,pso,dpso)
			
			call tip_dshape(ido,eita,dtipo)
            call tip_dshape(idi,zeta,dtipi)
			call reg_shape(idi,zeta,psi,dpsi)


			
			call findB1(phi,B1,xlo,xli,pso,dpso,psi,dpsi)
			call findE1(phi,nu,E1,xlo,xli,pso,dpso,psi,dpsi)

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
        
                         !...compute r
                         do k=1,2
              	             temp(k) = field(k) - source(k)
                         enddo
                         r = dsqrt(temp(1)*temp(1) + temp(2)*temp(2))

            !...calculate val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                    do k=1,NDOF
                    	do l=1,NDOF
                            if ((ido.eq.CTIP1).and.(idi.eq.CTIP2)) then
                        	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + (B1(l,k)*dtipo(i)*dtipi(j) + E1(l,k)*tipo(i)*tipi(j))*(eitapr+1.0d0)*(1.0d0-zetapr)*dlog(r)*wii(li)

								
			    elseif ((ido.eq.CTIP2).and.(idi.eq.CTIP1)) then
                              	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + (B1(l,k)*dtipo(i)*dtipi(j) + E1(l,k)*tipo(i)*tipi(j))*(1.0d0-eitapr)*(1.0d0+zetapr)*dlog(r)*wii(li)

								
			    else
                              	print*,'ekd_adt_tip_tip: impossible to have 2 tip elements with same type'
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
                    	if ((adt_type.eq.2).or.(adt_type.eq.4)) then
                        	ekbar2(NDOF*(i-1)+k,NDOF*(j-1)+l) = ekbar2(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		((1.0d0+rho)**3.0d0)*val_inner_integral(l,k,j,i)*wio(lo)
                        elseif ((adt_type.eq.1).or.(adt_type.eq.3)) then
                    		ekbar2(NDOF*(i-1)+k,NDOF*(j-1)+l) = ekbar2(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
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
      	ekbar2 = 0.25d0*ekbar2
    elseif ((adt_type.eq.1).or.(adt_type.eq.3)) then
    	ekbar2 = 0.5d0*ekbar2
    else
      	print*,'adt_type out of range'
        stop
    endif

	!--------------------------------------------------
    !...compute matrix ekbar1: singular if adt_type=1
	ekbar1=0.0d0
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
          	print*,'error in ekd_adt_tip_tip: adt_type out of range'
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
            !...distinguish adjacent types between 2/4.0d0 and 1/3
            if ((adt_type.eq.3).or.(adt_type.eq.4)) then
              	eitabar = alpha + beta
                eitapr = -1.0d0*eitabar
            else
            	eitapr = alpha + beta
            endif
            zetapr = alpha - beta
            !...transform to get rid of singularity of tip shape function
            if (ido.eq.CTIP1) then
            	eita = 0.5d0*(eitapr+1.0d0)*(eitapr+1.0d0)-1.0d0
            elseif (ido.eq.CTIP2) then
            	eita = 1.0d0-0.5d0*(1.0d0-eitapr)*(1.0d0-eitapr)
            else
              	print*,'ekd_adt_tip_tip: id of outer element is out of range'
                stop
            endif
            if (idi.eq.CTIP1) then
                zeta = 0.5d0*(zetapr+1.0d0)*(zetapr+1.0d0)-1.0d0
            elseif (idi.eq.CTIP2) then
                zeta = 1.0d0-0.5d0*(1.0d0-zetapr)*(1.0d0-zetapr)
            else
              	print*,'ekd_adt_tip_tip: id of inner element is out of range'
                stop
            endif
            !...obtain value of shape function and its derivatives of 
            call tip_shape(ido,eita,tipo)
            call tip_shape(idi,zeta,tipi)
			call reg_shape(ido,eita,pso,dpso)

            call tip_dshape(ido,eita,dtipo)
            call tip_dshape(idi,zeta,dtipi)
			call reg_shape(idi,zeta,psi,dpsi)

			
			call findB1(phi,B1,xlo,xli,pso,dpso,psi,dpsi)
			call findE1(phi,nu,E1,xlo,xli,pso,dpso,psi,dpsi)

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
        
                         !...compute r
                         do k=1,2
              	             temp(k) = field(k) - source(k)
                         enddo
                         r = dsqrt(temp(1)*temp(1) + temp(2)*temp(2))

            !...calculate val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                    do k=1,NDOF
                    	do l=1,NDOF
                            if ((ido.eq.CTIP1).and.(idi.eq.CTIP2)) then
                                val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + (B1(l,k)*dtipo(i)*dtipi(j) + E1(l,k)*tipo(i)*tipi(j))*(eitapr+1.0d0)*(1.0d0-zetapr)*dlog(r)*wii(li)

								
			    elseif ((ido.eq.CTIP2).and.(idi.eq.CTIP1)) then
                              	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + (B1(l,k)*dtipo(i)*dtipi(j) + E1(l,k)*tipo(i)*tipi(j))*(1.0d0-eitapr)*(1.0d0+zetapr)*dlog(r)*wii(li)

								
			    else
                              	print*,'ekd_adt_tip_tip: impossible to have 2 tip elements with same type'
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
                        	ekbar1(NDOF*(i-1)+k,NDOF*(j-1)+l) = ekbar1(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		((1.0d0-rho)**3.0d0)*val_inner_integral(l,k,j,i)*wio(lo)
                        elseif ((adt_type.eq.2).or.(adt_type.eq.4)) then
                    		ekbar1(NDOF*(i-1)+k,NDOF*(j-1)+l) = ekbar1(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
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
      	ekbar1 = 0.25d0*ekbar1
    elseif ((adt_type.eq.2).or.(adt_type.eq.4)) then
    	ekbar1 = 0.5d0*ekbar1
    else
      	print*,'adt_type out of range'
        stop
    endif
    !--------------------------------------------------
    !...evaluate ek2bar (regular integral)
    ek2bar = 0.0d0
    do lo=1,nto
      	!...coordinate of Gauss point
        eitapr = xio(lo)
        !...transform to get rid of singularity at tip shape functions
        if (ido.eq.CTIP1) then
        	eita = 0.5d0*(eitapr+1.0d0)*(eitapr+1.0d0)-1.0d0
        elseif (ido.eq.CTIP2) then
           	eita = 1.0d0-0.5d0*(1.0d0-eitapr)*(1.0d0-eitapr)
        else
           	print*,'ekd_adt_tip_tip: id of outer element is out of range'
            stop
        endif
        !...shape function & its derivatives at Gauss point
        call tip_shape(ido,eita,tipo)
        call tip_dshape(ido,eita,dtipo)
		call reg_shape(ido,eita,pso,dpso)


		
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
           		print*,'ekd_adt_tip_tip: id of inner element is out of range'
            	stop
        	endif
            !...shape function & its derivatives at Gauss point
        	call tip_shape(idi,zeta,tipi)
            call tip_dshape(idi,zeta,dtipi)
        	call reg_shape(idi,zeta,psi,dpsi)
			


						
            call findB2(phi,B2,xlo,xli,pso,dpso,psi,dpsi)
			call findE2(phi,nu,E2,xlo,xli,pso,dpso,psi,dpsi)
			!...compute val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                	do k=1,NDOF
                    	do l=1,NDOF
                            if ((ido.eq.CTIP1).and.(idi.eq.CTIP2)) then
                            	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + (B2(l,k)*dtipo(i)*dtipi(j) + E2(l,k)*tipo(i)*tipi(j))*(1.0d0+eitapr)*(1.0d0-zetapr)*wii(li)
		
								 
								
                            elseif ((ido.eq.CTIP2).and.(idi.eq.CTIP1)) then
                              	val_inner_integral(l,k,j,i) = val_inner_integral(l,k,j,i) + (B2(l,k)*dtipo(i)*dtipi(j) +  E2(l,k)*tipo(i)*tipi(j))*(1.0d0-eitapr)*(1.0d0+zetapr)*wii(li)

								 
			    else
                              	print*,'adt_tip_tip: impossible to have 2 tip elements with same type'
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
                        ek2bar(NDOF*(i-1)+k,NDOF*(j-1)+l) = ek2bar(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        val_inner_integral(l,k,j,i)*wio(lo)
                    enddo
                enddo
            enddo
        enddo
    enddo !over loop of outer integral
    
    !...element stiffness matrix
    ek = ekbar1 + ekbar2 + ek2bar
	End Subroutine ekd_adt_tip_tip
	!--------------------------------------------------
	Subroutine ekd_sed_reg_reg(region,ido,idi,ek,nto,nti,xio,xii,wio,wii,xlo,xli,phi,ht,nu)
	
    IMPLICIT NONE

    INTEGER, INTENT(IN)		:: ido,idi,nto,nti,region
    REAL(KIND=DBL), INTENT(OUT)	:: ek(:,:)
    REAL(KIND=DBL), INTENT(IN)	:: xio(:),xii(:),wio(:),wii(:)
    REAL(KIND=DBL), INTENT(IN)	:: xlo(:,:),xli(:,:)
	REAL(KIND=DBL), INTENT(IN)	:: phi,ht,nu
	!...local variables
    INTEGER		:: neo,nei,lo,li
    REAL(KIND=DBL)	:: ekbar(3*NDOF,3*NDOF),ek2bar(3*NDOF,3*NDOF)
    REAL(KIND=DBL)	:: eita,zeta
    REAL(KIND=DBL)	:: pso(3),dpso(3),psi(3),dpsi(3)
    INTEGER			:: i,j,k,l
    REAL(KIND=DBL)	:: val_inner_integral1(NDOF,NDOF,3,3),val_inner_integral2(NDOF,NDOF,3,3)
    REAL(KIND=DBL)	:: B1(NDOF,NDOF),B2(NDOF,NDOF)
    REAL(KIND=DBL)	:: E1(NDOF,NDOF),E2(NDOF,NDOF)
	REAL(KIND=DBL)	:: field(2),source(2),r,temp(2)
    
	
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


	    		
		!...initialize
        val_inner_integral1 = 0.0d0
        val_inner_integral2 = 0.0d0
		!...loop over inner integral (Gauss integration)
		do li =1,nti
			!...coordinate of Gauss point
        	zeta = xii(li)
            !...obtain value of shape function and its derivatives of 
            call reg_shape(idi,zeta,psi,dpsi)

	
			call findB1(phi,B1,xlo,xli,pso,dpso,psi,dpsi)
			
			call findB2(phi,B2,xlo,xli,pso,dpso,psi,dpsi)
			call findE1(phi,nu,E1,xlo,xli,pso,dpso,psi,dpsi)
                        call findE2(phi,nu,E2,xlo,xli,pso,dpso,psi,dpsi)

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
        
                         !...compute r
                         do k=1,2
              	             temp(k) = field(k) - source(k)
                         enddo
                         r = dsqrt(temp(1)*temp(1) + temp(2)*temp(2))

			!...calculate val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                    do k=1,NDOF
                    	do l=1,NDOF
                            val_inner_integral1(l,k,j,i) = val_inner_integral1(l,k,j,i) + (B1(l,k)*dpso(i)*dpsi(j) + E1(l,k)*pso(i)*psi(j))*dlog(r)*wii(li)
                            val_inner_integral2(l,k,j,i) = val_inner_integral2(l,k,j,i) + (B2(l,k)*dpso(i)*dpsi(j) + E2(l,k)*pso(i)*psi(j))*wii(li)

							
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
	End Subroutine ekd_sed_reg_reg
	!--------------------------------------------------
	Subroutine ekd_sed_reg_tip(region,ido,idi,ek,nto,nti,xio,xii,wio,wii,xlo,xli,phi,ht,nu)
	
    IMPLICIT NONE

    INTEGER, INTENT(IN)		:: ido,idi,nto,nti,region
    REAL(KIND=DBL), INTENT(OUT)	:: ek(:,:)
    REAL(KIND=DBL), INTENT(IN)	:: xio(:),xii(:),wio(:),wii(:)
    REAL(KIND=DBL), INTENT(IN)	:: xlo(:,:),xli(:,:)
	REAL(KIND=DBL), INTENT(IN)	:: phi,ht,nu
	!...local variables
    INTEGER		:: neo,nei,lo,li
    REAL(KIND=DBL)	:: ekbar(3*NDOF,3*NDOF),ek2bar(3*NDOF,3*NDOF)
    REAL(KIND=DBL)	:: eita,zeta,zetapr
    REAL(KIND=DBL)	:: pso(3),dpso(3),tipi(3),dtipi(3),psi(3),dpsi(3)
    INTEGER			:: i,j,k,l
    REAL(KIND=DBL)	:: val_inner_integral1(NDOF,NDOF,3,3),val_inner_integral2(NDOF,NDOF,3,3)
    REAL(KIND=DBL)	:: B1(NDOF,NDOF),B2(NDOF,NDOF)
	REAL(KIND=DBL)	:: E1(NDOF,NDOF),E2(NDOF,NDOF)
REAL(KIND=DBL)	:: field(2),source(2),r,temp(2)
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
              	print*,'ekd_sed_reg_tip: id of inner element is out of range'
                stop
            endif
            !...obtain value of shape function and its derivatives of 
            call tip_shape(idi,zeta,tipi)
            call tip_dshape(idi,zeta,dtipi)
            call reg_shape(idi,zeta,psi,dpsi)
			
	
		        call findB1(phi,B1,xlo,xli,pso,dpso,psi,dpsi)
			
			call findB2(phi,B2,xlo,xli,pso,dpso,psi,dpsi)
			call findE1(phi,nu,E1,xlo,xli,pso,dpso,psi,dpsi)
                        call findE2(phi,nu,E2,xlo,xli,pso,dpso,psi,dpsi)

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
        
                         !...compute r
                         do k=1,2
              	             temp(k) = field(k) - source(k)
                         enddo
                         r = dsqrt(temp(1)*temp(1) + temp(2)*temp(2))

			!...calculate val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                    do k=1,NDOF
                    	do l=1,NDOF
                            if (idi.eq.CTIP1) then
                        	val_inner_integral1(l,k,j,i) = val_inner_integral1(l,k,j,i) + (B1(l,k)*dpso(i)*dtipi(j) + E1(l,k)*pso(i)*tipi(j))*(1.0d0+zetapr)*dlog(r)*wii(li)
                            	val_inner_integral2(l,k,j,i) = val_inner_integral2(l,k,j,i) + (B2(l,k)*dpso(i)*dtipi(j) + E2(l,k)*pso(i)*tipi(j))*(1.0d0+zetapr)*wii(li)

								
			    elseif (idi.eq.CTIP2) then
                            	val_inner_integral1(l,k,j,i) = val_inner_integral1(l,k,j,i) + (B1(l,k)*dpso(i)*dtipi(j) + E1(l,k)*pso(i)*tipi(j))*(1.0d0-zetapr)*dlog(r)*wii(li)
                            	val_inner_integral2(l,k,j,i) = val_inner_integral2(l,k,j,i) + (B2(l,k)*dpso(i)*dtipi(j) + E2(l,k)*pso(i)*tipi(j))*(1.0d0-zetapr)*wii(li)

								
			    else
                              	print*,'sed_reg_tip: id of inner element is out of range'
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
	End Subroutine ekd_sed_reg_tip
	!--------------------------------------------------
	Subroutine ekd_sed_tip_reg(region,ido,idi,ek,nto,nti,xio,xii,wio,wii,xlo,xli,phi,ht,nu)
	
    IMPLICIT NONE

    INTEGER,INTENT(IN)		    :: ido,idi,nto,nti,region
    REAL(KIND=DBL),INTENT(OUT)	:: ek(:,:)
    REAL(KIND=DBL),INTENT(IN)	:: xio(:),xii(:),wio(:),wii(:)
    REAL(KIND=DBL),INTENT(IN)	:: xlo(:,:),xli(:,:)
	REAL(KIND=DBL), INTENT(IN)	:: phi,ht,nu
	!...local variables
    INTEGER		:: neo,nei,lo,li
    REAL(KIND=DBL)	:: ekbar(3*NDOF,3*NDOF),ek2bar(3*NDOF,3*NDOF)
    REAL(KIND=DBL)	:: eitapr,eita,zeta
    REAL(KIND=DBL)	:: tipo(3),dtipo(3),psi(3),dpsi(3),pso(3),dpso(3)
    INTEGER			:: i,j,k,l
    REAL(KIND=DBL)	:: val_inner_integral1(NDOF,NDOF,3,3),val_inner_integral2(NDOF,NDOF,3,3)
    REAL(KIND=DBL)	:: B1(NDOF,NDOF),B2(NDOF,NDOF)
	REAL(KIND=DBL)	:: E1(NDOF,NDOF),E2(NDOF,NDOF)
REAL(KIND=DBL)	:: field(2),source(2),r,temp(2)
	neo = NODE(ido)
    nei = NODE(idi)
	!...compute ekbar (regular integral)
	ekbar = 0.0d0
    ek2bar = 0.0d0
	!...begins outer integral loop
	do lo = 1,nto
    	!...coordinate of integration point
		eitapr = xio(lo)
        !...transform to remove singularity of tip shape functions
        if (ido.eq.CTIP1) then
        	eita = 0.5d0*(eitapr+1.0d0)*(eitapr+1.0d0)-1.0d0
        elseif (ido.eq.CTIP2) then
        	eita = 1.0d0-0.5d0*(1.0d0-eitapr)*(1.0d0-eitapr)
        else
          	print*,'ekd_sed_tip_reg: id of outer element is out of range'
            stop
        endif
        !...obtain value of shape function and its derivatives of 
        call tip_shape(ido,eita,tipo)
        call tip_dshape(ido,eita,dtipo)
		call reg_shape(ido,eita,pso,dpso)


		
		!...initialize
        val_inner_integral1 = 0.0d0
        val_inner_integral2 = 0.0d0
		!...loop over inner integral (Gauss integration)
		do li =1,nti
			!...coordinate of Gauss point
        	zeta = xii(li)
            !...obtain value of shape function and its derivatives of 
            call reg_shape(idi,zeta,psi,dpsi)



                        call findB1(phi,B1,xlo,xli,pso,dpso,psi,dpsi)
			
			call findB2(phi,B2,xlo,xli,pso,dpso,psi,dpsi)
			call findE1(phi,nu,E1,xlo,xli,pso,dpso,psi,dpsi)
                        call findE2(phi,nu,E2,xlo,xli,pso,dpso,psi,dpsi)

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
        
                         !...compute r
                         do k=1,2
              	             temp(k) = field(k) - source(k)
                         enddo
                         r = dsqrt(temp(1)*temp(1) + temp(2)*temp(2))

			!...calculate val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                    do k=1,NDOF
                    	do l=1,NDOF
                            val_inner_integral1(l,k,j,i) = val_inner_integral1(l,k,j,i) + (B1(l,k)*dtipo(i)*dpsi(j) + E1(l,k)*tipo(i)*psi(j))*dlog(r)*wii(li)
                            val_inner_integral2(l,k,j,i) = val_inner_integral2(l,k,j,i) + (B2(l,k)*dtipo(i)*dpsi(j) + E2(l,k)*tipo(i)*psi(j))*wii(li)

							
			enddo
                    enddo
                enddo
            enddo
        enddo !of loop over inner integral
        do i=1,neo
          	do j=1,nei
            	do k=1,NDOF
                	do l=1,NDOF
                    	if (ido.eq.CTIP1) then
                    		ekbar(NDOF*(i-1)+k,NDOF*(j-1)+l) = ekbar(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		val_inner_integral1(l,k,j,i)*(1.0d0+eitapr)*wio(lo)
                        	ek2bar(NDOF*(i-1)+k,NDOF*(j-1)+l) = ek2bar(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		val_inner_integral2(l,k,j,i)*(1.0d0+eitapr)*wio(lo)
                        elseif (ido.eq.CTIP2) then
                        	ekbar(NDOF*(i-1)+k,NDOF*(j-1)+l) = ekbar(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		val_inner_integral1(l,k,j,i)*(1.0d0-eitapr)*wio(lo)
                        	ek2bar(NDOF*(i-1)+k,NDOF*(j-1)+l) = ek2bar(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		val_inner_integral2(l,k,j,i)*(1.0d0-eitapr)*wio(lo)
                        else
                          	print*,'ekd_sed_tip_reg: id of outer element is out of range'
                            stop
                        endif
                    enddo
                enddo
            enddo
        enddo
    enddo !of loop over outer integral
    
    !...element stiffness matrix
    ek = ekbar + ek2bar
	End Subroutine ekd_sed_tip_reg
	!--------------------------------------------------
	Subroutine ekd_sed_tip_tip(region,ido,idi,ek,nto,nti,xio,xii,wio,wii,xlo,xli,phi,ht,nu)
	
    IMPLICIT NONE

    INTEGER,INTENT(IN)		    :: ido,idi,nto,nti,region
    REAL(KIND=DBL),INTENT(OUT)	:: ek(:,:)
    REAL(KIND=DBL),INTENT(IN)	:: xio(:),xii(:),wio(:),wii(:)
    REAL(KIND=DBL),INTENT(IN)	:: xlo(:,:),xli(:,:)
	REAL(KIND=DBL), INTENT(IN)	:: phi,ht,nu
	!...local variables
    INTEGER		:: neo,nei,lo,li
    REAL(KIND=DBL)	:: ekbar(3*NDOF,3*NDOF),ek2bar(3*NDOF,3*NDOF)
    REAL(KIND=DBL)	:: eitapr,eita,zetapr,zeta
    REAL(KIND=DBL)	:: tipo(3),tipi(3),dtipo(3),dtipi(3),pso(3),dpso(3),psi(3),dpsi(3)
    INTEGER			:: i,j,k,l
	REAL(KIND=DBL)	:: val_inner_integral1(NDOF,NDOF,3,3),val_inner_integral2(NDOF,NDOF,3,3)
	REAL(KIND=DBL)	:: B1(NDOF,NDOF),B2(NDOF,NDOF)
	REAL(KIND=DBL)	:: E1(NDOF,NDOF),E2(NDOF,NDOF)
REAL(KIND=DBL)	:: field(2),source(2),r,temp(2)
	neo = NODE(ido)
    nei = NODE(idi)
	!...compute ekbar (regular integral)
	ekbar = 0.0d0
    ek2bar = 0.0d0

	!...begins outer integral loop
	do lo = 1,nto
    	!...coordinate of integration point
		eitapr = xio(lo)
        !...transform to remove singularity of tip shape functions
        if (ido.eq.CTIP1) then
        	eita = 0.5d0*(eitapr+1.0d0)*(eitapr+1.0d0)-1.0d0
        elseif (ido.eq.CTIP2) then
        	eita = 1.0d0-0.5d0*(1.0d0-eitapr)*(1.0d0-eitapr)
        else
          	print*,'ekd_sed_tip_tip: id of outer element is out of range'
            stop
        endif
        !...obtain value of shape function and its derivatives of 
        call tip_shape(ido,eita,tipo)
        call tip_dshape(ido,eita,dtipo)
		call reg_shape(ido,eita,pso,dpso)


		
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
              	print*,'ekd_sed_tip_tip: id of inner element is out of range'
                stop
            endif
            !...obtain value of shape function and its derivatives of 
            call tip_shape(idi,zeta,tipi)
            call tip_dshape(idi,zeta,dtipi)
			call reg_shape(idi,zeta,psi,dpsi)


			
			call findB1(phi,B1,xlo,xli,pso,dpso,psi,dpsi)
			
			call findB2(phi,B2,xlo,xli,pso,dpso,psi,dpsi)
			call findE1(phi,nu,E1,xlo,xli,pso,dpso,psi,dpsi)
                        call findE2(phi,nu,E2,xlo,xli,pso,dpso,psi,dpsi)

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
        
                         !...compute r
                         do k=1,2
              	             temp(k) = field(k) - source(k)
                         enddo
                         r = dsqrt(temp(1)*temp(1) + temp(2)*temp(2))

			!...calculate val_inner_integral(l,k,j,i)
            do i=1,neo
              	do j=1,nei
                    do k=1,NDOF
                    	do l=1,NDOF
                            if (idi.eq.CTIP1) then
                        	val_inner_integral1(l,k,j,i) = val_inner_integral1(l,k,j,i) + (B1(l,k)*dtipo(i)*dtipi(j) + E1(l,k)*tipo(i)*tipi(j))*(1.0d0+zetapr)*dlog(r)*wii(li)
                            	val_inner_integral2(l,k,j,i) = val_inner_integral2(l,k,j,i) + (B2(l,k)*dtipo(i)*dtipi(j) + E2(l,k)*tipo(i)*tipi(j))*(1.0d0+zetapr)*wii(li)

			    elseif (idi.eq.CTIP2) then
                            	val_inner_integral1(l,k,j,i) = val_inner_integral1(l,k,j,i) + (B1(l,k)*dtipo(i)*dtipi(j) + E1(l,k)*tipo(i)*tipi(j))*(1.0d0-zetapr)*dlog(r)*wii(li)
                            	val_inner_integral2(l,k,j,i) = val_inner_integral2(l,k,j,i) + (B2(l,k)*dtipo(i)*dtipi(j) + E2(l,k)*tipo(i)*tipi(j))*(1.0d0-zetapr)*wii(li)

								 
			    else
                              	print*,'ekd_sed_tip_tip: id of inner element is out of range'
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
                    	if (ido.eq.CTIP1) then
                    		ekbar(NDOF*(i-1)+k,NDOF*(j-1)+l) = ekbar(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		val_inner_integral1(l,k,j,i)*(1.0d0+eitapr)*wio(lo)
                        	ek2bar(NDOF*(i-1)+k,NDOF*(j-1)+l) = ek2bar(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		val_inner_integral2(l,k,j,i)*(1.0d0+eitapr)*wio(lo)
                        elseif (ido.eq.CTIP2) then
                        	ekbar(NDOF*(i-1)+k,NDOF*(j-1)+l) = ekbar(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		val_inner_integral1(l,k,j,i)*(1.0d0-eitapr)*wio(lo)
                        	ek2bar(NDOF*(i-1)+k,NDOF*(j-1)+l) = ek2bar(NDOF*(i-1)+k,NDOF*(j-1)+l) + &
                        		val_inner_integral2(l,k,j,i)*(1.0d0-eitapr)*wio(lo)
                        else
                          	print*,'ekd_sed_tip_tip: id of outer element is out of range'
                            stop
                        endif
                    enddo
                enddo
            enddo
        enddo
    enddo !of loop over outer integral
    
    !...element stiffness matrix
    ek = ekbar + ek2bar
	End Subroutine ekd_sed_tip_tip
End Subroutine ek_d

!----------------------
Subroutine Normal_Position(xl,ps,dps,position,normal,dj)
   	!returns source point or field point,and normal
    REAL(KIND=DBL), INTENT(IN)	:: xl(2,3)
!...coordinates of nodes in real space
    REAL(KIND=DBL), INTENT(IN):: ps(3),dps(3)
!...shape function and derivative of shape function
    
    REAL(KIND=DBL), INTENT(OUT)	:: position(2),normal(3),dj
!...postion of source or field point i.e.position of point on real space of the corresponding gauss point,normal at the point on real space 
!...corresponding to gauss point,jacobian of mapping at the corresponding gauss point
    !...local variables
    INTEGER					:: i,j,k
    REAL(KIND=DBL)	:: dxds(2),tangent(2)
	!--------------
	dxds=0.0d0
	normal = 0.0d0
	dj = 0.0d0
!...initialising all local variables
    do j = 1,3
         do k = 1,2
              dxds(k) = dxds(k) + xl(k,j)*dps(j)
         enddo
    enddo
    dj = dsqrt(dxds(1)*dxds(1) + dxds(2)*dxds(2))
    !...tangent vector
    do j = 1,2
         tangent(j) = dxds(j)/dj
    enddo
    !...normal = tangent x e3
    normal(1) = tangent(2)
    normal(2) = -1.0d0*tangent(1)
!...this chunk of code is same as that in subroutine nodenormal in prep.f90
	!-----------
	position=0.0d0
    do i=1,3
        do k=1,2
            position(k) = position(k) + xl(k,i)*ps(i)
        enddo
    enddo
	!-----------
    
End Subroutine Normal_Position
!---------------------------------

Subroutine findr(field,source,r,temp)

REAL(KIND=DBL), INTENT(IN)	:: field(2),source(2)
!...coordinate of source and field point corresponding to the gauss point pair under consideration
   real(kind=dbl), intent(out) :: r,temp(2)
!...temp gives components of distance, r simply gives the distance
   !...local variables

   !real(kind=dbl)::temp(2)
   
   integer :: k
   r=0.0d0
   temp=0.0d0 
   !...compute the vector r = xi - y
   do k=1,2
       temp(k) = field(k) - source(k)
   enddo
   r = dsqrt(temp(1)*temp(1) + temp(2)*temp(2))

end subroutine findr

!--------------------------------------------
SUBROUTINE findB1(phi,B1,xlo,xli,pso,dpso,psi,dpsi)

   REAL(KIND=DBL), INTENT(IN):: phi,xlo(2,3),xli(2,3),pso(3),dpso(3),psi(3),dpsi(3)
   
   real(kind=dbl), intent(out)::B1(3,3) 
!...B1 is the kernel component that multiplies with ln (r) and derivative of shape function
   !...local variables
   real(kind=dbl):: x1,w1,a,field(2),source(2),r,temp(2),normal_y(3),jacobo,normal_xi(3),jacobi
   !...x1 is part of x that multiplies with ln r
   !...w1 is part of w that multiplies with ln r. 
   !...a is r^2
   integer:: l,k,kdelta(2,2)
!...kdelta is simple kronecker delta
   
   call Normal_Position(xlo,pso,dpso,source,normal_y,jacobo)

   call Normal_Position(xli,psi,dpsi,field,normal_xi,jacobi)   

   call findr(field,source,r,temp)   

   a = r*r
   w1 = 2.0d0*45045.0d0*a*ht**7.0d0*(495.0d0*a**2.0d0 - 1540.0d0*a*ht**2.0d0 + 504.0d0*ht**4.0d0)
   w1 = w1/(6936930.0d0*a*ht**6.0d0)
	
   kdelta(1,1)=1
   kdelta(2,2)=1
   kdelta(1,2)=0
   kdelta(2,1)=0

   x1 = -2.0d0*45045.0d0*ht**7.0d0*(2145.0d0*a**3.0d0 - 10010.0d0*a**2.0d0*ht**2.0d0 + 6552.0d0*a*ht**4.0d0 + 6336.0d0*ht**6.0d0)
   x1 = x1/(180360180.0d0*ht**6.0d0)

   B1=0.0d0

   do l = 1,2
        do k = 1,2
		     B1(l,k)=phi/(ht**6.0d0)*(w1*temp(l)*temp(k)+x1*kdelta(k,l))
        enddo
   enddo
B1(3,3)=1.0d0
END SUBROUTINE findB1
!---------------------------

SUBROUTINE findB2(phi,B2,xlo,xli,pso,dpso,psi,dpsi)

   REAL(KIND=DBL), INTENT(IN):: phi,xlo(2,3),xli(2,3),pso(3),dpso(3),psi(3),dpsi(3)
   
   real(kind=dbl), intent(out)::B2(3,3) 
!...B2 is the kernel component that multiplies with derivative of shape function but does not multiply with ln (r)
   !...local variables
   real(kind=dbl):: x2,w2,a,field(2),source(2),r,temp(2),normal_y(3),jacobo,normal_xi(3),jacobi
   !...x2 is part of x that does not multiply with ln r
   !...w2 is part of w that does not multiply with ln r
   !...a = r^2

   integer:: l,k,kdelta(2,2)
!...kdelta is kronecker delta	
   kdelta(1,1)=1
   kdelta(2,2)=1
   kdelta(1,2)=0
   kdelta(2,1)=0

   call Normal_Position(xlo,pso,dpso,source,normal_y,jacobo)

   call Normal_Position(xli,psi,dpsi,field,normal_xi,jacobi)     

   call findr(field,source,r,temp)   

   a = r*r

   w2 = -1280.0d0*r**13.0d0 - 3171168.0d0*r**7.0d0*ht**6.0d0 + 41621580.0d0*r**5.0d0*ht**8.0d0 - &
                41621580.0d0*r**3.0d0*ht**10.0d0 + 1280.0d0*a**6.0d0*dsqrt(a + 4.0d0*ht**2.0d0) - 2560.0d0*a**5.0d0*ht**2.0d0*&
	        dsqrt(a + 4.0d0*ht**2.0d0) + 7680.0d0*a**4.0d0*ht**4.0d0*dsqrt(a + 4.0d0*ht**2.0d0) + &
	        3145568.0d0*a**3.0d0*ht**6.0d0*dsqrt(a + 4.0d0*ht**2.0d0) - 25577041.0d0*a**2.0d0*ht**8.0d0*&
	        dsqrt(a + 4.0d0*ht**2.0d0) + 14740488.0d0*a*ht**10.d0*dsqrt(a + 4.0d0*ht**2.0d0) + &
                2744280.0d0*ht**12.0d0*dsqrt(a + 4.0d0*ht**2.0d0)           
			   
   w2 = w2*4.0d0
   
   w2 = w2 - 2.0d0*45045.0d0*a*ht**7.0d0*(495.0d0*a**2.0d0 - 1540.0d0*a*ht**2.0d0 + 504.0d0*ht**4.0d0)*&
                dlog(2.0d0*ht + dsqrt(a + 4.0d0*ht**2.0))
   
   w2 = w2/(6936930.0d0*a*ht**6.0d0)

   !...evaluated w2			   			      

   x2 = -2560.0d0*r**13.0d0 - 11778624.0d0*r**7.0d0*ht**6.0d0 + 216432216.0d0*r**5.0d0*ht**8.0d0 - &
                360720360.0d0*r**3.0d0*ht**10.0d0 + 2560.0d0*a**6.0d0*dsqrt(a + 4.0d0*ht**2.0d0) - 5120.0d0*a**5.0d0*ht**2.0d0*&
	        dsqrt(a + 4.0d0*ht**2.0d0) + 15360.0d0*a**4.0d0*ht**4.0d0*dsqrt(a + 4.0d0*ht**2.0d0) + &
	        11727424.0d0*a**3.0d0*ht**6.0d0*dsqrt(a + 4.0d0*ht**2.0d0) - 143188739.0d0*a**2.0d0*ht**8.0d0*&
	        dsqrt(a + 4.0d0*ht**2.0d0) + 155053566.0d0*a*ht**10.d0*dsqrt(a + 4.0d0*ht**2.0d0) + &
                68798664.0d0*ht**12.0d0*dsqrt(a + 4.0d0*ht**2.0d0)
			   
   x2 = x2*-4.0d0
   
   x2 = x2 + 2.0d0*45045.0d0*ht**7.0d0*(2145.0d0*a**3.0d0 - 10010.0d0*a**2.0d0*ht**2.0d0 + 6552.0d0*a*ht**4.0d0 + 6336.0d0*ht**6.0d0)*&
               dlog(2.0d0*ht + dsqrt(a + 4.0d0*ht**2.0d0))

   x2 = x2/(180360180.0d0*ht**6.0d0)


   !...evaluated x2

   B2 = 0.0d0

   do l = 1,2
        do k = 1,2
             B2(l,k) = phi/(ht**6.0d0)*(w2*temp(l)*temp(k)+x2*kdelta(k,l))
        enddo
   enddo

END SUBROUTINE findB2
!---------------------

SUBROUTINE findE1(phi,nu,E1,xlo,xli,pso,dpso,psi,dpsi)

REAL(KIND=DBL), INTENT(IN):: phi,nu,xlo(2,3),xli(2,3),pso(3),dpso(3),psi(3),dpsi(3)

real(kind=dbl), intent(out)::E1(3,3) 
!...E1 is the kernel component that multiplies with ln (r) and shape function

!...local variables
real(kind=dbl):: z1,y1,temporary1,temporary2,temporary3,temporary4,temporary5,a,field(2),source(2),r,temp(2),normal_y(3),jacobo,normal_xi(3),jacobi
!...z1 is part of z that multiplies with ln r
!...y1 is part of y that multiplies with ln r
!...a =r^2

integer:: l,k,kdelta(2,2)
!...kdelta is kronecker delta

call Normal_Position(xlo,pso,dpso,source,normal_y,jacobo)

call Normal_Position(xli,psi,dpsi,field,normal_xi,jacobi)

call findr(field,source,r,temp)

a = r*r

E1 = 0.0d0

z1 = -2.0d0*ht**7.0d0*(1485.0d0*a**2.0d0 - 3080.0d0*a*ht**2.0d0 + 504.0d0*ht**4.0d0)
z1 = z1/2772.0d0

y1 = 2.0d0*3465.0d0*a*ht**7.0d0*(27.0d0*a - 28.0d0*ht**2.0d0)
y1 = y1/(43659.0d0*a)

temporary1 = 36.0d0*z1*phi/(ht**12.0d0)

temporary2 = 36.0d0*y1*phi/(ht**12.0d0) 
	
kdelta(1,1)=1
kdelta(2,2)=1
kdelta(1,2)=0
kdelta(2,1)=0

do l = 1,2
     do k = 1,2
	  temporary3 = (1.0d0-nu)*kdelta(1,k)*kdelta(2,l)*normal_xi(1) + 2.0d0*nu*kdelta(2,k)*kdelta(1,l)*normal_xi(1) + kdelta(k,l)*normal_xi(2)&
                       -(1.0d0+nu)*kdelta(1,k)*kdelta(1,l)*normal_xi(2)
          temporary4 = (1.0d0-nu)*kdelta(2,k)*kdelta(1,l)*normal_xi(2) + 2.0d0*nu*kdelta(1,k)*kdelta(2,l)*normal_xi(2) + kdelta(k,l)*normal_xi(1)&
                       -(1.0d0+nu)*kdelta(2,k)*kdelta(2,l)*normal_xi(1)
          temporary5 = (normal_y(1)*normal_xi(1) + normal_y(2)*normal_xi(2))*temp(l)*temp(k)
          E1(l,k) = temporary1*(normal_y(2)*temporary3 + normal_y(1)*temporary4) + temporary2*temporary5	  
          E1(l,k) = E1(l,k)*jacobo*jacobi
     enddo
enddo
E1(3,3)=1.0d0
END SUBROUTINE findE1
!---------------------

SUBROUTINE findE2(phi,nu,E2,xlo,xli,pso,dpso,psi,dpsi)

REAL(KIND=DBL), INTENT(IN):: phi,nu,xlo(2,3),xli(2,3),pso(3),dpso(3),psi(3),dpsi(3)

real(kind=dbl), intent(out)::E2(3,3) 
!...E2 is kernel component that multiplies with shape function but not with ln (r)

!...local variables
real(kind=dbl):: z2,y2,temporary1,temporary2,temporary3,temporary4,temporary5,a,field(2),source(2),r,temp(2),normal_y(3),jacobo,normal_xi(3),jacobi
!...z2 is part of z that does not multiply with ln r
!...y2 is part of y that does not multiply with ln r
!...a = r^2

integer:: l,k,kdelta(2,2)
!...kdelta is kronecker delta

call Normal_Position(xlo,pso,dpso,source,normal_y,jacobo)

call Normal_Position(xli,psi,dpsi,field,normal_xi,jacobi)

call findr(field,source,r,temp)

a = r*r

E2=0.0d0
	
kdelta(1,1)=1
kdelta(2,2)=1
kdelta(1,2)=0
kdelta(2,1)=0

z2 = 640.0d0*r**11.0d0 + 853776.0d0*r**5.0d0*ht**6.0d0 - 8004150.0d0*r**3.0d0*ht**8.0d0 + 4802490.0d0*r*ht**10.0d0 - &
     640.0d0*a**5.0d0*dsqrt(a + 4.0d0*ht**2.0d0) + 1280.0d0*a**4.0d0*ht**2.0d0*dsqrt(a + 4.0d0*ht**2.0d0) - &
     3840.0d0*a**3.0d0*ht**4.0d0*dsqrt(a + 4.0d0*ht**2.0d0) - 840976.0d0*a**2.0d0*ht**6.0d0*dsqrt(a + 4.0d0*ht**2.0d0) + &
     4521377.0d0*a*ht**8.0d0*dsqrt(a + 4.0d0*ht**2.0d0) - 1378566.0d0*ht**10.0d0*dsqrt(a + 4.0d0*ht**2.0d0)

z2 = z2/2401245.0d0

z2 = z2 + 2.0d0/2772.0d0*ht**7.0d0*(1485.0d0*a**2.0d0 - 3080.0d0*a*ht**2.0d0 + 504.0d0*ht**4.0d0)*dlog(2.0d0*ht + dsqrt(a + 4.0d0*ht**2.0)) 

y2 = -64.0d0*r**11.0d0 - 38808.0d0*r**5.0d0*ht**6.0d0 + 218295.0d0*r**3.0d0*ht**8.0d0 - 43659.0d0*r*ht**10.0d0 + &
     64.0d0*a**5.0d0*dsqrt(a + 4.0d0*ht**2.0d0) - 128.0d0*a**4.0d0*ht**2.0d0*dsqrt(a + 4.0d0*ht**2.0d0) + &
     384.0d0*a**3.0d0*ht**4.0d0*dsqrt(a + 4.0d0*ht**2.0d0) + 37528.0d0*a**2.0d0*ht**6.0d0*dsqrt(a + 4.0d0*ht**2.0d0) - &
     104321.0d0*a*ht**8.0d0*dsqrt(a + 4.0d0*ht**2.0d0) + 3969.0d0*ht**10.0d0*dsqrt(a + 4.0d0*ht**2.0d0)

y2 = y2*2.0d0

y2 = y2 - 2.0d0*3465.0d0*a*ht**7.0d0*(27.0d0*a - 28.0d0*ht**2.0d0)*dlog(2.0d0*ht + dsqrt(a + 4.0d0*ht**2.0)) 

y2 = y2/(43659.0d0*a)

temporary1 = 36.0d0*z2*phi/(ht**12.0d0)

temporary2 = 36.0d0*y2*phi/(ht**12.0d0)


do l = 1,2
     do k = 1,2
          temporary3 = (1.0d0-nu)*kdelta(1,k)*kdelta(2,l)*normal_xi(1) + 2.0d0*nu*kdelta(2,k)*kdelta(1,l)*normal_xi(1) + kdelta(k,l)*normal_xi(2)&
                       -(1.0d0+nu)*kdelta(1,k)*kdelta(1,l)*normal_xi(2)
          temporary4 = (1.0d0-nu)*kdelta(2,k)*kdelta(1,l)*normal_xi(2) + 2.0d0*nu*kdelta(1,k)*kdelta(2,l)*normal_xi(2) + kdelta(k,l)*normal_xi(1)&
                       -(1.0d0+nu)*kdelta(2,k)*kdelta(2,l)*normal_xi(1)
          temporary5 = (normal_y(1)*normal_xi(1)+normal_y(2)*normal_xi(2))*temp(l)*temp(k)
          E2(l,k) = temporary1*(normal_y(2)*temporary3 + normal_y(1)*temporary4) + temporary2*temporary5
          E2(l,k) = E2(l,k)*jacobo*jacobi    
     enddo
enddo

END SUBROUTINE findE2
!------------------------
