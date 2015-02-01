Subroutine Tstress
	USE DefinitionConstant
    USE GlobalData
    USE ElementMatrix
    IMPLICIT NONE
	!--------------------------------------------------
    !...local variables
	REAL(KIND=DBL)	:: xlo(2,3),xli(2,3),pso(3),dpso(3)
    REAL(KIND=DBL)	:: el(3,3),el_inv(3,3),rhs(9),eb(9,9),eat(9,9),dxds(2)
    REAL(KIND=DBL)	:: xt,wt,temp,detj,detl,jacobo
    INTEGER			:: i,j,k,l,m,n,k0
    INTEGER			:: ido,idi,nno,nni
    REAL(KIND=DBL)	:: sigma_u(3,3),dsigma_u(3),tangent(2)
    REAL(KIND=DBL)	:: epsilon_s
    INTEGER			:: ieo,iei,alpha,beta,I_sign
    !--------------------------------------------------
    !...this portion of code work for the version of T-stress with tip element only
    print*,'--------------------------------------------------'
    print*,'T-stress calculated by 1ST APPROACH, integrate on TIP ONLY'
    I_sign = 1
    !...loop over elements, but only work with tip element for outer integral
	k0 = 0
    do k = 1,total_region
    	do i=1,total_region_elem(k)
        	!...global element # of outer element
        	ieo = k0 + i
        	!...only continue the procedure if element ie is a tip element
    		if ((elemid(ieo).ne.CTIP1).and.(elemid(ieo).ne.CTIP2)) cycle
        	!...initialize solution sigma_u of tip element ie
        	sigma_u = 0.0d0
            !...get id of the tip element
            ido = elemid(ieo)
        	!...get number of nodes of outer element
        	nno = NODE(ido)
        	!...get nodal coordinates of outer element
        	do m = 1,nno
        		do n =1,2
            		xlo(n,m)=node_coor(n,elnode(m,ieo))
            	enddo
        	end do
        	!--------------------------------------------------
        	!...compute matrix el
			el = 0.0d0
			!...begins the integration loop
			do l = 1,ng_ef
        		!...get Gauss point and weight from global data xi(:,:) and wi(:,:)
            	xt = xi(l,ng_ef)
            	wt = wi(l,ng_ef)
            	!...use regular shape functions for both sigma_u and test function
				call reg_shape(ido,xt,pso,dpso)
				!...calculate derivative of position vector: dxds(2)
    			dxds = 0.0d0
				do j = 1,nno
	    			do n = 1,2
            			dxds(n) = dxds(n)+xlo(n,j)*dpso(j)
            		enddo
				enddo
				!...compute jacobian
				detj= dsqrt(dxds(1)*dxds(1) + dxds(2)*dxds(2))
				!...calculate el(3,3)
        		temp=detj*wt
        		do j = 1,nno
            		do n = 1,nno
                		!...Note: temporarily delete 0.5 to get u directly (not sum of u)
                		!el(i,n) = el(i,n) + 0.5d0*pso(i)*pso(n)*temp
                		el(j,n) = el(j,n) + pso(j)*pso(n)*temp
                	enddo
            	enddo
			enddo !of loop over integration points
        	!--------------------------------------------------
        	!...compute inverse of matrix el
        	detl = el(1,1)*el(2,2)*el(3,3) + el(1,2)*el(2,3)*el(3,1) + el(1,3)*el(2,1)*el(3,2) - &
				   el(3,1)*el(2,2)*el(1,3)-	el(3,2)*el(2,3)*el(1,1) - el(3,3)*el(2,1)*el(1,2)
        	if (abs(detl).lt.SMALL_NUM) then
        		print*,'tstress.f90 error: determinant of el is zero!'
        		stop
        	endif
			temp = 1.0d0/detl
			el_inv(1,1)=(el(2,2)*el(3,3)-el(3,2)*el(2,3))*temp
			el_inv(2,1)=(el(3,1)*el(2,3)-el(2,1)*el(3,3))*temp
			el_inv(3,1)=(el(2,1)*el(3,2)-el(3,1)*el(2,2))*temp
			el_inv(2,2)=(el(1,1)*el(3,3)-el(3,1)*el(1,3))*temp
			el_inv(3,2)=(el(3,1)*el(1,2)-el(1,1)*el(3,2))*temp
			el_inv(3,3)=(el(1,1)*el(2,2)-el(1,2)*el(2,1))*temp
			el_inv(1,2)=el_inv(2,1)
        	el_inv(1,3)=el_inv(3,1)
			el_inv(2,3)=el_inv(3,2)
        	!--------------------------------------------------
        	!...compute the RHS
            rhs = 0.0d0
       		do j=1,total_region_elem(k)
                !...global element # of inner element
                iei = k0 + j
          		!...get id of inner element
        		idi = elemid(iei)
        		!...get number of nodes of inner element
        		nni = NODE(idi)
        		!...get nodal coordinates of inner element
        		do m = 1,nni
        			do n =1,2
            			xli(n,m)=node_coor(n,elnode(m,iei))
            		enddo
        		end do
                !...calculate eb and eat
                if ((ELTYPE(idi).eq.CTIP).or.(ELTYPE(idi).eq.CREGULAR)) then
                  	!...crack elements always have traction prescribed
                	call ek_at(k,ieo,iei,xlo,xli,eat,.true.)
                else
                    call ek_b(k,ieo,iei,xlo,xli,eb,.true.)
                    call ek_at(k,ieo,iei,xlo,xli,eat,.true.)
                endif
                !...loop over nodes of outer element (tip element)
                do m = 1,nno
                  	!...loop over directions of outer element
                	do alpha = 1,3
                        !...loop over nodes of inner element
                        do n = 1,nni
                          	!...loop over directions of inner element
                          	do beta = 1,3
                                SELECT CASE (elldid(beta,iei))
                				CASE (BDISP,BTRAC,BTRFREE,NOLOAD)
                                	if ((elem_region(iei).ne.k).and.(elemid(iei).eq.IQUAD)) then
                                    	I_sign = -1
                                    else
                                      	I_sign = 1
                                    endif
                                	rhs(3*(m-1)+alpha) = rhs(3*(m-1)+alpha) + &
                                    	I_sign*eb(3*(m-1)+alpha,3*(n-1)+beta)*sol_trac(beta,n,iei) - &
                                        I_sign*eat(3*(m-1)+alpha,3*(n-1)+beta)*sol_disp(beta,n,iei)
                                CASE (CTRAC,CTRFREE)
                                	rhs(3*(m-1)+alpha) = rhs(3*(m-1)+alpha) - &
                                        eat(3*(m-1)+alpha,3*(n-1)+beta)*sol_disp(beta,n,iei)
                                CASE default
                                	print*,'tstress.f90: current version has not support type of loading: ',elldid(beta,iei)
                                    stop
                                END SELECT
                            enddo
                        enddo
                    enddo
                enddo
            enddo   !...of j over total element of region k
        	!...solve for sigma_u of outer element i
        	sigma_u = 0.0d0
        	!...loop over 3 components of sigma_u
        	do l=1,3
            	do n=1,nno
              		do m = 1,nno
              			sigma_u(l,n) = sigma_u(l,n) + el_inv(n,m)*rhs(3*(m-1)+l)
                	enddo
            	enddo
        	enddo
!        	print'(a45,i5)','Sum of nodal displacements of tip element: ',elem_sys2user(ieo)
!        	print*,'node#     u1     u2     u3'
!        	do n = 1,nno
!          		print'(i5,3(f15.10,2x))',node_sys2user(elnode(n,ieo)),(sigma_u(l,n),l=1,3)
!        	enddo
        	!--------------------------------------------------
        	!...calculate epsilon_s = du_s/ds
        	!...natural coordinate of tip
        	if (ido.eq.CTIP1) then
          		xt = -1.0d0
        	elseif (ido.eq.CTIP2) then
        		xt = 1.0d0
        	endif
        	!...value of shape functions and its derivatives at xt
        	call reg_shape(ido,xt,pso,dpso)
        	!...derivative of position vector
        	dxds = 0.0d0
			do j = 1,nno
	    		do m = 1,2
           			dxds(m) = dxds(m)+xlo(m,j)*dpso(j)
           		enddo
			enddo
        	!...jacobian
        	jacobo = dsqrt(dxds(1)*dxds(1) + dxds(2)*dxds(2))
        	if (jacobo.lt.SMALL_NUM) then
          		print*,'jacobian is almost zero'
            	stop
        	endif
        	!...tangential vector
        	do m = 1,2
          		tangent(m) = dxds(m)/jacobo
        	enddo
        	!...derivative of (global) sum of displ wrt to arc length s
        	dsigma_u = 0.0d0
        	do n = 1,nno
          		do m = 1,3
            		dsigma_u(m) = dsigma_u(m) + sigma_u(m,n)*dpso(n)/jacobo
            	enddo
        	enddo
        	!...component of dsigma_u on tangential direction
        	epsilon_s = 0.0d0
        	do m = 1,2
          		epsilon_s = epsilon_s + dsigma_u(m)*tangent(m)
        	enddo
!        	print*
!        	print'(a50)','d(sum_uk)/ds at the tip: '
!        	print'(3(f15.10,2x))',(dsigma_u(m),m=1,3)
!        	print*
!			print'(a50)','Tangent vector at the tip node:'
!        	print'(2(f15.10,2x))',(tangent(m),m=1,2)
!        	print*
			print*,'Tip element: ',elem_sys2user(ieo)
        	print'(a50)','Dot product of d(sum_uk)/ds with tangent vector:'
        	print'(f15.10)',epsilon_s
    	enddo   !...of i over total elements of region k
        k0 = k0 + total_region_elem(k)
    enddo   !...of k over region
End Subroutine Tstress
!--------------------------------------------------