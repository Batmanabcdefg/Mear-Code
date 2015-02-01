Subroutine Tstress3
!...subroutine to calculate T-stress by 2nd approach: integral on tip element only
!...integrate only on tip element
	USE DefinitionConstant
    USE GlobalData
    USE ElementMatrix
    USE LinearSolver
    IMPLICIT NONE
	!--------------------------------------------------
    !...local variables
	REAL(KIND=DBL)	:: xlo(2,3),xli(2,3),pso(3),dpso(3)
    REAL(KIND=DBL)	:: dxds(2)
    REAL(KIND=DBL)	:: el(3*NDOF,3*NDOF),rhs(3*NDOF),el_vector(3*NDOF*(3*NDOF+1)/2)
    REAL(KIND=DBL)	:: eb(3*NDOF,3*NDOF),eat(3*NDOF,3*NDOF),eb3t(NDOF,3*NDOF),eb3o(NDOF,3*NDOF)
    REAL(KIND=DBL)	:: eat3t(NDOF,3*NDOF),eat3o(NDOF,3*NDOF)
    REAL(KIND=DBL)	:: xt,wt,temp,detj
    INTEGER			:: i,j,k,l,m,n,k0,jj,nn
    INTEGER			:: ido,idi,nno,nni,isrc,IGsrc,jfld,JGfld
    INTEGER			:: alpha,beta
    LOGICAL			:: success
    !--------------------------------------------------
	k0 = 0
    do k = 1,total_region
    	do isrc = 1,total_region_elem(k)
        	!...global element # of outer element
        	IGsrc = isrc + k0
        	!...only continue the procedure if element ie is a TIP element
    		if (eltype(elemid(IGsrc)).ne.CTIP) cycle
            !...get id of the outer element
            ido = elemid(IGsrc)
        	!...get number of nodes of outer element
        	nno = NODE(ido)
        	!...get nodal coordinates of outer element
        	do m = 1,nno
        		do n = 1,2
            		xlo(n,m) = node_coor(n,elnode(m,IGsrc))
            	enddo
        	end do
        	!--------------------------------------------------
        	!...compute matrix el
			el = 0.0d0
            el_vector = 0.0d0
			!...begins the integration loop: use #Gauss-point same as load vector
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
				!...calculate el(3*NDOF,3*NDOF)
        		temp=detj*wt
        		do j = 1,nno
                	jj = NDOF*(j-1) + 1
            		do n = 1,nno
                    	nn = NDOF*(n-1) + 1
                		el(jj,nn) = el(jj,nn) + pso(j)*pso(n)*temp
                	enddo
            	enddo
			enddo !of loop over integration points
            !...fill in other positions of el
            do j = 1,nno
            	jj = NDOF*(j-1) + 1
            	do n = 1,nno
                	nn = NDOF*(n-1) + 1
                    el(jj+1,nn+1) = el(jj,nn)
                    el(jj+2,nn+2) = el(jj,nn)
                enddo
            enddo
            !...put matrix el in vector form
            do i = 1,3*NDOF
              	do j = i,3*NDOF
                	n = j*(j-1)/2 + i
                    el_vector(n) = el(i,j)
                enddo
            enddo
            !--------------------------------------------------
            !...loop over inner element to calculate rhs
            rhs = 0.0d0
            do jfld = 1,total_region_elem(k)
              	JGfld = jfld + k0
                idi = elemid(JGfld)
                nni = NODE(idi)
                do m = 1,nni
                  	do n = 1,2
                    	xli(n,m) = node_coor(n,elnode(m,JGfld))
                    enddo
                enddo
                !--------------------------------------------------
                !...this algorithm still hasn't worked for the case of mixed boundary conditions (eg. symmetric BCs)
            	!...compute ebtip and eg/ehtip to get value of sigma_u at boundaries of the element
                !...depend on where field element is, calculate needed matrices
           		if ((eltype(idi).eq.CTIP).or.(eltype(idi).eq.CREGULAR)) then
                	!...for tip node
                  	call eat3tip(k,IGsrc,JGfld,xlo,xli,eat3t)
                    !...for node on the other side of tip element
                    call eat3oth(k,IGsrc,JGfld,xlo,xli,eat3o)
               	else
                  	!...for tip node
               		call eat3tip(k,IGsrc,JGfld,xlo,xli,eat3t)
                  	call eb3tip(k,IGsrc,JGfld,xlo,xli,eb3t)
                    !...for node on the other side of tip element
                    call eat3oth(k,IGsrc,JGfld,xlo,xli,eat3o)
                  	call eb3oth(k,IGsrc,JGfld,xlo,xli,eb3o)
                endif
                !...calculate eb and eat
                if ((eltype(idi).eq.CTIP).or.(eltype(idi).eq.CREGULAR)) then
                  	call eatts(k,IGsrc,JGfld,xlo,xli,eat)
                else
                  	call eatts(k,IGsrc,JGfld,xlo,xli,eat)
                    call ebts(k,IGsrc,JGfld,xlo,xli,eb)
                endif
                !...loop over nodes of outer element (tip element)
                do m = 1,nno
                  	!...loop over directions of outer element
                	do alpha = 1,NDOF
                        !...loop over nodes of inner element
                        do n = 1,nni
                          	!...loop over directions of inner element
                          	do beta = 1,NDOF
                                SELECT CASE (elldid(beta,JGfld))
                				CASE (BDISP,BTRAC,BTRFREE)
                                	rhs(NDOF*(m-1)+alpha) = rhs(NDOF*(m-1)+alpha) &
                                    	- eb(NDOF*(m-1)+alpha,NDOF*(n-1)+beta)*sol_trac(beta,n,JGfld) &
                                        + eat(NDOF*(m-1)+alpha,NDOF*(n-1)+beta)*sol_disp(beta,n,JGfld)
                                    if (ido.eq.CTIP1) then
                                      	if (m.eq.1) then
                                      		!...add the value of sigma_u at node 1
                                      		rhs(alpha) = rhs(alpha) &
                                        		- eb3t(alpha,NDOF*(n-1)+beta)*sol_trac(beta,n,JGfld) &
                                        		+ eat3t(alpha,NDOF*(n-1)+beta)*sol_disp(beta,n,JGfld)
                                        elseif (m.eq.2) then
                                        	!...for the other boundary node (node 2)
                                        	rhs(NDOF+alpha) = rhs(NDOF+alpha) &
                                        		+ eb3o(alpha,NDOF*(n-1)+beta)*sol_trac(beta,n,JGfld) &
                                        		- eat3o(alpha,NDOF*(n-1)+beta)*sol_disp(beta,n,JGfld)
                                        endif
                                    elseif (ido.eq.CTIP2) then
                                    	if (m.eq.2) then
                                    		!...add the value of sigma_u at node 2
                                      		rhs(NDOF+alpha) = rhs(NDOF+alpha) &
                                        		+ eb3t(alpha,NDOF*(n-1)+beta)*sol_trac(beta,n,JGfld) &
                                        		- eat3t(alpha,NDOF*(n-1)+beta)*sol_disp(beta,n,JGfld)
                                        elseif (m.eq.1) then
                                        	!...for the other boundary node (node 1)
                                        	rhs(alpha) = rhs(alpha) &
                                        		- eb3o(alpha,NDOF*(n-1)+beta)*sol_trac(beta,n,JGfld) &
                                        		+ eat3o(alpha,NDOF*(n-1)+beta)*sol_disp(beta,n,JGfld)
                                        endif
                                    endif
                                CASE (CTRAC,CTRFREE)
                                	rhs(NDOF*(m-1)+alpha) = rhs(NDOF*(m-1)+alpha) + &
                                        eat(NDOF*(m-1)+alpha,NDOF*(n-1)+beta)*sol_disp(beta,n,JGfld)
                                    if (ido.eq.CTIP1) then
                                      	if (m.eq.1) then
                                      		!...for tip node
                                      		rhs(alpha) = rhs(alpha) &
                                        		+ eat3t(alpha,NDOF*(n-1)+beta)*sol_disp(beta,n,JGfld)
                                        elseif (m.eq.2) then
                                        	!...for other boundary node (node 2)
                                        	rhs(NDOF+alpha) = rhs(NDOF+alpha) &
                                        		- eat3o(alpha,NDOF*(n-1)+beta)*sol_disp(beta,n,JGfld)
                                        endif
                                    elseif (ido.eq.CTIP2) then
                                    	if (m.eq.2) then
                                    		!...for tip node (node 2)
                                      		rhs(NDOF+alpha) = rhs(NDOF+alpha) &
                                        		- eat3t(alpha,NDOF*(n-1)+beta)*sol_disp(beta,n,JGfld)
                                        elseif (m.eq.1) then
                                        	!...for other boundary node (node 1)
                                        	rhs(alpha) = rhs(alpha) &
                                        		+ eat3o(alpha,NDOF*(n-1)+beta)*sol_disp(beta,n,JGfld)
                                        endif
                                    endif
                                CASE default
                                	print*,'tstress3.f90: current version has not support type of loading: ',elldid(beta,JGfld)
                                    stop
                                END SELECT
                            enddo
                        enddo
                    enddo
                enddo
            enddo   !...of jfld over total element of region k
            !debug
            !print*,'source element:',IGsrc
            !print*,'rhs='
            !do i = 1,9
            !  print*,rhs(i)
            !enddo
            !print*,'el='
            !do i = 1,9
            !  	print'(9(f15.7,1x))',(el(j,i),j=1,9)
            !enddo
            !...solve for Pi(k) of IGsrc
            success = dPCGSolver(3*NDOF,el_vector,rhs)
            !...print out value of pi(k) at tip
            print*,'Value of Pi:'
            print*,'Node      Pi1      Pi2      Pi3'
           	if (ido.eq.CTIP1) then
               	print'(i5,3(f15.7,2x))',node_sys2user(elnode(1,IGsrc)),(rhs(i),i=1,3)
            elseif (ido.eq.CTIP2) then
            	print'(i5,3(f15.7,2x))',node_sys2user(elnode(2,IGsrc)),(rhs(NDOF+i),i=1,3)
            endif
        enddo !of isrc over total elements of region k
        k0 = k0 + total_region_elem(k)
    enddo !of k
End Subroutine Tstress3
!--------------------------------------------------