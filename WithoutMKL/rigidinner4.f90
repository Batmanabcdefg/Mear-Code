Subroutine RigidInner4
!...calculate rigid disp at TWO points inside inner domain
	USE DefinitionConstant
    USE GlobalData
    USE ElementMatrix
    IMPLICIT NONE
	!--------------------------------------------------
    !...local variables
	REAL(KIND=DBL)	:: xli(2,3)
    REAL(KIND=DBL)	:: rigid_disp1(NDOF),rigid_disp2(NDOF),rigid_disp(2),u0(2)
    REAL(KIND=DBL)	:: eb(NDOF,3*NDOF),eat(NDOF,3*NDOF)
    INTEGER			:: i,j,k,m,n,k0
    INTEGER			:: idi,nni,isrc,jfld,IGsrc,JGfld
    INTEGER			:: alpha,beta,I_sign
    REAL(KIND=DBL)	:: omega,inside_point(2)
    !...12/8/09: tolerance for calculating of rigid displacements of pure-traction inner bound.
    REAL(KIND=DBL),PARAMETER	:: tol=1.0d-6

    !--------------------------------------------------
    u0 = 0.d0    !rigid translation
    omega = 0.d0 !rigid rotation
    k0 = 0       !total elements for each region
    I_sign = 1   !sign for interface elements
    do k = 1,total_region
      	!...check if this region contains any inner boundary with pure traction
        do i = 1,total_inner_bound
          	if (inner_region(i).eq.k) then
            !...inner boundary i is in region k, need to calculation rigid displacements
        		!...compute rigid displacements (translation + rotation) of point1
                do alpha = 1,2
                  	inside_point(alpha) = inner_inside_point(i,1,alpha)
                enddo
        		rigid_disp1 = 0.0d0
       			!...loop over inner elements
       			do jfld = 1,total_region_elem(k)
       				JGfld = jfld + k0
       				idi = elemid(JGfld)
       				nni = NODE(idi)
       				do m = 1,nni
       					do n = 1,2
       						xli(n,m) = node_coor(n,elnode(m,JGfld))
       					enddo
      				enddo
           			!...all field elements need eatt
           			call eat_rigid2(k,JGfld,inside_point,xli,eat)
            		!...if field element is NOT on crack, then need ebb
            		if ((eltype(idi).ne.CTIP).and.(eltype(idi).ne.CREGULAR)) then
            			call eb_rigid2(k,JGfld,inside_point,xli,eb)
            		endif
            		!...loop over directions of IGsrc
            		do alpha = 1,NDOF
            			!...loop over nodes of inner element
                		do n = 1,nni
                			!...loop over directions of inner element
                    		do beta = 1,NDOF
                    			SELECT CASE (elldid(beta,JGfld))
                				CASE (BDISP,BTRAC,BTRFREE,NOLOAD)
                                	!...if field element is on interface and not a "master" element
                                    !then normal vector is switched direction
                                	if ((elem_region(JGfld).ne.k).and.(elemid(JGfld).eq.IQUAD)) then
                                    	I_sign = -1
                                    else
                                      	I_sign = 1
                                    endif
                        			rigid_disp1(alpha) = rigid_disp1(alpha) + &
                                    	I_sign*eb(alpha,NDOF*(n-1)+beta)*sol_trac(beta,n,JGfld) - &
                                        I_sign*eat(alpha,NDOF*(n-1)+beta)*sol_disp(beta,n,JGfld)
                        		CASE (CTRAC,CTRFREE)
                            		rigid_disp1(alpha) = rigid_disp1(alpha) - &
                                        eat(alpha,NDOF*(n-1)+beta)*sol_disp(beta,n,JGfld)
                        		CASE default
                        			print*,'rigidinner4.f90: current version has not support type of loading: ',elldid(beta,JGfld)
                            		stop
                        		END SELECT
                    		enddo
                		enddo
            		enddo
        		enddo
        		!...print out computed displacement
        		print*,'Coordinates of node where rigid displacement was calculated:'
        		print'(a3,f10.2)','x= ',inner_inside_point(i,1,1)
        		print'(a3,f10.2)','y= ',inner_inside_point(i,1,2)
        		print*,'Rigid displacement:'
        		do alpha = 1,NDOF
        			print'(a1,i1,a1,f15.5)','u',alpha,'=',rigid_disp1(alpha)
        		enddo
        		!--------------------------------------------------
				!...compute rigid displacements (translation + rotation) of point2
                do alpha = 1,2
                  	inside_point(alpha) = inner_inside_point(i,2,alpha)
                enddo
        		rigid_disp2 = 0.0d0
       			!...loop over inner elements
       			do jfld = 1,total_region_elem(k)
       				JGfld = jfld + k0
       				idi = elemid(JGfld)
       				nni = NODE(idi)
       				do m = 1,nni
       					do n = 1,2
       						xli(n,m) = node_coor(n,elnode(m,JGfld))
       					enddo
      				enddo
           			!...all field elements need eatt
           			call eat_rigid2(k,JGfld,inside_point,xli,eat)
            		!...if field element is on ordinary boundary, then need ebb
            		if ((eltype(idi).ne.CTIP).and.(eltype(idi).ne.CREGULAR)) then
            			call eb_rigid2(k,JGfld,inside_point,xli,eb)
            		endif
            		!...loop over directions of IGsrc
            		do alpha = 1,NDOF
            			!...loop over nodes of inner element
                		do n = 1,nni
                			!...loop over directions of inner element
                    		do beta = 1,NDOF
                    			SELECT CASE (elldid(beta,JGfld))
                				CASE (BDISP,BTRAC,BTRFREE,NOLOAD)
                                	!...if field element is on interface and not a "master" element
                                    !then normal vector is switched direction
                                	if ((elem_region(JGfld).ne.k).and.(elemid(JGfld).eq.IQUAD)) then
                                    	I_sign = -1
                                    else
                                      	I_sign = 1
                                    endif
                        			rigid_disp2(alpha) = rigid_disp2(alpha) + &
                                    	I_sign*eb(alpha,NDOF*(n-1)+beta)*sol_trac(beta,n,JGfld) - &
                                        I_sign*eat(alpha,NDOF*(n-1)+beta)*sol_disp(beta,n,JGfld)
                        		CASE (CTRAC,CTRFREE)
                            		rigid_disp2(alpha) = rigid_disp2(alpha) - &
                                        eat(alpha,NDOF*(n-1)+beta)*sol_disp(beta,n,JGfld)
                        		CASE default
                        			print*,'rigidinner4.f90: current version has not support type of loading: ',elldid(beta,JGfld)
                            		stop
                        		END SELECT
                    		enddo
                		enddo
            		enddo
        		enddo
        		!...print out rigid displacement
        		print*,'Coordinates of node where rigid displacement was calculated:'
        		print'(a3,f10.2)','x= ',inner_inside_point(i,2,1)
        		print'(a3,f10.2)','y= ',inner_inside_point(i,2,2)
        		print*,'Rigid displacement:'
        		do alpha = 1,NDOF
        			print'(a1,i1,a1,f15.5)','u',alpha,'=',rigid_disp2(alpha)
        		enddo
        		!...calculate rigid rotation
        		if (dabs(inner_inside_point(i,2,2)-inner_inside_point(i,1,2)).le.tol) then
          			!...point 1 and point 2 have same x2-coordinates: using x1 to compute omega
          			omega = (rigid_disp2(2)-rigid_disp1(2))/(inner_inside_point(i,2,1)-inner_inside_point(i,1,1))
        		elseif (dabs(inner_inside_point(i,2,1)-inner_inside_point(i,1,1)).le.tol) then
          			!...point 1 and point 2 have same x1-coordinates: using x2 to compute omega
          			omega = (rigid_disp2(1)-rigid_disp1(1))/(inner_inside_point(i,1,2)-inner_inside_point(i,2,2))
        		else
          			!...point1 and point2 are on incline line, then we can use either x1 or x2 to compute omega
            		omega = (rigid_disp2(1)-rigid_disp1(1))/(inner_inside_point(i,1,2)-inner_inside_point(i,2,2))
        		endif
                print'(a24,f8.3)','Rigid rotation omega = ',omega
        		!...calculate rigid translation
        		u0(1) = rigid_disp1(1) + omega*(inner_inside_point(i,1,2)-node_coor(2,inner_node(i,1)))
        		u0(2) = rigid_disp1(2) - omega*(inner_inside_point(i,1,1)-node_coor(1,inner_node(i,1)))
        		print*,'Rigid translation:'
        		print'(a5,f8.3)','u1 = ',u0(1)
        		print'(a5,f8.3)','u2 = ',u0(2)
        		!...print out real displacement of all nodes on inner boundary
        		print'(a36,i2)','Real Displacement of Inner Boundary ',i
        		print*,'    Element  Node    Displacement-X   Displacement-Y   Displacement-Z'
        		do isrc = 1,inner_total_elem(i)
                	!...get USER # of element on inner boundary i
                    j = inner_elem(i,isrc)
                    !...loop over all elements to get corresponding system# of element j
                    IGsrc = 0
                    do m = 1,total_region_elem(k)
                      	!...global system#
                      	n = m + k0
                        if (j.eq.elem_sys2user(n)) then
                          	IGsrc = n
                            exit
                        endif
                    enddo
                    if (IGsrc.eq.0) then
                      	print*,'element ',j,'of inner boundary ',i,'is not in correct region'
                        stop
                    endif
            		rigid_disp = 0.0d0
            		do m = 1,NODE(elemid(IGsrc))
                		rigid_disp(1) = u0(1) - omega*(node_coor(2,elnode(m,IGsrc))-node_coor(2,inner_node(i,1)))
                    	rigid_disp(2) = u0(2) + omega*(node_coor(1,elnode(m,IGsrc))-node_coor(1,inner_node(i,1)))
                        !...re-compute real displacements
                        do alpha = 1,2
                          	sol_disp(alpha,m,IGsrc) = sol_disp(alpha,m,IGsrc) + rigid_disp(alpha)
                        enddo
            			!print'(i5,1x,i5,1x,2(f15.5,1x))',elem_sys2user(IGsrc),node_sys2user(elnode(m,IGsrc)),&
                    	!	(sol_disp(alpha,m,IGsrc)+rigid_disp(alpha),alpha=1,2)
                        if (NDOF.eq.1) then
                        	PRINT "(1x,i7,2x,i7,6x,e13.6,3x)",                           &
                   			elem_sys2user(IGsrc), node_sys2user(elnode(m,IGsrc)),                   &
                   			(sol_disp(j,m,IGsrc), j = 1,NDOF)
                    	elseif (NDOF.eq.2) then
                        	PRINT "(1x,i7,2x,i7,6x,2(e13.6,3x))",                           &
                   			elem_sys2user(IGsrc), node_sys2user(elnode(m,IGsrc)),                   &
                   			(sol_disp(j,m,IGsrc), j = 1,NDOF)
                    	elseif (NDOF.eq.3) then
              				PRINT "(1x,i7,2x,i7,6x,3(e13.6,3x))",                           &
                   			elem_sys2user(IGsrc), node_sys2user(elnode(m,IGsrc)),                   &
                   			(sol_disp(j,m,IGsrc), j = 1,NDOF)
                        elseif (NDOF.eq.4) then
                        	PRINT "(1x,i7,2x,i7,6x,4(e13.6,3x))",                           &
                   			elem_sys2user(IGsrc), node_sys2user(elnode(m,IGsrc)),                   &
                   			(sol_disp(j,m,IGsrc), j = 1,NDOF)
                        elseif (NDOF.eq.5) then
                        	PRINT "(1x,i7,2x,i7,6x,5(e13.6,3x))",                           &
                   			elem_sys2user(IGsrc), node_sys2user(elnode(m,IGsrc)),                   &
                   			(sol_disp(j,m,IGsrc), j = 1,NDOF)
                        else
                          	print*,'rigidinner4: not print displ for N=',NDOF
                            stop
                        endif
                	enddo
                    print*
        		enddo
            endif !of inner_region(i).eq.k
        enddo !...of i=1,total_inner_bound
        
        k0 = k0 + total_region_elem(k) !...this is needed for JGfld
    enddo !of k = 1,total_region
End subroutine RigidInner4
