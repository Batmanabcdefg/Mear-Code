Subroutine RigidInner3(rotcenter)
!...calculate rigid motion at 2 points ON INNER BOUNDARY to get rigid displ and rotation
	USE DefinitionConstant
    USE GlobalData
    USE ElementMatrix
    USE LinearSolver
    IMPLICIT NONE
    
	!...argument: center of rigid rotation
    !(this is the node where we fix all components to lock rigid translation)
    INTEGER	:: rotcenter
	!--------------------------------------------------
    !...local variables
	REAL(KIND=DBL)	:: xlo(2,3),xli(2,3)
    REAL(KIND=DBL)	:: rhs(NDOF,2),u_rigid(NDOF,2),temp(2)
    REAL(KIND=DBL)	:: eb(NDOF,3*NDOF),eat(NDOF,3*NDOF)
    INTEGER			:: i,k,m,n,k0
    INTEGER			:: ido,idi,nno,nni,isrc,jfld,IGsrc,JGfld
    INTEGER			:: alpha,beta
    REAL(KIND=DBL)	:: omega,u0(2)
    !--------------------------------------------------
    print*,'Results from rigidinner 3'
	k0 = 0
    do k = 1,total_region
    	do isrc = 1,total_region_elem(k)
        	!...global element #
        	IGsrc = isrc + k0
        	!...look for element on inner boundary that have two corner nodes are not SBLP
            !then these 2 nodes are use to calculate rigid displacement
            if ((elemid(IGsrc).eq.INNERBOUND).and.(node_id(elnode(1,IGsrc)).ne.SBLP)&
              								.and.(node_id(elnode(2,IGsrc)).ne.SBLP)) then
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
              	!...compute the rhs for node1 and node2 of IGsrc
                rhs = 0.0d0
                do i = 1,2
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
                    	call eat_rigid(k,IGsrc,JGfld,xlo,xli,i,eat)
                        !...if field element is on ordinary boundary, then need ebb
              			if ((eltype(idi).ne.CTIP).and.(eltype(idi).ne.CREGULAR)) then
                  			call eb_rigid(k,IGsrc,JGfld,xlo,xli,i,eb)
                        endif
                        !...loop over directions of IGsrc
                		do alpha = 1,NDOF
                       		!...loop over nodes of inner element
                       		do n = 1,nni
                       			!...loop over directions of inner element
                       			do beta = 1,NDOF
                             		SELECT CASE (elldid(beta,JGfld))
                					CASE (BDISP,BTRAC,BTRFREE)
                               			rhs(alpha,i) = rhs(alpha,i) + &
                                   				eb(alpha,NDOF*(n-1)+beta)*sol_trac(beta,n,JGfld) - &
                                       			eat(alpha,NDOF*(n-1)+beta)*sol_disp(beta,n,JGfld)
                               		CASE (CTRAC,CTRFREE)
                               			rhs(alpha,i) = rhs(alpha,i) - &
                                       			eat(alpha,NDOF*(n-1)+beta)*sol_disp(beta,n,JGfld)
                               		CASE default
                               			print*,'rigidinner1.f90: current version has not support type of loading: ',elldid(beta,JGfld)
                                   		stop
                               		END SELECT
                           		enddo
                       		enddo
                    	enddo
                    enddo
                enddo
                !...computed rhs is 0.5*(u_rigid+u) -> need to multiply by 2 to get (u_rigid+u)
                rhs = rhs*2.0d0
                print*,'computed displacement'
                print*,'node:',node_sys2user(elnode(1,IGsrc))
                print'(3(f15.5,1x))',(rhs(alpha,1),alpha=1,3)
                print*,'node:',node_sys2user(elnode(2,IGsrc))
                print'(3(f15.5,1x))',(rhs(alpha,2),alpha=1,3)
                !...get rigid displacement
                do i = 1,2
                	do alpha = 1,NDOF
                       	u_rigid(alpha,i) = 0.5d0*(rhs(alpha,i)-sol_disp(alpha,i,IGsrc))
                    enddo
                enddo
                print*,'rigid displacement'
                print*,'node:',node_sys2user(elnode(1,IGsrc))
                print'(3(f15.5,1x))',(u_rigid(alpha,1),alpha=1,3)
                print*,'node:',node_sys2user(elnode(2,IGsrc))
                print'(3(f15.5,1x))',(u_rigid(alpha,2),alpha=1,3)
                !...calculate rigid rotation
                if (dabs(node_coor(2,elnode(1,IGsrc))-node_coor(2,elnode(2,IGsrc))).le.SMALL_NUM) then
                  	!...node 1 and node 2 have same x2-coordinates: using x1 to compute omega
                  	omega = (u_rigid(2,2)-u_rigid(2,1))/(node_coor(1,elnode(2,IGsrc))-node_coor(1,elnode(1,IGsrc)))
                    print*,'case 1, omega=',omega
                elseif (dabs(node_coor(1,elnode(2,IGsrc))-node_coor(1,elnode(1,IGsrc))).le.SMALL_NUM) then
                	!...node 1 and node 2 have same x1-coordinates: using x2 to compute omega
                  	omega = (u_rigid(1,2)-u_rigid(1,1))/(node_coor(2,elnode(1,IGsrc))-node_coor(2,elnode(2,IGsrc)))
                    print*,'case 2, omega=',omega
                else
                  	!...node 1 and node 2 are on incline line, then we can use either x1 or x2 to compute omega
                    omega = (u_rigid(1,2)-u_rigid(1,1))/(node_coor(2,elnode(1,IGsrc))-node_coor(2,elnode(2,IGsrc)))
                    print*,'case 3, omega=',omega
                endif
                !...calculate rigid translation
                u0(1) = u_rigid(1,1) + omega*(node_coor(2,elnode(1,IGsrc))-node_coor(2,rotcenter))
                u0(2) = u_rigid(2,1) - omega*(node_coor(1,elnode(1,IGsrc))-node_coor(1,rotcenter))
                exit !of isrc = 1,total_region_elem(k)
            endif    
        enddo !of isrc = 1,total_region_elem(k)
        
        !...print out real displacement of all nodes on inner boundary
        print*,'Real displacement'
        print*,'Element  Node   u1   u2   u3'
        do isrc = 1,total_region_elem(k)
          	IGsrc = isrc + k0
          	if (elemid(IGsrc).eq.INNERBOUND) then
            	u_rigid = 0.0d0
            	do m = 1,NODE(elemid(IGsrc))
                	temp(1) = u0(1) - omega*(node_coor(2,elnode(m,IGsrc))-node_coor(2,rotcenter))
                    temp(2) = u0(2) + omega*(node_coor(1,elnode(m,IGsrc))-node_coor(1,rotcenter))
            		print'(i5,1x,i5,1x,2(f15.5,1x))',elem_sys2user(IGsrc),node_sys2user(elnode(m,IGsrc)),&
                    (sol_disp(alpha,m,IGsrc)+temp(alpha),alpha=1,2)
                enddo
            endif
        enddo
        k0 = k0 + total_region_elem(k)
    enddo !of k = 1,total_region
End subroutine rigidinner3
