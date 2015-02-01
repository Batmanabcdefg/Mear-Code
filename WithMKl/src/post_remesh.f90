SUBROUTINE Post
!...5/19/09: this subroutine is renamed from post.f90 after adding the subroutine CrackGrowth
!and deleting the temporary subroutine GrowthAngle (The rest are the same for both post.f90 and post_remesh.f90)
!...5/24/09: add the bilinear law of propagation (compute all growth data at once, then propagate each tip)
!...1/5/10: rotate material constants to geometry coordinates
	USE DefinitionConstant
	USE GlobalData
    !...4/29/09
    USE Remesh
	IMPLICIT NONE
	! local variables
	REAL(KIND=DBL)   	:: pre_ndisp(NDOF,total_node)
	REAL(KIND=DBL)   	:: ey1,pv1,mu,E(3,NDOF,NDOF,3),xele(2,3), uk_pressure_temp(NDOF),SIF_pressure_temp(NDOF)
        REAL(KIND=DBL)	    :: normal_i(2,3)!, tip_data
        REAL(KIND=DBL)      :: min_scaling(2) !,location of minumun pressure scaling in crack scaling matrix
	INTEGER             :: i,j,k,kk,l,p,q,m,n,ie,alpha,inode,ierr, ck, pr
    !...5/12/09: add variables of contracted moduli for easy control the material matrx
    REAL(KIND=DBL)	:: C11,C12,C13,C14,C15,C16,&
                           C22,C23,C24,C25,C26,C33,C34,C35,C36,&
                           C44,C45,C46,C55,C56,C66
	REAL(KIND=DBL)	:: e11,e12,e13,e14,e15,e16,&
    				   e21,e22,e23,e24,e25,e26,&
                       e31,e32,e33,e34,e35,e36
	REAL(KIND=DBL)	:: k11,k12,k13,k22,k23,k33
    REAL(KIND=DBL)	:: h16,&
    				   h21,h22
	REAL(KIND=DBL)	:: b11,b22
    REAL(KIND=DBL)	:: g11,g22
    INTEGER			:: matl_no,matl_type
    REAL(KIND=DBL)	:: E_temp(3,NDOF,NDOF,3),dir_cosine(NDOF,NDOF)
    

	!------------------------------------------------------------
    !...deallocate if growth is simulated
    if ((growth_flag.eq.1).and.(nstep.ge.2)) then
    	deallocate(sol_trac,sol_trac_traction,sol_trac_pressure,STAT=ierr)
    	if (ierr.ne.0) then
      		print*,'sol_trac...: deallocation request denied!'
        	stop
    	endif
	deallocate(sol_disp,opening_disp,opening_disp_traction, &
                   opening_disp_pressure,sol_disp_traction, sol_disp_pressure,STAT=ierr)
    	if (ierr.ne.0) then
      		print*,'sol_disp...: deallocation request denied!'
        	stop
    	endif
	deallocate(uk_traction,uk_pressure,STAT=ierr)
    	if (ierr.ne.0) then
      		print*,'uk...: deallocation request denied!'
        	stop
    	endif
	deallocate(SIF_traction,SIF_pressure,&
                   nodal_SIF,nodal_SIF_total,nodal_SIF_traction,nodal_SIF_pressure,SIF_effective,STAT=ierr)
    	if (ierr.ne.0) then
      		print*,'SIF...: deallocation request denied!'
        	stop
    	endif
	deallocate(volume_ratio,arc_length,crack_volume_total,crack_volume_tractions,crack_volume_pressure,STAT=ierr)
    	if (ierr.ne.0) then
      		print*,'volume and arc length...: deallocation request denied!'
        	stop
    	endif
        	deallocate(crack_scaling, crack_scaling_temp,tipnode_sys2user,STAT=ierr)
    	if (ierr.ne.0) then
      		print*,'crack_scaling...: deallocation request denied!'
        	stop
    	endif



     endif
    !------------------------------------------------------------
	!...allocate global variables to store nodal data (prescribed/result)
    !modified from local variables (4/6/09)
    allocate(sol_trac(NDOF,3,total_elem),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'sol_trac: allocation request denied'
        stop
    endif
    allocate(sol_trac_traction(NDOF,3,total_elem),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'sol_trac_traction: allocation request denied'
        stop
    endif
    allocate(sol_trac_pressure(total_cracks,NDOF,3,total_elem),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'sol_trac_pressure: allocation request denied'
        stop
    endif
    allocate(sol_disp(NDOF,3,total_elem),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'sol_disp: allocation request denied'
        stop
    endif
    allocate(opening_disp(total_elem,3),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'opening_disp: allocation request denied'
        stop
    endif
    allocate(opening_disp_traction(total_elem,3),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'opening_disp_traction: allocation request denied'
        stop
    endif
    allocate(opening_disp_pressure(total_cracks,total_elem,3),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'opening_disp_pressure: allocation request denied'
        stop
    endif
    allocate(sol_disp_traction(NDOF,3,total_elem),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'sol_disp_traction: allocation request denied'
        stop
    endif
    allocate(sol_disp_pressure(total_cracks,NDOF,3,total_elem),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'sol_disp_pressure: allocation request denied'
        stop
    endif
    allocate(arc_length(total_cracks),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'arc_length: allocation request denied'
        stop
    endif
    allocate(volume_ratio(total_cracks),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'volume_ratio: allocation request denied'
        stop
    endif
    allocate(crack_volume_total(total_cracks),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'crack_volume_total: allocation request denied'
        stop
    endif
    allocate(crack_volume_tractions(total_cracks),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'crack_volume_tractions: allocation request denied'
        stop
    endif
    allocate(crack_volume_pressure(total_cracks,total_cracks),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'crack_volume_pressure: allocation request denied'
        stop
    endif
    !...5/10/09: allocate global variables to store SIF data
    total_tip_node = 0   !total number of tip nodes of ALL regions
    do i = 1,total_elem
      	if ((elemid(i).eq.CTIP1).or.(elemid(i).eq.CTIP2)) then
        	total_tip_node = total_tip_node + 1
        endif
    enddo
    allocate(nodal_SIF(NDOF,total_tip_node),tipnode_sys2user(total_tip_node),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'nodal_SIF and tipnode_sys2user: allocation request denied'
        stop
    endif
     allocate(nodal_SIF_total(NDOF,total_tip_node),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'nodal_SIF_total: allocation request denied'
        stop
    endif
    allocate(nodal_SIF_traction(NDOF,total_tip_node),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'nodal_SIF_traction: allocation request denied'
        stop
    endif
    allocate(nodal_SIF_pressure(total_cracks,NDOF,total_tip_node),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'nodal_SIF_pressure: allocation request denied'
        stop
     endif
    allocate(SIF_effective(total_tip_node),STAT=ierr)
    if (ierr .ne. 0) then
       print*,'SIF_effective: allocation request denied'
       stop
    endif
    allocate(uk_traction(NDOF),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'uk_traction: allocation request denied'
        stop
     endif    
     allocate(uk_pressure(total_cracks,NDOF),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'uk_pressure: allocation request denied'
        stop
     endif
    allocate(SIF_traction(NDOF),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'SIF_traction: allocation request denied'
        stop
     endif
    allocate(SIF_pressure(total_cracks,NDOF),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'SIF_pressure: allocation request denied'
        stop
     endif
    allocate(crack_scaling(total_cracks),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'crack_scaling: allocation request denied'
        stop
     endif

    allocate(crack_scaling_temp(total_cracks,total_tip_node),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'crack_scaling_temp: allocation request denied'
        stop
     endif




    !------------------------------------------------------------
	! Obtain the prescribed nodal value
	DO i = 1, total_elem   
    	DO inode = 1, NODE(elemid(i))
      		DO alpha = 1,NDOF
        		IF (unknown(alpha,elnode(inode,i)) == TRACTION) THEN
          			pre_ndisp(alpha,elnode(inode,i)) = elnode_val(alpha,inode,i)
        		END IF
      		END DO
    	END DO
  	END DO
	!------------------------------------------------------------
  	! Get tractions and displacements
  	ie = 0
  	DO k = 1, total_region
    	!...7/20/09: deactivate this unnecessary portion
    	!IF (material_type(region_mat_no(k))==1) THEN 
       	!	ey1  = material(1, region_mat_no(k))
       	!	pv1 = material(2, region_mat_no(k))
       	!	mu = 0.5d0 * ey1 / (1.d0 + pv1)  
		!ELSE
	   	!	DO kk=1,21
	    !		Matprop(kk)=material(kk, region_mat_no(k))
	   	!	END DO
		!END IF
        
    	DO i = 1, total_region_elem(k)
      		ie = ie + 1  ! global element no.
      		DO alpha = 1,NDOF
        		SELECT CASE (elldid(alpha,ie))
                	!----------
          			CASE (BTRAC,BTRFREE,CTRAC,CTRFREE)
          			! known traction
          			DO inode = 1, NODE(elemid(ie))
            			!sol_trac(alpha,inode,ie) = elnode_val(alpha,inode,ie)
                                sol_trac_traction(alpha,inode,ie) = elnode_val_trac(alpha,inode,ie)
                                
                                do ck = 1,total_cracks
                                   sol_trac_pressure(ck,alpha,inode,ie) = elnode_val_pressure(ck,alpha,inode,ie)
                                enddo
                                
            			IF (unknown(alpha,elnode(inode,ie)) == TRACTION) THEN
              				sol_disp(alpha,inode,ie) = pre_ndisp(alpha,elnode(inode,ie))
            			ELSE
                        	! unknown is DISPLACEMENT
              				if ((node_id(elnode(inode,ie)).eq.SBLP).and.(ELTYPE(elemid(ie)).ne.BREGULAR)) then
                            	!...SBL node of crack element (associated with upper node) -> need to get delta_u
                                !print*,'correct SBLB node: ',node_sys2user(elnode(inode,ie)),'component: ',alpha
                                !sol_disp(alpha,inode,ie) = dF(equno(inode,ie)+alpha-1) - dF(equno(inode,ie)+NDOF+alpha-1)
                                sol_disp_traction(alpha,inode,ie) = dF_traction(equno(inode,ie)+alpha-1) - dF_traction(equno(inode,ie)+NDOF+alpha-1)

                                do ck = 1,total_cracks
                                   sol_disp_pressure(ck,alpha,inode,ie) = dF_pressure(ck,equno(inode,ie)+alpha-1) - dF_pressure(ck,equno(inode,ie)+NDOF+alpha-1)
                                enddo

                            else
              					!sol_disp(alpha,inode,ie) = dF(equno(inode,ie)+alpha-1)
              					sol_disp_traction(alpha,inode,ie) = dF_traction(equno(inode,ie)+alpha-1)
                                                do ck = 1,total_cracks
                                                   sol_disp_pressure(ck,alpha,inode,ie) = dF_pressure(ck,equno(inode,ie)+alpha-1)
              					enddo

                            endif
            			END IF
          			END DO
					!----------
          			CASE (BDISP)
          			! known displacement
          			DO inode = 1, NODE(elemid(ie))
            			sol_disp(alpha,inode,ie) = elnode_val(alpha,inode,ie)
            			! restore the solved traction
            			sol_trac(alpha,inode,ie) = dF(equno(inode,ie)+alpha-1)
          			END DO
					!----------
          			CASE (NOLOAD)
          			! interface element, both TRACTION and DISPLACEMENT are unknown
          			DO inode = 1, NODE(elemid(ie))
            			sol_disp(alpha, inode, ie) = dF(equno(inode,ie)+alpha-1)
            			sol_trac(alpha, inode, ie) = dF(equno(inode,ie)+alpha+NDOF-1)
          			END DO
        		END SELECT
      		END DO ! loop alpha
    	END DO ! loop i
	END DO ! loop k(region)


	!------------------------------------------------------------
    !...5/10/09: put the calculation of SIFs in this subroutine
    ie = 0
    tip_node_count = 0
    DO k = 1, total_region
      	matl_no = region_mat_no(k)
        matl_type = material_type(matl_no)
        E = 0.d0
      	if (matl_type.eq.1) then
        !...isotropy
        	if ((NDOF.eq.3).or.(NDOF.eq.2)) then
       			ey1 = material(1, region_mat_no(k))
       			pv1 = material(2, region_mat_no(k))
       			mu = 0.5d0 * ey1 / (1.0d0 + pv1)
            else
              	print*,'post.f90: not yet isotropy for media:',NDOF
                stop
            endif
		elseif (matl_type.eq.2) then
        !...cubic material
	   		if (NDOF.eq.3) then
            	C11 = material(1,matl_no)
                C12 = material(2,matl_no)
                C44 = material(3,matl_no)
	     		E(1,1,1,1) = C11
		 		E(1,1,2,2) = C12
		 		E(1,1,3,3) = C12
	     		E(2,2,1,1) = C12
		 		E(2,2,2,2) = C11
	     		E(2,2,3,3) = C12
		 		E(3,3,1,1) = C12
		 		E(3,3,2,2) = C12
	     		E(3,3,3,3) = C11
	     		E(1,2,1,2) = C44
	     		E(1,2,2,1) = C44
	     		E(2,1,1,2) = C44
	     		E(2,1,2,1) = C44
	     		E(1,3,1,3) = C44
	     		E(1,3,3,1) = C44
	     		E(3,1,1,3) = C44
	     		E(3,1,3,1) = C44
	     		E(2,3,2,3) = C44
	     		E(2,3,3,2) = C44
	     		E(3,2,2,3) = C44
	     		E(3,2,3,2) = C44
            else
              	print*,'not yet cubic for other media'
            	stop
        	endif
		elseif (matl_type.eq.3) then
		!...Transversely isotropy
        	if (NDOF.eq.2) then
            !...plane strain of elastic media
            	!...2 is elastic symmetry: order of input C11,C12,C13,C22,C44
				C11 = material(1,matl_no)
        		C12 = material(2,matl_no)
            	!...Note: C13 is not need for E_ijkl but it is needed for calculating E_3113
            	C13 = material(3,matl_no)
        		C22 = material(4,matl_no)
        		C44 = material(5,matl_no)
				E(1,1,1,1) = C11
        		E(1,1,2,2) = C12
        		E(2,2,1,1) = C12
        		E(2,2,2,2) = C22
        		E(3,2,2,3) = C44
        		E(3,1,1,3) = 0.5d0*(C11-C13)
        		E(1,2,1,2) = C44
        		E(1,2,2,1) = C44
        		E(2,1,1,2) = C44
        		E(2,1,2,1) = C44
    		elseif (NDOF.eq.3) then
        	!...elastic media
        	!-----------------------------------------------------------------------------------
			!This code for the case of 3 is elastic symmetry: order of input C11,C12,C13,C33,C55
			!C11 = material(1,matl_no)
        	!C12 = material(2,matl_no)
        	!C13 = material(3,matl_no)
        	!C33 = material(4,matl_no)
        	!C55 = material(5,matl_no)
			!E(1,1,1,1)=C11
			!E(1,1,2,2)=C12
			!E(1,1,3,3)=C13
			!E(2,2,1,1)=C12
			!E(2,2,2,2)=C11
			!E(2,2,3,3)=C13
			!E(3,3,1,1)=C13
			!E(3,3,2,2)=C13
			!E(3,3,3,3)=C33
			!E(1,2,1,2)=0.50d00*(C11-C12)
			!E(1,2,2,1)=E(1,2,1,2)
			!E(2,1,1,2)=E(1,2,1,2)
			!E(2,1,2,1)=E(1,2,1,2)
			!E(1,3,1,3)=C55
			!E(1,3,3,1)=C55
			!E(3,1,1,3)=C55
			!E(3,1,3,1)=C55
			!E(2,3,2,3)=C55
			!E(2,3,3,2)=C55
			!E(3,2,2,3)=C55
			!E(3,2,3,2)=C55
			!-----------------------------------------------------------------------------------
				!...2 is elastic symmetry: order of input C11,C12,C13,C22,C44
				C11 = material(1,matl_no)
        		C12 = material(2,matl_no)
        		C13 = material(3,matl_no)
        		C22 = material(4,matl_no)
        		C44 = material(5,matl_no)
				E(1,1,1,1) = C11
        		E(1,1,2,2) = C12
        		E(1,1,3,3) = C13
        		E(2,2,1,1) = C12
        		E(2,2,2,2) = C22
        		E(2,2,3,3) = C12
        		E(3,3,1,1) = C13
        		E(3,3,2,2) = C12
        		E(3,3,3,3) = C11
        		E(2,3,2,3) = C44
        		E(2,3,3,2) = C44
        		E(3,2,2,3) = C44
        		E(3,2,3,2) = C44
        		E(3,1,3,1) = 0.5d0*(C11-C13)
        		E(3,1,1,3) = E(3,1,3,1)
        		E(1,3,3,1) = E(3,1,3,1)
        		E(1,3,1,3) = E(3,1,3,1)
        		E(1,2,1,2) = C44
        		E(1,2,2,1) = C44
        		E(2,1,1,2) = C44
        		E(2,1,2,1) = C44
        	elseif (NDOF.eq.4) then
        	!...piezoelectric media: 2 is elastic symmetry
        		C11 = material(1,matl_no)
        		C12 = material(2,matl_no)
        		C13 = material(3,matl_no)
        		C22 = material(4,matl_no)
        		C44 = material(5,matl_no)
				E(1,1,1,1) = C11
        		E(1,1,2,2) = C12
        		E(1,1,3,3) = C13
        		E(2,2,1,1) = C12
        		E(2,2,2,2) = C22
        		E(2,2,3,3) = C12
        		E(3,3,1,1) = C13
        		E(3,3,2,2) = C12
        		E(3,3,3,3) = C11
        		E(2,3,2,3) = C44
        		E(2,3,3,2) = C44
        		E(3,2,2,3) = C44
        		E(3,2,3,2) = C44
        		E(3,1,3,1) = 0.5d0*(C11-C13)
        		E(3,1,1,3) = E(3,1,3,1)
        		E(1,3,3,1) = E(3,1,3,1)
        		E(1,3,1,3) = E(3,1,3,1)
        		E(1,2,1,2) = C44
        		E(1,2,2,1) = C44
        		E(2,1,1,2) = C44
        		E(2,1,2,1) = C44
            	e16 = material(6,matl_no)
            	e21 = material(7,matl_no)
            	e22 = material(8,matl_no)
            	k11 = material(9,matl_no)
            	k22 = material(10,matl_no)
            	E(1,2,4,1) = e16
            	E(2,1,4,1) = e16
            	E(1,4,2,1) = e16
            	E(1,4,1,2) = e16
            	E(1,1,4,2) = e21
            	E(2,4,1,1) = e21
            	E(2,2,4,2) = e22
            	E(2,4,2,2) = e22
            	E(3,3,4,2) = e21
            	E(2,4,3,3) = e21
            	E(2,3,4,3) = e16
            	E(3,2,4,3) = e16
            	E(3,4,3,2) = e16
            	E(3,4,2,3) = e16
            	E(1,4,4,1) = -k11
            	E(2,4,4,2) = -k22
            	E(3,4,4,3) = -k11
            elseif (NDOF.eq.5) then
        		!...this code will be for the case of 2 is elastic symmetry
            	C11 = material(1,matl_no)
        		C12 = material(2,matl_no)
        		C13 = material(3,matl_no)
        		C22 = material(4,matl_no)
        		C44 = material(5,matl_no)
            	e16 = material(6,matl_no)
            	e21 = material(7,matl_no)
            	e22 = material(8,matl_no)
            	k11 = material(9,matl_no)
            	k22 = material(10,matl_no)
            	h16 = material(11,matl_no)
            	h21 = material(12,matl_no)
            	h22 = material(13,matl_no)
            	b11 = material(14,matl_no)
            	b22 = material(15,matl_no)
            	g11 = material(16,matl_no)
            	g22 = material(17,matl_no)
            
				E(1,1,1,1) = C11
        		E(1,1,2,2) = C12
        		E(1,1,3,3) = C13
        		E(2,2,1,1) = C12
        		E(2,2,2,2) = C22
        		E(2,2,3,3) = C12
        		E(3,3,1,1) = C13
        		E(3,3,2,2) = C12
        		E(3,3,3,3) = C11
        		E(2,3,2,3) = C44
        		E(2,3,3,2) = C44
        		E(3,2,2,3) = C44
        		E(3,2,3,2) = C44
        		E(3,1,3,1) = 0.5d0*(C11-C13)
        		E(3,1,1,3) = E(3,1,3,1)
        		E(1,3,3,1) = 0.5d0*(C11-C13)
        		E(1,3,1,3) = E(3,1,3,1)
        		E(1,2,1,2) = C44
        		E(1,2,2,1) = C44
        		E(2,1,1,2) = C44
        		E(2,1,2,1) = C44
            
            	E(1,2,4,1) = e16
            	E(2,1,4,1) = e16
            	E(1,4,2,1) = e16
            	E(1,4,1,2) = e16
            	E(1,1,4,2) = e21
            	E(2,4,1,1) = e21
            	E(2,2,4,2) = e22
            	E(2,4,2,2) = e22
            	E(3,3,4,2) = e21
            	E(2,4,3,3) = e21
            	E(2,3,4,3) = e16
            	E(3,2,4,3) = e16
            	E(3,4,3,2) = e16
            	E(3,4,2,3) = e16
            
            	E(1,4,4,1) = -k11
            	E(2,4,4,2) = -k22
            	E(3,4,4,3) = -k11

            	E(1,2,5,1) = h16
            	E(2,1,5,1) = h16
            	E(1,5,2,1) = h16
            	E(1,5,1,2) = h16
            	E(1,1,5,2) = h21
            	E(2,5,1,1) = h21
            	E(2,2,5,2) = h22
            	E(2,5,2,2) = h22
            	E(3,3,5,2) = h21
            	E(2,5,3,3) = h21
            	E(2,3,5,3) = h16
            	E(3,2,5,3) = h16
            	E(3,5,3,2) = h16
            	E(3,5,2,3) = h16

            	E(1,4,5,1) = -b11
            	E(1,5,4,1) = -b11
            	E(2,4,5,2) = -b22
            	E(2,5,4,2) = -b22
            	E(3,4,5,3) = -b11
            	E(3,5,4,3) = -b11
            
            	E(1,5,5,1) = -g11
            	E(2,5,5,2) = -g22
            	E(3,5,5,3) = -g11
        	else
        		print*,'not coding yet forvtransversely isotropy of media:',NDOF
        		stop
        	endif
		elseif (matl_type.eq.4) then
    	!...orthotropy
        	if (NDOF.eq.2) then
        	!...plane strain of elastic media
        		C11 = material(1,matl_no)
        		C12 = material(2,matl_no)
        		C22 = material(3,matl_no)
        		C44 = material(4,matl_no)
        		C55 = material(5,matl_no)
       			C66 = material(6,matl_no)
				E(1,1,1,1)=C11
	    		E(1,1,2,2)=C12
	    		E(2,2,1,1)=C12
	    		E(2,2,2,2)=C22
                E(3,2,2,3)=C44
                E(3,1,1,3)=C55
	    		E(1,2,1,2)=C66
	    		E(1,2,2,1)=C66
	    		E(2,1,1,2)=C66
	    		E(2,1,2,1)=C66
        	elseif (NDOF.eq.3) then
			!...elastic media, order of input C11,C12,C13,C22,C23,C33,C44,C55,C66
        	!...6/26/09: modify 23->4, 31->5, 12->6
    			C11 = material(1,matl_no)
        		C12 = material(2,matl_no)
        		C13 = material(3,matl_no)
        		C22 = material(4,matl_no)
        		C23 = material(5,matl_no)
        		C33 = material(6,matl_no)
        		C44 = material(7,matl_no)
        		C55 = material(8,matl_no)
       			C66 = material(9,matl_no)
				E(1,1,1,1)=C11
	    		E(1,1,2,2)=C12
	    		E(1,1,3,3)=C13
	    		E(2,2,1,1)=C12
	    		E(2,2,2,2)=C22
	    		E(2,2,3,3)=C23
	    		E(3,3,1,1)=C13
	    		E(3,3,2,2)=C23
	    		E(3,3,3,3)=C33
                E(2,3,2,3)=C44
	    		E(2,3,3,2)=C44
	    		E(3,2,2,3)=C44
	    		E(3,2,3,2)=C44
                E(1,3,1,3)=C55
	    		E(1,3,3,1)=C55
	    		E(3,1,1,3)=C55
	    		E(3,1,3,1)=C55
	    		E(1,2,1,2)=C66
	    		E(1,2,2,1)=C66
	    		E(2,1,1,2)=C66
	    		E(2,1,2,1)=C66
            !...1/21/10 add orthotropy for NDOF=4, only works for the case of polling direction is 3-axis
        	!(this part is developed to compare with Denda's paper, not sure in general case)
        	elseif (NDOF.eq.4) then
        		!...elastic moduli
            	C11 = material(1,matl_no)
        		C12 = material(2,matl_no)
            	!...No need to have C13, C23 and C33 (need to modify prep.f90 also)
	       		C13 = material(3,matl_no)
        		C22 = material(4,matl_no)
        		C23 = material(5,matl_no)
	       		C33 = material(6,matl_no)
        		C44 = material(7,matl_no)
        		C55 = material(8,matl_no)
       			C66 = material(9,matl_no)
				E(1,1,1,1)=C11
	    		E(1,1,2,2)=C12
	    		E(1,1,3,3)=C13
	    		E(2,2,1,1)=C12
	    		E(2,2,2,2)=C22
	    		E(2,2,3,3)=C23
	    		E(3,3,1,1)=C13
	    		E(3,3,2,2)=C23
	    		E(3,3,3,3)=C33
	    		E(2,3,2,3)=C44
	    		E(2,3,3,2)=C44
	    		E(3,2,2,3)=C44
	    		E(3,2,3,2)=C44
	    		E(1,3,1,3)=C55
	    		E(1,3,3,1)=C55
	    		E(3,1,1,3)=C55
	    		E(3,1,3,1)=C55
            	E(1,2,1,2)=C66
	    		E(1,2,2,1)=C66
	    		E(2,1,1,2)=C66
	    		E(2,1,2,1)=C66
            	!...piezoelectric constants
            	e15 = material(10,matl_no)
            	e24 = material(11,matl_no)
            	e31 = material(12,matl_no)
            	e32 = material(13,matl_no)
            	e33 = material(14,matl_no)
            	E(3,1,4,1) = e15
            	E(1,3,4,1) = e15
            	E(1,4,3,1) = e15
            	E(1,4,1,3) = e15
            	E(2,3,4,2) = e24
            	E(3,2,4,2) = e24
            	E(2,4,2,3) = e24
            	E(2,4,3,2) = e24
            	E(1,1,4,3) = e31
            	E(3,4,1,1) = e31
            	E(2,2,4,3) = e32
            	E(3,4,2,2) = e32
            	E(3,3,4,3) = e33
            	E(3,4,3,3) = e33
            	!...dielectric permittivity
            	k11 = material(15,matl_no)
            	k22 = material(16,matl_no)
            	k33 = material(17,matl_no)
            	E(1,4,4,1) = -k11
            	E(2,4,4,2) = -k22
            	E(3,4,4,3) = -k33
        	else
          		print*,'not yet orthotropy for media:',NDOF
            	stop
        	endif
		elseif (matl_type.eq.5) then
    	!...monoclinic material
        	if (NDOF.eq.2) then
			!...elastic media, plane strain: z=0 is plane of symmetry
    			!...7/16/09: modify 12->6, 23->4, 31->5
    			C11 = material(1,matl_no)
        		C12 = material(2,matl_no)
        		C16 = material(3,matl_no)        
        		C22 = material(4,matl_no)
        		C26 = material(5,matl_no)
        		C44 = material(6,matl_no)
        		C45 = material(7,matl_no)
        		C55 = material(8,matl_no)
        		C66 = material(9,matl_no)
				E(1,1,1,1)=C11
	    		E(1,1,2,2)=C12
				E(1,1,1,2)=C16
				E(1,1,2,1)=C16
				E(2,2,1,1)=C12
				E(2,2,2,2)=C22
				E(2,2,1,2)=C26
				E(2,2,2,1)=C26
            	E(3,2,2,3)=C44
            	E(3,2,1,3)=C45
            	E(3,1,2,3)=C45
            	E(3,1,1,3)=C55
            	E(1,2,1,1)=C16
            	E(1,2,2,2)=C26
            	E(1,2,1,2)=C66
				E(1,2,2,1)=C66
            	E(2,1,1,1)=C16
            	E(2,1,2,2)=C26
				E(2,1,1,2)=C66
				E(2,1,2,1)=C66
    		elseif (NDOF.eq.3) then
			!...elastic media: z=0 is plane of symmetry
    		!...7/16/09: modify 12->6, 23->4, 31->5
    			C11 = material(1,matl_no)
        		C12 = material(2,matl_no)
        		C13 = material(3,matl_no)
        		C16 = material(4,matl_no)        
        		C22 = material(5,matl_no)
        		C23 = material(6,matl_no)
        		C26 = material(7,matl_no)
        		C33 = material(8,matl_no)
        		C36 = material(9,matl_no)
        		C44 = material(10,matl_no)
        		C45 = material(11,matl_no)
        		C55 = material(12,matl_no)
        		C66 = material(13,matl_no)
				E(1,1,1,1)=C11
	    		E(1,1,2,2)=C12
	    		E(1,1,3,3)=C13
				E(1,1,1,2)=C16
				E(1,1,2,1)=C16
				E(1,2,1,1)=C16
				E(2,1,1,1)=C16
				E(2,2,1,1)=C12
				E(2,2,2,2)=C22
				E(2,2,3,3)=C23
				E(2,2,1,2)=C26
				E(2,2,2,1)=C26
				E(1,2,2,2)=C26
				E(2,1,2,2)=C26
				E(3,3,1,1)=C13
				E(3,3,2,2)=C23
				E(3,3,3,3)=C33
				E(3,3,1,2)=C36
				E(3,3,2,1)=C36
				E(1,2,3,3)=C36
				E(2,1,3,3)=C36
				E(2,3,2,3)=C44
				E(2,3,3,2)=C44
				E(3,2,2,3)=C44
				E(3,2,3,2)=C44
				E(2,3,3,1)=C45
				E(2,3,1,3)=C45
				E(3,2,3,1)=C45
				E(3,2,1,3)=C45
				E(3,1,2,3)=C45
				E(3,1,3,2)=C45
				E(1,3,2,3)=C45
				E(1,3,3,2)=C45
				E(3,1,3,1)=C55
				E(3,1,1,3)=C55
				E(1,3,3,1)=C55
				E(1,3,1,3)=C55
				E(1,2,1,2)=C66
				E(1,2,2,1)=C66
				E(2,1,1,2)=C66
				E(2,1,2,1)=C66
        	else
          		print*,'not yet monoclinic for media:',NDOF
            	stop
        	endif
		elseif (matl_type.eq.6) then
    	!...general anisotropy
        	if (NDOF.eq.1) then
        		C11 = material(1,matl_no)
            	C12 = material(2,matl_no)
            	C13 = material(3,matl_no)
            	C22 = material(4,matl_no)
            	C23 = material(5,matl_no)
            	C33 = material(6,matl_no)
            	E(1,1,1,1) = C11
            	E(1,1,1,2) = C12
            	E(1,1,1,3) = C13
            	E(2,1,1,1) = C12
            	E(2,1,1,2) = C22
            	E(2,1,1,3) = C23
            	E(3,1,1,1) = C13
            	E(3,1,1,2) = C23
            	E(3,1,1,3) = C33
        	elseif (NDOF.eq.2) then
        		!...plane strain of elastic media
            	C11 = material(1,matl_no)
        		C12 = material(2,matl_no)
        		C14 = material(3,matl_no)        
        		C15 = material(4,matl_no)
        		C16 = material(5,matl_no)
        		C22 = material(6,matl_no)
        		C24 = material(7,matl_no)
        		C25 = material(8,matl_no)
        		C26 = material(9,matl_no)
        		C44 = material(10,matl_no)
        		C45 = material(11,matl_no)        
        		C46 = material(12,matl_no)
        		C55 = material(13,matl_no)
        		C56 = material(14,matl_no)
        		C66 = material(15,matl_no)
				E(1,1,1,1)=C11
	    		E(1,1,2,2)=C12
            	E(1,1,2,3)=C14
            	E(1,1,1,3)=C15
				E(1,1,1,2)=C16
				E(1,1,2,1)=C16
            	E(2,2,1,1)=C12
				E(2,2,2,2)=C22
            	E(2,2,2,3)=C24
            	E(2,2,1,3)=C25
            	E(2,2,1,2)=C26
				E(2,2,2,1)=C26
            	E(3,2,1,1)=C14
            	E(3,2,2,2)=C24
            	E(3,2,2,3)=C44
            	E(3,2,1,3)=C45
            	E(3,2,1,2)=C46
				E(3,2,2,1)=C46
            	E(3,1,1,1)=C15
            	E(3,1,2,2)=C25
            	E(3,1,2,3)=C45
            	E(3,1,1,3)=C55
            	E(3,1,1,2)=C56
				E(3,1,2,1)=C56
				E(1,2,1,1)=C16
            	E(1,2,2,2)=C26
            	E(1,2,2,3)=C46
            	E(1,2,1,3)=C56
            	E(1,2,1,2)=C66
				E(1,2,2,1)=C66
				E(2,1,1,1)=C16
            	E(2,1,2,2)=C26
				E(2,1,2,3)=C46
            	E(2,1,1,3)=C56
            	E(2,1,1,2)=C66
				E(2,1,2,1)=C66
    		elseif ((NDOF.eq.3).or.(NDOF.eq.4)) then
			!...elastic media, convention: 23->4,31->5,12->6
        	!...Order of input: see following notation
    			C11 = material(1,matl_no)
        		C12 = material(2,matl_no)
        		C13 = material(3,matl_no)
        		C14 = material(4,matl_no)        
        		C15 = material(5,matl_no)
        		C16 = material(6,matl_no)
        		C22 = material(7,matl_no)
        		C23 = material(8,matl_no)
        		C24 = material(9,matl_no)
        		C25 = material(10,matl_no)
        		C26 = material(11,matl_no)
        		C33 = material(12,matl_no)
        		C34 = material(13,matl_no)
        		C35 = material(14,matl_no)
        		C36 = material(15,matl_no)
        		C44 = material(16,matl_no)
        		C45 = material(17,matl_no)        
        		C46 = material(18,matl_no)
        		C55 = material(19,matl_no)
        		C56 = material(20,matl_no)
        		C66 = material(21,matl_no)
				E(1,1,1,1)=C11
	    		E(1,1,2,2)=C12
	    		E(1,1,3,3)=C13
            	E(1,1,2,3)=C14
				E(1,1,3,2)=C14
            	E(1,1,1,3)=C15
				E(1,1,3,1)=C15
				E(1,1,1,2)=C16
				E(1,1,2,1)=C16
            	E(2,2,1,1)=C12
				E(2,2,2,2)=C22
				E(2,2,3,3)=C23
            	E(2,2,2,3)=C24
				E(2,2,3,2)=C24
            	E(2,2,1,3)=C25
				E(2,2,3,1)=C25
            	E(2,2,1,2)=C26
				E(2,2,2,1)=C26
            	E(3,3,1,1)=C13
				E(3,3,2,2)=C23
				E(3,3,3,3)=C33
            	E(3,3,2,3)=C34
				E(3,3,3,2)=C34
            	E(3,3,1,3)=C35
				E(3,3,3,1)=C35
            	E(3,3,1,2)=C36
				E(3,3,2,1)=C36
            	E(2,3,1,1)=C14
            	E(2,3,2,2)=C24
            	E(2,3,3,3)=C34
            	E(2,3,2,3)=C44
				E(2,3,3,2)=C44
            	E(2,3,1,3)=C45
				E(2,3,3,1)=C45
            	E(2,3,1,2)=C46
				E(2,3,2,1)=C46
            	E(3,2,1,1)=C14
            	E(3,2,2,2)=C24
            	E(3,2,3,3)=C34
            	E(3,2,2,3)=C44
				E(3,2,3,2)=C44
            	E(3,2,1,3)=C45
				E(3,2,3,1)=C45
            	E(3,2,1,2)=C46
				E(3,2,2,1)=C46
            	E(1,3,1,1)=C15
            	E(1,3,2,2)=C25
            	E(1,3,3,3)=C35
            	E(1,3,2,3)=C45
				E(1,3,3,2)=C45
            	E(1,3,1,3)=C55
				E(1,3,3,1)=C55
            	E(1,3,1,2)=C56
				E(1,3,2,1)=C56
            	E(3,1,1,1)=C15
            	E(3,1,2,2)=C25
            	E(3,1,3,3)=C35
            	E(3,1,2,3)=C45
				E(3,1,3,2)=C45
            	E(3,1,1,3)=C55
				E(3,1,3,1)=C55
            	E(3,1,1,2)=C56
				E(3,1,2,1)=C56
				E(1,2,1,1)=C16
            	E(1,2,2,2)=C26
            	E(1,2,3,3)=C36
            	E(1,2,2,3)=C46
				E(1,2,3,2)=C46
            	E(1,2,1,3)=C56
				E(1,2,3,1)=C56
            	E(1,2,1,2)=C66
				E(1,2,2,1)=C66
				E(2,1,1,1)=C16
				E(2,1,2,2)=C26
				E(2,1,3,3)=C36
            	E(2,1,2,3)=C46
				E(2,1,3,2)=C46
            	E(2,1,1,3)=C56
				E(2,1,3,1)=C56
				E(2,1,1,2)=C66
				E(2,1,2,1)=C66
        		if (NDOF.eq.4) then
            	!...piezoelectric media
            		e11 = material(22,matl_no)
            		e12 = material(23,matl_no)
            		e13 = material(24,matl_no)
            		e14 = material(25,matl_no)
            		e15 = material(26,matl_no)
            		e16 = material(27,matl_no)
            		e21 = material(28,matl_no)
            		e22 = material(29,matl_no)
            		e23 = material(30,matl_no)
            		e24 = material(31,matl_no)
            		e25 = material(32,matl_no)
            		e26 = material(33,matl_no)
            		e31 = material(34,matl_no)
            		e32 = material(35,matl_no)
            		e33 = material(36,matl_no)
            		e34 = material(37,matl_no)
            		e35 = material(38,matl_no)
            		e36 = material(39,matl_no)
            		k11 = material(40,matl_no)
            		k12 = material(41,matl_no)
            		k13 = material(42,matl_no)
            		k22 = material(43,matl_no)
            		k23 = material(44,matl_no)
            		k33 = material(45,matl_no)
                	E(1,1,4,1) = e11
                	E(1,4,1,1) = e11
                	E(2,2,4,1) = e12
                	E(1,4,2,2) = e12
                	E(3,3,4,1) = e13
                	E(1,4,3,3) = e13
                	E(2,3,4,1) = e14
                	E(3,2,4,1) = e14
                	E(1,4,2,3) = e14
                	E(1,4,3,2) = e14
                	E(3,1,4,1) = e15
                	E(1,3,4,1) = e15
                	E(1,4,3,1) = e15
                	E(1,4,1,3) = e15
                	E(1,2,4,1) = e16
                	E(2,1,4,1) = e16
                	E(1,4,2,1) = e16
                	E(1,4,1,2) = e16

                	E(1,1,4,2) = e21
                	E(2,4,1,1) = e21
                	E(2,2,4,2) = e22
                	E(2,4,2,2) = e22
                	E(3,3,4,2) = e23
                	E(2,4,3,3) = e23
                	E(2,3,4,2) = e24
                	E(3,2,4,2) = e24
                	E(2,4,2,3) = e24
                	E(2,4,3,2) = e24
                	E(3,1,4,2) = e25
                	E(1,3,4,2) = e25
                	E(2,4,3,1) = e25
                	E(2,4,1,3) = e25
                	E(1,2,4,2) = e26
                	E(2,1,4,2) = e26
                	E(2,4,2,1) = e26
                	E(2,4,1,2) = e26

                	E(1,1,4,3) = e31
                	E(3,4,1,1) = e31
                	E(2,2,4,3) = e32
                	E(3,4,2,2) = e32
                	E(3,3,4,3) = e33
                	E(3,4,3,3) = e33
                	E(2,3,4,3) = e34
                	E(3,2,4,3) = e34
                	E(3,4,2,3) = e34
                	E(3,4,3,2) = e34
                	E(3,1,4,3) = e35
                	E(1,3,4,3) = e35
                	E(3,4,3,1) = e35
                	E(3,4,1,3) = e35
                	E(1,2,4,3) = e36
                	E(2,1,4,3) = e36
                	E(3,4,2,1) = e36
                	E(3,4,1,2) = e36

                	E(1,4,4,1) = -k11
                	E(1,4,4,2) = -k12
                	E(2,4,4,1) = -k12
                	E(1,4,4,3) = -k13
                	E(3,4,4,1) = -k13
                	E(2,4,4,2) = -k22
                	E(2,4,4,3) = -k23
                	E(3,4,4,2) = -k23
                	E(3,4,4,3) = -k33
            	endif   !...of if (NDOF.eq.4)
        	else
          		print*,'general anisotropy not done for media:',NDOF
            	stop
        	endif
    	else
      		print*,'kerne1.f90: material type is out of range:',matl_type
        	stop
    	endif   !of matl_type
        !--------------------------------------------------------------------------------
        !...following portion of the code is copied from kernel1c.f90
        !...11/24/2010: delete class_problem
        !if ((class_problem.eq.2).and.(rotate_flag).eq.1) then
        if (rotate_flag.eq.1) then
      		!...form direction cosin matrix a_ij
        	dir_cosine = 0.0d0
        	if (NDOF.eq.1) then
          		dir_cosine = 1.0d0
        	else
          		dir_cosine(1,1) = dcos(rotate_angle)
            	dir_cosine(1,2) = dsin(rotate_angle)
            	dir_cosine(2,1) = -dsin(rotate_angle)
            	dir_cosine(2,2) = dcos(rotate_angle)
        	endif
        	if (NDOF.ge.3) then
          		!...for the 1st and 2nd row
            	do kk = 3,NDOF
              		dir_cosine(1,kk) = 0.d0
                	dir_cosine(2,kk) = 0.d0
            	enddo
            	!...for the 3th, 4th,...NDOFth row
            	do kk = 3,NDOF
          			do n = 1,NDOF
            			if (kk.eq.n) then
                			dir_cosine(kk,n) = 1.0d0
                		else
                  			dir_cosine(kk,n) = 0.0d0
                		endif
            		enddo
        		enddo
        	endif
      		!...anisotropy and material coordinates are different with geometry coordinates
            E_temp = 0.0d0
        	do p = 1,3
          		do q = 1,NDOF
            		do m = 1,NDOF
                		do n = 1,3
                        	do i = 1,3
                          		do j = 1,NDOF
                            		do kk = 1,NDOF
                                		do l = 1,3
                                    		E_temp(p,q,m,n) = E_temp(p,q,m,n) + &
                                        	dir_cosine(i,p)*dir_cosine(j,q)*dir_cosine(kk,m)*dir_cosine(l,n)*E(i,j,kk,l)
                                    	enddo
                                	enddo
                            	enddo
                        	enddo
                    	enddo
                	enddo
            	enddo
        	enddo
        	!...updata moduli matrix after transformation
        	E = E_temp
    	endif
		!--------------------------------------------------------------------------------
       	DO i = 1, total_region_elem(k)
       		ie = ie + 1
       		IF (ELTYPE(elemid(ie)) /= CTIP) CYCLE
       		tip_node_count = tip_node_count + 1
            !...get coordinates of tip element
       		xele = node_coor(:,elnode(:,ie))
            !...get tip "displacement derivatives" from solution
            if (elemid(ie).eq.CTIP1) then
          		!uk(:) = sol_disp(:,1,ie)
          		uk_traction(:) = sol_disp_traction(:,1,ie)
                        do ck = 1,total_cracks
                           uk_pressure(ck,:) = sol_disp_pressure(ck,:,1,ie)
          		enddo

            elseif (elemid(ie).eq.CTIP2) then
               	!uk(:) = sol_disp(:,2,ie)
               	uk_traction(:) = sol_disp_traction(:,2,ie)

                do ck = 1,total_cracks
                   uk_pressure(ck,:) = sol_disp_pressure(ck,:,2,ie)
                enddo

            else
               	print*,'id error of tip element: ',elem_sys2user(ie)
                stop
            endif
       		! compute the SIF
	  		IF (material_type(region_mat_no(k))==1) THEN
           		!CALL CompSIF(elemid(ie),xele,uk,SIF,ey1,pv1)
           		CALL CompSIF(elemid(ie),xele,uk_traction,SIF_traction,ey1,pv1)
                        do ck = 1,total_cracks

                           uk_pressure_temp = uk_pressure(ck,:)

                           CALL CompSIF(elemid(ie),xele,uk_pressure_temp,SIF_pressure_temp,ey1,pv1)

                           SIF_pressure(ck,:)=SIF_pressure_temp

                        enddo
           		

	  		ELSE
           		!CALL CompSIFaniso(NDOF,elemid(ie),xele,uk,SIF,E)
           		CALL CompSIFaniso(NDOF,elemid(ie),xele,uk_traction,SIF_traction,E)
                        do ck = 1,total_cracks
                           CALL CompSIFaniso(NDOF,elemid(ie),xele,uk_pressure(ck,:),SIF_pressure(ck,:),E)
                        enddo

	  		END IF
       		!...put SIF of all tip elements together and calculate pressure scaling...3/7/12...aje567
       		nodal_SIF_traction(:,tip_node_count) = SIF_traction

                do ck = 1,total_cracks
                   nodal_SIF_pressure(ck,:,tip_node_count) = SIF_pressure(ck,:)
                enddo
       	

            !...user # of tip node
       		if (elemid(ie).eq.CTIP1) then
       			tipnode_sys2user(tip_node_count) = node_sys2user(elnode(1,ie))
          	elseif (elemid(ie).eq.CTIP2) then
       			tipnode_sys2user(tip_node_count) = node_sys2user(elnode(2,ie))
            endif
       	END DO ! loop over elements in each region
    END DO ! loop over regions


    !--------------------------------------------------------------------------------------
    !!!... CALCULATE VOLUME OF CRACK DUE TO REMOTE TRACTIONS AND PRESSURES ... added 3/11/12 ... aje567 ...!!!
    !--------------------------------------------------------------------------------------


    !.....update solution displacements to satisfy crack closure condition
    DO ie = 1,total_elem
       DO alpha = 1, NDOF
          DO inode = 1, NODE(elemid(ie))
             if((elemid(ie).eq.CTIP1).and.(inode.eq.1))then
                sol_disp_traction(alpha,inode,ie)=0.0d0
                do ck = 1,total_cracks
                   sol_disp_pressure(ck,alpha,inode,ie)=0.0d0
                enddo
             elseif((elemid(ie).eq.CTIP2).and.(inode.eq.2))then
                sol_disp_traction(alpha,inode,ie)=0.0d0
                do ck = 1,total_cracks
                   sol_disp_pressure(ck,alpha,inode,ie)=0.0d0
                enddo
             endif
          END DO
       END DO
    END DO

    !.....calculate opening displacements due to remote tractions and pressures

         do i = 1,total_elem
            !...use the nodal coordinate of each node on element i
            xele = node_coor(:,elnode(:,i))

            !...compute normals at nodes of element i
            call NodeNormal(elemid(i),xele,normal_i)
               
            !...compute opening displacement at node j of element i
            do j = 1, NODE(elemid(i))

               opening_disp_traction(i,j) = -sol_disp_traction(1,j,i)*normal_i(1,j)-sol_disp_traction(2,j,i)*normal_i(2,j)
               do ck = 1,total_cracks
                  opening_disp_pressure(ck,i,j) = -sol_disp_pressure(ck,1,j,i)*normal_i(1,j)-sol_disp_pressure(ck,2,j,i)*normal_i(2,j)
               enddo

            enddo !..j = 1, NODE(elemid(i))
         enddo !..iee = 1,total_elem

      !...compute the volumes of each crack due to remote tractions and pressures.
         do ck = 1,total_cracks
            call CrackVolume(ck,opening_disp_traction,arc_length(ck),crack_volume_tractions(ck)) !crack volume due to remote tractions
         
            do pr = 1,total_cracks
               call CrackVolume(ck,opening_disp_pressure(pr,:,:),arc_length(ck),crack_volume_pressure(ck,pr)) !crack ck volume due to pressure pr
            enddo
         enddo

         !..print volumes to output file
         do ck = 1,total_cracks
            print*
            print*,"Crack", ck, "length is: ", arc_length(ck)
            print*
            print*,"Crack",ck, "volume due to remote tractions is: ",crack_volume_tractions(ck)
            do pr = 1,total_cracks
               print*,"Crack", ck, "volume due to pressure", pr, " is: ",crack_volume_pressure(ck,pr)
               print*
            enddo
         enddo



!---------------------------------------------------------------------
!------- Calculate Pressure Scaling
!---------------------------------------------------------------------


      if (growth_flag.eq.0)then
         globalK(1)=1.0d0
         globalK(2)=1.0d0
         globalK(3)=1.0d0
      endif
      
      volume_ratio = 1.0d0


      tip_node_count = 0   !total number of tip nodes of ALL regions

      do ie=1,total_elem  !loop over all elements
         if ((elemid(ie).eq.CTIP1).or.(elemid(ie).eq.CTIP2))then   !only choose elements that are tip elements
            
            tip_node_count = tip_node_count + 1 !keep track of tip nodes
            


            ! calculate pressure necessary to grow each tip...OLD LOGIC FOR TWO CRACKS ONLY
!            if (total_cracks .eq. 2.0d0) then
!
!               if (volume_ratio(2) .eq. 1.0d0) then
!
!                  crack_scaling_temp(1,tip_node_count) = (globalK(1)-nodal_SIF_traction(1,tip_node_count))/(nodal_SIF_pressure(1,1,tip_node_count)&
!                       +nodal_SIF_pressure(2,1,tip_node_count))
!
!                  crack_scaling_temp(2,tip_node_count) = crack_scaling_temp(1,tip_node_count)
!
!               elseif (volume_ratio(2) .ne. 1.0d0) then
!
!                  crack_scaling_temp(1,tip_node_count) = (nodal_SIF_pressure(2,1,tip_node_count)*(volume_ratio(2)*crack_volume_tractions(2)-crack_volume_tractions(1))&
!                    +(nodal_SIF_traction(1,tip_node_count)-globalK(1))*(crack_volume_pressure(1,2)-volume_ratio(2)*crack_volume_pressure(2,2)))&
!                    /(nodal_SIF_pressure(2,1,tip_node_count)*crack_volume_pressure(1,1)-nodal_SIF_pressure(1,1,tip_node_count)*crack_volume_pressure(1,2)-volume_ratio(2)&
!                    *nodal_SIF_pressure(2,1,tip_node_count)*crack_volume_pressure(2,1)+volume_ratio(2)*nodal_SIF_pressure(1,1,tip_node_count)*crack_volume_pressure(2,2))
!               
!                  crack_scaling_temp(2,tip_node_count) = (nodal_SIF_pressure(1,1,tip_node_count)*(-volume_ratio(2)*crack_volume_tractions(2)+crack_volume_tractions(1))&
!                    -(nodal_SIF_traction(1,tip_node_count)-globalK(1))*(crack_volume_pressure(1,1)-volume_ratio(2)*crack_volume_pressure(2,1)))&
!                    /(nodal_SIF_pressure(2,1,tip_node_count)*crack_volume_pressure(1,1)-nodal_SIF_pressure(1,1,tip_node_count)*crack_volume_pressure(1,2)-volume_ratio(2)&
!                    *nodal_SIF_pressure(2,1,tip_node_count)*crack_volume_pressure(2,1)+volume_ratio(2)*nodal_SIF_pressure(1,1,tip_node_count)*crack_volume_pressure(2,2))
!               endif   !..volume_ratio
!            else
               !...solve for pressure for multiple cracks

            if (pressure_flag .eq. 1) then

               call Volume_Solver(total_cracks,nodal_SIF_traction,nodal_SIF_pressure,crack_volume_tractions,crack_volume_pressure,volume_ratio,crack_scaling_temp)
               
            elseif (pressure_flag .eq. 2) then
               
               call Pressure_Solver(total_cracks,nodal_SIF_traction,nodal_SIF_pressure,crack_scaling_temp)

            endif

!            endif !..total_cracks

         endif   !..crack_tip_elem
      enddo !...loop over all elements

      print*
      print*,'The temporary pressure scaling values are: '
      print*,crack_scaling_temp


      !min_scaling = minloc(crack_scaling_temp,DIM=2, MASK = crack_scaling_temp .ge. 0.0d0)
      min_scaling = minloc(crack_scaling_temp,DIM=2)
      do ck = 1,total_cracks
         crack_scaling(ck) = crack_scaling_temp(ck,min_scaling(1))
      enddo


!---------------------------------------------------------------------
!----------- Update Solution Stress Intensity Factors
!---------------------------------------------------------------------

      nodal_SIF = nodal_SIF_traction 

      do pr = 1,total_cracks
         nodal_SIF = nodal_SIF+ crack_scaling(pr)*nodal_SIF_pressure(pr,:,:)
      enddo


!---------------------------------------------------------------------
!----------- Update Solution Displacements
!---------------------------------------------------------------------

    !...calculate solution displacements

    do i=1,total_elem
       do j=1,NDOF
          do k=1,NODE(elemid(i))

             sol_disp(j,k,i) = -sol_disp_traction(j,k,i)

             do pr=1,total_cracks
                sol_disp(j,k,i)=sol_disp(j,k,i)-crack_scaling(pr)*sol_disp_pressure(pr,j,k,i)
             enddo

             if((elemid(i).eq.CTIP1).and.(k.eq.1))then
                sol_disp(j,k,i)=0.0d0
             elseif((elemid(i).eq.CTIP2).and.(k.eq.2))then
                sol_disp(j,k,i)=0.0d0
             endif
          enddo !...k=1,NODE
       enddo  !...j=1,NDOF
    enddo  !...i=1,total_elem

    !.....calculate opening displacements

    do i = 1,total_elem
       !...use the nodal coordinate of each node on element i
       xele = node_coor(:,elnode(:,i))
       
       !...compute normals at nodes of element i
       call NodeNormal(elemid(i),xele,normal_i)
       
       !...compute opening displacement at node j of element i
       do j = 1, NODE(elemid(i))
          opening_disp(i,j) = sol_disp(1,j,i)*normal_i(1,j)+sol_disp(2,j,i)*normal_i(2,j)
       enddo !..j = 1, NODE(elemid(i))
    enddo !..iee = 1,total_elem



!---------------------------------------------------------------------
!----------- Update Solution Crack Volumes
!---------------------------------------------------------------------

    do ck = 1,total_cracks
       
       crack_volume_total(ck)=crack_volume_tractions(ck)
 
       do pr = 1,total_cracks
          crack_volume_total(ck) = crack_volume_total(ck)+crack_scaling(pr)*crack_volume_pressure(ck,pr)
       enddo

    enddo



!-----------------------------------------------------------------------
	CALL OutputResult
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!--------- Crack Growth
!-----------------------------------------------------------------------

    !...5/10/09: crack growth added
    if (growth_flag.eq.1) call CrackGrowth





!------------------------------------------------------------
CONTAINS
	SUBROUTINE OutputResult
    IMPLICIT NONE
    
    INTEGER 				:: i, j, k, m, ie, pr, ck
    REAL(KIND=DBL) 	:: eforce(NDOF) !...u(3,3): 1st index is direction, 2nd is node

    !...5/25/09: print step
    print*,'--------------------------------------------------'
    print*,'Result of step: ',nstep
                                   
    PRINT *; PRINT *
    PRINT *, "O U T P U T  O F  R E S U L T S"
    PRINT *
    IF (boundary_flag == YES) THEN
    	PRINT *
      	PRINT *, "TOTAL DISPLACEMENT/JUMP ON BOUNDARIES/CRACKS"
      	PRINT *
      	PRINT *,"Boundary/  Element     Node    Displacement-X  Displacement-Y  ",         &
              "Displacement-Z"
      	PRINT *, "  Crack"

      	! Output displacemets/jump on boundaries/cracks
      	ie = 0
      	DO k = 1, total_region
        	PRINT *, "REGION:",k
                print*,'Total number of elements in region k: ', total_region_elem(k)
        	DO i = 1, total_region_elem(k)
          		ie = ie + 1  ! global element no.
          		IF (i <= total_region_bd_elem(k)) THEN  ! on boundaries
            		DO m = 1, NODE(elemid(ie))
                    	if (NDOF.eq.1) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,e13.6,3x)",                           &
                   			'B',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                    	elseif (NDOF.eq.2) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,2(e13.6,3x))",                           &
                   			'B',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                    	elseif (NDOF.eq.3) then
              				PRINT "(4x,a1,7x,i7,2x,i7,6x,3(e13.6,3x))",                           &
                   			'B',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)

                        elseif (NDOF.eq.4) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,4(e13.6,3x))",                           &
                   			'B',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                        elseif (NDOF.eq.5) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,5(e13.6,3x))",                           &
                   			'B',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                        else
                          	print*,'post.f90: not print displ for N=',NDOF
                            stop
                        endif
            		END DO
          		ELSE IF (i <= total_region_bd_elem(k) + total_region_ck_elem(k)) THEN  
            		! on crack
            		DO m = 1, NODE(elemid(ie))
                    	if (NDOF.eq.1) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,e13.6,3x)",                           &
                   			'C',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                    	elseif (NDOF.eq.2) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,2(e13.6,3x))",                           &
                   			'C',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                    	elseif (NDOF.eq.3) then
              				PRINT "(4x,a1,7x,i7,2x,i7,6x,3(e13.6,3x))",                           &
                   			'C',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)

                        elseif (NDOF.eq.4) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,4(e13.6,3x))",                           &
                   			'C',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                        elseif (NDOF.eq.5) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,5(e13.6,3x))",                           &
                   			'C',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
						endif  
            		END DO
          		ELSE
            		DO m = 1, NODE(elemid(ie))
                    	if (NDOF.eq.1) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,e13.6,3x)",                           &
                   			'I',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                    	elseif (NDOF.eq.2) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,2(e13.6,3x))",                           &
                   			'I',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                    	elseif (NDOF.eq.3) then
              				PRINT "(4x,a1,7x,i7,2x,i7,6x,3(e13.6,3x))",                           &
                   			'I',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)

                        elseif (NDOF.eq.4) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,4(e13.6,3x))",                           &
                   			'I',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                        elseif (NDOF.eq.5) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,5(e13.6,3x))",                           &
                   			'I',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                        endif
            		END DO
          		END IF
          		PRINT *
        	END DO
			PRINT *
      	END DO
    	PRINT *
      	PRINT *, "DISPLACEMENT/JUMP ON BOUNDARIES/CRACKS DUE TO TRACTIONS"
      	PRINT *
      	PRINT *,"Boundary/  Element     Node    Displacement-X  Displacement-Y  ",         &
              "Displacement-Z"
      	PRINT *, "  Crack"

      	! Output displacemets/jump on boundaries/cracks
      	ie = 0
      	DO k = 1, total_region
        	PRINT *, "REGION:",k
                print*,'Total number of elements in region k: ', total_region_elem(k)
        	DO i = 1, total_region_elem(k)
          		ie = ie + 1  ! global element no.
          		IF (i <= total_region_bd_elem(k)) THEN  ! on boundaries
            		DO m = 1, NODE(elemid(ie))
                    	if (NDOF.eq.1) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,e13.6,3x)",                           &
                   			'B',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                    	elseif (NDOF.eq.2) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,2(e13.6,3x))",                           &
                   			'B',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                    	elseif (NDOF.eq.3) then
              				PRINT "(4x,a1,7x,i7,2x,i7,6x,3(e13.6,3x))",                           &
                   			'B',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp_traction(j,m,ie), j = 1,NDOF)

                        elseif (NDOF.eq.4) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,4(e13.6,3x))",                           &
                   			'B',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                        elseif (NDOF.eq.5) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,5(e13.6,3x))",                           &
                   			'B',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                        else
                          	print*,'post.f90: not print displ for N=',NDOF
                            stop
                        endif
            		END DO
          		ELSE IF (i <= total_region_bd_elem(k) + total_region_ck_elem(k)) THEN  
            		! on crack
            		DO m = 1, NODE(elemid(ie))
                    	if (NDOF.eq.1) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,e13.6,3x)",                           &
                   			'C',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                    	elseif (NDOF.eq.2) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,2(e13.6,3x))",                           &
                   			'C',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                    	elseif (NDOF.eq.3) then
              				PRINT "(4x,a1,7x,i7,2x,i7,6x,3(e13.6,3x))",                           &
                   			'C',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp_traction(j,m,ie), j = 1,NDOF)

                        elseif (NDOF.eq.4) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,4(e13.6,3x))",                           &
                   			'C',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                        elseif (NDOF.eq.5) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,5(e13.6,3x))",                           &
                   			'C',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
						endif  
            		END DO
          		ELSE
            		DO m = 1, NODE(elemid(ie))
                    	if (NDOF.eq.1) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,e13.6,3x)",                           &
                   			'I',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                    	elseif (NDOF.eq.2) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,2(e13.6,3x))",                           &
                   			'I',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                    	elseif (NDOF.eq.3) then
              				PRINT "(4x,a1,7x,i7,2x,i7,6x,3(e13.6,3x))",                           &
                   			'I',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp_traction(j,m,ie), j = 1,NDOF)

                        elseif (NDOF.eq.4) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,4(e13.6,3x))",                           &
                   			'I',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                        elseif (NDOF.eq.5) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,5(e13.6,3x))",                           &
                   			'I',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                        endif
            		END DO
          		END IF
          		PRINT *
        	END DO
			PRINT *
      	END DO
    	PRINT *
      	PRINT *, "DISPLACEMENT/JUMP ON BOUNDARIES/CRACKS DUE TO PRESSURE IN CRACK 1"
      	PRINT *
      	PRINT *,"Boundary/  Element     Node    Displacement-X  Displacement-Y  ",         &
              "Displacement-Z"
      	PRINT *, "  Crack"

      	! Output displacemets/jump on boundaries/cracks
      	ie = 0
      	DO k = 1, total_region
        	PRINT *, "REGION:",k
                print*,'Total number of elements in region k: ', total_region_elem(k)
        	DO i = 1, total_region_elem(k)
          		ie = ie + 1  ! global element no.
          		IF (i <= total_region_bd_elem(k)) THEN  ! on boundaries
            		DO m = 1, NODE(elemid(ie))
                    	if (NDOF.eq.1) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,e13.6,3x)",                           &
                   			'B',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                    	elseif (NDOF.eq.2) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,2(e13.6,3x))",                           &
                   			'B',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                    	elseif (NDOF.eq.3) then

              				PRINT "(4x,a1,7x,i7,2x,i7,6x,3(e13.6,3x))",                           &
                   			'I',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp_pressure(1,j,m,ie), j = 1,NDOF)

                        elseif (NDOF.eq.4) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,4(e13.6,3x))",                           &
                   			'B',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                        elseif (NDOF.eq.5) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,5(e13.6,3x))",                           &
                   			'B',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                        else
                          	print*,'post.f90: not print displ for N=',NDOF
                            stop
                        endif
            		END DO
          		ELSE IF (i <= total_region_bd_elem(k) + total_region_ck_elem(k)) THEN  
            		! on crack
            		DO m = 1, NODE(elemid(ie))
                    	if (NDOF.eq.1) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,e13.6,3x)",                           &
                   			'C',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                    	elseif (NDOF.eq.2) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,2(e13.6,3x))",                           &
                   			'C',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                    	elseif (NDOF.eq.3) then

              				PRINT "(4x,a1,7x,i7,2x,i7,6x,3(e13.6,3x))",                           &
                   			'I',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp_pressure(1,j,m,ie), j = 1,NDOF)

                        elseif (NDOF.eq.4) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,4(e13.6,3x))",                           &
                   			'C',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                        elseif (NDOF.eq.5) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,5(e13.6,3x))",                           &
                   			'C',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
						endif  
            		END DO
          		ELSE
            		DO m = 1, NODE(elemid(ie))
                    	if (NDOF.eq.1) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,e13.6,3x)",                           &
                   			'I',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                    	elseif (NDOF.eq.2) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,2(e13.6,3x))",                           &
                   			'I',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                    	elseif (NDOF.eq.3) then

              				PRINT "(4x,a1,7x,i7,2x,i7,6x,3(e13.6,3x))",                           &
                   			'I',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp_pressure(1,j,m,ie), j = 1,NDOF)

                        elseif (NDOF.eq.4) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,4(e13.6,3x))",                           &
                   			'I',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                        elseif (NDOF.eq.5) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,5(e13.6,3x))",                           &
                   			'I',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                        endif
            		END DO
          		END IF
          		PRINT *
        	END DO
			PRINT *
      	END DO
    	PRINT *
      	PRINT *, "DISPLACEMENT/JUMP ON BOUNDARIES/CRACKS DUE TO PRESSURE IN CRACK 2"
      	PRINT *
      	PRINT *,"Boundary/  Element     Node    Displacement-X  Displacement-Y  ",         &
              "Displacement-Z"
      	PRINT *, "  Crack"

      	! Output displacemets/jump on boundaries/cracks
      	ie = 0
      	DO k = 1, total_region
        	PRINT *, "REGION:",k
                print*,'Total number of elements in region k: ', total_region_elem(k)
        	DO i = 1, total_region_elem(k)
          		ie = ie + 1  ! global element no.
          		IF (i <= total_region_bd_elem(k)) THEN  ! on boundaries
            		DO m = 1, NODE(elemid(ie))
                    	if (NDOF.eq.1) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,e13.6,3x)",                           &
                   			'B',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                    	elseif (NDOF.eq.2) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,2(e13.6,3x))",                           &
                   			'B',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                    	elseif (NDOF.eq.3) then
              				PRINT "(4x,a1,7x,i7,2x,i7,6x,3(e13.6,3x))",                           &
                   			'I',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp_pressure(2,j,m,ie), j = 1,NDOF)
                        elseif (NDOF.eq.4) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,4(e13.6,3x))",                           &
                   			'B',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                        elseif (NDOF.eq.5) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,5(e13.6,3x))",                           &
                   			'B',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                        else
                          	print*,'post.f90: not print displ for N=',NDOF
                            stop
                        endif
            		END DO
          		ELSE IF (i <= total_region_bd_elem(k) + total_region_ck_elem(k)) THEN  
            		! on crack
            		DO m = 1, NODE(elemid(ie))
                    	if (NDOF.eq.1) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,e13.6,3x)",                           &
                   			'C',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                    	elseif (NDOF.eq.2) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,2(e13.6,3x))",                           &
                   			'C',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                    	elseif (NDOF.eq.3) then
              				!PRINT "(4x,a1,7x,i7,2x,i7,6x,3(e13.6,3x))",                           &
                   			!'I',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			!(sol_disp_pressure(2,j,m,ie), j = 1,NDOF)
                        elseif (NDOF.eq.4) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,4(e13.6,3x))",                           &
                   			'C',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                        elseif (NDOF.eq.5) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,5(e13.6,3x))",                           &
                   			'C',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
						endif  
            		END DO
          		ELSE
            		DO m = 1, NODE(elemid(ie))
                    	if (NDOF.eq.1) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,e13.6,3x)",                           &
                   			'I',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                    	elseif (NDOF.eq.2) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,2(e13.6,3x))",                           &
                   			'I',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                    	elseif (NDOF.eq.3) then
              				PRINT "(4x,a1,7x,i7,2x,i7,6x,3(e13.6,3x))",                           &
                   			'I',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp_pressure(2,j,m,ie), j = 1,NDOF)
                        elseif (NDOF.eq.4) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,4(e13.6,3x))",                           &
                   			'I',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                        elseif (NDOF.eq.5) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,5(e13.6,3x))",                           &
                   			'I',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                   			(sol_disp(j,m,ie), j = 1,NDOF)
                        endif
            		END DO
          		END IF
          		PRINT *
        	END DO
			PRINT *
      	END DO

      	! Output tractions on boundaries/cracks
      	PRINT *
      	PRINT *
      	PRINT *, "TRACTIONS ON BOUNDARIES/CRACKS DUE TO REMOTE STRESS"
      	PRINT *
      	PRINT *,"Boundary/  Element     Node        Traction-X      Traction-Y     ",      &
              " Traction-Z"
      	PRINT *,"  Crack"
      	! values on boundaries
      	ie = 0
      	DO k = 1, total_region
        	!...7/20/09: deactivate this unnecessary portion
 	    	!IF (material_type(region_mat_no(k))==1) THEN 
           	!	ey1  = material(1, region_mat_no(k))
           	!	pv1 = material(2, region_mat_no(k))
           	!	mu = 0.5d0 * ey1 / (1.d0 + pv1)  
	    	!ELSE
	       	!	DO kk=1,21
	        !  		Matprop(kk)=material(kk, region_mat_no(k))
	       	!	END DO
	    	!END IF
        	PRINT *, "REGION:",k
                print*,'Total number of elements in region k: ', total_region_elem(k)

        	DO i = 1, total_region_elem(k)
          		ie = ie + 1  ! global element no.
                        !print*, 'Element number: ' , ie
          		IF (elldid(1,ie) == BTRFREE .or. elldid(1,ie) == CTRFREE) cycle
          		IF (i <= total_region_bd_elem(k)) THEN  ! on boundaries
            		DO m = 1, NODE(elemid(ie))
                    	if (NDOF.eq.1) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,e13.6,3x)",                          &
                    		'B',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                  &
                    		(sol_trac(j,m,ie), j = 1,NDOF)
                    	elseif (NDOF.eq.2) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,2(e13.6,3x))",                          &
                    		'B',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                  &
                    		(sol_trac(j,m,ie), j = 1,NDOF)
                    	elseif (NDOF.eq.3) then
              				PRINT "(4x,a1,7x,i7,2x,i7,6x,3(e13.6,3x))",                          &
                    		'B',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                  &
                    		(sol_trac(j,m,ie), j = 1,NDOF)
                        elseif (NDOF.eq.4) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,4(e13.6,3x))",                          &
                    		'B',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                  &
                    		(sol_trac(j,m,ie), j = 1,NDOF)
                        elseif (NDOF.eq.5) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,5(e13.6,3x))",                          &
                    		'B',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                  &
                    		(sol_trac(j,m,ie), j = 1,NDOF)
                        endif
            		END DO
            		CALL net_force(ie, sol_trac, eforce)
                    if (NDOF.eq.1) then
                      	PRINT "(12x, i7, 15x, e13.6,3x)",                                     &
                  		elem_sys2user(i), (eforce(j), j = 1,NDOF)
                    elseif (NDOF.eq.2) then
                      	PRINT "(12x, i7, 15x, 2(e13.6,3x))",                                     &
                  		elem_sys2user(i), (eforce(j), j = 1,NDOF)
                    elseif (NDOF.eq.3) then
            			PRINT "(12x, i7, 15x, 3(e13.6,3x))",                                     &
                  		elem_sys2user(i), (eforce(j), j = 1,NDOF)
                    elseif (NDOF.eq.4) then
                    	PRINT "(12x, i7, 15x, 4(e13.6,3x))",                                     &
                  		elem_sys2user(i), (eforce(j), j = 1,NDOF)
                    elseif (NDOF.eq.5) then
                    	PRINT "(12x, i7, 15x, 5(e13.6,3x))",                                     &
                  		elem_sys2user(i), (eforce(j), j = 1,NDOF)
                    endif
          		ELSE IF (i <= total_region_bd_elem(k) + total_region_ck_elem(k)) THEN  
            		! on cracks 
            		DO m = 1, NODE(elemid(ie))
                    	if (NDOF.eq.1) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,e13.6,3x)",                              &
                    		'C',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                  &
                    		(sol_trac(j,m,ie), j = 1,NDOF)
                    	elseif (NDOF.eq.2) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,2(e13.6,3x))",                              &
                    		'C',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                  &
                    		(sol_trac(j,m,ie), j = 1,NDOF)
                    	elseif (NDOF.eq.3) then
              				PRINT "(4x,a1,7x,i7,2x,i7,6x,3(e13.6,3x))",                              &
                    		'C',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                  &
                    		(sol_trac_traction(j,m,ie), j = 1,NDOF)
                        elseif (NDOF.eq.4) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,4(e13.6,3x))",                              &
                    		'C',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                  &
                    		(sol_trac(j,m,ie), j = 1,NDOF)
                        elseif (NDOF.eq.5) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,5(e13.6,3x))",                              &
                    		'C',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                  &
                    		(sol_trac(j,m,ie), j = 1,NDOF)
                        endif
            		END DO
          		ELSE  ! on interface element
            		DO m = 1, NODE(elemid(ie))
                    	if (NDOF.eq.1) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,e13.6,3x)",                                &
                    			'I',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                  &
                    			(sol_trac(j,m,ie), j = 1,NDOF)
                    	elseif (NDOF.eq.2) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,2(e13.6,3x))",                                &
                    			'I',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                  &
                    			(sol_trac(j,m,ie), j = 1,NDOF)
                    	elseif (NDOF.eq.3) then
              				PRINT "(4x,a1,7x,i7,2x,i7,6x,3(e13.6,3x))",                                &
                    			'I',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                  &
                    			(sol_trac(j,m,ie), j = 1,NDOF)
                        elseif (NDOF.eq.4) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,4(e13.6,3x))",                              &
                    		'C',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                  &
                    		(sol_trac(j,m,ie), j = 1,NDOF)
                        elseif (NDOF.eq.5) then
                        	PRINT "(4x,a1,7x,i7,2x,i7,6x,5(e13.6,3x))",                              &
                    		'C',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                  &
                    		(sol_trac(j,m,ie), j = 1,NDOF)
                        endif
            		END DO
          		END IF
          		PRINT *
        	END DO
        	PRINT *
             END DO

      	PRINT *
      	PRINT *
      	PRINT *, "TRACTIONS ON BOUNDARIES/CRACKS DUE TO PRESSURE IN CRACK 1"
      	PRINT *
      	PRINT *,"Boundary/  Element     Node        Traction-X      Traction-Y     ",      &
              " Traction-Z"
      	PRINT *,"  Crack"
        ie=0
        DO k = 1, total_region
        	DO i = 1, total_region_elem(k)
          		ie = ie + 1  ! global element no.
                        !print*, 'Element number: ' , ie
            		DO m = 1, NODE(elemid(ie))
             				PRINT "(4x,a1,7x,i7,2x,i7,6x,3(e13.6,3x))",                              &
                    		'C',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                  &
                    		(sol_trac_pressure(1,j,m,ie), j = 1,NDOF)
                                     END DO
                                  END DO
                               END DO
      	PRINT *
      	PRINT *
      	PRINT *, "TRACTIONS ON BOUNDARIES/CRACKS DUE TO PRESSURE IN CRACK 2"
      	PRINT *
      	PRINT *,"Boundary/  Element     Node        Traction-X      Traction-Y     ",      &
              " Traction-Z"
      	PRINT *,"  Crack"
        ie=0

      	DO k = 1, total_region
        	DO i = 1, total_region_elem(k)
          		ie = ie + 1  ! global element no.
                        !print*, 'Element number: ' , ie
            		DO m = 1, NODE(elemid(ie))
             			!	PRINT "(4x,a1,7x,i7,2x,i7,6x,3(e13.6,3x))",                              &
                    		!'C',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                  &
                    		!(sol_trac_pressure(2,j,m,ie), j = 1,NDOF)
                                     END DO
                                  END DO
                               END DO

	END IF   !...of (boundary_flag.eq.yes)

	IF (sif_flag == YES) THEN
    	!...output SIF
      	IF (total_ckelem == 0) RETURN
      	PRINT *; PRINT *
      	PRINT *, "  S T R E S S   I N T E N S I T Y   F A C T O R S"
      	PRINT *
        if (NDOF.eq.1) then
          	PRINT *, "   Node          KI"
      		PRINT *
      		DO i = 1, tip_node_count
        		PRINT "(1x, i7, 3x, f15.7,1x)", &
               	tipnode_sys2user(i), (nodal_SIF(k,i),k=1,NDOF)
      		END DO
        elseif (NDOF.eq.2) then
          	PRINT *, "   Node          KI              KII"
      		PRINT *
      		DO i = 1, tip_node_count
        		PRINT "(1x, i7, 3x, 2(f15.7,1x))", &
               	tipnode_sys2user(i), (nodal_SIF(k,i),k=1,NDOF)
      		END DO
        elseif (NDOF.eq.3) then
                PRINT *, "      From Total Tractions   "
      		PRINT *, "   Node          KI              KII             KIII"
      		PRINT *
      		DO i = 1, tip_node_count
        		PRINT "(1x, i7, 3x, 3(f15.7,1x))", &
               	tipnode_sys2user(i), (nodal_SIF(k,i),k=1,NDOF)
      		END DO
                PRINT *, "      From Remote Tractions   "
    		PRINT *, "   Node          KI              KII             KIII"
      		PRINT *
      		DO i = 1, tip_node_count
        		PRINT "(1x, i7, 3x, 3(f15.7,1x))", &
               	tipnode_sys2user(i), (nodal_SIF_traction(k,i),k=1,NDOF)
      		END DO
                DO pr=1,total_cracks
                   PRINT *, "      From Internal Pressure ",pr
                   PRINT *, "   Node          KI              KII             KIII"
                   PRINT *
                   DO i = 1, tip_node_count
                      PRINT "(1x, i7, 3x, 3(f15.7,1x))", &
                           tipnode_sys2user(i), (nodal_SIF_pressure(pr,k,i),k=1,NDOF)
                   END DO
                ENDDO
            
        elseif (NDOF.eq.4) then
        	PRINT *, "   Node          KI              KII             KIII             KIV"
      		PRINT *
      		DO i = 1, tip_node_count
        		PRINT "(1x, i7, 3x, 4(es15.7,1x))", &
               	tipnode_sys2user(i), (nodal_SIF(k,i),k=1,NDOF)
      		END DO
        elseif (NDOF.eq.5) then
        	PRINT '(a85)', "   Node          KI              KII             KIII             KIV              KV"
      		PRINT *
      		DO i = 1, tip_node_count
        		PRINT "(1x, i7, 3x, 5(es15.7,1x))", &
               	tipnode_sys2user(i), (nodal_SIF(k,i),k=1,NDOF)
      		END DO
             endif
 

     	PRINT *; PRINT *
      	PRINT *, "  C A L C U L A T E D       P R E S S U R E    "
        PRINT *, "---------------------------------------------------"
        PRINT *, "The calculated pressures for load step",nstep," are: "
        DO ck=1,total_cracks
           PRINT*, "Pressure in crack ", ck, " is: ", crack_scaling(ck)*pressure_value
        ENDDO
      	PRINT *; PRINT *
     	PRINT *; PRINT *
      	PRINT *, "  C A L C U L A T E D       V O L U M E   "
        PRINT *, "---------------------------------------------------"
        PRINT *, "The calculated volumes for load step",nstep," are: "
        DO ck = 1,total_cracks
           PRINT *, "Volume of crack ", ck, " is: ", crack_volume_total(ck)
        ENDDO
      	PRINT *; PRINT *
        PRINT *, " C A L C U L A T E D   T O T A L   O P E N I N G      D I S P L A C E M E N T S "
        PRINT *, "---------------------------------------------------"
     	PRINT *,"Boundary/  Element     Node    Opening Displacement  "
      	PRINT *, "  Crack"

      	! Output opening displacemets on cracks
      	ie = 0
      	do k = 1, total_region
           PRINT *, "REGION:",k
           print*,'Total number of elements in region k: ', total_region_elem(k)
           do i = 1, total_region_elem(k)
              ie = ie + 1  ! global element no.
              !do m = 1, NODE(elemid(ie))
                 if (NDOF.eq.3) then
                    PRINT "(4x,a1,7x,i7,2x,i7,6x,e13.6,3x)",                           &
                         'C',elem_sys2user(ie), node_sys2user(elnode(1,ie)),                   &
                         opening_disp(ie,1)
                    PRINT "(4x,a1,7x,i7,2x,i7,6x,e13.6,3x)",                           &
                         'C',elem_sys2user(ie), node_sys2user(elnode(3,ie)),                   &
                         opening_disp(ie,3)
                 endif
              !enddo
           enddo
           
        enddo

        !PRINT *, " C A L C U L A T E D    O P E N I N G      D I S P L A C E M E N T S  F R O M   T R A C T I O N S"
        !PRINT *, "---------------------------------------------------"
     	!PRINT *,"Boundary/  Element     Node    Opening Displacement  "
      	!PRINT *, "  Crack"

      	! Output opening displacemets on cracks
      	!ie = 0
      	!do k = 1, total_region
           !PRINT *, "REGION:",k
           !print*,'Total number of elements in region k: ', total_region_elem(k)
           !do i = 1, total_region_elem(k)
              !ie = ie + 1  ! global element no.
              !do m = 1, NODE(elemid(ie))
                 !if (NDOF.eq.3) then
                    !PRINT "(4x,a1,7x,i7,2x,i7,6x,3(e13.6,3x))",                           &
                         !'C',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                         !opening_disp_traction(ie,m)
                 !endif
              !enddo
           !enddo
        !enddo
        !DO ck = 1,total_cracks
           !PRINT *, " C A L C U L A T E D    O P E N I N G      D I S P L A C E M E N T S  F R O M   P R E S S U R E ", ck
           !PRINT *, "---------------------------------------------------"
           !PRINT *,"Boundary/  Element     Node    Opening Displacement  "
           !PRINT *, "  Crack"

           ! Output opening displacemets on cracks
           ie = 0
           !do k = 1, total_region
              !PRINT *, "REGION:",k
              !print*,'Total number of elements in region k: ', total_region_elem(k)
              !do i = 1, total_region_elem(k)
                 !ie = ie + 1  ! global element no.
                 !do m = 1, NODE(elemid(ie))
                    !if (NDOF.eq.3) then
                       !PRINT "(4x,a1,7x,i7,2x,i7,6x,3(e13.6,3x))",                           &
                            !'C',elem_sys2user(ie), node_sys2user(elnode(m,ie)),                   &
                            !opening_disp_pressure(ck,ie,m)
                    !endif
                 !enddo
              !enddo
           !enddo
        !ENDDO


        !....Output of Pressure and Displacement Data for Creating Pressure vs. Volume Graphs

        !..Pressure_Profile.dat
        if (total_cracks .eq. 2) then
           write(115,'(i5,2x,f10.5,2x,f10.5,2x,f10.5,2x,f10.5,2x,f10.5,2x,f10.5)') nstep, arc_length(1),&
             crack_volume_total(1), crack_scaling(1)*pressure_value,arc_length(2), crack_volume_total(2),&
             crack_scaling(2)*pressure_value

        elseif (total_cracks .eq. 3) then

           write(115,'(i5,2x,f10.5,2x,f10.5,2x,f10.5,2x,f10.5,2x,f10.5,2x,f10.5,2x,f10.5,2x,f10.5,2x,f10.5)') nstep, arc_length(1),&
            crack_volume_total(1), crack_scaling(1)*pressure_value,arc_length(2), crack_volume_total(2),&
             crack_scaling(2)*pressure_value,arc_length(3),&
             crack_volume_total(3), crack_scaling(3)*pressure_value

        elseif (total_cracks .eq. 4) then

           write(115,'(i5,2x,f10.5,2x,f10.5,2x,f10.5,2x,f10.5,2x,f10.5,2x,f10.5,2x,f10.5,2x,f10.5,2x,f10.5,&
                2x,f10.5,2x,f10.5,2x,f10.5)') nstep, arc_length(1), crack_volume_total(1), &
                crack_scaling(1)*pressure_value,arc_length(2), crack_volume_total(2),&
                crack_scaling(2)*pressure_value,arc_length(3),&
                crack_volume_total(3), crack_scaling(3)*pressure_value,arc_length(4),&
                crack_volume_total(4), crack_scaling(4)*pressure_value


        endif


          ENDIF
	END SUBROUTINE OutputResult
	!--------------------------------------------------
	SUBROUTINE net_force(ie, sol_trac, force)
    ! 
    ! compute the net force at the boundary(for each element)
    !
    ! arguments
    ! ie        -- the element number(IN)
    ! sol_trac  -- solution of traction(IN)
    ! force     -- net force on the element(OUT)
	!
    IMPLICIT NONE
    ! arguments
    INTEGER, INTENT(IN)          :: ie
    REAL(DBL), INTENT(IN)  :: sol_trac(:,:,:)
    REAL(DBL), INTENT(OUT) :: force(:)
    ! locall variables
    INTEGER         :: nodes, i, j, k, ni
    REAL(DBL) :: x, xyz(2,3), psi(3), dpsi(3), Jac, dxds(2), Jwi, temp

    nodes = NODE(elemid(ie))  ! number of nodes of the element
    xyz =  node_coor(:, elnode(:,ie))
    !...integration order
    ni = 30
    force = 0.0d0
    !...loop over integration points
    DO i = 1,ni
      	x = xi(i,ni)
      	! call shaphe function w.r.t x
      	CALL reg_shape(elemid(ie), x, psi, dpsi)
		!...derivative of position vector
      	dxds = 0.0d0
      	DO j = 1, nodes
        	DO k = 1,2
				dxds(k) = dxds(k) + xyz(k,j)*dpsi(j)
          	END DO
        END DO
      	! compute Jacobian(x)
      	Jac = dsqrt(dxds(1)*dxds(1) + dxds(2)*dxds(2))
        !...jacobian*weight
      	Jwi = Jac * wi(i,ni)
      	DO k = 1,NDOF
        	temp = 0.0d0
        	DO j = 1, nodes
          		temp = temp + sol_trac(k, j, ie) * psi(j) * Jwi
        	END DO
        	force(k) = force(k) + temp
      	END DO
    END DO
  	END SUBROUTINE net_force

    !------------------------------------------------------------
    !...5/10/09: subroutine to simulate growth of crack
    Subroutine CrackGrowth
    IMPLICIT NONE
	!...local variables
    INTEGER					:: ie,i,j,k,m,index,tip_type,ielem_adjacent,ierr,tnode,iequ, ck
    INTEGER					:: equ_i(3),tcount(3),dcount(3)
    REAL(KIND=DBL)	:: xt1(3),xt2(3),xt3(3),xt4(3),xt5(3),coor_new(3)
    REAL(KIND=DBL)	:: node4_elemi(3),node8_elemi(3),node8_elemadj(3),node4_elemnew(3),node8_elemnew(3)
    REAL(KIND=DBL)	:: dxds(2),dj,xitip,ps(3),dps(3),xele(2,3)
    REAL(KIND=DBL)	:: eita(2),fa(2),e1(3),e2(3),e3(3),C(6,6),S(6,6),Slocal(6,6)
    REAL(KIND=DBL)	:: sif1,sif2,sif3,sif1_traction,sif2_traction,sif3_traction,prop_angle
    REAL(KIND=DBL), ALLOCATABLE   :: sif1_pressure(:),sif2_pressure(:),sif3_pressure(:)
    REAL(KIND=DBL)	:: ratio_hoopstressKc,ratio_hoopstressKc_max,ratio
    REAL(KIND=DBL)	:: d_check1,d_check2,ad_step_size,e_prop(2),adv_vector(2)
    REAL(KIND=DBL)	:: temp
    INTEGER					:: total_region_tip
    !...5/27/09: add variables to store nodal normal vectors of 3 modified/new elements
    REAL(KIND=DBL)	:: normal_iee(2,3),normal_adj(2,3),normal_avg(2)
    INTEGER					:: ielem_adj
    
    INTEGER, ALLOCATABLE				:: node_sys2user_temp(:),node_id_temp(:)
    INTEGER, ALLOCATABLE				:: elem_sys2user_temp(:),elemid_temp(:),elem_crack_id_temp(:),elem_region_temp(:),&
    									   elnode_temp(:,:),elldid_temp(:,:)
	REAL(KIND=DBL), ALLOCATABLE	:: node_coor_temp(:,:),node_disp_val_temp(:,:),elnode_val_temp(:,:,:),&
                                           elnode_val_trac_temp(:,:,:),sol_trac_temp(:,:,:)
        !....add temp variable to keep track of traction..2/29/12...aje567
    
	!...variables to store new nodes and elements: max new element is 50, can be modified here if needed
    INTEGER, parameter		:: max_nnew_node = 100, max_nnew_elem = 50
    INTEGER					:: node_sys2user_new(max_nnew_node)
    REAL(KIND=DBL)	:: node_coor_new(2,max_nnew_node)
    INTEGER					:: elem_sys2user_new(max_nnew_elem),elnode_new(3,max_nnew_elem)
    INTEGER					:: max_usernode,max_userelem,nnew_node,nnew_elem
    INTEGER					:: iee,first_inter_elem
    INTEGER					:: i_newelem
    INTEGER					:: zero_elem !...this is used to save to (global) position of the first elem of current region

    !...5/24/09: add variable to store growth data for all crack tip elements (in one region)
    INTEGER					:: growth_data1(max_nnew_elem) !store node_sys2user of tip node
    REAL(KIND=DBL)	:: growth_data2(max_nnew_elem,4) !store local frame e1 and e2
    REAL(KIND=DBL)	:: growth_data3(max_nnew_elem,2) !store growth angle, ratio of hoop stress/Kc

	!------------------------------------------------------------
	!...output current mesh to techplot file
    write(112,'(a31,i5,a4,i5)'),'ZONE F=FEPOINT, ET=LINESEG, N= ',total_node,',E= ',total_elem
    do i = 1,total_node
      	write(112,'(i5,2x,f15.8,2x,f15.8,2x,a4)'),node_sys2user(i),node_coor(1,i),node_coor(2,i),'0'
    enddo
    ie = 0
    do k = 1,total_region
      	do i = 1,total_region_elem(k)
        	ie = ie + 1
            write(112,'(i5,2x,i5)'),elnode(1,ie),elnode(2,ie)
        enddo
    enddo
    !------------------------------------------------------------
    !...5/25/09: for debugging, create input file for current step result
    write(113,*)'Step: ',nstep
    write(113,'(4(i1,2x))')model_flag,boundary_flag,sif_flag
    !write(113,'(i1)')class_problem !11/24/2010: delete class_problem since it is irrelevant
    !if (class_problem.eq.2) then
    !  	!...anisotropy: current version rotate specimen -> no rotate any more
    !  	write(113,'(a1)') '0'
    !endif
    !...total material
    write(113,'(i2)') total_material
    do i = 1,total_material
    	write(113,*) material_name(i)
    	write(113,'(i1)') material_type(i)
        select case (material_type(i))
           	case (1)
               	!...isotropic: read value of E and nu
                write(113,'(2(f15.5,1x))')material(1,i),material(2,i) !...E and nu
           	case (2)
              	!...cubic
                write(113,'(3(f15.5,1x))')material(1,i),material(2,i),material(3,i)
            case (3)
             	!...transversely isotropic
                write(113,'(5(f15.5,1x))')material(1,i),material(2,i),material(3,i),&
                		material(4,i),material(5,i)
    		case (4)
              	!...orthotropic
                write(113,'(9(f15.5,1x))')material(1,i),material(2,i),material(3,i),&
                   		material(4,i),material(5,i),material(6,i),material(7,i),&
                        material(8,i),material(9,i)
            case (5)
               	!...monoclinic
                write(113,'(13(f15.5,1x))')material(1,i),material(2,i),material(3,i),&
                  		material(4,i),material(5,i),material(6,i),material(7,i),&
                        material(8,i),material(9,i),material(10,i),material(11,i),&
                        material(12,i),material(13,i)
            case (6)
               	!...general anisotropic
                write(113,'(21(f15.5,1x))')material(1,i),material(2,i),material(3,i),&
                  		material(4,i),material(5,i),material(6,i),material(7,i),&
                        material(8,i),material(9,i),material(10,i),material(11,i),&
                        material(12,i),material(13,i),material(14,i),material(15,i),&
                        material(16,i),material(17,i),material(18,i),material(19,i),&
                        material(20,i),material(21,i)
            case default
               	print*,'material type is out of range for anisotropy (1-6):',material_type(i)
                stop
        end select
    enddo
    write(113,'(i2)')total_region
    do i = 1,total_region
    	write(113,*)region_name(i)
    	write(113,'(i2)')region_mat_no(i)
    enddo
    do i = 1,total_region
    	write(113,'(3(i5,1x))')total_region_elem(1),total_region_bd_elem(1),total_region_ck_elem(1)
    	write(113,'(3(i5,1x))')total_region_node(1),total_region_bd_node(1),total_region_ck_node(1)
    enddo
    write(113,'(a1)')'1'
    write(113,*)'This is the input after remeshing of cycle ',nstep
    do i = 1,total_node
      	write(113,'(i5,1x,2(f10.3,1x),i5)')node_sys2user(i),(node_coor(j,i),j=1,2),node_id(i)
    enddo
    do i = 1,total_elem
      	write(113,'(10(i5,1x))')elem_sys2user(i),elemid(i),elem_crack_id(i),elem_region(i),(elnode(j,i),j=1,3),(elldid(j,i),j=1,3)
    enddo
    do i=1,total_elem
      	if (elldid(1,i).ne.BTRFREE.and.elldid(1,i).ne.CTRFREE.and.elldid(1,i).ne.NOLOAD) then
        	!...element i is prescribed either BDISP or BTRAC or CTRAC
            do j = 1,NODE(elemid(i))
              	write(113,'(3(f15.5,1x))') (elnode_val(k,j,i),k=1,3)
            enddo
        else
        	!...do nothing
        endif
    enddo
    write(113,'(i2)')total_rigid_constrain_point
    do i = 1,total_rigid_constrain_point
      	write(113,'(4(i5,1x))')rigid_constrain_point(i),(rigid_constrain_dir(i,j),j=1,3)
    enddo
    write(113,'(a1)')'0' !for integration flag
    !...5/11/09: parameters for remeshing
    write(113,'(a1)')growth_flag
    if (growth_flag.eq.YES) then
      	write(113,'(i5)')total_step
      	write(113,'(f8.2,1x,f8.2)')size_cracktip,grow_factor
        write(113,'(f8.2,1x,f8.2)')dcheck1_factor,dcheck2_factor
        write(113,'(f8.2,1x,f8.2,1x,f8.2)')globalK(1),globalK(2),globalK(3)
        write(113,'(f8.4,1x,f8.4,1x,f8.4)')parameter_alpha,parameter_beta,parameter_epsilon
        write(113,'(i1)')pressure_flag
       	if (pressure_flag.ne.0) write(113,'(f12.5)')pressure_value
    endif
    !------------------------------------------------------------
    !...output current mesh to Abaqus file
    !...print out node data
    write(114,*)'Step: ',nstep
    write(114,*)'*Node'
    do i = 1,total_node
      	write(114,'(i10,3(a1,f15.5))'),node_sys2user(i),',',node_coor(1,i),',',node_coor(2,i)
    enddo
    !...print out 3-node elements
    write(114,*)'*Element,type=B22'
    do i = 1,total_elem
      	if ((elemid(i).eq.BQUAD).or.(elemid(i).eq.CQUAD).or.(elemid(i).eq.CTIP1)&
        .or.(elemid(i).eq.CTIP2).or.(elemid(i).eq.IQUAD)) then
        	write(114,'(i10,3(a1,i10))'),elem_sys2user(i),',',node_sys2user(elnode(1,i)),',',node_sys2user(elnode(3,i)),&
            ',',node_sys2user(elnode(2,i))
        endif
    enddo
    !...print out 2-node elements
    write(114,*)'*Element,type=B21'
    do i = 1,total_elem
      	if ((elemid(i).eq.BLIN).or.(elemid(i).eq.CLIN).or.(elemid(i).eq.ILIN)) then
        	write(114,'(i10,2(a1,i10))'),elem_sys2user(i),(',',node_sys2user(elnode(j,i)),j=1,2)
        endif
    enddo
    !...create set of nodes (exclude mid-side nodes) for displaying the mesh
    write(114,*)'*Nset,nset=edge_node'
    do i = 1,total_elem
      	write(114,*)node_sys2user(elnode(1,i))
        write(114,*)node_sys2user(elnode(2,i))
    enddo
    

    !--------------------------------------------------------------------------------

    !...S T A R T   R E M E S H I N G

    !--------------------------------------------------------------------------------

    ie = 0
    zero_elem = 0
    !...loop over regions
    do k = 1,total_region
       if (total_region_ck_elem(k).eq.0) then
        	!...no crack in region k, then just need to get the global element ie
            ie = ie + total_region_elem(k)
            cycle
        elseif (total_region_ck_elem(k).gt.0) then
        	!...5/27/09: get the global position for the first elem of region k
            zero_elem = ie
        	!...get compliance matrix needed for calculation of growth angle
      	    if (material_type(region_mat_no(k)).eq.3) then
           		!...transversely isotropic
	    		C(1,1)=material(1,region_mat_no(k))
	    		C(1,2)=material(2,region_mat_no(k))
	    		C(1,3)=material(3,region_mat_no(k))
	    		C(1,4)=0.0d0
	    		C(1,5)=0.0d0
	    		C(1,6)=0.0d0
	    		C(2,2)=material(4,region_mat_no(k))
	    		C(2,3)=C(1,2)
	    		C(2,4)=0.0d0
	    		C(2,5)=0.0d0
	    		C(2,6)=0.0d0
	    		C(3,3)=C(1,1)
	    		C(3,4)=0.0d0
	    		C(3,5)=0.0d0
	    		C(3,6)=0.0d0
	    		C(4,4)=material(5,region_mat_no(k))
	    		C(4,5)=0.0d0
	    		C(4,6)=0.0d0
	    		C(5,5)=0.5*(C(1,1)-C(1,3))
	    		C(5,6)=0.0d0
	    		C(6,6)=C(4,4)
            	do i=2,6
    	   			do j=1,i-1
            			C(i,j)=C(j,i)
    	   			enddo
	    		enddo
	    		!...calculate compliance matrix by inverting moduli matrix
	    		call matrix6_inverse(C,S)
        	elseif (material_type(region_mat_no(k)).eq.1) then
            	!...isotropy -> do nothing
            else
          		print*,'CrackGrowth: current version does not support material type: ',material_type(region_mat_no(k))
            	stop
        	endif
            !--------------------------------------------------
        	!...obtain the first element # after all crack elements of region k
            first_inter_elem = 0
        	iee = ie
        	do i = 1,total_region_elem(k)
            	!...element number of region k
          		iee = iee + 1
          		if ((elemid(iee).eq.IQUAD).or.(elemid(iee).eq.ILIN)) then
                	!...loop hit the first interface of this region
                	first_inter_elem = iee
					exit
                endif
        	enddo
            !...defensive programming: this will work only when interface elems defined after crack elems
            if ((first_inter_elem.eq.0).and.(k.lt.total_region)) then
              	!...region k has to have interface elements
                print*,'CrackGrowth: error with numbering of interface elements of region: ',k
                stop
            endif
            !--------------------------------------------------
        	!...initialize number of new nodes and elements
            nnew_node = 0
    		nnew_elem = 0
            !--------------------------------------------------
            !...get current maximum user node # and element #
        	max_usernode = 0
        	max_userelem = 0
        	do j = 1,total_node
           		if (node_sys2user(j).gt.max_usernode) max_usernode = node_sys2user(j)
        	enddo
        	do j = 1,total_elem
           		if (elem_sys2user(j).gt.max_userelem) max_userelem = elem_sys2user(j)
        	enddo
            !--------------------------------------------------
            !...initialize variables to store new mesh of region k
            nnew_node = 0
            nnew_elem = 0
            node_sys2user_new = 0
            node_coor_new = 0.0
            elem_sys2user_new = 0
            elnode_new = 0
            !--------------------------------------------------
            !...allocate variables for temporarily save all mesh data
    		allocate(node_sys2user_temp(total_node),node_coor_temp(2,total_node),&
        		node_id_temp(total_node),node_disp_val_temp(3,total_node),STAT=ierr)
    		if (ierr.ne.0) then
      			print*,'temporary node data: allocation request denied'
        		stop
    		endif
    		allocate(elem_sys2user_temp(total_elem),elemid_temp(total_elem),&
                        elem_crack_id_temp(total_elem),&
        		elem_region_temp(total_elem),STAT=ierr)
    		if (ierr.ne.0) then
      			print*,'temporary elem data: allocation request denied'
        		stop
    		endif
    		allocate(elnode_temp(3,total_elem),elldid_temp(3,total_elem),&
               	elnode_val_temp(3,3,total_elem),elnode_val_trac_temp(3,3,total_elem),&
                sol_trac_temp(3,3,total_elem),STAT=ierr)
    		if (ierr.ne.0) then
      			print*,'temporary elnode and elldid: allocation request denied'
        		stop
    		endif
    		allocate(sif1_pressure(total_cracks),sif2_pressure(total_cracks),sif3_pressure(total_cracks),STAT=ierr)
    		if (ierr.ne.0) then
      			print*,'sif_pressure: allocation request denied'
        		stop
    		endif

        	!--------------------------------------------------
        	!...transfer all old data to temporary variables
        	do i = 1,total_node
          		node_sys2user_temp(i) = node_sys2user(i)
            	node_id_temp(i) = node_id(i)
            	do j = 1,2
            		node_coor_temp(j,i) = node_coor(j,i)
            	enddo
                do j = 1,3
                   	node_disp_val_temp(j,i) = node_disp_val(j,i)
                enddo
        	enddo
        	do i = 1,total_elem
          		elem_sys2user_temp(i) = elem_sys2user(i)
            	elemid_temp(i) = elemid(i)
                elem_crack_id_temp(i) = elem_crack_id(i)
            	elem_region_temp(i) = elem_region(i)
            	do j = 1,3
            		elnode_temp(j,i) = elnode(j,i)
                	elldid_temp(j,i) = elldid(j,i)
                    do m = 1,3
                       	elnode_val_temp(m,j,i) = elnode_val(m,j,i)
                        elnode_val_trac_temp(m,j,i)=elnode_val_trac(m,j,i) !..added 2/29/12..aje567
                        sol_trac_temp(m,j,i)=sol_trac(m,j,i);
                    enddo
            	enddo
        	enddo
            !--------------------------------------------------
            !...compute growth data for all tips of region k
            iee = ie
            total_region_tip = 0
            growth_data1 = 0
            growth_data2 = 0.0
            growth_data3 = 0.0
            ratio_hoopstressKc_max = 0.0
            do i = 1,total_region_elem(k) !ends at line 1836
              	iee = iee + 1
                !...skip process if element iee is not a tip element
            	if (eltype(elemid(iee)).ne.CTIP) cycle  !ends at line 1836
                !--------------------------------------------------
                !...count for number of tips in region k
              	total_region_tip = total_region_tip + 1
                !...get type of tip element iee
            	if (elemid(iee).eq.CTIP1) then
              		tip_type = 1
            	elseif (elemid(iee).eq.CTIP2) then
            		tip_type = 2
            	else
              		print*,'CrackGrowth: wrong id of tip element: ',elem_sys2user(iee)
                	stop
                     endif
            	!--------------------------------------------------
                !...compute growth angle (local)
            	xele = node_coor(:,elnode(:,iee))
            	!...local coordinates and user node # at tip
    			if (tip_type.eq.1) then
					xitip = -1.0d0
                	tnode = node_sys2user(elnode(1,iee))
    			elseif (tip_type.eq.2) then
    				xitip = 1.0d0
                	tnode = node_sys2user(elnode(2,iee))
                     endif
            	!...get SIF's
            	do j = 1,total_tip_node
              		if (tnode.eq.tipnode_sys2user(j)) then


                	sif1_traction = nodal_SIF_traction(1,j)
                    	sif2_traction = nodal_SIF_traction(2,j)
                    	sif3_traction = nodal_SIF_traction(3,j)

                        do ck = 1,total_cracks
                           sif1_pressure(ck) = nodal_SIF_pressure(ck,1,j)
                           sif2_pressure(ck) = nodal_SIF_pressure(ck,2,j)
                           sif3_pressure(ck) = nodal_SIF_pressure(ck,3,j)
                        enddo

                      	sif1 = nodal_SIF(1,j)
                    	sif2 = nodal_SIF(2,j)
                    	sif3 = nodal_SIF(3,j)
                        


                    	exit
                	endif
                	!...control the loop
                	if (j.eq.total_tip_node) then
                  		print*,'CrackGrowth: something wrong with total number of tip nodes'
                    	stop
                	endif
                     enddo
            	!...compute shape functions at node 1 (for type 1 of tip element) or node 2 (for type 2)
    			CALL reg_shape(elemid(iee),xitip,ps,dps)
				!...compute derivative of position vector r
                dxds = 0.0
    			do m = 1,NODE(elemid(iee))
      				do j=1,2
        				dxds(j) = dxds(j) + xele(j,m)*dps(m)
        			enddo
    			enddo
				dj = dsqrt(dxds(1)*dxds(1)+dxds(2)*dxds(2))
				!...compute the UNIT tangential vector eita() (1-axis of local coordinates)
    			DO m = 1,2
      				if (tip_type.eq.1) then
						eita(m) = -dxds(m)/dj
        			else
          				eita(m) = dxds(m)/dj
        			endif
                             END DO
    			!...compute the UNIT normal vector fa
    			!fa is the 2-axis of local coordinates: and is opposite to the normal vector of the Sc+
    			if (tip_type.eq.1) then
					fa(1) = eita(2)
					fa(2) = -eita(1)
    			elseif (tip_type.eq.2) then
    				fa(1) = -eita(2)
					fa(2) = eita(1)
    			endif
                !...store user node #
                growth_data1(total_region_tip) = tnode
                !...store local frame to growth data
                growth_data2(total_region_tip,1) = eita(1)
                growth_data2(total_region_tip,2) = eita(2)
                growth_data2(total_region_tip,3) = fa(1)
                growth_data2(total_region_tip,4) = fa(2)
            	!...form the full local coordinate system: this is needed since we're going to use
                !the subroutines comtransform and growth_direction_sym from 3D code
            	do m = 1,2
              		e1(m) = eita(m)
                	e2(m) = fa(m)
            	enddo
            	e1(3) = 0.0
            	e2(3) = 0.0
            	e3(1) = 0.0
            	e3(2) = 0.0
            	e3(3) = e1(1)*e2(2) - e1(2)*e2(1)


            	if (material_type(region_mat_no(k)).eq.1) then
              		!...isotropic material, using criteria of max hoop stress
                	if (abs(sif2).le.SMALL_NUM) then
                  		prop_angle = 0.0
                	else
                  		!prop_angle = 2.0d0*datan((1.0d0-dsqrt(1.0d0+8.0*(sif2/sif1)*(sif2/sif1)))/(4.0*sif2/sif1))
                        prop_angle = 2.0d0*datan((sif1-dsqrt(sif1*sif1+8.0*sif2*sif2))/(4.0*sif2))

                	endif

                    !...isotropy: ratio of hoop stress/Kc is the value of hoop stress itself

                    temp = 0.5d0*prop_angle
                    ratio_hoopstressKc = sif1*(0.75d0*dcos(temp)+0.25d0*dcos(3.0d0*temp)) - &
                    					sif2*(0.75d0*dsin(temp)+0.75d0*dsin(3.0d0*temp))
                    !...convert angle of growth to degree
                    prop_angle = prop_angle*180.0/PI
            	elseif (material_type(region_mat_no(k)).eq.3) then
            		!...transversely isotropic
					!...get compliance in local coord
					call comtransform(S,e1,e2,e3,Slocal)
            		!...obtain direction of propagation and the ratio of Ke/Kec: no need of e3 and sif3
            		call growth_direction_sym(sif1,sif2,e1,e2,Slocal,globalK,prop_angle,ratio_hoopstressKc)
                     endif
                !...store to growth data
                growth_data3(total_region_tip,1) = prop_angle
                growth_data3(total_region_tip,2) = ratio_hoopstressKc
                !...get the maximum ratio of hoop stress/Kc
                if (ratio_hoopstressKc.ge.ratio_hoopstressKc_max) ratio_hoopstressKc_max = ratio_hoopstressKc

             enddo ! end loop over all tip elements from line 1715


            !--------------------------------------------------
            !...output growth data for checking
            print*,'Node#         e1            e2          Growth angle    Ratio of hoop stress/Kc'
            do i = 1,total_region_tip
              	print'(i5,2x,2(f6.3,1x),5x,2(f6.3,1x),f8.3,5x,f10.5)',growth_data1(i),(growth_data2(i,j),j=1,4),growth_data3(i,1),growth_data3(i,2)
            enddo
            !--------------------------------------------------
            print *,nstep
            !...begin propagation process
            do i = 1,total_region_elem(k)
        		!...global (system) element#
            	ie = ie + 1
                !...skip process if element ie is not a tip element
            	if (eltype(elemid(ie)).ne.CTIP) cycle
                !--------------------------------------------------
                !...get type of tip element iee
            	if (elemid(ie).eq.CTIP1) then
              		tip_type = 1
            	elseif (elemid(ie).eq.CTIP2) then
            		tip_type = 2
            	else
              		print*,'CrackGrowth: wrong id of tip element: ',elem_sys2user(iee)
                	stop
            	endif
            	!--------------------------------------------------
                !...get # of adjacent element ielem_adjacent
            	ielem_adjacent = 0
            	do j = 1,total_elem
              		if ((tip_type.eq.1).and.(elnode(2,ie).eq.elnode(1,j))) then
                		ielem_adjacent = j
                    	exit
                	elseif ((tip_type.eq.2).and.(elnode(1,ie).eq.elnode(2,j))) then
                		ielem_adjacent = j
                    	exit
                	endif
            	enddo
            	if (ielem_adjacent.eq.0) then
              		print*,'CrackGrowth: no adjacent element to tip elem: ',elem_sys2user(ie)
                	stop
            	endif
                
            	!--------------------------------------------------
                !...get coordinates of 5 points needed for remeshing
                xt1 = 0.0 !initialize to get 0 for xti(3)
                xt2 = 0.0
                xt3 = 0.0
                xt4 = 0.0
                xt5 = 0.0
            	do j = 1,2
                	if (tip_type.eq.1) then
            			xt1(j) = node_coor(j,elnode(1,ie))
                    	xt2(j) = node_coor(j,elnode(3,ie))
                    	xt3(j) = node_coor(j,elnode(2,ie))
                    	xt4(j) = node_coor(j,elnode(3,ielem_adjacent))
                    	xt5(j) = node_coor(j,elnode(2,ielem_adjacent))
                	elseif (tip_type.eq.2) then
                		xt1(j) = node_coor(j,elnode(2,ie))
                    	xt2(j) = node_coor(j,elnode(3,ie))
                    	xt3(j) = node_coor(j,elnode(1,ie))
                    	xt4(j) = node_coor(j,elnode(3,ielem_adjacent))
                    	xt5(j) = node_coor(j,elnode(1,ielem_adjacent))
                	endif
            	enddo
            	!--------------------------------------------------
            	!...compute adaptive parameters to form new elements
				d_check1 = dcheck1_factor*size_cracktip
				d_check2 = dcheck2_factor*size_cracktip
            	!...set maximum advance step size: this is going to be advancing of crack front at the point of max K
            	ad_step_size = grow_factor*size_cracktip
                !...obtain the advance of crack tip using bilinear law
                !...local coordinates and user node # at tip
    			if (tip_type.eq.1) then
				   	tnode = node_sys2user(elnode(1,ie))
    			elseif (tip_type.eq.2) then
                	tnode = node_sys2user(elnode(2,ie))
    			endif
                do j = 1,total_region_tip
                  	if (growth_data1(j).eq.tnode) then
                    	prop_angle = growth_data3(j,1)
                        ratio = growth_data3(j,2)/ratio_hoopstressKc_max
                        do m = 1,2
                          	e1(m) = growth_data2(j,m)
                            e2(m) = growth_data2(j,m+2)
                        enddo
                        if (ratio.lt.(1.0d0-parameter_beta)) then
                          	ad_step_size = 0.0
                        elseif (ratio.lt.(1.0d0-parameter_alpha)) then
                        	ad_step_size = (ratio - 1.0d0 + parameter_beta)*(1.0d0 - parameter_alpha)/&
                            				(parameter_beta - parameter_alpha)*ad_step_size
                        !...add parameter epsilon for the case of same ratio_hoopstressKc
                        elseif (ratio.le.(1.0d0-parameter_epsilon)) then
                        	ad_step_size = ratio*ad_step_size
                        elseif (ratio.le.1.0d0) then
                        	!...do nothing, ad_step_size is the maximum ad_step_size associated with max ration_hoopstressKc
                        elseif ((ratio.lt.0.0).or.(ratio.gt.1.0d0)) then
                          	print*,'CrackGrowth: error, ratio_hoop stressKc/ratio_hoopstressKc_max out of range'
                            stop
                        endif
                        exit
                    endif
                enddo
            	!--------------------------------------------------
            	!...calculate propagation direction
				do j=1,2
	    			e_prop(j)=dcos(prop_angle*PI/180.0)*e1(j)+dsin(prop_angle*PI/180.0)*e2(j)
				enddo
				!...vector of advance
				do j=1,2
	    			adv_vector(j)=ad_step_size*e_prop(j)
				enddo
            	!--------------------------------------------------
            	!...obtain coordinates of the new node
                coor_new = 0.0
				do j=1,2
            		if (tip_type.eq.1) then
						coor_new(j)=node_coor(j,elnode(1,ie))+adv_vector(j)
                	elseif (tip_type.eq.2) then
                		coor_new(j)=node_coor(j,elnode(2,ie))+adv_vector(j)
                	endif
				end do
            	!--------------------------------------------------
            	call edge_node_propagate(xt1,xt2,xt3,xt4,xt5,coor_new,size_cracktip,d_check1,d_check2,node4_elemi,node8_elemi,&
					& node8_elemadj,i_newelem,node4_elemnew,node8_elemnew)
				!...update new coordinates of node1/2,node2/1,node3 of element ie and node3 of adjacent element
                !(this is needed when there is no new element created)
				do j = 1,2
                	if (tip_type.eq.1) then
						node_coor(j,elnode(1,ie)) = coor_new(j)
						node_coor(j,elnode(2,ie)) = node4_elemi(j)
                    elseif (tip_type.eq.2) then
                    	node_coor(j,elnode(2,ie)) = coor_new(j)
						node_coor(j,elnode(1,ie)) = node4_elemi(j)
                    endif
					node_coor(j,elnode(3,ie)) = node8_elemi(j)
					node_coor(j,elnode(3,ielem_adjacent)) = node8_elemadj(j)
				enddo
                !...also, need to update data in temporary for the case there are new elements, then all old data are deallocated
                do j = 1,2
                   if (tip_type.eq.1) then
                      node_coor_temp(j,elnode(1,ie)) = coor_new(j)
                      node_coor_temp(j,elnode(2,ie)) = node4_elemi(j)
                   elseif (tip_type.eq.2) then
                      node_coor_temp(j,elnode(2,ie)) = coor_new(j)
                      node_coor_temp(j,elnode(1,ie)) = node4_elemi(j)
                   endif
                   node_coor_temp(j,elnode(3,ie)) = node8_elemi(j)
                   node_coor_temp(j,elnode(3,ielem_adjacent)) = node8_elemadj(j)
                enddo
                !--------------------------------------------------
                !...new element and nodes created:
                if (i_newelem.eq.1) then
                  
                   !...first new node
                   nnew_node = nnew_node + 1

                   if (nnew_node.gt.max_nnew_node) then
                      print*,'number of new nodes exceed max_nnew_node, increase it in post_remesh.f90'
                      stop
                   endif

                   node_sys2user_new(nnew_node) = max_usernode + nnew_node
                   
                   do j = 1,2
                      node_coor_new(j,nnew_node) = node4_elemnew(j)
                   enddo
                   
                   !...second new node
                   nnew_node = nnew_node + 1

                   if (nnew_node.gt.max_nnew_node) then

                      print*,'number of new nodes exceed max_nnew_node, increase it in post_remesh.f90'
                      stop

                   endif

                   node_sys2user_new(nnew_node) = max_usernode + nnew_node

                   do j = 1,2
                      node_coor_new(j,nnew_node) = node8_elemnew(j)
                   enddo
                   
                   !...update connectivity of temporary variables for adjacent element
                   if (tip_type.eq.1) then

                      elnode_temp(1,ielem_adjacent) = total_node + nnew_node - 1

                   elseif (tip_type.eq.2) then

                      elnode_temp(2,ielem_adjacent) = total_node + nnew_node - 1

                   endif

                   !...new element
                   nnew_elem = nnew_elem + 1
                   
                   if (nnew_elem.gt.max_nnew_elem) then
                      print*,'number of new elements exceed max_nnew_elem, increase it in post_remesh.f90'
                      stop
                   endif

                   elem_sys2user_new(nnew_elem) = max_userelem + nnew_elem

                   if (tip_type.eq.1) then

                      elnode_new(1,nnew_elem) = elnode_temp(2,ie)
                      elnode_new(2,nnew_elem) = elnode_temp(1,ielem_adjacent)
                      elnode_new(3,nnew_elem) = total_node + nnew_node

                   elseif (tip_type.eq.2) then

                      elnode_new(1,nnew_elem) = elnode_temp(2,ielem_adjacent)
                      elnode_new(2,nnew_elem) = elnode_temp(1,ie)
                      elnode_new(3,nnew_elem) = total_node + nnew_node

                   endif

                   !...elemid: always CQUAD, elem_region = k, elldid: same as adjacent elem
                endif !of i_newelem
             enddo !...of i=1,total_region_elem(k)

            !----------------------------------------------------------------------------
            !...FINISH REMESHING of all crack tip elements of region k
            !----------------------------------------------------------------------------

            if (nnew_node.gt.0) then
          		!...there are new nodes (and certainly new elements) created during remeshing of region k
        		!--------------------------------------------------
                !...deallocate global variables so that we can re-allocate
                deallocate(node_sys2user,node_coor,node_id,node_disp_val,unknown,node_type,STAT=ierr)
                if (ierr.ne.0) then
                  	print*,'node_sys2user... deallocation request denied'
                    stop
                endif
        		!...re-allocate global node variables
    			allocate(node_sys2user(total_node+nnew_node),node_coor(2,total_node+nnew_node),&
        			node_id(total_node+nnew_node),STAT=ierr)
    			if (ierr.ne.0) then
      				print*,'node data: allocation request denied'
        			stop
    			endif
				!...Note: unknown and node_type will be redo from beginning, thus no need to save
                allocate(node_disp_val(NDOF,total_node+nnew_node),unknown(NDOF,total_node+nnew_node),&
                	node_type(NDOF,total_node+nnew_node),STAT=ierr)
                if (ierr.ne.0) then
                  	print*,'node_disp_val: allocation request denied'
                    stop
                endif
                !--------------------------------------------------
                !...deallocate global variables so that we can re-allocate
                deallocate(elem_sys2user,elemid,elem_crack_id,elem_region,elnode,elldid,elnode_val,&
                     elnode_val_trac,elnode_val_pressure,sol_trac,equno,STAT=ierr)
                if (ierr.ne.0) then
                  	print*,'elem_sys2user... deallocation request denied'
                    stop
                endif
                !...re-allocate global element variables
    			allocate(elem_sys2user(total_elem+nnew_elem),elemid(total_elem+nnew_elem),&
                                elem_crack_id(total_elem+nnew_elem),&
        			elem_region(total_elem+nnew_elem),STAT=ierr)
    			if (ierr.ne.0) then
      				print*,'elem data: allocation request denied'
        			stop
    			endif
                        
                        elem_crack_id = 0
                        
    			allocate(elnode(3,total_elem+nnew_elem),elldid(3,total_elem+nnew_elem),&
                	elnode_val(3,3,total_elem+nnew_elem),&
                        elnode_val_trac(3,3,total_elem+nnew_elem),&
                        elnode_val_pressure(total_cracks,3,3,total_elem+nnew_elem),&
                        sol_trac(3,3,total_elem+nnew_elem),STAT=ierr)
    			if (ierr.ne.0) then
      				print*,'elnode and elldid: allocation request denied'
        			stop
    			endif
                !...Note: equno will be redo from beginning, no need to save
                allocate(equno(3,total_elem+nnew_elem),STAT=ierr)
    			if (ierr.ne.0) then
      				print*,'equno: allocation request denied'
        			stop
    			endif
            	!--------------------------------------------------
            	!...transfer back data from temporary varibles to global data
            	do i = 1,total_node
          			node_sys2user(i) = node_sys2user_temp(i)
            		node_id(i) = node_id_temp(i)
            		do j = 1,2
              			node_coor(j,i) = node_coor_temp(j,i)
            		enddo
                    !...5/26/09: debug this error (this may not need since we will set its value later depend on ellid and elnode_val)
                    do j = 1,3
                      	node_disp_val(j,i) = node_disp_val_temp(j,i)
                    enddo
        		enddo
        		do i = 1,total_elem
          			elem_sys2user(i) = elem_sys2user_temp(i)
            		elemid(i) = elemid_temp(i)
                        elem_crack_id(i) = elem_crack_id_temp(i)
            		elem_region(i) = elem_region_temp(i)
            		do j = 1,3
            			elnode(j,i) = elnode_temp(j,i)
                		elldid(j,i) = elldid_temp(j,i)
                        do m = 1,3
                          	elnode_val(m,j,i) = elnode_val_trac_temp(m,j,i)
                          	elnode_val_trac(m,j,i) = elnode_val_trac_temp(m,j,i)
                                sol_trac(m,j,i) = sol_trac_temp(m,j,i)

                        enddo
            		enddo
        		enddo
            	!--------------------------------------------------
            	!...update data of new nodes
                do i = 1,nnew_node
                  	node_sys2user(total_node + i) = node_sys2user_new(i)
                    node_id(total_node + i) = NORMAL
                    do j = 1,2
                      	node_coor(j,total_node + i) = node_coor_new(j,i)
                    enddo
                enddo
                !...Note: node_disp_val is just for nodes with displacement prescribed -> no need to update
                !...Note: no need to update rigid constraints since the old system node # doesn't change
                !--------------------------------------------------
                !...update data of new elements: this process need to rearrange all elements
                !after the last crack element of region k

                if (first_inter_elem.gt.0) then
                	!...shift all elements after the last crack element up to nnew_elem

!                	do i = first_inter_elem + nnew_elem,total_elem + nnew_elem
                	do i = total_elem + 1,total_elem + nnew_elem !...this is a temporary solution for the problem.3/2/12..!!!
!                	do i = total_elem + nnew_elem,first_inter_elem + nnew_elem !!!...changed direction of loop?..2/29/12
                  		elem_sys2user(i) = elem_sys2user(i - nnew_elem)
                    	elemid(i) = elemid(i - nnew_elem)
                    	elem_region(i) = elem_region(i - nnew_elem)
                    	do j = 1,3
                      		elnode(j,i) = elnode(j,i - nnew_elem)
                        	elldid(j,i) = elldid(j,i - nnew_elem)
                        	do m = 1,3
                          		elnode_val_trac(m,j,i) = elnode_val_trac(m,j,i - nnew_elem)
                                        elnode_val(m,j,i)=elnode_val_trac(m,j,i)
                                        sol_trac(m,j,i)=elnode_val(m,j,i)
                        	enddo
                    	enddo
                	enddo
                endif

                if (nnew_elem.gt.0) then !!!...temporary fix for the problem...3/2/12...!!!
                   do i=total_elem+1,total_elem+nnew_elem
                      do j = 1,3
                         elnode(j,i) = elnode(j,i - nnew_elem)
                         elldid(j,i) = elldid(j,i - nnew_elem)
                         do m = 1,3
                            elnode_val_trac(m,j,i) = elnode_val_trac(m,j,i - nnew_elem)
                            elnode_val(m,j,i)=elnode_val_trac(m,j,i)
                            sol_trac(m,j,i)=elnode_val(m,j,i) !...sol_trac is correct here!...3/2/12..!!!
                         enddo
                      enddo
                   enddo
                endif
                

                !...transfer data from temporary variables to global variables
                if (first_inter_elem.gt.0) then
                	index = first_inter_elem - 1
                elseif (first_inter_elem.eq.0) then
                	!...region k does not have interface elements, thus last ck elem is the last elem
                	index = total_elem
                endif
                
                do i = 1,nnew_elem
                  	elem_sys2user(index + i) = elem_sys2user_new(i)
                    elemid(index + i) = CQUAD
                    elem_region(index + i) = k
                    do j = 1,3
                      	elnode(j,index + i) = elnode_new(j,i)
                        !...5/29/09: add the case of traction applied on crack
                        elldid(j,index + i) = CTRAC
                !        if (pressure_flag.eq.NO) then
                !        	elldid(j,index + i) = CTRFREE
                !        	do m = 1,3
                !          		elnode_val(m,j,index+i) = 0.0
                !        	enddo
                !        elseif (pressure_flag.eq.YES) then
                !        	elldid(j,index + i) = CTRAC
                !            !defer updating elnode_val later (after recompute normals for all ck elem)
                !        else
                !          	print*,'pressure_flag out of range!'
                !            stop
                !        endif
                    enddo

                    !----------------------------------------------------------------
                    !---- Update elem_crack_id

                    do j=1,total_elem+nnew_elem
                       if ((elnode(1,index+i).eq.elnode(2,j)) .and. (elem_crack_id(index+i) .ne. 1) .and. (elem_crack_id(index+i) .ne. 2) .and. (elem_crack_id(index+i) .ne. 3)) then
                          elem_crack_id(index+i) = elem_crack_id(j)

                          print*
                          print*, 'The added element and crack id associated with node ', index+i,'is:'
                          print*, elem_sys2user(index+i), elem_crack_id(index+i)
                          print*, 'This crack id is taken from the adjacent element and crack id: '
                          print*, elem_sys2user(j), elem_crack_id(j)
                          print*

                       elseif((elnode(2,index+i).eq.elnode(1,j)) .and. (elem_crack_id(index+i) .ne. 1) .and. (elem_crack_id(index+i) .ne. 2) .and. (elem_crack_id(index+i) .ne. 3)) then
                          elem_crack_id(index+i) = elem_crack_id(j)

                          print*
                          print*, 'The added element and crack id associated with node ', index+i,'is:'
                          print*, elem_sys2user(index+i), elem_crack_id(index+i)
                          print*, 'This crack id is taken from the adjacent element and crack id: '
                          print*, elem_sys2user(j), elem_crack_id(j)
                          print*

                       endif
                    enddo
                 enddo
                !--------------------------------------------------
            	!...update total number of nodes and elements
            	total_node = total_node + nnew_node
            	total_elem = total_elem + nnew_elem
                total_cknode = total_cknode + nnew_node
                total_ckelem = total_ckelem + nnew_elem
                !--------------------------------------------------
                !...update region data
                total_region_node(k) = total_region_node(k) + nnew_node
                total_region_ck_node(k) = total_region_ck_node(k) + nnew_node
                total_region_elem(k) = total_region_elem(k) + nnew_elem
                total_region_ck_elem(k) = total_region_ck_elem(k) + nnew_elem
             endif !of nnew_node.gt.0



            !--------------------------------------------------
            !...deallocate temporary variables so that we can allocate it again for next region
            deallocate(node_sys2user_temp,node_coor_temp,node_id_temp,node_disp_val_temp,STAT=ierr)
            if (ierr.ne.0) then
               	print*,'node_sys2user_temp...: deallocation request denied'
                stop
            endif
            deallocate(elem_sys2user_temp,elemid_temp,elem_region_temp,elnode_temp,&
               	elldid_temp,elnode_val_temp,elnode_val_trac_temp,sol_trac_temp,STAT=ierr)
            if (ierr.ne.0) then
               	print*,'elem_sys2user_temp...: deallocation request denied'
                stop
            endif
            !--------------------------------------------------
            !...re-calculate traction on crack elements
            iee = zero_elem
            if (pressure_flag.ne.0) then
              	do i = 1,total_region_elem(k)
              		!...global element #
                	iee = iee + 1
                	!...skip process if iee is not a crack element
                	if ((eltype(elemid(iee)).ne.CTIP).and.(eltype(elemid(iee)).ne.CREGULAR)) cycle
                	!...get nodal coordinates of element iee
                	xele = node_coor(:,elnode(:,iee))
                	!...compute normals at nodes of element iee
                    call NodeNormal(elemid(iee),xele,normal_iee)

                    !delete
                    !print*,'elem: ',elem_sys2user(iee)
                    !do j = 1,3
                    !  	temp = dsqrt(normal_iee(1,j)*normal_iee(1,j)+normal_iee(2,j)*normal_iee(2,j))
                    !  	print'(2(f15.8,1x),f15.8)',(normal_iee(m,j),m=1,2),temp
                    !enddo

                    !...search for adjacent element on the side of node1 of iee
                    ielem_adj = 0
                    do j = 1,total_elem
                      	if ((eltype(elemid(j)).eq.CTIP).or.(eltype(elemid(j)).eq.CREGULAR)) then
                      		if (elnode(1,iee).eq.elnode(2,j)) then
                        		ielem_adj = j
                            	exit
                        	endif
                        endif
                    enddo
                    !...averaging normal at node1 of iee
                    if (ielem_adj.ne.0) then
                    	xele = node_coor(:,elnode(:,ielem_adj))
                    	call NodeNormal(elemid(ielem_adj),xele,normal_adj)
                        !delete
                        !print*,'adj elem on node 1:',elem_sys2user(ielem_adj)
                        !do j = 1,3
                        !  	temp = dsqrt(normal_adj(1,j)*normal_adj(1,j)+normal_adj(2,j)*normal_adj(2,j))
                      	!	print'(2(f15.8,1x),f15.8)',(normal_adj(m,j),m=1,2),temp
                    	!enddo
                        !...average normal at node1 of iee
                        do j = 1,2
                          	normal_avg(j) = normal_iee(j,1) + normal_adj(j,2)
                        enddo
                        dj = dsqrt(normal_avg(1)*normal_avg(1) + normal_avg(2)*normal_avg(2))
                        normal_iee(1,1) = normal_avg(1)/dj
                        normal_iee(2,1) = normal_avg(2)/dj
                    endif
                    !--------------------------------------------------
                    !...search for adjacent element on the side of node2 of iee
                    ielem_adj = 0
                    do j = 1,total_elem
                      	if ((eltype(elemid(j)).eq.CTIP).or.(eltype(elemid(j)).eq.CREGULAR)) then
                      		if (elnode(2,iee).eq.elnode(1,j)) then
                        		ielem_adj = j
                            	exit
                        	endif
                        endif
                    enddo
                    !...averaging normal at node2 of iee
                    if (ielem_adj.ne.0) then
                    	xele = node_coor(:,elnode(:,ielem_adj))
                    	call NodeNormal(elemid(ielem_adj),xele,normal_adj)
                        !delete
                        !print*,'adj elem on node 2:',elem_sys2user(ielem_adj)
                        !do j = 1,3
                        !  	temp = dsqrt(normal_adj(1,j)*normal_adj(1,j)+normal_adj(2,j)*normal_adj(2,j))
                      	!	print'(2(f15.8,1x),f15.8)',(normal_adj(m,j),m=1,2),temp
                    	!enddo
                        !...average normal at node1 of iee
                        do j = 1,2
                          	normal_avg(j) = normal_iee(j,2) + normal_adj(j,1)
                        enddo
                        dj = dsqrt(normal_avg(1)*normal_avg(1) + normal_avg(2)*normal_avg(2))
                        normal_iee(1,2) = normal_avg(1)/dj
                        normal_iee(2,2) = normal_avg(2)/dj
                    endif
                    !--------------------------------------------------
                    !...re-calculate traction applied on crack elements
                    do j = 1,NODE(elemid(iee))

                       elnode_val_trac(1,j,iee) = -1.0d0*sigma_x*normal_iee(1,j)
                          
                       elnode_val_trac(2,j,iee) = -1.0d0*sigma_y*normal_iee(2,j)

                       elnode_val_trac(3,j,iee) = 0.0d0

                       do ck = 1,total_cracks !..loop over all cracks
                          do m = 1,2 !..loop over first two  NDOF
                   
                             if(ck .eq. elem_crack_id(iee))then   !..if current crack is the crack that elem iee is on
                                elnode_val_pressure(ck,m,j,iee) = -1.0d0*pressure_value*normal_iee(m,j)  !input pressure used here..aje567
                             else
                                elnode_val_pressure(ck,m,j,iee) = 0.0d0
                             endif
                          enddo
                          elnode_val_pressure(ck,3,j,iee) = 0.0d0

                       enddo


                    enddo ! of j=1,NODE(elemid(iee))
                 enddo !of i=1,total_region_elem(k)
            endif !of (pressure_flag.eq.YES)
            !--------------------------------------------------
            !...print out traction for checking
            !iee = zero_elem
            !do i = 1,total_region_elem(k)
            !  	iee = iee + 1
            !    if ((eltype(elemid(iee)).ne.CTIP).and.(eltype(elemid(iee)).ne.CREGULAR)) cycle
            !    print*,'Element: ',elem_sys2user(iee)
            !    do j = 1,NODE(elemid(iee))
            !      	temp = dsqrt(elnode_val(1,j,iee)**2.0d0+elnode_val(2,j,iee)**2.0d0+elnode_val(3,j,iee)**2.0d0)
            !       	print'(3(f15.8,1x),f15.8)',(elnode_val(m,j,iee),m=1,3),temp
            !    enddo
            !enddo
            !--------------------------------------------------
         endif !of total_region_ck_elem(k).gt.0
      enddo !of k=1,total_region

    !--------------------------------------------------------------------------------------
    !!!............END OF CRACK PROPOGATION.........................!!!
    !--------------------------------------------------------------------------------------



    !--------------------------------------------------------------------------------------
    !!!...REPEAT PREPROCESSING TO SET UP NEXT LOAD STEP.............!!!
    !--------------------------------------------------------------------------------------


    !...5/26/09: found a bug here
    !...nodal prescribed displacements. Note: up to here, node_disp_val is not initiated yet!
    do i = 1,total_elem
      	do k = 1,3
        	if (elldid(k,i).eq.BDISP) then
            	do j = 1,NODE(elemid(i))
                	node_disp_val(k,elnode(j,i)) = elnode_val(k,j,i)
                enddo
            endif
        enddo
    enddo
    !--------------------------------------------------
    !...preprocessing for unknown,node_type,equno
    !...compute the equation number corresponding to each node of element
    iequ = 1
    do i = 1,total_node
      	if (node_id(i).eq.SBLP) then
        	!...surface breaking line node
            do k = 1,total_elem
              	do j = 1,NODE(elemid(k))
                	if (elnode(j,k).eq.i) then
                    	equno(j,k) = iequ !...initial no is coressponding to C+ node
                    endif
                enddo
            enddo
            !...reserve extra dof here for later assigning to C- node
            iequ = iequ + 2*NDOF
        elseif (node_id(i).eq.INTER) then
        	!...interface node: right now put here, but not figure out how it works in 2D with interface nodes
            do k = 1,total_elem
              	do j = 1,NODE(elemid(k))
                	if (elnode(j,k).eq.i) then
                    	equno(j,k) = iequ
                    endif
                enddo
            enddo
            !...reserve extra dof here for later assigning
            iequ = iequ + 2*NDOF
        elseif (node_id(i).eq.NORMAL) then
        	!...normal node
            do k = 1,total_elem
              	do j = 1,NODE(elemid(k))
                	if (elnode(j,k).eq.i) then
                    	!...5/13/09: check error
                        if (eltype(elemid(k)).eq.INTERFE) then
                          	print*,'error!, normal node: ',node_sys2user(i),'is on interface element: ',elem_sys2user(k)
                            stop
                        endif
                    	equno(j,k) = iequ
                    endif
                enddo
            enddo
            iequ = iequ + NDOF
        else
          	print*,'node_id is out of range, node:',node_sys2user(i)
            stop
         endif
    enddo
    !------------------------------
    !...determine equation number for extra nodes on surface breaking line
    do i = 1,total_elem
      	!...just look on regular elements on boundary
      	if (ELTYPE(elemid(i)).ne.BREGULAR) cycle
        equ_i = equno(:,i)
        do j = 1,NODE(elemid(i))
          	if (node_id(elnode(j,i)).eq.SBLP) then
            	!...search crack element sharing this node j
                do k = 1,total_elem
                  	if (ELTYPE(elemid(k)).eq.BREGULAR) cycle
                    do m = 1,NODE(elemid(k))
                      	if ((elnode(j,i).eq.elnode(m,k)).and.(j.eq.m)) then
                        	!boundary element i sharing SBLP node j with crack element k
                           	!...element k is associated with C-
                            equno(j,i) = equ_i(j) + NDOF
                        endif
                    enddo
                enddo
            endif
        enddo
    enddo
    !------------------------------
    !...compute the rigid point eqn
    do i = 1,total_rigid_constrain_point
      	do j = 1,total_elem
        	do k = 1,NODE(elemid(j))
            	if (elnode(k,j).eq.rigid_constrain_point(i)) then
                	rigid_constrain_point_eqn(i) = equno(k,j)
                endif
            enddo
        enddo
    enddo
    !------------------------------
    !...total dof
    !...5/3/09: fix the bug when last node is interface/SBLP node
    if ((node_id(total_node).eq.INTER).or.(node_id(total_node).eq.SBLP)) then
      	total_dof = MAXVAL(equno) + 2*NDOF -1
    elseif (node_id(total_node).eq.NORMAL) then
    	total_dof = MAXVAL(equno) + 2
    else
      	print*,'error from prep.f90: id of last node is out of range'
        stop
    endif
    !------------------------------
    !...determine the unknowns of each node
    do i = 1,total_node
      	if (node_id(i).eq.INTER) then
        	unknown(:,i) = DISP_TRAC
            cycle
        endif
        tcount = 0 !...counter of prescribed traction
        dcount = 0 !...counter of prescribed displacement
        !...search for element containing node i
        do j = 1,total_elem
          	do k = 1,NODE(elemid(j))
            	if (elnode(k,j).eq.i) then
                	do alpha = 1,NDOF
                    	select case (elldid(alpha,j))
                        	case (BTRAC,BTRFREE,CTRAC,CTRFREE)
                            	!...traction is prescribed
                                tcount(alpha) = tcount(alpha) + 1
                            case (BDISP)
                            	!...displacement is prescribed
                                dcount(alpha) = dcount(alpha) + 1
                            case default
                            	print*,'id of prescribed data (ellid) on element:',elem_sys2user (j),'is out of range'
                                stop
                        end select
                    enddo
                endif
            enddo
		enddo
        do alpha = 1,NDOF
          	if ((tcount(alpha).eq.0).and.(dcount(alpha).gt.0)) then
            	!...all elements sharing node i have prescribed displacement in alpha-dir
                unknown(alpha,i) = TRACTION
                node_type(alpha,i) = 0      !...not an edge node
            elseif ((dcount(alpha).eq.0).and.(tcount(alpha).gt.0)) then
            	!...all elements sharing node i have prescribed traction (or traction free) in alpha-dir
                unknown(alpha,i) = DISPLACEMENT
                node_type(alpha,i) = 0      !...not an edge node
            elseif ((dcount(alpha).gt.0).and.(tcount(alpha).gt.0)) then
                !...elements sharing node i some are traction prescribed, some are displacement prescribed
                unknown(alpha,i) = TRACTION
                node_type(alpha,i) = EDGE_NODE
                !...Note: even with an real edge node but if all elements sharing that node have either traction
                !or displacement prescribed, then it is still considered "not an edge node"
            else
              	print*,'something wrong with prescribed data of elements sharing node:',node_sys2user(i),'in direction:',alpha
                stop
            endif
        enddo
    enddo   !of i=1,total_node
	!------------------------------
    !...deallocate
    deallocate(dF,dF_traction,dF_pressure,dAK,dAK_traction,dAK_pressure,STAT=ierr)
    if (ierr.ne.0) then
      	print*,'dF,dAK: deallocation request denied!'
        stop
    endif
	!...allocate space for AK and F
    allocate(dF(total_dof),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'dF: allocation request denied!'
        stop
    endif
    allocate(dF_traction(total_dof),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'dF: allocation request denied!'
        stop
    endif
    allocate(dF_pressure(total_cracks,total_dof),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'dF_pressure: allocation request denied!'
        stop
    endif
    allocate(dAK(total_dof*(1+total_dof)/2),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'dAK: allocation request denied!'
        stop
    endif
    allocate(dAK_traction(total_dof*(1+total_dof)/2),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'dAK_traction: allocation request denied!'
        stop
    endif
    allocate(dAK_pressure(total_cracks,total_dof*(1+total_dof)/2),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'dAK_pressure: allocation request denied!'
        stop
    endif


	!------------------------------
    !...initialize dF and dAK
    dF = 0.0d0
    dF_traction = 0.0d0
    dF_pressure = 0.0d0

    dAK = 0.0d0
    dAK_traction = 0.0d0
    dAK_pressure = 0.0d0

	End Subroutine CrackGrowth

    !--------------------------------------------------


    Subroutine NodeNormal(id,xel,normal)
    !...subroutine to calculate normals at nodes of element
    INTEGER, INTENT(IN)		:: id
    REAL(KIND=DBL), INTENT(IN)	:: xel(2,3)
    REAL(KIND=DBL), INTENT(OUT)	:: normal(2,3)
    !...local variables
    INTEGER					:: i,j,k
    REAL(KIND=DBL)	:: xi,ps(3),dps(3),dxds(2),dj,tangent(2)
	do i = 1,3
		!...coordinates of master element
        select case (i)
        case (1)
           	xi = -1.0d0
        case (2)
           	xi = 1.0d0
        case (3)
           	xi = 0.0
        end select
        call reg_shape(id,xi,ps,dps)
        dxds = 0.0d0
        do j = 1,NODE(id)
           	do k = 1,2
               	dxds(k) = dxds(k) + xel(k,j)*dps(j)
            enddo
        enddo
        dj = dsqrt(dxds(1)*dxds(1) + dxds(2)*dxds(2))
        !...tangent vector
        do j = 1,2
           	tangent(j) = dxds(j)/dj
        enddo
        !...normal = tangent x e3
        normal(1,i) = tangent(2)
        normal(2,i) = -1.0d0*tangent(1)
    enddo
    End Subroutine NodeNormal


    !--------------------------------------------------
    Subroutine CrackVolume(crack,solution_disp,arc_length,crack_volume)
      !subroutine to calculate the arc length and volume of a crack

      !this subroutine is completely new by Andrew Erickson
      ! added 3/15/12
      !..updated 4/4/12 to handle multiple cracks

      INTEGER, INTENT(IN) :: crack
      REAL (KIND=DBL), INTENT(IN) :: solution_disp(:,:)
      REAL (KIND=DBL), INTENT(OUT) :: arc_length
      REAL (KIND=DBL), INTENT(OUT) :: crack_volume

      !...local variables

      INTEGER             :: nint, i, ii, j, k
      REAL(KIND=DBL), ALLOCATABLE	:: point(:),weight(:)
      REAL(KIND=DBL)	:: ps(3), dps(3),ps_u(3), dps_u(3), elem_length, elem_volume,&
                           opening_displacement, w_analytic, dS_squared, xi_coord




      !...set number of gaussian integration points
      nint=60

      !...allocate space for integration points and weights

      allocate(point(nint),weight(nint),STAT=ierr)
      if (ierr.ne.0) then
         print*,'point,weight: allocation request denied'
         stop
      endif

      !...get gaussian integration points and weights
      call GaussIntPoint(nint,point,weight)

      !...calculate arc length

      arc_length = 0       !...initialize arc length
      crack_volume = 0     !...initialize crack volume
      

      do ii=1,total_elem !..loop over all elements

         if(elem_crack_id(ii)==crack)then !only use elements that are on the specific crack

            elem_length = 0 !...initialize element arc length
            elem_volume = 0 !...initialize element volume

            do k=1,nint     !..loop over all integration points

               call reg_shape(elemid(ii),point(k),ps,dps)

               if(elemid(ii).eq.CQUAD)then !..if normal crack element use quadratic shape functions
                  call reg_shape(elemid(ii),point(k),ps_u,dps_u)
               elseif((elemid(ii).eq.CTIP1).or.(elemid(ii).eq.CTIP2))then  !...if tip element use square root shape functions
                  call tip_shape(elemid(ii),point(k),ps_u)
                  call tip_dshape(elemid(ii),point(k),dps_u)
               endif
               
               dS_squared = 0
               opening_displacement = 0
               xi_coord = 0
               w_analytic = 0

               do i = 1,NODE(elemid(ii))   !..loop over nodes of element ii

                  do j=1,NODE(elemid(ii))  !..loop over nodes of element ii

                     dS_squared = dS_squared + (node_coor(1,elnode(i,ii))*node_coor(1,elnode(j,ii))+&
                          node_coor(2,elnode(i,ii))*node_coor(2,elnode(j,ii)))*dps(i)*dps(j)
                  enddo
               

                  xi_coord = xi_coord + node_coor(1,elnode(i,ii))*ps(i) !...calculate x-coordinate of the gauss point

                  opening_displacement = opening_displacement + solution_disp(ii,i)*ps_u(i)  !...calculate opening displacement at the gauss point

               enddo
                        
               elem_length = elem_length + weight(k)*dsqrt(dS_squared)
               elem_volume = elem_volume + weight(k)*opening_displacement*dsqrt(dS_squared)



            enddo

            arc_length = arc_length + elem_length  !...calculate the total arc length of the crack
            crack_volume = crack_volume + elem_volume  !...calculate the total crack volume
         endif
      enddo




    	deallocate(point, weight, STAT=ierr)
    	if (ierr.ne.0) then
      		print*,'integration weight and points...: deallocation request denied!'
        	stop
    	endif


    End Subroutine CrackVolume
END SUBROUTINE Post
