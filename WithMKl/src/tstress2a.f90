Subroutine Tstress2a
!...subroutine to calculate T-stress by 2nd approach: weak form integrates on whole crack
!...difference with tstress2: sigma_u at SBL node is directly from solution (u+ + u-)
!...MUST have exact # of ck nodes: SBLP nodes are considered as crack nodes
	USE DefinitionConstant
    USE GlobalData
    USE ElementMatrix
    USE LinearSolver
    IMPLICIT NONE
    INTERFACE
        Subroutine assts_gl(isrc,ifld,sequ,fequ,el,gl)
        	INTEGER, INTENT(IN)			:: isrc,ifld
    		INTEGER, INTENT(IN)			:: sequ(:),fequ(:)
    		REAL(KIND=SELECTED_REAL_KIND(15)), INTENT(IN)	:: el(:,:)
    		REAL(KIND=SELECTED_REAL_KIND(15)), INTENT(OUT)	:: gl(:)
        End Subroutine assts_gl
        Subroutine assts_grhs(isrc,sequ,rhs,grhs)
        	INTEGER, INTENT(IN)			:: isrc
    		INTEGER, INTENT(IN)			:: sequ(:)
    		REAL(KIND=SELECTED_REAL_KIND(15)), INTENT(IN)	:: rhs(:)
    		REAL(KIND=SELECTED_REAL_KIND(15)), INTENT(OUT)	:: grhs(:)
        End Subroutine assts_grhs
    END INTERFACE

    !...local variables
	REAL(KIND=DBL)	:: xlo(2,3),xli(2,3),pso(3),dpso(3)
    REAL(KIND=DBL)	:: dxds(2)
    REAL(KIND=DBL)	:: el(3*NDOF,3*NDOF),rhs(3*NDOF)
    REAL(KIND=DBL)	:: eb(3*NDOF,3*NDOF),ebt(NDOF,3*NDOF)
    REAL(KIND=DBL)	:: eat(3*NDOF,3*NDOF),eatt(NDOF,3*NDOF)
    REAL(KIND=DBL)	:: xt,wt,detj
    INTEGER			:: i,j,k,l,m,n,ierr,total_ckdof,k0,jj,nn
    INTEGER			:: ido,idi,nno,nni,isrc,IGsrc,jfld,JGfld
    INTEGER			:: alpha,beta,SBLP_node,local_SBLP
    REAL(KIND=DBL), ALLOCATABLE	:: gl(:),grhs(:)
    LOGICAL			:: success,ckregion_allocated
    INTEGER, ALLOCATABLE		:: ckelem_sys2user(:),ckelem_id(:),ckelem_connect(:,:)
    INTEGER, ALLOCATABLE		:: cknode_sys2user(:),ckequno(:,:)
	REAL(KIND=DBL), ALLOCATABLE	:: cknode_coor(:,:)
    INTEGER						:: nelem,nnode,idone,itemp,ck_IGsrc
    REAL(KIND=DBL)	:: pii(NDOF)
    REAL(KIND=DBL)	:: sum_uk_SBLP(NDOF)
	INTEGER			:: I_sign
    REAL(KIND=DBL)	:: tangent(2),epsilon_s
    !--------------------------------------------------
    print*,'--------------------------------------------------'
	print*,'Tstress calculation by 2ND APPROACH, integrate on WHOLE CRACK'
    	print*,'(displacement of SBLP directly from solution)'
    !...initiate number of element of current region
	k0 = 0
    !...indicator of allocation for variable to store ck data for current region
    ckregion_allocated = .false.

	I_sign = 1
    do k = 1,total_region
      	!--------------------------------------------------------------------------------
        !...GET CRACK DATA INTO CK VARIABLES
        nelem = 0
        nnode = 0
		!...total ck elements of region k
        nelem = total_region_ck_elem(k)
        !...total ck nodes of region k
        nnode = total_region_ck_node(k)
        !----------------------------------------------------------------------
        !...exit whole process if region k does not have crack
        if ((nelem.eq.0).or.(nnode.eq.0)) then
          	k0 = k0 + total_region_elem(k)
            cycle
        endif
        !----------------------------------------------------------------------
        !...allocate space for gl base on total number of crack node
    	total_ckdof = NDOF*nnode
        !...allocate space to store crack data
        if (ckregion_allocated) then
          	deallocate(gl,grhs,ckelem_sys2user,ckelem_id,ckelem_connect,&
            	cknode_sys2user,cknode_coor,ckequno,STAT=ierr)
            if (ierr.ne.0) then
              	print*,'tstress2a.f90: deallocation requested denied'
                stop
            endif
        endif
    	allocate (gl(total_ckdof*(total_ckdof+1)/2),STAT=ierr)
    	if (ierr.ne.0) then
      		print*,'tstress2a.f90: gl allocation request denied'
        	stop
    	endif
    	!...allocate space for grhs
    	allocate (grhs(total_ckdof),STAT=ierr)
    	if (ierr.ne.0) then
      		print*,'tstress2a.f90: grhs allocation request denied'
        	stop
    	endif
        allocate(ckelem_sys2user(nelem),ckelem_id(nelem),&
        	ckelem_connect(3,nelem),STAT=ierr)
        if (ierr.ne.0) then
          	print*,'tstress2a.f90: ckelem allocation request denied'
            stop
        endif
        allocate(cknode_sys2user(nnode),cknode_coor(2,nnode),STAT=ierr)
        if (ierr.ne.0) then
          	print*,'tstress2a.f90: cknode allocation request denied'
            stop
        endif
        allocate(ckequno(3,nelem),STAT=ierr)
        if (ierr.ne.0) then
          	print*,'tstress2a.f90: ckequno allocation request denied'
            stop
        endif
        !...indicate that variables to store crack data is already allocated,
        !so that next region will need to deallocate them before re-allocate
        ckregion_allocated = .true.
        !----------------------------------------------------------------------
        !...initialize
        gl = 0.0d0
    	grhs = 0.0d0
        I_sign = 1
        !...reset the counter for number of crack elements/nodes
		nelem = 0 !counter for number of crack elements
		nnode = 0 !counter for number of crack nodes
		do isrc = 1,total_region_elem(k)
        	!...global element number
        	i = isrc + k0
			if ((eltype(elemid(i)).ne.CTIP).and.(eltype(elemid(i)).ne.CREGULAR)) cycle
			nelem=nelem+1
			ckelem_sys2user(nelem) = elem_sys2user(i)
			ckelem_id(nelem) = elemid(i)
			do j=1,node(elemid(i))
				ckelem_connect(j,nelem) = elnode(j,i)
				idone = 0
                !...loop over stored node data to search if it was stored yet
				do n = 1,nnode
					if (node_sys2user(elnode(j,i)).eq.cknode_sys2user(n)) then
						idone = 1
						exit
					endif
				enddo
				if (idone.eq.0) then
					!...elnode(j,i) has not stored to nodal data of crack yet
					nnode = nnode + 1
					cknode_sys2user(nnode) = node_sys2user(elnode(j,i))
                    do n = 1,2
						cknode_coor(n,nnode) = node_coor(n,elnode(j,i))
					enddo
				endif
			enddo !...of j
		enddo !...of isrc
        !...defensive programming
        if (nelem.ne.total_region_ck_elem(k)) then
          	print*,'tstress2a.f90 error: number of crack elements is wrong for region:',k
            stop
        endif
        if (nnode.ne.total_region_ck_node(k)) then
          	print*,'tstress2a.f90 error: number of crack nodes is wrong for region:',k
            stop
        endif
        !...translate system# of ck_elnode(:,:) which is currently refer to node_coor(:,:) to system# of ck_node_coor(:,:)
		do i = 1,nelem !...loop over crack elements
			do j = 1,node(ckelem_id(i))
				!...get the user# of node j (currently still corresponding to system# of global data)
				itemp = node_sys2user(ckelem_connect(j,i))
				!...search the position (system#) of itemp in cknode_sys2user(:)
                idone = 0
				do n = 1,nnode
					if (itemp.eq.cknode_sys2user(n)) then
                    	idone = 1
                        exit
                    endif
				enddo
				!...replace system# to new system# which is corresponding to crack data
                if (idone.eq.1) then
					ckelem_connect(j,i) = n
                else
                  	print*,'tstress2a.f90: something wrong with crack data, user node #:',itemp
                    stop
                endif
			enddo
		enddo
        !...determine eqn number of crack data
        itemp = 1
        !...loop over crack nodes
        do i = 1,nnode
            do j = 1,nelem
              	do n = 1,node(ckelem_id(j))
                	if (ckelem_connect(n,j).eq.i) then
                    	ckequno(n,j) = itemp
                    endif
                enddo
            enddo
            !...all nodes on crack have NDOF for each node, nothing special such as SBLP...
            itemp = itemp + NDOF
        enddo
        !--------------------------------------------------------------------------------
        !...START CALCULATING T-STRESS
    	do isrc = 1,total_region_elem(k)
        	!...global element # of outer element
        	IGsrc = isrc + k0
        	!...only continue the process if element is a crack element
    		if ((eltype(elemid(IGsrc)).ne.CTIP).and.(eltype(elemid(IGsrc)).ne.CREGULAR)) cycle
            !...id of the outer element (crack element)
            ido = elemid(IGsrc)
        	!...number of nodes of outer element
        	nno = NODE(ido)
        	!...nodal coordinates of outer element
        	do m = 1,nno
        		do n = 1,2
            		xlo(n,m) = node_coor(n,elnode(m,IGsrc))
            	enddo
        	end do
            !...check if IGsrc has SBL node
            SBLP_node = 0
            local_SBLP = 0
            do m = 1,nno
              	if (node_id(elnode(m,IGsrc)).eq.SBLP) then
                	SBLP_node = elnode(m,IGsrc)
                    local_SBLP = m
                    exit
                endif
            enddo
            sum_uk_SBLP = 0.0d0
            if (SBLP_node.ne.0) then
              	!...initialize sum_uk_SBLP
                sum_uk_SBLP = 0.d0
              	!...search element sharing SBLP_node
                do m = 1,total_region_elem(k)
                  	i = m + k0
                    do n = 1,node(elemid(i))
                      	if ((elnode(n,i).eq.SBLP_node).and.(i.ne.IGsrc)) then
                        	do alpha = 1,NDOF
                            	sum_uk_SBLP(alpha) = sum_uk_SBLP(alpha) + sol_disp(alpha,n,i)
                            enddo
                            exit
                        endif
                    enddo
                enddo
                !debug
                !print*,'value of sum of displacement at SBLP node (from solution):',node_sys2user(SBLP_node)
                !do alpha=1,3
                !  	print*,sum_uk_SBLP(alpha)
                !enddo
            endif
        	!--------------------------------------------------
        	!...compute matrix el of element IGsrc
			el = 0.0d0
			!...begins the integration loop: use number of Gauss-point same as load vector
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
        		do j = 1,nno
                	jj = NDOF*(j-1) + 1
            		do n = 1,nno
                    	nn = NDOF*(n-1) + 1
                		el(jj,nn) = el(jj,nn) + pso(j)*pso(n)*detj*wt
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
            !...get element number associated with IGsrc in ck data
            ck_IGsrc = 0
            itemp = elem_sys2user(IGsrc)
            do i = 1,nelem
              	if (itemp.eq.ckelem_sys2user(i)) then
                	ck_IGsrc = i
                    exit
                endif
            enddo
            if (ck_IGsrc.eq.0) then
              	print*,'tstress2a.f90: something wrong with ck elem data'
                stop
            endif
            !...assemble to global gl using ck data
            call assts_gl(nno,nno,ckequno(:,ck_IGsrc),ckequno(:,ck_IGsrc),el,gl)
            !--------------------------------------------------
            !...loop over inner element to form the RHS of element IGsrc
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
            	!...compute ebtip and eg/ehtip if source element is a tip element
            	if (eltype(ido).eq.CTIP) then
                    !...depend on where field element is, calculate needed matrices
              		if ((eltype(idi).eq.CTIP).or.(eltype(idi).eq.CREGULAR)) then
                  		call eattip(k,IGsrc,JGfld,xlo,xli,eatt)
                	else
                  		call eattip(k,IGsrc,JGfld,xlo,xli,eatt)
                    	call ebtip(k,IGsrc,JGfld,xlo,xli,ebt)
                	endif
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
                				CASE (BDISP,BTRAC,BTRFREE,NOLOAD)
                                	if ((elem_region(JGfld).ne.k).and.(elemid(JGfld).eq.IQUAD)) then
                                    	I_sign = -1
                                    else
                                      	I_sign = 1
                                    endif
                                	rhs(NDOF*(m-1)+alpha) = rhs(NDOF*(m-1)+alpha) - &
                                    	I_sign*eb(NDOF*(m-1)+alpha,NDOF*(n-1)+beta)*sol_trac(beta,n,JGfld) + &
                                        I_sign*eat(NDOF*(m-1)+alpha,NDOF*(n-1)+beta)*sol_disp(beta,n,JGfld)
                                    !...if outer element is a tip element, then need to add eb/eat at tip node
                                    if ((ido.eq.CTIP1).and.(m.eq.1)) then
                                      	rhs(alpha) = rhs(alpha) &
                                        	- I_sign*ebt(alpha,NDOF*(n-1)+beta)*sol_trac(beta,n,JGfld) &
                                        	+ I_sign*eatt(alpha,NDOF*(n-1)+beta)*sol_disp(beta,n,JGfld)
                                    elseif ((ido.eq.CTIP2).and.(m.eq.2)) then
                                      	rhs(NDOF+alpha) = rhs(NDOF+alpha) &
                                        	+ I_sign*ebt(alpha,NDOF*(n-1)+beta)*sol_trac(beta,n,JGfld) &
                                        	- I_sign*eatt(alpha,NDOF*(n-1)+beta)*sol_disp(beta,n,JGfld)
                                    endif
                                CASE (CTRAC,CTRFREE)
                                	rhs(NDOF*(m-1)+alpha) = rhs(NDOF*(m-1)+alpha) + &
                                        eat(NDOF*(m-1)+alpha,NDOF*(n-1)+beta)*sol_disp(beta,n,JGfld)
                                    !...if outer element is a tip element, then need to add eb/eat at tip node
                                    if ((ido.eq.CTIP1).and.(m.eq.1)) then
                                      	rhs(alpha) = rhs(alpha) &
                                        	+ eatt(alpha,NDOF*(n-1)+beta)*sol_disp(beta,n,JGfld)
                                    elseif ((ido.eq.CTIP2).and.(m.eq.2)) then
                                      	rhs(NDOF+alpha) = rhs(NDOF+alpha) &
                                        	- eatt(alpha,NDOF*(n-1)+beta)*sol_disp(beta,n,JGfld)
                                    endif
                                CASE default
                                	print*,'tstress2a.f90: current version does not support type of loading: ',elldid(beta,JGfld)
                                    stop
                                END SELECT
                            enddo
                        enddo
                    enddo
                enddo
            enddo   !...of jfld over total element of region k
            !debug
            !print*,'rhs of element before adding nodal value',elem_sys2user(IGsrc)
            !do i = 1,9
            !  	print*,rhs(i)
            !enddo
            !...add sigma_u when IGsrc has SBLP node
            if (SBLP_node.ne.0) then
              	do alpha = 1,NDOF
              		if (local_SBLP.eq.1) then
                    	rhs(alpha) = rhs(alpha) - 0.5d0*sum_uk_SBLP(alpha)
                    elseif (local_SBLP.eq.2) then
                    	rhs(NDOF+alpha) = rhs(NDOF+alpha) + 0.5d0*sum_uk_SBLP(alpha)
                    else
                      	print*,'tstress2a: wrong local node of SBLP'
                        stop
                    endif
                enddo
            endif
            !debug
            !print*,'rhs of element',elem_sys2user(IGsrc)
            !do i = 1,9
            !  	print*,rhs(i)
            !enddo
            !...assemble eb and eat to global RHS
            call assts_grhs(nno,ckequno(:,ck_IGsrc),rhs,grhs)
        enddo !of isrc over total elements of region k
        !--------------------------------------------------------------------------------
    	!...solve using preconditioner solver
    	success = dPCGSolver(total_ckdof,gl,grhs)
    	!...get value of pi at tip
    	pii = 0.d0
    	print*,'Value of Pii at tip node:'
    	do i = 1,nelem
      		if ((ckelem_id(i).ne.CTIP1).and.(ckelem_id(i).ne.CTIP2)) cycle
            !...get du_i/ds from solution
        	if (ckelem_id(i).eq.CTIP1) then
          		do j = 1,NDOF
          			pii(j) = grhs(ckequno(1,i)+j-1)
          		enddo
            	!print*,'Node        Pii1        Pii2          Pii3'
            	!print'(i5,3(f15.10,2x))',cknode_sys2user(ckelem_connect(1,i)),(pii(j),j=1,NDOF)
        	elseif (ckelem_id(i).eq.CTIP2) then
          		do j = 1,NDOF
          			pii(j) = grhs(ckequno(2,i)+j-1)
          		enddo
            	!print*,'Node        Pii1        Pii2          Pii3'
            	!print'(i5,3(f15.10,2x))',cknode_sys2user(ckelem_connect(2,i)),(pii(j),j=1,NDOF)
        	endif
            !...compute tangent vector
            nno = NODE(ckelem_id(i))
            ido = ckelem_id(i)
            do m = 1,nno
        		do n =1,2
            		xlo(n,m)=cknode_coor(n,ckelem_connect(m,i))
            	enddo
        	end do
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
        	detj = dsqrt(dxds(1)*dxds(1) + dxds(2)*dxds(2))
        	if (detj.lt.SMALL_NUM) then
          		print*,'tstress2a.f90: jacobian is almost zero'
            	stop
        	endif
        	!...tangential vector
        	do m = 1,2
          		tangent(m) = dxds(m)/detj
        	enddo
            !...component of pi() on tangential direction
            epsilon_s = 0.0d0
            do m = 1,2
              	epsilon_s = epsilon_s  + pii(m)*tangent(m)
            enddo
            !...print out
            if (ido.eq.CTIP1) then
              	itemp = cknode_sys2user(ckelem_connect(1,i))
            elseif (ido.eq.CTIP2) then
            	itemp = cknode_sys2user(ckelem_connect(2,i))
            endif
            print*,'Strain epsilon_11 in local (tip) coordinate'
            print'(a7,i5,a5,i5,1x,f15.10)','Element',ckelem_sys2user(i),'Node',itemp,epsilon_s
      	enddo
        k0 = k0 + total_region_elem(k)
    enddo !of k
    
End Subroutine Tstress2a
!--------------------------------------------------