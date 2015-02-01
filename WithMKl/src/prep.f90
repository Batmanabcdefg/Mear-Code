Subroutine DataPreparation
!...My version !!!!	
!...in this subroutine, you basically read the data from input file and echo the input data
USE DefinitionConstant		
!...defined in global.f90
USE GlobalData				
!...defined in global.f90
USE Remesh 
!...for matrix inverse
IMPLICIT NONE
!...in this subroutine, there are calls to 5 subroutines, namely, PressureTractions,InitGaussIntegration,InputDataEcho,NodeNormal & reg_shape
!...PressureTractions,InputDataEcho,NodeNormal are defined in this file whereas reg_shape is defined in shape.f90

	!...local variables
    INTEGER		:: iono = 5 !...unit number of input
	INTEGER		:: i,j,k,m,ierr,alpha,equ_i(3),ntitle,kk
    INTEGER		:: integration_flag,iequ
    INTEGER, ALLOCATABLE	:: tcount(:),dcount(:)
    INTEGER		:: itemp(3),icheck
    
	!...variables for rotating geometry and loading
!    REAL(KIND=DBL)	:: rotate_angle,coor_temp(2)
    REAL(KIND=DBL),ALLOCATABLE	:: elnode_val_temp(:)	  
	!...what goes inside is NDOF
    REAL(KIND=DBL),ALLOCATABLE	:: elnode_val_trac_temp(:)	  
	!...what goes inside is again NDOF
    REAL(KIND=DBL),ALLOCATABLE	:: elnode_val_pressure_temp(:,:)    
	!...what goes inside is total_cracks and NDOF

!    INTEGER			:: rotate_flag
    !...12/30/09: add option to input engineering constants for orthotropy
!    INTEGER			:: mat_flag
!    REAL(KIND=DBL)	:: E1,E2,E3,nu23,nu31,nu12,G23,G31,G12  !engineering constants for orthotropic material
!    REAL(KIND=DBL)	:: compliance(6,6),moduli(6,6),moduli_reduce(6,6)
	!5/27/09: add variables for output files
    CHARACTER(32)	:: techplot_file,abaqus_file,input_file, pressure_file    
	!...all these files take 32 characters in 1 line
    !...7/13/09: to store number of independent constants of material
    INTEGER		:: matconst													  
	!...stores the number of independent material constants
    !...12/8/09: tolerance for locking of rigid displacements of pure-traction inner bound.
    REAL(KIND=DBL),PARAMETER	:: tol=1.0d-6
	
	!...open the input file and associate it with unit iono = 5
	!open (iono,file='example1.inp',status='unknown')
    
	!----------------------------------------------------------------------
    !...degree of freedom per node, also indicator of media type
    read(UNIT=iono,FMT=*)NDOF												  
	!...number of dofs per node

    !----------------------------------------------------------------------
    !...flag for displaying data/result
    read(UNIT=iono,FMT=*) model_flag,boundary_flag,sif_flag
	!...model flag =1 if you want to print out input data in output file
	!...boundary flag=1 if you want to print out results of displacements and stresses on boundary
	!...sif flag=1 if you want to print out results of sifs at crack tips
    !...class of problem
    !read(UNIT=iono,FMT=*) class_problem
    !if (class_problem.eq.2) then
      	!...anisotropy: material coordinates are different with geometry coordinates ?(1-yes;0-no)
      	read(UNIT=iono,FMT=*) rotate_flag		
		!...this flag again explained in page 7 of FADD2D manual
        if (rotate_flag.eq.1) then
       	   	!...angle from 1-axis of geometry to 1-axis of material, positive ccw, [-180,180] deg.
          	read(UNIT=iono,FMT=*)rotate_angle
            !...convert to radian
            rotate_angle = rotate_angle*PI/180.0d0		
			!...if material coordinates are different from geometry coordinates, use rotate_angle
        endif
    !...total material
    read(UNIT=iono,FMT=*) total_material                      
	!...total_material is the number of materials
    read(UNIT=iono,FMT=*) total_cracks 
	!...added total number of cracks to input file 4/3/12..aje567
	read(UNIT=iono,FMT=*) ht			   
	!...added height of cracks, Saumik
    
	!------------------------------------------------------------
    !...calculate number of independent material constants
    if ((NDOF.ge.3).and.(NDOF.le.5)) then
    	matconst = (9*NDOF*NDOF - 15*NDOF)/2 + 3
    elseif (NDOF.eq.2) then
    	matconst = 15
    elseif (NDOF.eq.1) then
    	matconst = 6					   
		!...number of material constants that need to be specified to make the constitutive matrix
    else
      	print*,'prep.f90: current version does not support NDOF=',NDOF
        stop
    endif
    !------------------------------------------------------------
    !...note: AllocateStatus is an integer variable. AllocateStatus takes the value 0 if allocation is successful or some other machine dependent 
	!...value if there is insufficient memory
    !------------------------------------------------------------
	!...matconst refers to total # of material constants and total_material refers to total # of materials
	
	allocate(material(matconst,total_material),STAT=ierr)            
    if (ierr.ne.0) then
       	print*,'material: allocation request denied'
        stop
    endif
    allocate(material_type(total_material),STAT=ierr)
    if (ierr.ne.0) then
       	print*,'material_type: allocation request denied'
        stop
    endif
    allocate(material_name(total_material),STAT=ierr)
    if (ierr.ne.0) then
       	print*,'material_name: allocation request denied'			 
		!...allocation request is denied if ierr is not equal to 0
        stop
    endif															 !material(:,:),material_type(:,:),material_name(:,:)
   
    !------------------------------------------------------------
    allocate(elnode_val_temp(NDOF),STAT=ierr)
    if (ierr.ne.0) then
       	print*,'elnode_val_temp: allocation request denied'
        stop
    endif
    !------------------------------------------------------------
    !...2/28/12..added temp variables for both tractions and pressure..aje567
    allocate(elnode_val_trac_temp(NDOF),STAT=ierr)
    if (ierr.ne.0) then
       	print*,'elnode_val_trac_temp: allocation request denied'
        stop
    endif
    !------------------------------------------------------------
    !..added functionality for multiple cracks...aje567
    allocate(elnode_val_pressure_temp(total_cracks,NDOF),STAT=ierr)
    if (ierr.ne.0) then
       	print*,'elnode_val_pressure_temp: allocation request denied'
        stop
    endif															
	!...elnode_val_temp(:),elnode_val_trac_temp(:),elnode_val_pressure_temp(:)

    !------------------------------------------------------------
    allocate(rigid_constrain_dir(max_rigid_constrain_point,NDOF),STAT=ierr)	
	!...first index is the constraint # and second index is DOF #
    if (ierr.ne.0) then
       	print*,'rigid_constrain_dir: allocation request denied'
        stop
    endif
    !------------------------------------------------------------
    allocate(BI1(NDOF,NDOF),BI2(NDOF,NDOF),STAT=ierr)
	!...changed this, Saumik	  
    if (ierr.ne.0) then
       	print*,'Kernels 1: allocation request denied'
        stop
    endif
    allocate(EI1(NDOF,NDOF),EI2(NDOF,NDOF),STAT=ierr)
	!...changed this, Saumik
    if (ierr.ne.0) then
       	print*,'Kernels 2: allocation request denied'
        stop
    endif														   
    !------------------------------------------------------------
    do i=1,total_material						
	!...this loop fills up the constitutive matrix, first index goes from 1 to matconst, second index goes from 1 to total_material
       	read(UNIT=iono,FMT='(a32)')material_name(i)				   
		!...input the material name
        read(UNIT=iono,FMT=*)material_type(i)					   
		!...input the material type
        select case (material_type(i))
        case (1)												   
        !...isotropy
        	if ((NDOF.eq.2).or.(NDOF.eq.3)) then
            	!...elastic plane strain/generalized plane strain: E & nu
                read(UNIT=iono,FMT=*)material(1,i),material(2,i)
            elseif (NDOF.eq.4) then
            	!...piezoelectric media
                read(UNIT=iono,FMT=*)material(1,i),material(2,i),material(3,i)
            else
                print*,'prep.f90: not develop isotropy of dof =',NDOF
                stop
            endif
        case (2)												   
		!...material type 2 for cubic
        !...cubic
        	if (NDOF.eq.3) then
            !...elastic media
                read(UNIT=iono,FMT=*)material(1,i),material(2,i),material(3,i)
            else
                print*,'prep.f90: not develop cubic material of dof =',NDOF
                stop
            endif
        case (3)
        !...transversely isotropy								   
		!...material type 3 for transversely isotropic
        	if (NDOF.eq.2) then
            	!...plane strain of elastic media
                !...still need C13 for calculating C55: wrong! no need to have C13 (1/1/2010)
            	read(UNIT=iono,FMT=*)material(1,i),material(2,i),material(3,i),&
                	material(4,i),material(5,i)
        	elseif (NDOF.eq.3) then
            !...elastic media
                read(UNIT=iono,FMT=*)material(1,i),material(2,i),material(3,i),&
                	material(4,i),material(5,i)
            elseif (NDOF.eq.4) then
            !...piezoelectric media
            	read(UNIT=iono,FMT=*)material(1,i),material(2,i),material(3,i),&
                	material(4,i),material(5,i),material(6,i),material(7,i),&
                    material(8,i),material(9,i),material(10,i)
            elseif (NDOF.eq.5) then
            !...magnetoelectroelastic media: need 17 constants
            	read(UNIT=iono,FMT=*)material(1,i),material(2,i),material(3,i),&
                	material(4,i),material(5,i),material(6,i),material(7,i),&
                    material(8,i),material(9,i),material(10,i),material(11,i),&
                    material(12,i),material(13,i),material(14,i),material(15,i),&
                    material(16,i),material(17,i)
            else
                print*,'prep.f90: not develop transversely isotropy of dof =',NDOF
                stop
            endif
    	case (4)												   
		!...material type 4 for orthotropic
                !...input moduli directly
        		if (NDOF.eq.2) then
            	!...plane strain of elastic media
            		read(UNIT=iono,FMT=*)material(1,i),material(2,i),material(3,i),&
                   	material(4,i),material(5,i),material(6,i)
        		elseif (NDOF.eq.3) then
            	!...elastic media
                	read(UNIT=iono,FMT=*)material(1,i),material(2,i),material(3,i),&
                   	material(4,i),material(5,i),material(6,i),material(7,i),&
                    material(8,i),material(9,i)
            	elseif (NDOF.eq.4) then
                	read(UNIT=iono,FMT=*)material(1,i),material(2,i),material(3,i),&
                   	material(4,i),material(5,i),material(6,i),material(7,i),&
                    material(8,i),material(9,i),material(10,i),material(11,i),&
                    material(12,i),material(13,i),material(14,i),material(15,i),&
                    material(16,i),material(17,i)
                else
                	print*,'prep.f90: not develop orthotropy of dof =',NDOF
                	stop
            	endif
        case (5)											
		!...material type 5 for monoclinic
        !...monoclinic
        	if (NDOF.eq.2) then
            !...plane strain of elastic media: 9 independent constants
            	read(UNIT=iono,FMT=*)material(1,i),material(2,i),material(3,i),&
                  	material(4,i),material(5,i),material(6,i),material(7,i),&
                    material(8,i),material(9,i)
        	elseif (NDOF.eq.3) then
            !...elastic media: 13 independent constants
                read(UNIT=iono,FMT=*)material(1,i),material(2,i),material(3,i),&
                  	material(4,i),material(5,i),material(6,i),material(7,i),&
                    material(8,i),material(9,i),material(10,i),material(11,i),&
                    material(12,i),material(13,i)
            else
                print*,'prep.f90: not develop monoclinic of dof =',NDOF
                stop
            endif
        case (6)											
		!...material type 6 for general anisotropic
       	!...general anisotropic
        	if (NDOF.eq.1) then
            	read(UNIT=iono,FMT=*)material(1,i),material(2,i),material(3,i),&
                  	material(4,i),material(5,i),material(6,i)
        	elseif (NDOF.eq.2) then
            !...plane strain of elastic media: 15 independent constants
            	read(UNIT=iono,FMT=*)material(1,i),material(2,i),material(3,i),&
                  	material(4,i),material(5,i),material(6,i),material(7,i),&
                    material(8,i),material(9,i),material(10,i),material(11,i),&
                    material(12,i),material(13,i),material(14,i),material(15,i)
        	elseif (NDOF.eq.3) then
            !...elastic media: 21 independent constants
                read(UNIT=iono,FMT=*)material(1,i),material(2,i),material(3,i),&
                  	material(4,i),material(5,i),material(6,i),material(7,i),&
                    material(8,i),material(9,i),material(10,i),material(11,i),&
                    material(12,i),material(13,i),material(14,i),material(15,i),&
                    material(16,i),material(17,i),material(18,i),material(19,i),&
                    material(20,i),material(21,i)
            elseif (NDOF.eq.4) then
            !...piezoelectric media: 45 independent constants
            	read(UNIT=iono,FMT=*)material(1,i),material(2,i),material(3,i),&
                  	material(4,i),material(5,i),material(6,i),material(7,i),&
                    material(8,i),material(9,i),material(10,i),material(11,i),&
                    material(12,i),material(13,i),material(14,i),material(15,i),&
                    material(16,i),material(17,i),material(18,i),material(19,i),&
                    material(20,i),material(21,i),material(22,i),material(23,i),&
                    material(24,i),material(25,i),material(26,i),material(27,i),&
                    material(28,i),material(29,i),material(30,i),material(31,i),&
                    material(32,i),material(33,i),material(34,i),material(35,i),&
                    material(36,i),material(37,i),material(38,i),material(39,i),&
                    material(40,i),material(41,i),material(42,i),material(43,i),&
                    material(44,i),material(45,i)
             elseif (NDOF.eq.5) then
             !...magnetoelectroelastic media
             	read(UNIT=iono,FMT=*)material(1,i),material(2,i),material(3,i),&
                  	material(4,i),material(5,i),material(6,i),material(7,i),&
                    material(8,i),material(9,i),material(10,i),material(11,i),&
                    material(12,i),material(13,i),material(14,i),material(15,i),&
                    material(16,i),material(17,i),material(18,i),material(19,i),&
                    material(20,i),material(21,i),material(22,i),material(23,i),&
                    material(24,i),material(25,i),material(26,i),material(27,i),&
                    material(28,i),material(29,i),material(30,i),material(31,i),&
                    material(32,i),material(33,i),material(34,i),material(35,i),&
                    material(36,i),material(37,i),material(38,i),material(39,i),&
                    material(40,i),material(41,i),material(42,i),material(43,i),&
                    material(44,i),material(45,i),material(46,i),material(47,i),&
                    material(48,i),material(49,i),material(50,i),material(51,i),&
                    material(52,i),material(53,i),material(54,i),material(55,i),&
                    material(56,i),material(57,i),material(58,i),material(59,i),&
                    material(60,i),material(61,i),material(62,i),material(63,i),&
                    material(64,i),material(65,i),material(66,i),material(67,i),&
                    material(68,i),material(69,i),material(70,i),material(71,i),&
                    material(72,i),material(73,i),material(74,i),material(75,i),&
                    material(76,i),material(77,i),material(78,i)
             else
                print*,'prep.f90: not develop general anisotropy of dof =',NDOF
                stop
            endif
        case default
           	print*,'not support material type:',material_type(i)
            stop
        end select
    enddo
    !------------------------------------------------------------
    read(UNIT=iono,FMT=*) total_region                
	!...total number of regions,important in multi region problem
    
    allocate(region_mat_no(total_region),region_name(total_region),total_region_node(total_region),&
    	total_region_elem(total_region),total_region_bd_elem(total_region),&
        total_region_ck_elem(total_region),total_region_bd_node(total_region),&
        total_region_ck_node(total_region),STAT=ierr) 
		!...all region related variables
    if (ierr.ne.0) then
      	print*,'region variables: allocation request denied'
        stop
    endif											  
	!...all the region related variables are allocated space here
    
	do i=1,total_region
      	!...read name of region i
      	read(UNIT=iono,FMT='(a32)'),region_name(i)	  
		!...read the region name
        !...read material # of region i
        read(UNIT=iono,FMT=*) region_mat_no(i)		  
		!...read the material number of region
    enddo

    !------------------------------------------------------------
    do i=1,total_region
      	read(UNIT=iono,FMT=*)total_region_elem(i),total_region_bd_elem(i),&
        	total_region_ck_elem(i)	
			!...total number of elements, total number of boundary elements and crack elements in each of the regions
    enddo
    
	do i=1,total_region
      	read(UNIT=iono,FMT=*)total_region_node(i),total_region_bd_node(i),&
        	total_region_ck_node(i)	
			!...total number of nodes, total number of boundary nodes and crack nodes in each of the regions
    enddo
    !------------------------------------------------------------
    !...compute total number of nodes and elements
    total_node = 0
    total_elem = 0                                    
	!...initialize total_node and total_elem
    do i=1,total_region
      	total_node = total_node + total_region_node(i)
        total_elem = total_elem + total_region_elem(i)
		!...total number of nodes and elements
    enddo
    !...compute total boundary nodes, crack nodes and crack elements
    total_bdnode = 0
    total_cknode = 0
    total_ckelem = 0
    do i=1,total_region
      	total_bdnode = total_bdnode + total_region_bd_node(i) 
		!...total number of boundary nodes
        total_cknode = total_cknode + total_region_ck_node(i) 
		!...total number of crack nodes
        total_ckelem = total_ckelem + total_region_ck_elem(i) 
		!...total number of crack elements
    enddo
    !------------------------------------------------------------
    read(UNIT=iono,FMT=*)ntitle
    do i=1,ntitle
      	read(unit=iono,FMT='(a80)')title(i)
    enddo													  
	!...number of titles

    !------------------------------------------------------------
    allocate(node_sys2user(total_node),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'node_sys2user: allocation request denied'
        stop
    endif
    allocate(node_coor(2,total_node),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'node_coor: allocation request denied'
        stop
    endif
    allocate(node_id(total_node),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'node_id: allocation request denied'
        stop
    endif
	!------------------------------------------------------------
    allocate(elem_sys2user(total_elem),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'elem_sys2user: allocation request denied'
        stop
    endif
    allocate(elemid(total_elem),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'elemid: allocation request denied'
        stop
    endif
    allocate(elem_crack_id(total_elem),STAT=ierr) !..added 4/3/12..aje567
    if (ierr.ne.0) then
      	print*,'elem_crack_id: allocation request denied'
        stop
    endif
    allocate(elem_region(total_elem),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'elem_region: allocation request denied'
        stop
    endif
    allocate(elnode(3,total_elem),STAT=ierr)
	!...this 3 refers to the max no. of nodes in each element, that still is 3
    if (ierr.ne.0) then
      	print*,'elnode: allocation request denied'
        stop
    endif
    !...7/14/09: change from 3->NDOF
    allocate(elldid(NDOF,total_elem),STAT=ierr)
    !...say if DOF #2 of elem # 3 is BTRAC, then the DOF # 2 of all nodes of elem # 3 is BTRAC... 
	if (ierr.ne.0) then
      	print*,'elldid: allocation request denied'
        stop
    endif
    !------------------------------------------------------------
    allocate(unknown(NDOF,total_node),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'unknown: allocation request denied'
        stop
    endif
    allocate(node_type(NDOF,total_node),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'node_type: allocation request denied'
        stop
    endif
    allocate(equno(3,total_elem),STAT=ierr)
	!...this 3 refers to the max no. of nodes in each element, that is still 3	           
    if (ierr.ne.0) then
      	print*,'equno: allocation request denied'
        stop
    endif
    allocate(elnode_val(NDOF,3,total_elem),STAT=ierr)
	!...again this 3 is still max no. of nodes in each element, that is still 3
    if (ierr.ne.0) then
      	print*,'elnode_val: allocation request denied'
        stop
    endif
    !...2/28/12..added allocate for traction and pressure..aje567
    allocate(elnode_val_trac(NDOF,3,total_elem),STAT=ierr)
	!...this 3 is max no. of nodes in each element
    if (ierr.ne.0) then
      	print*,'elnode_val_trac: allocation request denied'
        stop
    endif
    !...added functionality for multiple cracks..4/24/12...aje567
    allocate(elnode_val_pressure(total_cracks,NDOF,3,total_elem),STAT=ierr)
	!...this 3 is max no. of nodes in each element
    if (ierr.ne.0) then
      	print*,'elnode_val_pressure: allocation request denied'
        stop
    endif

    allocate(node_disp_val(NDOF,total_node),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'node_disp_val: allocation request denied'
        stop
    endif											  
	!...all these declared in global.f90
	!------------------------------------------------------------
    !...nodal data
    !...1/5/10: deactive the rotation of geometry and rotation of applied loading 
    !...since we updated material rotation in kernel1c.f90 and post_remesh.f90
    do i = 1,total_node
!      	if ((class_problem.eq.1).or.&
!        	((class_problem.eq.2).and.(rotate_flag.ne.1))) then
        	!...isotropy, or anisotropy but not rotate structure
      		read(UNIT=iono,FMT=*)node_sys2user(i),(node_coor(j,i),j=1,2),node_id(i)	
			!...node numbers, nodal coordinates and node IDs
    enddo

    !------------------------------------------------------------
    !...element data
    !...Dec9/09: user will provide user# of nodes for element connectivity
    do i = 1,total_elem
      	!read(UNIT=iono,FMT=*)elem_sys2user(i),elemid(i),elem_region(i), &
		!	(elnode(j,i),j=1,3),(elldid(j,i),j=1,NDOF)
        read(UNIT=iono,FMT=*)elem_sys2user(i),elemid(i),elem_crack_id(i),elem_region(i), &
			(itemp(j),j=1,3),(elldid(j,i),j=1,NDOF)
			!...itemp(3) is the middle node in 3 node element
        !...get system# node nodes of element i
        do j = 1,NODE(elemid(i))			
		!...elemid(i) can be anything from 1 to 8, depends on the type of element in question
          	icheck = 0
          	do k = 1,total_node
            	if (itemp(j).eq.node_sys2user(k)) then
                	elnode(j,i) = k	 !elnode(j,i) is the jth node of ith element
                                      
                    icheck = 1
                    exit
                endif						
				!...node_sys2user(k) is the global node number of the kth node
            enddo                           
			!...itemp(j) is the global node number of jth node of ith element
            if (icheck.eq.0) then
              	print*,'error, element ',elem_sys2user(i),'contains undeclared node:',itemp(j)
                stop
            endif
        enddo
    enddo

    !------------------------------------------------------------
    do i=1,total_elem
	    do j=1,NODE(elemid(i))
		    do k=1,NDOF
			    elnode_val_trac(k,j,i)=0.0d0
			enddo
		enddo
	enddo
	!...since we have the elnode_val_trac in formulation and we don't specify tractions on crack nodes, we prescribe zero values, Saumik
	
	!------------------------------------------------------------
    !...rigid body motion constraints
    !...Note: constrain point is the SYSTEM # of node (to be consistent with fadd3d)
    !...Dec9/09: change to USER # (by changing the code below in the session "compute the rigid point eqn")
    read(UNIT=iono,FMT=*)total_rigid_constrain_point  !total number of nodes used for constraints of rigid-body motions
    do i = 1,total_rigid_constrain_point
      	!...read(UNIT=iono,FMT=*)rigid_constrain_point(i),(rigid_constrain_dir(i,j),j=1,NDOF)
        read(UNIT=iono,FMT=*)itemp(1),(rigid_constrain_dir(i,j),j=1,NDOF)  
		!...rigid_constrain_dir(i,j)=1 means that you are fixing the ith
        !...get sys# of itemp(1)										   
		!...rigid constraint point in the jth direction
        j = 0
        do j = 1,total_node
          	if (node_sys2user(j).eq.itemp(1)) exit
        enddo
        if (j.eq.0) then
          	print*,'rigid node ',itemp(1),'has not declared yet'
            stop
        endif
        rigid_constrain_point(i) = j
    enddo

	!------------------------------------------------------------
    !...integration order: read flag to indicate if user need to provide int order
    read(UNIT=iono,FMT=*)integration_flag
    if (integration_flag.eq.YES) then
      	!...order of outer/ln/inner integral when computing stiffness matrix
      	read(UNIT=iono,FMT=*)ng_out,ng_ln,ng_in
        !...integration order when computing load vector
        read(UNIT=iono,FMT=*)ng_ef
        !...order of out/ln/in integral when computing matrix of T-stress
!        read(UNIT=iono,FMT=*)ngt_out,ngt_ln,ngt_in
        !...order when computing L matrix of T-stress
!        read(UNIT=iono,FMT=*)ngt_ef
    else
      	!...use default integration order
        ng_out = 20
                !...number of integration points for outer integral
        ng_ln = 16
                !...16 is the maximum # of log gauss points you can have here
		!...number of integration points for load integral
        ng_in = 30
                !...number of integration points for inner integral
        ng_ef = 20			
		!...number of integration points for the force vector integral
	endif
    !------------------------------------------------------------
    !...check bounds of integration order
    if ((ng_out.gt.max_nint).or.(ng_in.gt.max_nint).or.(ng_ef.gt.max_nint)) then
      	print*,'Gauss integration order is bigger than max order'
        stop
    endif
    if (ng_ln.gt.max_logint) then
      	print*,'logarith integration order is bigger than max order!'
        stop
    endif                   
	!...simply check bounds of integration order
    !------------------------------------------------------------
    !...flag of pure traction problem with inner boundary
    read(UNIT=iono,FMT=*)inner_flag	   
	!...flag for problems with inner boundaries subjected to pure traction
    if (inner_flag.eq.YES) then
        !...read total number of inner boundaries that have pure traction applied
        read(UNIT=iono,FMT=*) total_inner_bound		 
		!...total number of inner boundaries with pure traction applied
        !...allocate variables to store data for inner boundaries
        allocate(inner_region(total_inner_bound),inner_total_elem(total_inner_bound),STAT=ierr)
        if (ierr.ne.0) then
          	print*,'inner_region, inner_total_elem: allocation request denied'
            stop
        endif
        allocate(inner_node(total_inner_bound,2),inner_inside_point(total_inner_bound,2,2),STAT=ierr)
        if (ierr.ne.0) then
          	print*,'inner_node,inner_total_elem: allocation request denied'
            stop
        endif
        !...read region#, total elements
        !...2 separate nodes (USER number)
        !...coordinates of 2 separate inside points for each inner boundary
        do i = 1,total_inner_bound
          	read(UNIT=iono,FMT=*)inner_region(i),inner_total_elem(i)
            !...read(UNIT=iono,FMT=*)inner_node(i,1),inner_node(i,2)
            read(UNIT=iono,FMT=*)itemp(1),itemp(2)   
			!...this is user# of nodes
            !...get sys# of these nodes
            j = 0 
			!...sys# of inner_node(i,1)
            k = 0 
			!...sys# of inner_node(i,2)
            do j = 1,total_node
              	if (node_sys2user(j).eq.itemp(1)) exit  
            enddo
            do k = 1,total_node
              	if (node_sys2user(k).eq.itemp(2)) exit	
            enddo
            if ((j.eq.0).or.(k.eq.0)) then
              	print*,'undeclared nodes for inner boundary ',i
                stop
            endif
            inner_node(i,1) = j		  
			!...inner_node(i,1)=itemp(1)
            inner_node(i,2) = k		  
			!...inner_node(i,2)=itemp(2)
            read(UNIT=iono,FMT=*)(inner_inside_point(i,1,j),j=1,2),(inner_inside_point(i,2,j),j=1,2)
        enddo						  
		!...coordinates of itemp(1) and itemp(2)
        !...read element numbers (USER #) on each inner boundary
        !...first, find max. total number of elements on inner boundaries for allocation
        k = 0   
		!...max number of elements on inner boundary
        do i = 1,total_inner_bound
          	if (k.le.inner_total_elem(i)) then
            	k = inner_total_elem(i)			
				!...k is the total number of inner elements in ith inner boundary
            endif
        enddo
        allocate(inner_elem(total_inner_bound,k),STAT=ierr)
        if (ierr.ne.0) then
          	print*,'inner_elem: allocation request denied'
            stop
        endif
        !...then read element number (USER #) of each inner boundary
        !...NOTE: later in subroutine rigidinner4.f90, we will convert these user#
        !...to sys# for re-computing the real displacements of these elements, AND on inner boundary
        !...no posibility to have two sys# associated with one user# as it is on interface, so the
        !...conversion always have unique image.
        !...(reason not to save sys# here because it will take a lot of memory for temporary variable)
        do i = 1,total_inner_bound
          	read(UNIT=iono,FMT=*) (inner_elem(i,j),j=1,inner_total_elem(i))
        enddo	  
		!...inner_elem(i,j) is the jth inner element of ith inner boundary
    endif
    
    !...fix nodes of each inner boundary to lock its rigid displacements
    if (inner_flag.eq.YES) then
      	!...check bound
        if ((total_rigid_constrain_point + 2*total_inner_bound).gt.max_rigid_constrain_point) then
      		print*,'need to increase max_rigid_constrain_point in global.f90'
        	stop
    	endif
      	do i = 1,total_inner_bound
        	!...fix the first node in all direction
            rigid_constrain_point(total_rigid_constrain_point + 1) = inner_node(i,1)
            do j = 1,NDOF
              	rigid_constrain_dir(total_rigid_constrain_point + 1,j) = 1
            enddo
            !...fix the second node to lock rigid rotation about 3-axis
            if (dabs(node_coor(1,inner_node(i,1))-node_coor(1,inner_node(i,2))).lt.tol) then
              	!...node 1 and 2 are on vertical straight line, fix node 2 in 1-direction
                rigid_constrain_point(total_rigid_constrain_point + 2) = inner_node(i,2)
            	do j = 1,NDOF
              		rigid_constrain_dir(total_rigid_constrain_point + 2,j) = 0
            	enddo
                rigid_constrain_dir(total_rigid_constrain_point + 2,1) = 1
            elseif (dabs(node_coor(2,inner_node(i,1))-node_coor(2,inner_node(i,2))).lt.tol) then
              	!...node 1 and 2 are on horizontal straight line, fix node 2 in 2-direction
                rigid_constrain_point(total_rigid_constrain_point + 2) = inner_node(i,2)
            	do j = 1,NDOF
              		rigid_constrain_dir(total_rigid_constrain_point + 2,j) = 0
            	enddo
                rigid_constrain_dir(total_rigid_constrain_point + 2,2) = 1
            else
              	!...other cases: node 1 and 2 are on inclined line, node 2 can be fixed either in 1 or 2-dir
                rigid_constrain_point(total_rigid_constrain_point + 2) = inner_node(i,2)
            	do j = 1,NDOF
              		rigid_constrain_dir(total_rigid_constrain_point + 2,j) = 0
            	enddo
                rigid_constrain_dir(total_rigid_constrain_point + 2,2) = 1
            endif
      		!...increase the total number of rigid constraints
        	total_rigid_constrain_point = total_rigid_constrain_point + 2
        enddo
    endif
    !------------------------------------------------------------
    read(UNIT=iono,FMT=*),tstress_flag
    !------------------------------------------------------------
    !...4/20/09: output initial mesh to ABAQUS file
    open (15,file='abaqus.inp',status='unknown')
    !...print out node data
    write(15,*)'*Node'
    do i = 1,total_node
      	write(15,'(i10,3(a1,f15.5))'),node_sys2user(i),',',node_coor(1,i),','&
        	,node_coor(2,i)
    enddo
    !...print out 3-node elements
    write(15,*)'*Element,type=B22'
    do i = 1,total_elem
      	if ((elemid(i).eq.BQUAD).or.(elemid(i).eq.CQUAD).or.(elemid(i).eq.CTIP1)&
        	.or.(elemid(i).eq.CTIP2).or.(elemid(i).eq.IQUAD)) then
        	write(15,'(i10,3(a1,i10))'),elem_sys2user(i),','&
            	,node_sys2user(elnode(1,i)),',',node_sys2user(elnode(3,i)),','&
            	,node_sys2user(elnode(2,i))
        endif
    enddo
    !...print out 2-node elements
    write(15,*)'*Element,type=B21'
    do i = 1,total_elem
      	if ((elemid(i).eq.BLIN).or.(elemid(i).eq.CLIN).or.(elemid(i).eq.ILIN)) then
        	write(15,'(i10,2(a1,i10))'),elem_sys2user(i),(','&
            	,node_sys2user(elnode(j,i)),j=1,2)
        endif
    enddo
    !...create set of nodes (exclude mid-side nodes) for displaying the mesh
    write(15,*)'*Nset,nset=edge_node'
    do i = 1,total_elem
      	write(15,*)node_sys2user(elnode(1,i))
        write(15,*)node_sys2user(elnode(2,i))
    enddo
    !--------------------------------------------------------------------------------
    !...2/28/12: Moved pressure flag outside of crack growth case ...aje567
    !!!!!!...NOW THERE ARE SIX FLAG OPTIONS INSTEAD OF FIVE....!!!!!
    !!!!!!...Added remote stresses with tractions calculated the same as internal pressure...!!!!
    !--------------------------------------------------------------------------------

        read(UNIT=iono,FMT=*)pressure_flag
        if ((pressure_flag.ne.0).and.(pressure_flag.ne.1).and.(pressure_flag.ne.2)) then
      		print*,'Error: pressure_flag must be either:' 
                print*,'0 - No Internal Pressure'
                print*,'1 - Volume Constraint'
                print*,'2 - Pressure Constraint'    !pressure_flag=0 means no internal pressure
        	stop
    	endif

        if (pressure_flag.ne.0) then
           read(UNIT=iono,FMT=*)sigma_x  
		   !...read in remote stress in x-direction
           read(UNIT=iono,FMT=*)sigma_y  
		   !...read in remote stress in y-direction
           read(UNIT=iono,FMT=*)pressure_value 
		   !...read in pressure value 
           call PressureTractions
		   !...don't bother about following module since we feed a non-zero pressure flag
        !---------------------------------------
		else
           do i=1,total_elem
              !if (elldid(1,i).ne.BTRFREE.and.elldid(1,i).ne.CTRFREE.and.&
                   !elldid(1,i).ne.NOLOAD) then
                 !...element i is prescribed either BDISP or BTRAC or CTRAC
                 do j = 1,NODE(elemid(i))
                    do k = 1,NDOF
                         elnode_val(k,j,i)=elnode_val_trac(k,j,i)+elnode_val_pressure(elem_crack_id(i),k,j,i)
                    enddo
                 enddo
              !else
                 !...element i is prescribed with traction free or no load
              !   elnode_val(:,:,i) = 0.0d0
              !endif
           enddo
        !---------------------------------------   
        endif

    !--------------------------------------------------------------------------------
    !...CRACK GROWTH FLAG AND GROWTH DATA
    !...5/11/09: parameters for remeshing
    read(UNIT=iono,FMT=*)growth_flag
    if ((growth_flag.ne.0).and.(growth_flag.ne.1)) then
      	print*,'error: growth_flag must be either 0 or 1'
        stop
    endif
    if (growth_flag.eq.YES) then
      	read(UNIT=iono,FMT=*)total_step
        if (total_step.lt.0) then
          	print*,'number of cycle for growth is less than 0!'
            stop
        endif
      	read(UNIT=iono,FMT=*)size_cracktip,grow_factor
        read(UNIT=iono,FMT=*)dcheck1_factor,dcheck2_factor
        read(UNIT=iono,FMT=*)globalK(1),globalK(2),globalK(3)
        !...parameter epsilon is for the case of if (K/Kmax)<(1-epsilon) 
        !...then da/da_max = 1
        read(UNIT=iono,FMT=*)parameter_alpha,parameter_beta,parameter_epsilon
        !-------------------------------------------------------------------------
        !...flag for problem of traction (pressure) applied on crack during growth
        !....MOVED THIS OUTSIDE OF CRACK GROWTH SUBROUTINE....FLAG IN INPUT FILE NOW
        !....BEFORE GROWTH FLAG!......aje567......2/28/12
        !  read(UNIT=iono,FMT=*)pressure_flag
        !  if ((pressure_flag.ne.0).and.(pressure_flag.ne.1)) then
        !		print*,'error: pressure_flag must be either 0 or 1'
        !  	stop
        !	endif
        !  if (pressure_flag.eq.YES) then
        !    	read(UNIT=iono,FMT=*)pressure_value
        !  endif
        !-------------------------------------------------------------------------
        !...5/27/09: read names for output files
        !...output for techplot display of growth (only corner nodes)
        read(UNIT=iono,FMT=*)techplot_file
        !...data at each step for input file to restart
        read(UNIT=iono,FMT=*)input_file
        !...mesh at each step for plotting in Abaqus (corner & midside nodes)
        read(UNIT=iono,FMT=*)abaqus_file
        !..pressure data at each load step
        read(UNIT=iono,FMT=*)pressure_file

        open(112,file=techplot_file,status='unknown')
		!...112 - tecplot file
        open(113,file=input_file,status='unknown')
		!...113 - input file
        open(114,file=abaqus_file,status='unknown')
		!...114 - abaqus file
        open(115,file=pressure_file,status='unknown')
        !...115 - pressure file
    endif
    !...E N D   O F   R E A D I N G   D A T A   F R O M   I N P U T   F I L E
	!--------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------
	!...set Gauss points and weights for 1D: 1 -> max_nint
	call InitGaussIntegration
    !--------------------------------------------------------------------------------
    !...compute the global (average) shear modulus
    !...(temporarily use data of region 1)
!    if (class_problem.eq.1) then
!        !...isotropy
!        gmu = 0.0_REAL_KIND
!        do k = 1,total_region
!            ey = material(1,region_mat_no(k))
!            nu = material(2,region_mat_no(k))
!            mu = 0.5_REAL_KIND*ey/(1.0_REAL_KIND + nu)
!            gmu = gmu + mu
!        enddo
!        gmu = gmu/total_region
!		gmu = 1.0_DBL
!    elseif (class_problem.eq.2) then
!        !...anisotropy: following Jaroon's code, use gmu = 1
!        gmu = 1.0_DBL
!    else
!        print*,'class problem is out of range (1-2): ',class_problem
!        stop
!    endif
    !------------------------------------------------------------
    !...compute the equation number corresponding to each node of element
    iequ = 1		
	!...you basically number the equations corresponding to each node of each element
    do i = 1,total_node
      	if (node_id(i).eq.SBLP) then			
		!...if node is on intersection of crack and boundary
        	!...surface breaking line node
            do k = 1,total_elem
              	do j = 1,NODE(elemid(k))
                	if (elnode(j,k).eq.i) then	
					!...if jth node of kth element is ith global node 
                    	equno(j,k) = iequ       
						!...initial no is corresponding to C+ node
                         endif
                enddo
            enddo
            !...reserve extra dof here for later assigning to C- node
            iequ = iequ + 2*NDOF			   

        elseif (node_id(i).eq.INTER) then      
		!...if node is on interface between 2 regions of multi-region problems
        	!...interface node: traction and displacement are unknown
            do k = 1,total_elem
              	do j = 1,NODE(elemid(k))
                	if (elnode(j,k).eq.i) then
                    	equno(j,k) = iequ
                    endif
                enddo
            enddo
            iequ = iequ + 2*NDOF			   

        elseif (node_id(i).eq.NORMAL) then
        	!...normal node
            do k = 1,total_elem
              	do j = 1,NODE(elemid(k))
                	if (elnode(j,k).eq.i) then
                    	!...5/13/09: control error of input file
                        if (eltype(elemid(k)).eq.INTERFE) then
                          	print*,'error!, normal node: ',node_sys2user(i),&
                            	'is on interface element: ',elem_sys2user(k)
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
    !------------------------------------------------------------
    !...determine equation number for extra nodes on surface breaking line
    do i = 1,total_elem
      	!...just look on regular elements on boundary
      	if (ELTYPE(elemid(i)).ne.BREGULAR) cycle  
		!...cycle statement tells the compiler not to execute any code below the statement and return to the start of the loop and continue the 
		!...loop, using the next value in the index
        !...since all elements are crack elements, this part is not executed due to cycle statement
		!...regular boundary element,regular crack element,crack tip element,element on interface of 2 regions of multi region problem
		equ_i = equno(:,i)						  
        do j = 1,NODE(elemid(i))
          	if (node_id(elnode(j,i)).eq.SBLP) then
            	!...search crack element sharing this node j
                do k = 1,total_elem
                  	if (ELTYPE(elemid(k)).eq.BREGULAR) cycle
                    do m = 1,NODE(elemid(k))
                      	if ((elnode(j,i).eq.elnode(m,k)).and.(j.eq.m)) then
                        	!...boundary element i sharing SBLP node j with crack element k
                           	!...element k is associated with C-
                            equno(j,i) = equ_i(j) + NDOF
                        endif
                    enddo
                enddo
            endif
        enddo
    enddo
    
    !------------------------------------------------------------
    !...compute the rigid point eqn
    !...Dec9/09: rigid_constrain_point(i) is USER # of node being constrainted
    do i = 1,total_rigid_constrain_point
      	do j = 1,total_elem
        	do k = 1,NODE(elemid(j))
                !...if (node_sys2user(elnode(k,j)).eq.rigid_constrain_point(i)) then
            	if (elnode(k,j).eq.rigid_constrain_point(i)) then
                	rigid_constrain_point_eqn(i) = equno(k,j)
                endif
            enddo
        enddo
    enddo
    !------------------------------------------------------------
    !...total dof
    !...5/3/09: fix the bug when last node is interface/SBLP node
    if ((node_id(total_node).eq.INTER).or.(node_id(total_node).eq.SBLP)) then
      	total_dof = MAXVAL(equno) + 2*NDOF -1
    elseif (node_id(total_node).eq.NORMAL) then
    	total_dof = MAXVAL(equno) + NDOF - 1
		!...as we checked,the maximum in equno is 25 (for 1 crack 4 element 1 region 3 DOF hydraulic fracture problem)
		!...add 3 to it, becomes 28, subtract 1 from that, it becomes 27 (total # of DOF)
    else
      	print*,'error from prep.f90: id of last node is out of range'
        stop
    endif
    !------------------------------------------------------------
    !...determine the unknowns of each node
    
	allocate(tcount(NDOF),dcount(NDOF),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'tcount,dcount: allocation request denied'
		print*,'tcount: allocation request denied'
        stop
    endif
    
	do i = 1,total_node
      	if (node_id(i).eq.INTER) then
		!...if node is on interface 
        	unknown(:,i) = DISP_TRAC
            cycle
        endif
		!...now we are considering nodes with either traction prescribed or displacement prescribed
        tcount = 0 
		!...counter of prescribed traction
        dcount = 0 
		!...counter of prescribed displacement
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
                            	print*,'id of prescribed data (ellid) on element:'&
                                	,elem_sys2user (j),'is out of range'
                                stop
                        end select
                    enddo
                endif
            enddo
		enddo
        do alpha = 1,NDOF
          	if ((tcount(alpha).eq.0).and.(dcount(alpha).gt.0)) then
            	!...all elements sharing node i have prescribed displ in alpha-dir
                unknown(alpha,i) = TRACTION
                node_type(alpha,i) = 0      !...not an edge node
            elseif ((dcount(alpha).eq.0).and.(tcount(alpha).gt.0)) then
            	!...all elements sharing node i have prescribed traction (or traction
                !...free) in alpha-direction
                unknown(alpha,i) = DISPLACEMENT
                node_type(alpha,i) = 0      !...not an edge node
            elseif ((dcount(alpha).gt.0).and.(tcount(alpha).gt.0)) then
                !...elements sharing node i some are traction prescribed, 
                !...some are displacement prescribed
                unknown(alpha,i) = TRACTION
                node_type(alpha,i) = EDGE_NODE
                !...Note: even with an real edge node but if all elements sharing that
                !...node have either traction or displacement prescribed, then it is
                !...still considered "not an edge node"
            else
              	print*,'sth wrong with prescribed data of elements sharing node:'&
                ,node_sys2user(i),'in direction:',alpha
                stop
            endif
        enddo
    enddo   !of i=1,total_node
	!------------------------------------------------------------
	!...allocate space for AK and F
    allocate(dF(total_dof),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'dF: allocation request denied!'
        stop
    endif

    allocate(dF_traction(total_dof),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'dF_traction: allocation request denied!'
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


	!------------------------------------------------------------
    !...initialize dF and dAK
    dF = 0.0d0
    dF_traction = 0.0d0
    dF_pressure = 0.0d0

    dAK = 0.0d0
    dAK_traction = 0.0d0
    dAK_pressure = 0.0d0

	!------------------------------------------------------------


    !...display basis information of problem
    print*,'TITLE'
    do i = 1,ntitle
    	print'(a80)',title(i)
    enddo
    if (NDOF.eq.1) then
      	print'(a30)','Darcy flow in porous media    '
    elseif (NDOF.eq.2) then
      	print'(a30)','Plane strain of elastic media '
    elseif (NDOF.eq.3) then
      	print'(a30)','Elastic Media                 '
    elseif (NDOF.eq.4) then
    	print'(a30)','Piezoelectric Media           '
    elseif (NDOF.eq.5) then
    	print'(a30)','MagnetoPiezo Elastic Media    '
    else
      	print*,'current version has not supported this media type;',NDOF
        stop
    endif
    print'(a30,i15)','Total number of cracks: ',total_cracks
    print'(a30,i15)','Total number of nodes: ',total_node
    print'(a30,i15)','Total number of elements: ',total_elem
    print'(a30,i15)','Total number of DOF: ',total_dof
    print'(a30,i15)','Length of stiffness vector: ',total_dof*(total_dof+1)/2
    print*,'INTEGRATION PARAMETERS'
    print'(a50,i15)','Number of int points of outer integral: ',ng_out
    print'(a50,i15)','Number of int points of inner integral: ',ng_in
    print'(a50,i15)','Number of int points of logarith integral: ',ng_ln
    print'(a50,i15)','Number of int points of force vector integral: ',ng_ef
	print*,'DEFINITION OF MATERIALS'
   	do i = 1,total_material
       	print*,'Material name: ',material_name(i)
        print*,'Material No. : ',i
        print*,'Following is material properties for elastic media case!'
        select case (material_type(i))
           	case (1)
               	print*,'Isotropic material'
                print*,"Young's modulus = ",material(1,i)
          		print*,"Poisson ratio = ",material(2,i)
           	case (2)
               	print*,'Cubic material'
                print'(a10,e15.3)','C11 = ',material(1,i)
                print'(a10,e15.3)','C12 = ',material(2,i)
                print'(a10,e15.3)','C44 = ',material(3,i)
            case (3)
               	print*,'Transversely isotropyic material: 2 is elastic symmetry'
                print'(a10,e15.3)','C11 = ',material(1,i)
                print'(a10,e15.3)','C12 = ',material(2,i)
                print'(a10,e15.3)','C13 = ',material(3,i)
                print'(a10,e15.3)','C22 = ',material(4,i)
                print'(a10,e15.3)','C44 = ',material(5,i)
                if (NDOF.eq.4) then
                  	print'(a10,e15.3)','e16 = ',material(6,i)
                    print'(a10,e15.3)','e21 = ',material(7,i)
                    print'(a10,e15.3)','e22 = ',material(8,i)
                    print'(a10,e15.3)','k11 = ',material(9,i)
                    print'(a10,e15.3)','k22 = ',material(10,i)
                endif
            case (4)
            	print*,'Orthotropic material'
            	if (NDOF.eq.2) then
                	print'(a10,e15.3)','C11 = ',material(1,i)
                    print'(a10,e15.3)','C12 = ',material(2,i)
                    print'(a10,e15.3)','C22 = ',material(3,i)
                    print'(a10,e15.3)','C44 = ',material(4,i)
                    print'(a10,e15.3)','C55 = ',material(5,i)
                    print'(a10,e15.3)','C66 = ',material(6,i)
               	elseif (NDOF.eq.3) then
                	print'(a10,e15.3)','C11 = ',material(1,i)
                	print'(a10,e15.3)','C12 = ',material(2,i)
                	print'(a10,e15.3)','C13 = ',material(3,i)
                	print'(a10,e15.3)','C22 = ',material(4,i)
                	print'(a10,e15.3)','C23 = ',material(5,i)
                	print'(a10,e15.3)','C33 = ',material(6,i)
                	print'(a10,e15.3)','C44 = ',material(7,i)
                	print'(a10,e15.3)','C55 = ',material(8,i)
                	print'(a10,e15.3)','C66 = ',material(9,i)
                endif
            case (5)
              	print*,'Monoclinic material'
                if (NDOF.eq.2) then
                  	print'(a10,e15.3)','C11 = ',material(1,i)
                	print'(a10,e15.3)','C12 = ',material(2,i)
                	print'(a10,e15.3)','C16 = ',material(3,i)
                	print'(a10,e15.3)','C22 = ',material(4,i)
                	print'(a10,e15.3)','C26 = ',material(5,i)
                	print'(a10,e15.3)','C44 = ',material(6,i)
                	print'(a10,e15.3)','C45 = ',material(7,i)
                	print'(a10,e15.3)','C55 = ',material(8,i)
                	print'(a10,e15.3)','C66 = ',material(9,i)
                elseif (NDOF.eq.3) then
                	print'(a10,e15.3)','C11 = ',material(1,i)
                	print'(a10,e15.3)','C12 = ',material(2,i)
                	print'(a10,e15.3)','C13 = ',material(3,i)
                	print'(a10,e15.3)','C16 = ',material(4,i)
                	print'(a10,e15.3)','C22 = ',material(5,i)
                	print'(a10,e15.3)','C23 = ',material(6,i)
                	print'(a10,e15.3)','C26 = ',material(7,i)
                	print'(a10,e15.3)','C33 = ',material(8,i)
                	print'(a10,e15.3)','C36 = ',material(9,i)
                	print'(a10,e15.3)','C44 = ',material(10,i)
                	print'(a10,e15.3)','C45 = ',material(11,i)
                	print'(a10,e15.3)','C55 = ',material(12,i)
                	print'(a10,e15.3)','C66 = ',material(13,i)
                endif
            case (6)
              	print*,'General anisotropic material'
                if (NDOF.eq.2) then
                  	print'(a10,e15.3)','C11 = ',material(1,i)
                	print'(a10,e15.3)','C12 = ',material(2,i)
                	print'(a10,e15.3)','C14 = ',material(3,i)
                	print'(a10,e15.3)','C15 = ',material(4,i)
                	print'(a10,e15.3)','C16 = ',material(5,i)
                	print'(a10,e15.3)','C22 = ',material(6,i)
                	print'(a10,e15.3)','C24 = ',material(7,i)
                	print'(a10,e15.3)','C25 = ',material(8,i)
                	print'(a10,e15.3)','C26 = ',material(9,i)
                	print'(a10,e15.3)','C44 = ',material(10,i)
                	print'(a10,e15.3)','C45 = ',material(11,i)
                	print'(a10,e15.3)','C46 = ',material(12,i)
                	print'(a10,e15.3)','C55 = ',material(13,i)
                	print'(a10,e15.3)','C56 = ',material(14,i)
                	print'(a10,e15.3)','C66 = ',material(15,i)
                elseif ((NDOF.eq.3).or.(NDOF.eq.4)) then
                	print'(a10,e15.3)','C11 = ',material(1,i)
                	print'(a10,e15.3)','C12 = ',material(2,i)
                	print'(a10,e15.3)','C13 = ',material(3,i)
                	print'(a10,e15.3)','C14 = ',material(4,i)
                	print'(a10,e15.3)','C15 = ',material(5,i)
                	print'(a10,e15.3)','C16 = ',material(6,i)
                	print'(a10,e15.3)','C22 = ',material(7,i)
                	print'(a10,e15.3)','C23 = ',material(8,i)
                	print'(a10,e15.3)','C24 = ',material(9,i)
                	print'(a10,e15.3)','C25 = ',material(10,i)
                	print'(a10,e15.3)','C26 = ',material(11,i)
                	print'(a10,e15.3)','C33 = ',material(12,i)
                	print'(a10,e15.3)','C34 = ',material(13,i)
                	print'(a10,e15.3)','C35 = ',material(14,i)
                	print'(a10,e15.3)','C36 = ',material(15,i)
                	print'(a10,e15.3)','C44 = ',material(16,i)
                	print'(a10,e15.3)','C45 = ',material(17,i)
                	print'(a10,e15.3)','C46 = ',material(18,i)
                	print'(a10,e15.3)','C55 = ',material(19,i)
                	print'(a10,e15.3)','C56 = ',material(20,i)
                	print'(a10,e15.3)','C66 = ',material(21,i)
                	if (NDOF.eq.4) then
                    	print'(a10,e15.3)','e11 = ',material(22,i)
                        print'(a10,e15.3)','e12 = ',material(23,i)
                        print'(a10,e15.3)','e13 = ',material(24,i)
                        print'(a10,e15.3)','e14 = ',material(25,i)
                        print'(a10,e15.3)','e15 = ',material(26,i)
                        print'(a10,e15.3)','e16 = ',material(27,i)
                        print'(a10,e15.3)','e21 = ',material(28,i)
                        print'(a10,e15.3)','e22 = ',material(29,i)
                        print'(a10,e15.3)','e23 = ',material(30,i)
                        print'(a10,e15.3)','e24 = ',material(31,i)
                        print'(a10,e15.3)','e25 = ',material(32,i)
                        print'(a10,e15.3)','e26 = ',material(33,i)
                        print'(a10,e15.3)','e31 = ',material(34,i)
                        print'(a10,e15.3)','e32 = ',material(35,i)
                        print'(a10,e15.3)','e33 = ',material(36,i)
                        print'(a10,e15.3)','e34 = ',material(37,i)
                        print'(a10,e15.3)','e35 = ',material(38,i)
                        print'(a10,e15.3)','e36 = ',material(39,i)
                        print'(a10,e15.3)','k11 = ',material(40,i)
                        print'(a10,e15.3)','k12 = ',material(41,i)
                        print'(a10,e15.3)','k13 = ',material(42,i)
                        print'(a10,e15.3)','k22 = ',material(43,i)
                        print'(a10,e15.3)','k23 = ',material(44,i)
                        print'(a10,e15.3)','k33 = ',material(45,i)
                    endif
                endif
            case default
              	print*,'this material type is out of range: ',material_type(i)
                stop
        end select
    enddo
    print*,'DEFINITION OF REGIONS'
    do i = 1,total_region
      	print*,'Region: ',region_name(i)
        print*,'Material: ',material_name(region_mat_no(i))
        print*,'Total elements in the region: ',total_region_elem(i)
    enddo
    if (model_flag.eq.YES) then
      	!...print out input data'
        call InputDataEcho
    endif

    !------------------------------------------------------------


    CONTAINS
    Subroutine InputDataEcho
		implicit none
        integer	:: i,j,k,m,ie
        print*,'DEFINITION OF NODAL POINTS'
        print*,'Node      x      y      node_id'
        do i = 1,total_node
          	print'(i5,2x,f10.5,2x,f10.5,2x,i5)',node_sys2user(i),&
            	(node_coor(k,i),k=1,2),node_id(i)
        enddo
        if (total_rigid_constrain_point.ne.0) then
          	print*,'CONSTRAINS TO REMOVE RIGID BODY MOTION'
            print*,'Node      x-dir      y-dir     z-dir...'
            do i = 1,total_rigid_constrain_point
              	print*,rigid_constrain_point(i),&
                	(rigid_constrain_dir(i,j),j=1,NDOF)
            enddo
        endif
        !...output element definition
        print*,'DEFINITION OF ELEMENTS'
        ie = 0
        do k = 1,total_region
          	print'(a8,1x,a32,1x,a1,i2,a1)','Region:',region_name(k),'(',k,')'
            do i = 1,total_region_elem(k)
              	ie = ie + 1
                print'(i8,3(i8,1x))',elem_sys2user(ie),&
                	(node_sys2user(elnode(j,ie)),j=1,3)
            enddo
        enddo
        !...output prescribed value
        print*,'PRESCRIBED BOUNDARY VALUES'
        ie = 0
        do k = 1,total_region
          	do i = 1,total_region_elem(k)
            	ie = ie + 1
                !if (elldid(1,ie).eq.BTRFREE.or.elldid(1,ie).eq.CTRFREE.or.&
                !  	elldid(1,ie).eq.NOLOAD) cycle
					!if element is either boundary traction free or crack is traction free or no load
                do m = 1,NODE(elemid(ie))
                  	print'(5x,i5,5x,a1,i1,3(f10.3,2x,i2))',elem_sys2user(ie),&
                    	'n',m,(elnode_val_trac(kk,m,ie),elldid(kk,ie),kk=1,NDOF)
                enddo
            enddo
        enddo
    End Subroutine InputDataEcho


!-------------------------------------------------------------------------------------------------------------------
    !...2/29/12: subroutine to calculate the tractions on the crack due to the user input pressure


    Subroutine PressureTractions
    IMPLICIT NONE
	!...local variables
    INTEGER		:: ie,iee,i,j,k,m, ck
    REAL(KIND=DBL)	:: normal_iee(2,3),normal_adj(2,3),normal_avg(2)
    REAL(KIND=DBL)	:: dj,xele(2,3)
    INTEGER		:: ielem_adj
    INTEGER		:: zero_elem !...this is used to save to (global) position of the first elem of current region

    !---------------------------------------------------------------------------------------------
    !...convert user input pressure value into traction due to pressure...added 2/29/12..aje567
    !---------------------------------------------------------------------------------------------

    !...the the global position for the first element of region k
    ie = 0
    zero_elem = 0
    do k = 1,total_region
       if (total_region_ck_elem(k).eq.0) then                             
	   !...I have changed the index from 1 to k        
          !...no crack in region k, then just need to get the global element ie
          ie = ie + total_region_elem(k)
          cycle
       elseif (total_region_ck_elem(k).gt.0) then
          !...5/27/09: get the global position for the first elem of region k
          zero_elem = ie              
		  !...zero_elem is the 'zeroth' element of the first region with a crack in it
       endif			                                     

    !...k takes region # of regions with cracks in them 
	
	!...if for instance the first 2 regions have no cracks in them, and third regions has cracks, and given that region 1 has 12 elements and 
	!...region 2 has 14 elements, then when k hits 3, zero_elem=12+14=26, and k=3
	
	!...calculate traction on crack elements
    iee = zero_elem						 
	!...iee is the zeroth element of first region (k) with a crack in it
    if (pressure_flag.ne.0) then         
	!...there is either a volume constraint (pressure_flag=1) or pressure constraint (pressure_flag=2)
       do i = 1,total_region_elem(k)	 
	   !...again here I have changed the index from 1 to k
          !...global element #
          iee = iee + 1					 
		  !...with the above scenario in mind, iee=27 and onwards
          !...skip process if iee is not a crack element
          if ((eltype(elemid(iee)).ne.CTIP).and.(eltype(elemid(iee)).ne.CREGULAR)) cycle
          !...get nodal coordinates of element iee
		  !...use the cycle statement to loop over the non-crack elements
		  !...so if first 4 elements are non-crack elements, then iee loops over 27,28,29,30 and takes value 31
		  !...then we check whether element # 31 is regular crack element or crack tip element and so on...
          xele = node_coor(:,elnode(:,iee))	 
		  !...nodal coordinates of either the crack tip element	or regular crack element, whichever comes first
          !...compute normals at nodes of element iee
          !...xele has been declared as xele(2,3)...So it contains the coordinates of nodes in element under consideration
		  call NodeNormal(elemid(iee),xele,normal_iee)	
		  !...elemid could be 1 through 8 anything
          !...iee is first crack element
		  !...search for adjacent element on the side of node1 of iee
          ielem_adj = 0
          do j = 1,total_elem
             if ((eltype(elemid(j)).eq.CTIP).or.(eltype(elemid(j)).eq.CREGULAR)) then
                if (elnode(1,iee).eq.elnode(2,j)) then
                   ielem_adj = j	
				   !...if first node of crack tip/regular crack element under consideration is same as 2nd node of jth element
                   exit
                endif
             endif
          enddo
          !...averaging normal at node1 of iee
          if (ielem_adj.ne.0) then
             xele = node_coor(:,elnode(:,ielem_adj))
             call NodeNormal(elemid(ielem_adj),xele,normal_adj)

             !...average normal at node1 of iee
             do j = 1,2
                normal_avg(j) = normal_iee(j,1) + normal_adj(j,3)  
				!...node normal average
             enddo
             dj = dsqrt(normal_avg(1)*normal_avg(1) + normal_avg(2)*normal_avg(2))
             normal_iee(1,1) = normal_avg(1)/dj
             normal_iee(2,1) = normal_avg(2)/dj	   
			 !...final node normal on node 1 of iee
          endif

          !--------------------------------------------------
          !...search for adjacent element on the side of node2 of iee
          ielem_adj = 0
          do j = 1,total_elem
             if ((eltype(elemid(j)).eq.CTIP).or.(eltype(elemid(j)).eq.CREGULAR)) then
			 !...if jth element is crack tip element or regular crack element
                if (elnode(2,iee).eq.elnode(1,j)) then
                   ielem_adj = j  
				   !...if second node of crack tip/regular crack element under consideration is same as 1st node of jth element
                   exit
                endif
             endif
          enddo

          !...averaging normal at node2 of iee
          if (ielem_adj.ne.0) then
             xele = node_coor(:,elnode(:,ielem_adj))
             call NodeNormal(elemid(ielem_adj),xele,normal_adj)

             !...average normal at node2 of iee
             do j = 1,2
                normal_avg(j) = normal_iee(j,3) + normal_adj(j,1)
             enddo
             dj = dsqrt(normal_avg(1)*normal_avg(1) + normal_avg(2)*normal_avg(2))
             normal_iee(1,2) = normal_avg(1)/dj
             normal_iee(2,2) = normal_avg(2)/dj	   
			 !...final node normal on node 2 of iee
          endif

		  !--------------------------------------------------------------------------------------------
          !...calculate traction applied on crack elements due to pressure and remote stress
          do j = 1,NODE(elemid(iee))

             elnode_val_trac(1,j,iee) = -1.0d0*sigma_x*normal_iee(1,j)
                          
             elnode_val_trac(2,j,iee) = -1.0d0*sigma_y*normal_iee(2,j)

             elnode_val_trac(3,j,iee) = 0.0d0
			 !..Changed it back to 3DOF	problem  

             do ck = 1,total_cracks !..loop over all cracks
                do m = 1,2 !..loop over first two  NDOF

                   if(ck .eq. elem_crack_id(iee))then   
				   !..if current crack is the crack that elem iee is on
                      elnode_val_pressure(ck,m,j,iee) = -1.0d0*pressure_value*normal_iee(m,j)   
					  !...elnode_val_pressure is contribution to traction on crack due to pressure on crack
                   else
                      elnode_val_pressure(ck,m,j,iee) = 0.0d0
					  !...This means that pressure in crack #ck has no contribution to traction in element iee which does not belong to crack
					  !...#ck
                   endif 

                enddo
                
                elnode_val_pressure(ck,3,j,iee) = 0.0d0
				!...changed it back to 3 dof proeblem
             enddo


          enddo ! of j=1,NODE(elemid(iee))
       enddo !of i=1,total_region_elem(k)
    endif !of (pressure_flag.eq.YES)

    !---------------------------------------------------------------------------------------------
    !...set traction value as sum of tractions and pressures
    !---------------------------------------------------------------------------------------------
 enddo
  End Subroutine PressureTractions

    !---------------------------------------------------------------------------------------------


  Subroutine NodeNormal(id,xel,normal)
    !...subroutine to calculate normals at nodes of element
	!...normal directed into the discontinuity
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
			!...the 3 coordinates of master element i.e.xi=-1,0,1
        end select
        call reg_shape(id,xi,ps,dps)    
		!...id is element ID, could be 1 through 8 depending on what type of element
        dxds = 0.0d0
        do j = 1,NODE(id)
           	do k = 1,2
               	dxds(k) = dxds(k) + xel(k,j)*dps(j)	  
				!...since isoparametric elements are used, the same shape functions are used to interpolate the nodal coordinates as well
            enddo
        enddo
        dj = dsqrt(dxds(1)*dxds(1) + dxds(2)*dxds(2))   
		!...dj is magnitude of tangent vector
        !...tangent vector
        do j = 1,2
           	tangent(j) = dxds(j)/dj
        enddo
        !...normal = tangent x e3
		!...n = t x e_3 = i(t_2) - j(t_1) = i(t_2) + j(-t_1)
		!...n_x = t_2, n_y = -t_1
        normal(1,i) = tangent(2)
        normal(2,i) = -1.0d0*tangent(1)   
		!...components of normal drawn to the element at xi=-1,0,1
    enddo
  End Subroutine NodeNormal


    !...include the source file of set up integration points and weights for 1D problem
    !...with this command, we need a specific command in makefile (same for elem.f90)
    INCLUDE "setint.f90"	             
End Subroutine DataPreparation
