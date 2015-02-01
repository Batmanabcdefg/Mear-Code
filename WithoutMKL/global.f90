MODULE DefinitionConstant   
!...this module sets numbers to flags (nodeID flag,type of element flag,prescribed data on element flag,unknowns 
!...on nodes flag,yes/no flag)   
!...Note: this version use double precision for real variables/constants
!...and all intrinsic function are associated with double precision too
	IMPLICIT NONE
    !...Node ID
    INTEGER,PARAMETER :: NORMAL      = 0, &    !...normal node i.e. node on crack or node on boundary
    					 SBLP        = -1,&	   !...node on intersection of crack and boundary
                         INTER       = -2      !...node on interface between 2 regions of multi-region problems
    !...Element ID
    !...Current version ONLY work for quadratic elements with 3 nodes
	INTEGER,PARAMETER :: CTIP1 = 1, & !...3-node tip element (tip-node1)
                         CTIP2 = 2, & !...3-node tip element (tip-node2)  
	!...CTIP1 is for one crack tip and CTIP2 is for other crack tip of the same crack
                         CQUAD = 3, & !...3-node quadratic crack element
                         CLIN  = 4, & !...2-node linear crack element
                         BQUAD = 5, & !...3-node quadratic boundary element
                         BLIN  = 6, & !...2-node linear boundary element
                         !...4/28/09: add 2 types of interface element
                         IQUAD = 7, & !...3-node quadratic interface element  
						 !interface elements are elements on interface for multi-region problems
                         ILIN  = 8    !...2-node linear interface element	  		   
						 															   
	!...Type of prescribed data on elements
    INTEGER,PARAMETER :: BDISP       = 0, &   !...element is on boundary with prescribed displacement
    					 BTRAC       = 1, &   !...element is on boundary with prescribed traction
                         CTRAC       = 2, &   !...element is on crack with prescribed traction
                         BTRFREE     = 3, &   !...element is on traction free boundary
                         CTRFREE     = 4, &   !...element is on traction free crack
                         NOLOAD      = 5	  !...element is on interface for multi-region problems
                          
	!...Unknown of nodes
    INTEGER,PARAMETER :: TRACTION     = 0, &  !...traction on nodes unknown
    					 DISPLACEMENT = 1, &  !...displacement on nodes unknown
                         DISP_TRAC    = 2	  !...traction and displacement on nodes unknown
	!...Edge node flag
    INTEGER,PARAMETER :: EDGE_NODE    = 1
    !...Type of elements
    INTEGER,PARAMETER :: BREGULAR     = 1, &  !...regular boundary element 
    					 CREGULAR     = 2, &  !...regular crack element
                         CTIP         = 3, &  !...crack tip element
                         INTERFE      = 4	  !...element on the interface of multi-region problems
	!...Flags
    INTEGER,PARAMETER :: YES          = 1, &
                         NO           = 0, &
		    	         LOW          = 1, &
                         HIGH         = 2
    !...double precision (15 digits after decimal)
	INTEGER,PARAMETER :: DBL=SELECTED_REAL_KIND(15)
    !...define biggest number and smallest number
	REAL(KIND=DBL) :: LARGE_NUM = 1.0d+16, SMALL_NUM = 1.0d-16
END MODULE DefinitionConstant
!----------------------------------------------------------------------
MODULE GlobalData               
!...this module is exclusively used to introduce/declare the variables in the code                    
	USE DefinitionConstant
	!...title of the program
	CHARACTER(80), SAVE :: title(30)
    !...8/6/09: NDOF is now also the indicator for type of media
    !...NDOF = 1 -> Darcy's flow in porous media
    !...NDOF = 2 -> Plane strain of elastic media
    !...NDOF = 3 -> Elastic media (plane strain + anti-plane shear)
    !...NDOF = 4 -> Piezoelectric media
    !...NDOF = 5 -> Magnetoelectroelastic media
    INTEGER, SAVE	:: NDOF									 
	!...number of degrees of freedom
    
    REAL(KIND=DBL), SAVE	:: ht                             
	!...height of crack
	!...added this

	!...indicator of problem with inner boundary under pure traction
    INTEGER, SAVE	:: inner_flag              
	!...if there is an inner boundary subjected to pure traction like plate with doubly cracked hole
    !...total number of inner boundary
    INTEGER, SAVE	:: total_inner_bound       
	!...total number of inner boundaries
    !...region for each inner boundary
    INTEGER, ALLOCATABLE, SAVE	:: inner_region(:)   
	!...what goes inside is total number of inner boundaries
    !...total elements for each inner boundary
    INTEGER, ALLOCATABLE, SAVE	:: inner_total_elem(:)	
	!...what goes inside is again total number of inner boundaries
    !...element (system) number of each inner boundary
    INTEGER, ALLOCATABLE, SAVE	:: inner_elem(:,:)		
	!...first index is inner boundary number and second index is element # in each of these inner boundaries
    !...two separate nodes on of each inner boundary
    !...(these nodes will be use to lock rigid displacemens)
    INTEGER, ALLOCATABLE, SAVE	:: inner_node(:,:)		
	!...first index is inner boundary number, second index is 1 or 2, for the endpoints of the inner boundary
    !...two separate points inside each inner boundary
    !...(these points will be used to calculate rigid translation and rotation)
    REAL(KIND=DBL), ALLOCATABLE, SAVE	:: inner_inside_point(:,:,:)  
	!...first and second index same as inner_node, third index is 1 or 2, feeding coordinates of the inner nodes
    
    !...T-stress calculation flag
    INTEGER, SAVE	:: tstress_flag									
	!...whether or not you want to calculate T stress
    !...flags
    INTEGER, SAVE	:: model_flag = 1, & 
    				   boundary_flag = 1, sif_flag = 1	
					   !...model flag is whether or not you want to echo input data in output file, boundary flag is whether or not
					   !...you want to print out the results of displacements and stresses on boundary, sif flag is whether or not
    				   !...you want to print out the results of stress intensity factors at crack tips
    !...global matrix and load vector
    REAL(KIND=DBL),ALLOCATABLE,SAVE	:: dAK(:),dAK_traction(:),dAK_pressure(:,:),dF(:),dF_traction(:),dF_pressure(:,:)
    !..added stiffness matrix and load vectors for pressure and tractions...3/5/12...aje567
    !..added stiffness matrix and load vectors for pressure for two cracks...4/3/12...aje567
    !..changed for pressure to be allocatable for number of cracks...4/24/12...aje567


    !...4/6/09: global data to serve tstress subroutine:
    REAL(KIND=DBL),ALLOCATABLE,SAVE	:: sol_trac(:,:,:), sol_disp(:,:,:), sol_trac_traction(:,:,:), sol_disp_traction(:,:,:),&
                                           sol_trac_pressure(:,:,:,:), sol_disp_pressure(:,:,:,:),opening_disp(:,:),&
                                           opening_disp_traction(:,:), opening_disp_pressure(:,:,:)
	!...5/10/09: global data:
    REAL(KIND=DBL),ALLOCATABLE,SAVE	:: nodal_SIF(:,:),nodal_SIF_total(:,:), nodal_SIF_traction(:,:),nodal_SIF_pressure(:,:,:),SIF_effective(:)

    REAL(KIND=DBL), ALLOCATABLE     :: uk_traction(:),uk_pressure(:,:),SIF_traction(:), SIF_pressure(:,:)


    REAL(KIND=DBL), ALLOCATABLE,SAVE:: crack_scaling(:), crack_scaling_temp(:,:)

    !...system to user node # of tip nodes
	INTEGER,ALLOCATABLE,SAVE		:: tipnode_sys2user(:)
    
    !...kernels
    REAL(KIND=DBL),ALLOCATABLE,SAVE :: BI1(:,:),EI1(:,:)
	!...added this
    REAL(KIND=DBL),ALLOCATABLE,SAVE :: BI2(:,:),EI2(:,:)
	!...added this
    !...all these come in the kernel, space for all this is allocated in prep.f90

    INTEGER, ALLOCATABLE, SAVE   	:: node_sys2user(:),elem_sys2user(:)  
	!...memory allocated in prep.f90
    
    !...number of nodes/elements
    INTEGER, SAVE		:: total_bdnode, total_cknode, total_bdelem, &	
    					   total_ckelem, total_node, total_elem			
						   !...bdnode and bdelem are boundary nodes and boundary elements, cknode is crack node
						   !...total number of nodes and elements, ckelem is crack element
    INTEGER, SAVE		:: total_material,total_region, total_cracks 
	!..added total number of cracks...4/3/12

    INTEGER				:: total_tip_node, tip_node_count				
	!...total number of tip nodes, tip node counter

    !...class problem
    !...1 - isotropy
    !...2 - anisotropy
    !...INTEGER, SAVE       :: class_problem !11/24/2010: delete this variable since it is irrelevant
    !...5/1/10: rotate_flag and rotate_angle becomes global so that we can
    !...transform material constants to geometry coordinates
    REAL(KIND=DBL)		:: rotate_angle	
	!...if coordinates of material and structure are different, rotate_angle is the angle between x1 axis of
    !...structure and x1 axis of material
	INTEGER				:: rotate_flag	
	
    INTEGER,ALLOCATABLE,SAVE		:: material_type(:)					  
	!...what goes inside is material # 
    CHARACTER(32),ALLOCATABLE,SAVE 	:: material_name(:), region_name(:)	  
	!...what goes inside material_name is material #, what goes inside region_name is region #
    !...data for each region:
    INTEGER,ALLOCATABLE,SAVE	:: region_mat_no(:), total_region_elem(:),&
    							   total_region_node(:), &
                                   total_region_bd_elem(:), &			  
                                   total_region_ck_elem(:), &
                                   total_region_bd_node(:), &
                                   total_region_ck_node(:)
								   !...bd is boundary and ck is crack				  
								   !...all these are allocated memory in prep.f90
	!...total number of dof of problem:
    INTEGER, SAVE		:: total_dof									  
	!...total number of dofs
    REAL(KIND=DBL), ALLOCATABLE, SAVE	:: node_coor(:,:)
    INTEGER, ALLOCATABLE, SAVE   		:: node_id(:),unknown(:,:), &	  
    									   node_type(:,:)
    !...node ID and node type, allocated in prep.f90
	!...element connectivity:
    INTEGER, ALLOCATABLE, SAVE   		:: elnode(:,:), elemid(:),elem_crack_id(:), &	   
    									   elldid(:,:), elem_region(:)	  
										   !...elldid is element load ID
										   !...element connectivity, allocated in prep.f90
    !...nodal value of element
    REAL(KIND=DBL), ALLOCATABLE, SAVE   :: elnode_val(:,:,:)			  
	!...allocated in prep.f90
    !...nodal value of element due to tractions...added 2/28/12...aje567
    REAL(KIND=DBL), ALLOCATABLE, SAVE   :: elnode_val_trac(:,:,:)		  
	!...allocated in prep.f90
    !...nodal value of element due to pressure....added 2/28/12...aje567
    REAL(KIND=DBL), ALLOCATABLE, SAVE   :: elnode_val_pressure(:,:,:,:) 
	!...added number of cracks, allocated in prep.f90


    !...prescribed nodal displacement value								  
    REAL(KIND=DBL), ALLOCATABLE, SAVE   :: node_disp_val(:,:)			  
	!...allocated in prep.f90
    !...equation number in global system
    INTEGER, ALLOCATABLE, SAVE   		:: equno(:,:)					  
	!...equation number in global system, memory allocated in prep.f90

    !...rigid constraint: max number of rigid constraints is 6
    !...can be expanded later
    !...also, change 3->NDOF for rigid_constrain_dir
    INTEGER, parameter	:: max_rigid_constrain_point = 10                 
	!...maximum number of rigid constraints
    INTEGER, SAVE		:: total_rigid_constrain_point = 0
    INTEGER, SAVE		:: rigid_constrain_point(max_rigid_constrain_point)
	!...node used for constraints of rigid body motions
    INTEGER, ALLOCATABLE, SAVE		:: rigid_constrain_dir(:,:)	          
	!...space for this is allocated in prep.f90
    INTEGER, SAVE		:: rigid_constrain_point_eqn(max_rigid_constrain_point)

    !...max integration points
    INTEGER, SAVE		:: max_nint = 100 
	!...for Gauss integration
    INTEGER, SAVE		:: max_logint = 16 
	!...for logarith integration
    !...integration points and weights
    REAL(KIND=DBL), ALLOCATABLE, SAVE	:: xi(:,:),wi(:,:)				 
	!...integration points and weights, gaussian
    REAL(KIND=DBL), ALLOCATABLE, SAVE	:: xi_log(:,:),wi_log(:,:)		 
	!...integration points and weights, log gaussian

    !...number of nodes associated with elemid
    INTEGER, SAVE		:: NODE(8) = (/3,3,3,2,3,2,3,2/)				 
	!...number of nodes associated with each of those element types
    !...type of element associated with elemid: !11/16/09: add element on inner boundary
    INTEGER, SAVE		:: ELTYPE(8) = (/CTIP,CTIP,CREGULAR,CREGULAR,BREGULAR,&
    									BREGULAR,INTERFE,INTERFE/)		 
										!...type of element associated with each of those element IDs
																		 
	!...material constants
	REAL(KIND=DBL), ALLOCATABLE, SAVE   	:: material(:,:)             
	!...space for this is allocated in prep.f90    

    !...integration order
    !...number of Gauss points for outer/inner/logarith integral:
    INTEGER, SAVE		:: ng_out,ng_in,ng_ln							 
	!...number of Gauss points for outer/inner/log intregral
    !...number of Gauss points for load vector integral:
    INTEGER, SAVE		:: ng_ef										 
	!...number of Gauss points for load vector integral
    !...number of Gauss points for outer/inner/logarith integral for T-stress calculation
    INTEGER, SAVE		:: ngt_out,ngt_in,ngt_ln						 
	!...number of Gauss points for outer/inner/log integral for T-stress calculation
    !...number of Gauss points for load vector integral for T-stress calculation
    INTEGER, SAVE		:: ngt_ef										 
	!...number of Gauss points for load vector integral for T-stress calculation

    !...global shear modulus
    !...(used as a constant make terms in coefficient matrix be in the same order of magnitude)
    !...(this constant is from Xiao's code, and he had to divide the kernels
    !...in isotropic case by gmu, then multiply back gmu with the solution.
    !...However, this constant won't be used in this code -> later make gmu = 1.0)
    !...REAL(KIND=DBL):: gmu
    
    !...Constant Pi (for double precison with total of 16 significant decimal digits)
	REAL(KIND=DBL), PARAMETER :: PI = 3.141592653589793d0

    !...5/11/09: parameters for remeshing
    INTEGER, SAVE		:: growth_flag,total_step=0,nstep
    REAL(KIND=DBL), SAVE:: size_cracktip   
	!...size of crack tip element
    REAL(KIND=DBL), SAVE:: dcheck1_factor,dcheck2_factor   
	!...parameters to control the dividing of elements behind crack-tip element during propogation
    REAL(KIND=DBL), SAVE:: grow_factor 
	!...a coefficient(usually less than 1.0)used to get max advance for the crack tip during growth i.e.max crack advance = growth_factor*size-cracktip
    REAL(KIND=DBL), SAVE:: globalK(3)
    REAL(KIND=DBL), SAVE:: parameter_alpha,parameter_beta,&
    					   parameter_epsilon  
						   !...these factors appear in the bilinear crack growth law
    INTEGER, SAVE		:: pressure_flag
    REAL(KIND=DBL), SAVE:: sigma_x, sigma_y, pressure_value

    REAL(KIND=DBL), ALLOCATABLE :: volume_ratio(:)

    REAL(KIND=DBL), ALLOCATABLE	:: arc_length(:),crack_volume_total(:), crack_volume_tractions(:), crack_volume_pressure(:,:)

    !----------------------------------------------------------------------
	CONTAINS
    Integer Function adjacent(selem,felem)
    !...Subroutine to determine whether selem and felem are adjacent elements)
        Implicit none
        Integer	:: selem,felem
        Integer	:: i,j

        adjacent = 0
        if (selem.eq.felem) return   
		!...exclude coincident case
        do i = 1,NODE(elemid(selem)) 
		!...elemid takes numbers from 1 to 8 ~ depending on the type of element, returns numbers from 1 to 8, and NODE returns number of nodes
          	do j = 1,NODE(elemid(felem))
            	if (elnode(i,selem).eq.elnode(j,felem)) then
                	adjacent = 1
                    return			 
					!...if any of those node numbers are same, then selem and felem are adjacent elements
                endif
            enddo
        enddo
    End Function adjacent
END MODULE GlobalData
