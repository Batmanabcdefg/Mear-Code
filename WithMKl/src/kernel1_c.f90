Subroutine Kernel1(region)
!...Subroutine to calculate integral D2 by dividing into 3 sub-domains and using
!...logarith rule for the 1st integration, Gauss for the 2nd, 
!...and Gauss with stretching for the 3nd integration
!...difference with kernel1_b: delete unneeded material constants E_3alphabeta3 (this is wrong,
!...because we need the full matrix so that transformation can be correctly applied => going back to kernel1.b, 
!...but we need to modify prep.f90 so that user only input needed constants for each type of material
!...transform material constants to geometry coordinates (no need to rotate geometry)
	USE DefinitionConstant
    USE GlobalData
    USE Remesh
    IMPLICIT NONE

	INTEGER, INTENT(IN)		:: region
    
	!...local variables
    INTEGER		:: n_theta,nit1,nit2,nln,ierr
    REAL(KIND=DBL)	:: f_theta,theta_step
    REAL(KIND=DBL), ALLOCATABLE	:: point(:),weight(:),pointln(:),weightln(:)
    INTEGER		:: i,j,k,l,p,q,m,n,i1,i2,i3,i4
    INTEGER		:: beta,gamma,eita,rho
    REAL(KIND=DBL)	:: CI1(NDOF,NDOF),UI1(NDOF,NDOF),GI1(NDOF,NDOF)
	REAL(KIND=DBL)	:: CI2_table(NDOF,NDOF,NDOF),UI2_table(NDOF,NDOF,NDOF),GI2_table(NDOF,NDOF,NDOF)
	REAL(KIND=DBL)	:: D1(2,2,NDOF,NDOF),D2(2,2,NDOF,NDOF)
    REAL(KIND=DBL)	:: D2a(2,2,NDOF,NDOF),D2b(2,2,NDOF,NDOF)
    REAL(KIND=DBL)	:: D3(2,2,NDOF,NDOF),D4(2,2,NDOF,NDOF),D5(2,2,NDOF,NDOF),D6(2,2,NDOF,NDOF)
    REAL(KIND=DBL)	:: D7(2,2,NDOF,NDOF)
    REAL(KIND=DBL)	:: f(2,2,NDOF,NDOF),h(2,2,NDOF,NDOF)
    REAL(KIND=DBL)	:: phi,theta,z(2)
    REAL(KIND=DBL)	:: zz(NDOF,NDOF),zzin(NDOF,NDOF)
    REAL(KIND=DBL)	:: temp,abound,bbound,c,t,tprime
    REAL(KIND=DBL)	:: A(NDOF,NDOF,2,2,NDOF,NDOF)
    INTEGER			:: matl_no,matl_type
    REAL(KIND=DBL)	:: E(3,NDOF,NDOF,3),E_temp(3,NDOF,NDOF,3),dir_cosine(NDOF,NDOF)
    REAL(KIND=DBL)	:: ey1,pv1,C11,C12,C13,C14,C15,C16,&
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

	!--------------------------------------------------------------------------------
    !...temporarily set values of theta for all cases
    n_theta = 181
    f_theta = PI
    theta_step = f_theta/dble(n_theta-1)								
	!...theta_step = pi/180, increment of 1 degree in radian
    !...number of Gauss points to compute integration of kernel
	nit1 = 60 !...for D1 calculation
    nit2 = 60 !...for D2 calculation
    nln = 16  !...for D2 calculation
    !...allocate Gauss points and weights for D1 calculation
	allocate(point(nit1),weight(nit1),STAT=ierr)
	if (ierr.ne.0) then
    	print*,'point,weight: allocation request denied'
        stop
	endif
    !CALL GaussIntPoint(nit1,point,weight)
    if ((nit1.gt.max_nint).or.(nit2.gt.max_nint)) then	                !max_nint is the maximum number of Gauss integration points
      	print*,'kernel1_a.f90: nit1 or nit2 greater than max_nint'
        stop
    endif
    do i = 1,nit1
      	point(i) = xi(i,nit1)											!xi's are the gauss points
        weight(i) = wi(i,nit1)											!wi's are the weights
    enddo
    !--------------------------------------------------------------------------------
	!...id # of material
    matl_no = region_mat_no(region)
    !...type of material
    matl_type = material_type(matl_no)
	!...set material properties and necessary information
    E = 0.0d0
!----------------------------------------------------------------------------------------------------
    !Complete moduli E_iJKl from independent constants of material						   
	!populating the E_{ijkl} matrix in this module
    if (matl_type.eq.1) then															   
	!...material type 1 means isotropy, 2 independent material constants
    !...isotropy
      	if ((NDOF.eq.3).or.(NDOF.eq.2)) then
        !...for isotropy, plane strain and antiplane are uncoupled, and the closed form of
        !kernels are independent
        	ey1 = material(1,matl_no)   !...Young's modulus
        	pv1 = material(2,matl_no)   !...Poisson ratio
        else
          	print*,'not develop yet isotropy for other media than elasticity'
            stop
        endif
    elseif (matl_type.eq.2) then														   
	!material type 2 means cubic material, 3 independent material constants
	!...Cubic material
    	if (NDOF.eq.3) then
        !...elastic media
        	C11 = material(1,matl_no)
        	C12 = material(2,matl_no)
        	C44 = material(3,matl_no)
			E(1,1,1,1)=C11
	    	E(1,1,2,2)=C12
	    	E(1,1,3,3)=C12
	    	E(2,2,1,1)=C12
	    	E(2,2,2,2)=C11
	    	E(2,2,3,3)=C12
	    	E(3,3,1,1)=C12
	    	E(3,3,2,2)=C12
	    	E(3,3,3,3)=C11
        	!...Jaroon worked for 12->4, but it doesn't matter for this case of cubic material
        	!...so this is also true for our convention 23->4
	    	E(1,2,1,2)=C44
	    	E(1,2,2,1)=C44
	    	E(2,1,1,2)=C44
	    	E(2,1,2,1)=C44
	    	E(1,3,1,3)=C44
	    	E(1,3,3,1)=C44
	    	E(3,1,1,3)=C44
	    	E(3,1,3,1)=C44
	    	E(2,3,2,3)=C44
	    	E(2,3,3,2)=C44
	    	E(3,2,2,3)=C44
	    	E(3,2,3,2)=C44
        else
          	print*,'not developed for cubic material of NDOF=',NDOF
            stop
        endif
	elseif (matl_type.eq.3) then                                                           
	!material type 3 means transversely isotropic material
	!...Transversely isotropy															   
	!NDOF=2 or 3, 5 independent material constants
    	if (NDOF.eq.2) then
            !...plane strain of elastic media
            	!...2 is elastic symmetry: order of input C11,C12,C13,C22,C44 
				!(C44=C66 when 2 is elastic symmetric axis)
				C11 = material(1,matl_no)
        		C12 = material(2,matl_no)
            	!...Note: No need to have C13 (Need to modify prep.f90 also!)
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
        	E(1,3,3,1) = 0.5d0*(C11-C13)
        	E(1,3,1,3) = E(3,1,3,1)
        	E(1,2,1,2) = C44
        	E(1,2,2,1) = C44
        	E(2,1,1,2) = C44
        	E(2,1,2,1) = C44
        elseif (NDOF.eq.4) then									
		!...NDOF=4, 10 independent material constants
        !...piezoelectric media: 2 is elastic symmetry
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
			!NDOF=5, 17 independent material constants
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
        	print*,'not coding yet for transversely isotropy of media:',NDOF
        	stop
        endif
	elseif (matl_type.eq.4) then                                               
	!material type means orthotropic material
    	!...orthotropy														   
		!NDOF=2, 6 independent material constants
        if (NDOF.eq.2) then
        	!...plane strain of elastic media
        	C11 = material(1,matl_no)
        	C12 = material(2,matl_no)
        	C22 = material(3,matl_no)
            !...No need to have C44 and C55 (need to modify prep.f90 also)
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
			!NDOF=3 or 4, 9 independent material constants
        elseif (NDOF.eq.3) then
		!...elastic media, order of input C11,C12,C13,C22,C23,C33,C44,C55,C66
        	!...6/26/09: modify 23->4, 31->5, 12->6
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
	!material type 5 means monoclinic material
    !...monoclinic material												
	!NDOF=2, 9 independent material constants
    	if (NDOF.eq.2) then
		!...elastic media, plane strain: z=0 is plane of symmetry
    		!...7/16/09: modify 12->6, 23->4, 31->5
    		C11 = material(1,matl_no)
        	C12 = material(2,matl_no)
        	C16 = material(3,matl_no)        
        	C22 = material(4,matl_no)
        	C26 = material(5,matl_no)
            !...No need to have C44,C45,C55
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
		!NDOF=3, 13 independent material constants
		!...elastic media: z=0 is plane of symmetry
    		!...7/16/09: modify 12->6, 23->4, 31->5
    		C11 = material(1,matl_no)
        	C12 = material(2,matl_no)
            !...No need to have C13,C23,C33,C36
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
	!material type 6 means general anisotropic material
    !...general anisotropy
    	if (NDOF.eq.1) then												   
		!NDOF=1, 6 independent material constants
        	C11 = material(1,matl_no)
            C12 = material(2,matl_no)
            !...No need to have C13,C23,C33
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
		!NDOF=2, 15 independent material constants
        	!...plane strain of elastic media
            C11 = material(1,matl_no)
        	C12 = material(2,matl_no)
            !...No need to have C14,C15,C24,C25,C44,C45,C46,C55,C56
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
    	!...aside: since the programming is too long to repeat Eijkl for piezoelectric
        !we combine both elastic and piezoelectric with common Eijkl here
        !then Ei4kl or Eij4l separately for piezoelectric
    	elseif ((NDOF.eq.3).or.(NDOF.eq.4)) then                                        
		!NDOF=3 or 4, 21 independent material constants
			!...elastic media, convention: 23->4,31->5,12->6
        	!...Order of input: see following notation
    		C11 = material(1,matl_no)
        	C12 = material(2,matl_no)
            !...no need to have C13,C23,C33,C34,C35,C36
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
                !...no need to have e13,e23,e31,e32,e33,e34,e35,e36,k13,k23,k33
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
!----------------------------------------------------------------------------------------------------
    !...1/5/09: transform material constants to geometry coordinates
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
            do k = 3,NDOF
              	dir_cosine(1,k) = 0.d0
                dir_cosine(2,k) = 0.d0
            enddo
            !...for the 3th, 4th,...NDOFth row
            do k = 3,NDOF
          		do n = 1,NDOF
            		if (k.eq.n) then
                		dir_cosine(k,n) = 1.0d0
                	else
                  		dir_cosine(k,n) = 0.0d0
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
                            	do k = 1,NDOF
                                	do l = 1,3
                                    	E_temp(p,q,m,n) = E_temp(p,q,m,n) + &
                                        dir_cosine(i,p)*dir_cosine(j,q)*dir_cosine(k,m)*dir_cosine(l,n)*E(i,j,k,l)
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
!print*,'E and E_temp'
!do p = 1,3
!  	do q = 1,NDOF
!    	do m = 1,NDOF
!        	do n = 1,3
!            	print'(f15.5,2x,f15.5)',E(p,q,m,n),E_temp(p,q,m,n)
!            enddo
!        enddo
!    enddo
!enddo

!----------------------------------------------------------------------------------------------------
	!...initialize kernels (that are global data, allocated in prep.f90)
	CI1 = 0.0d0
    CI2_table = 0.0d0
    GI1 = 0.0d0
    GI2_table = 0.0d0
    UI1 = 0.0d0
    UI2_table = 0.0d0
!---------------------------------------------------------------------------------------------------
    if (matl_type.eq.1) then                                              
	!for isotropic material, we directly get CI1,GI1,UI1
      	!...isotropy,
        
		if(NDOF.eq.2) then
		    !...elastic media
      		temp = 0.25d0*ey1/(PI*(1.0d0-pv1*pv1))	                      
			!ey1 is Young's Modulus and pv1 is Poisson's Ratio
      		!...compute CI1 (see note for explanation)
      		CI1(1,1) = temp
        	CI1(2,2) = temp												  !E/(4*pi*(1-\nu^2))
        	!CI1(3,3) = 0.25d0*ey1/(PI*(1.0d0+pv1))						  !E/(4*pi*(1+\nu))
        	!...compute GI1
        	temp = 0.25d0/(PI*(1.0d0-pv1))
        	GI1(1,2) = temp*(1.0d0-2.0d0*pv1)							  !E*(1-2\nu)/(4*pi*(1-\nu))
        	GI1(2,1) = -GI1(1,2)
        	!...compute UI1
        	temp = 0.25d0*(1.0d0+pv1)/(PI*ey1*(1.0d0-pv1))
        	UI1(1,1) = temp*(4.0d0*pv1-3.0d0)							  
			!(4*\nu-3)*(1+\nu)/(4*pi*E*(1-\nu))
        	UI1(2,2) = temp*(4.0d0*pv1-3.0d0)
        	!UI1(3,3) = -(1.0d0+pv1)/(PI*ey1)							  !-(1+\nu)/(pi*E)
        	!...NOTE: for isotropy, no need to calculate CI2_table,GI2_table,UI2_table
        	!since they are all in closed form, no need to interpolate 
			!(get directly from subroutine kernel2)
		endif
		
		if (NDOF.eq.3) then
          	!...elastic media
      		temp = 0.25d0*ey1/(PI*(1.0d0-pv1*pv1))	                      
			!ey1 is Young's Modulus and pv1 is Poisson's Ratio
      		!...compute CI1 (see note for explanation)
      		CI1(1,1) = temp
        	CI1(2,2) = temp												  !E/(4*pi*(1-\nu^2))
        	CI1(3,3) = 0.25d0*ey1/(PI*(1.0d0+pv1))						  !E/(4*pi*(1+\nu))
        	!...compute GI1
        	temp = 0.25d0/(PI*(1.0d0-pv1))
        	GI1(1,2) = temp*(1.0d0-2.0d0*pv1)							  !E*(1-2\nu)/(4*pi*(1-\nu))
        	GI1(2,1) = -GI1(1,2)
        	!...compute UI1
        	temp = 0.25d0*(1.0d0+pv1)/(PI*ey1*(1.0d0-pv1))
        	UI1(1,1) = temp*(4.0d0*pv1-3.0d0)							  
			!(4*\nu-3)*(1+\nu)/(4*pi*E*(1-\nu))
        	UI1(2,2) = temp*(4.0d0*pv1-3.0d0)
        	UI1(3,3) = -(1.0d0+pv1)/(PI*ey1)							  !-(1+\nu)/(pi*E)
        	!...NOTE: for isotropy, no need to calculate CI2_table,GI2_table,UI2_table
        	!since they are all in closed form, no need to interpolate 
			!(get directly from subroutine kernel2)
         !else
         !  	print*,'not yet develop for isotropy of piezo/magneto elastic'
      !      stop
        endif
    else
      	!...anisotropy													  !for anything else,we need to do the loops
		!...calculate D1(2,2,NDOF,NDOF)
		D1 = 0.0d0														  !initialising D1
    	!...loop over Gauss points
    	do n=1,nit1
      		!...angle phi = (PI/2)*zeta									  !pi/2*point()
        	phi = 0.5d0*point(n)*PI                                       !point() are the gauss points 
			!...components of z-vector
			z(1)= dcos(phi)
			z(2)= dsin(phi)
			!...compute (z,z) tensor
            zz = 0.d0
			do i1=1,NDOF
				do i2=i1,NDOF
					do i3=1,2
						do i4=1,2										   !zz is computed this way, need to understand Radon Transform for this
                            zz(i1,i2)=zz(i1,i2)+z(i3)*E(i3,i1,i2,i4)*z(i4)        
						enddo
					enddo
				enddo
			enddo
        	!...use symmetry of tensor (z,z)
            do i1 = 2,NDOF
              	do i2 = 1,(i1-1)
                	zz(i1,i2) = zz(i2,i1)							       !zz tensor is symmetric
                enddo
            enddo

			!...calculate inverse of (z,z)
            if (NDOF.eq.1) then
              	zzin = 1.0d0/zz
            elseif (NDOF.eq.2) then
              	call matrix2_inverse(zz,zzin)
            elseif (NDOF.eq.3) then
              	call matrix3_inverse(zz,zzin)
            elseif (NDOF.eq.4) then
            	call matrix4_inverse(zz,zzin)
            elseif (NDOF.eq.5) then
            	call matrix5_inverse(zz,zzin)							   !inverse of matrix zz
            else
              	print*,'kernel1.f90: cannot compute inverse of zz for NDOF=',NDOF
                stop
            endif

			do i1=1,2
				do i2=i1,2
					do i3=1,NDOF
						do i4=i3,NDOF									   !D1 is equation (43) in the computational paper, the part that multiplies equation (42)
							D1(i1,i2,i3,i4)=D1(i1,i2,i3,i4)+z(i1)*z(i2)*zzin(i3,i4)*weight(n)	  
						enddo
					enddo
				enddo
			enddo
		enddo ! end loop contour integral over phi
		!...use symmetry of D1
    	do i1=1,2
      		do i2=1,2
        		if (i2.ge.i1) then
            	!...case: already have data for D1(i1,i2,i3,i4) with i4>=i3
            		do i3=1,NDOF
                		do i4=1,NDOF
                    		if (i4.lt.i3) then
                        		D1(i1,i2,i3,i4) = D1(i1,i2,i4,i3)          !D1 is symmetric
                        	endif
                    	enddo
                	enddo
            	else
            	!...case: not yet have data for D1(i1,i2,i3,i4) at all
            		do i3=1,NDOF
                		do i4=1,NDOF
                    		if (i4.lt.i3) then
                        		D1(i1,i2,i3,i4) = D1(i2,i1,i4,i3)
                        	else
                        		D1(i1,i2,i3,i4) = D1(i2,i1,i3,i4)          !again, D1 is symmetric
                        	endif
                    	enddo
                	enddo
            	endif
        	enddo
    	enddo
    	!...multiply D1 with constant (-1/(4*PI)) which is included jacobian PI/2
        !...(-1/(4*Pi*Pi)*2*(Pi/2) = -1/(4*Pi)											  
    	D1 = (-0.25d0/PI)*D1
    	!...calculate CI1
        A = 0.0d0
    	do j=1,NDOF
        	do k=1,NDOF
            	do beta=1,2
               		do gamma=1,2
                   		do p=1,NDOF
                       		do m=1,NDOF									   !this is equation (17) in computational paper
                           		do eita=1,2
                               		A(j,k,beta,gamma,p,m) = A(j,k,beta,gamma,p,m) + &
                                    	E(eita,k,p,beta)*E(eita,j,m,gamma) - &
                                    	E(eita,j,k,eita)*E(beta,p,m,gamma)/dble(NDOF)  !or real(NDOF)   
                            	enddo
                           		CI1(j,k) = CI1(j,k) + A(j,k,beta,gamma,p,m)*D1(beta,gamma,p,m)			
                        	enddo
                    	enddo											   !this is equation 23(3) in computational paper
                	enddo
            	enddo
        	enddo
    	enddo
        !--------------------------------------------------------------------------------
        !...calculate GI1
    	do j=1,NDOF
        	do p=1,NDOF
            	do k=1,NDOF
               		do rho=1,2
                   		GI1(j,p) = GI1(j,p) + E(1,j,k,rho)*D1(2,rho,p,k) - E(2,j,k,rho)*D1(1,rho,p,k)	
                	enddo												  !this is equation 23(2) in computational paper
            	enddo
        	enddo
    	enddo
        !--------------------------------------------------------------------------------
        !...calculate UI1
    	do j=1,NDOF
        	do p=1,NDOF
            	UI1(j,p) = D1(1,1,p,j) + D1(2,2,p,j)					  !this is equation 23(1) in computational paper
        	enddo
    	enddo
		!--------------------------------------------------------------------------------
        !...7/29/09: reallocate point and weight for calculating D2
        if (nit1.ne.nit2) then
			deallocate(point,weight,STAT=ierr)
        	if (ierr.ne.0) then
          		print*,'kernel1.f90: point and weight deallocation request denied'
            	stop
        	endif
        	allocate(point(nit2),weight(nit2),STAT=ierr)
			if (ierr.ne.0) then	
    			print*,'point,weight: allocation request denied'
        		stop
			endif
            do i = 1,nit2
              	point(i) = xi(i,nit2)									  !coordinates of gauss points
                weight(i) = wi(i,nit2)									  !weights on the gauss points
            enddo
        endif
        allocate(pointln(nln),weightln(nln),STAT=ierr)
		if (ierr.ne.0) then	
    		print*,'pointln,weightln: allocation request denied'
        	stop
		endif
        do i = 1,nln
      		pointln(i) = xi_log(i,nln)									  !coordinates of log points
        	weightln(i) = wi_log(i,nln)									  !weights on the log points
    	enddo
        !...value of abound
        abound = 0.05d0
        bbound = 0.95d0													  !this is a in the computational paper
        c = 1.0d0 - dsqrt(2.0d0*(1.0d0 - bbound))					      !this is equation 54(2) in computational paper, a is chosen to be 0.95
    	!CALL GaussIntPoint(nit2,point,weight)
        !call LogIntPoint(nln,pointln,weightln)
		!...loop on theta-direction: calculate integral D2, then C2,G2,U2 for table of values of theta
		do i=1,n_theta
			!...angle theta: this angle will make table of n_theta values for interpolation of K2 later
			theta=theta_step*dble(i-1)
			D2 = 0.0d0
            !...calculate D2a
            D2a = 0.0d0
            D3 = 0.0d0
            D4 = 0.0d0
            D5 = 0.0d0
            D6 = 0.0d0
            D7 = 0.0d0
            !------------------------------------------------------------
            !...loop over regular Gauss points to calculate D5
			do n=1,nit2																   !we use nit1 for UI1,GI1,CI1 and nit2 for UI2,CI2,UI2
            	!...t' in term of t''
                tprime = 0.5d0*(point(n)+1.0d0)                                        
                !...t in term of t'
                t = abound*tprime
                !------------------------------------------------------------
                !...compute z and inverse of (z,z)
                temp = dsin(dacos(t))
                z(1) = t*dcos(theta) + temp*dsin(theta)
                z(2) = -temp*dcos(theta) + t*dsin(theta)
                zz = 0.0d0
				do i1 = 1,NDOF
					do i2 = i1,NDOF
						do i3 = 1,2
							do i4 = 1,2
								zz(i1,i2) = zz(i1,i2) + z(i3)*E(i3,i1,i2,i4)*z(i4)	   !zz is computed this way, need to understand Radon Transform
							enddo
						enddo
					enddo
				enddo
            	!...use symmetry of tensor (z,z)
            	do i1 = 2,NDOF
              		do i2 = 1,(i1-1)
                		zz(i1,i2) = zz(i2,i1)										   !zz is symmetric
                	enddo
            	enddo
				!...calculate inverse of (z,z)
                if (NDOF.eq.1) then
                  	zzin = 1.0d0/zz
                elseif (NDOF.eq.2) then
                  	call matrix2_inverse(zz,zzin)
                elseif (NDOF.eq.3) then
              		call matrix3_inverse(zz,zzin)
            	elseif (NDOF.eq.4) then
            		call matrix4_inverse(zz,zzin)
            	elseif (NDOF.eq.5) then
            		call matrix5_inverse(zz,zzin)									   !zzin is the inverse of zz
            	else
              		print*,'current version cannot compute inverse of zz for NDOF=',NDOF
                	stop
            	endif
            	!------------------------------------------------------------
				!...calculate f
				do i1 = 1,2
					do i2 = i1,2
						do i3 = 1,NDOF
							do i4 = i3,NDOF
								f(i1,i2,i3,i4) = z(i1)*z(i2)*zzin(i3,i4)			   !this is equation (47) in computational paper
							enddo
						enddo
					enddo
				enddo
                !...calculate h(t)
                if ((dsqrt(1.0d0+t).le.SMALL_NUM).or.(dsqrt(1.0d0-t).le.SMALL_NUM)) then
                  	print*,'kernel1.f90: error, dsqrt(1-t*t) is too small'
                    stop
                endif
                do i1 = 1,2
                  	do i2 = i1,2
                    	do i3 = 1,NDOF
                        	do i4 = i3,NDOF
                            	h(i1,i2,i3,i4) = f(i1,i2,i3,i4)/dsqrt(1.0d0-t*t)	   !h=f/sqrt(1-t^2)
                            enddo
                        enddo
                    enddo
                enddo
                !...calculate D5
                do i1 = 1,2
                  	do i2 = i1,2
                    	do i3 = 1,NDOF
                        	do i4 = i3,NDOF
                            	D5(i1,i2,i3,i4) = D5(i1,i2,i3,i4) + h(i1,i2,i3,i4)*weight(n)   !multiply with weights to get D5
                            enddo
                        enddo
                    enddo
                enddo
			enddo !...of n=1,nit2
            !------------------------------------------------------------
        	!...complete D5 by using symmetry of D5
			do i1=1,2
				do i2=1,2
					if (i2.ge.i1) then
                    !...case: already have data for D5(i1,i2,i3,i4) with i4>=i3
						do i3=1,NDOF
							do i4=1,NDOF
								if (i4.lt.i3) then
						    		D5(i1,i2,i3,i4) = D5(i1,i2,i4,i3)
						  		endif
							enddo
						enddo
					else
                    !...case: not yet have data for D5(i1,i2,i3,i4) at all
						do i3=1,NDOF
							do i4=1,NDOF
								if (i4.lt.i3) then
									D5(i1,i2,i3,i4)=D5(i2,i1,i4,i3)
								else
							    	D5(i1,i2,i3,i4)=D5(i2,i1,i3,i4)							   !D5 is symmetric, D5 is f_{j\beta}^{I\alpha}/sqrt(1-t^2)
								endif
							enddo
					 	enddo
					endif
				enddo
			enddo
        	!...multiply D5 with constant											  !this is how integral of f_{j\beta}^{I\alpha}*ln(t)/sqrt(1-t^2) is computed,but didn't
        	D5 = 0.5d0*abound*dlog(abound)*D5										  !understand how 
            !------------------------------------------------------------
            !...loop over logarith points to calculate D6							  !we loop over regular gauss points to calculate D5, we loop over logarithmic points to 
            do n = 1,nln															  !calculate D6
              	!...t in term of t'
                t = abound*pointln(n)
                !------------------------------------------------------------
                !...compute z and inverse of (z,z)
                temp = dsin(dacos(t))
                z(1) = t*dcos(theta) + temp*dsin(theta)
                z(2) = -temp*dcos(theta) + t*dsin(theta)
                zz = 0.0d0
				do i1 = 1,NDOF
					do i2 = i1,NDOF
						do i3 = 1,2
							do i4 = 1,2
								zz(i1,i2) = zz(i1,i2) + z(i3)*E(i3,i1,i2,i4)*z(i4)			 !this is the way zz is computed, need to understand Radon Transform
							enddo
						enddo
					enddo
				enddo
            	!...use symmetry of tensor (z,z)
            	do i1 = 2,NDOF
              		do i2 = 1,(i1-1)
                		zz(i1,i2) = zz(i2,i1)                                          !zz is symmetric
                	enddo
            	enddo
				!...calculate inverse of (z,z)
                if (NDOF.eq.1) then
                  	zzin = 1.0d0/zz
                elseif (NDOF.eq.2) then
                  	call matrix2_inverse(zz,zzin)
                elseif (NDOF.eq.3) then
              		call matrix3_inverse(zz,zzin)
            	elseif (NDOF.eq.4) then
            		call matrix4_inverse(zz,zzin)
            	elseif (NDOF.eq.5) then												  !zzin is inverse of zz
            		call matrix5_inverse(zz,zzin)
            	else
              		print*,'current version cannot compute inverse of zz for NDOF=',NDOF
                	stop
            	endif
            	!--------------------------------------------------------------------------------
				!...calculate f
				do i1 = 1,2
					do i2 = i1,2
						do i3 = 1,NDOF
							do i4 = i3,NDOF
								f(i1,i2,i3,i4) = z(i1)*z(i2)*zzin(i3,i4)			  !this is equation (47) in computational paper
							enddo
						enddo
					enddo
				enddo
                !...calculate h(t)
                if ((dsqrt(1.0d0+t).le.SMALL_NUM).or.(dsqrt(1.0d0-t).le.SMALL_NUM)) then
                  	print*,'kernel1.f90: error, dsqrt(1-t*t) is too small'
                    stop
                endif
                do i1 = 1,2
                  	do i2 = i1,2
                    	do i3 = 1,NDOF
                        	do i4 = i3,NDOF
                            	h(i1,i2,i3,i4) = f(i1,i2,i3,i4)/dsqrt(1.0d0-t*t)	  !similar to D5
                            enddo
                        enddo
                    enddo
                enddo
                !...calculate D6
                do i1 = 1,2
                  	do i2 = i1,2
                    	do i3 = 1,NDOF
                        	do i4 = i3,NDOF
                            	D6(i1,i2,i3,i4) = D6(i1,i2,i3,i4) + h(i1,i2,i3,i4)*weightln(n)		!similar to D5
                            enddo
                        enddo
                    enddo
                enddo
			enddo !...of n=1,nln
            !------------------------------------------------------------
        	!...use symmetry of D6
			do i1=1,2
				do i2=1,2
					if (i2.ge.i1) then
                    !...case: already have data for D6(i1,i2,i3,i4) with i4>=i3
						do i3=1,NDOF
							do i4=1,NDOF
								if (i4.lt.i3) then
						    		D6(i1,i2,i3,i4) = D6(i1,i2,i4,i3)
						  		endif
							enddo
						enddo
					else
                    !...case: not yet have data for D6(i1,i2,i3,i4) at all
						do i3=1,NDOF
							do i4=1,NDOF
								if (i4.lt.i3) then
									D6(i1,i2,i3,i4)=D6(i2,i1,i4,i3)
								else
							    	D6(i1,i2,i3,i4)=D6(i2,i1,i3,i4)									!similar to D5
								endif
							enddo
					 	enddo
					endif
				enddo
			enddo
        	!...multiply D2a12 with constant
        	D6 = -1.0d0*abound*D6																	!what is this constant?
            !------------------------------------------------------------
            D3 = D5 + D6
            !------------------------------------------------------------
            !...loop over regular Gauss points to calculate D4
            do n = 1,nit2
              	!...t' in term of t''																!
                tprime = 0.5d0*((1.0d0 - c)*point(n) + 1.0d0 + c)
                !...t in term of t'
                t = 1.0d0 - 0.5d0*(1.0d0 - tprime)*(1.0d0 - tprime)									!this is equation (52) in computational paper
                !------------------------------------------------------------
                !...compute z and inverse of (z,z)
                temp = dsin(dacos(t))
                z(1) = t*dcos(theta) + temp*dsin(theta)
                z(2) = -temp*dcos(theta) + t*dsin(theta)
                zz = 0.0d0
				do i1 = 1,NDOF
					do i2 = i1,NDOF
						do i3 = 1,2
							do i4 = 1,2
								zz(i1,i2) = zz(i1,i2) + z(i3)*E(i3,i1,i2,i4)*z(i4)				   !again, the way zz is computed, need to understand Radon transform
							enddo
						enddo
					enddo
				enddo
            	!...use symmetry of tensor (z,z)
            	do i1 = 2,NDOF
              		do i2 = 1,(i1-1)
                		zz(i1,i2) = zz(i2,i1)													   !symmetry of zz
                	enddo
            	enddo
				!...calculate inverse of (z,z)
                if (NDOF.eq.1) then
                  	zzin = 1.0d0/zz
                elseif (NDOF.eq.2) then
                  	call matrix2_inverse(zz,zzin)
                elseif (NDOF.eq.3) then
              		call matrix3_inverse(zz,zzin)
            	elseif (NDOF.eq.4) then
            		call matrix4_inverse(zz,zzin)
            	elseif (NDOF.eq.5) then
            		call matrix5_inverse(zz,zzin)												   !zzin is the inverse of zz
            	else
              		print*,'current version cannot compute inverse of zz for NDOF=',NDOF
                	stop
            	endif
            	!--------------------------------------------------------------------------------
				!...calculate f
				do i1 = 1,2
					do i2 = i1,2
						do i3 = 1,NDOF
							do i4 = i3,NDOF
								f(i1,i2,i3,i4) = z(i1)*z(i2)*zzin(i3,i4)						  !f is equation (47) in computational paper
							enddo
						enddo
					enddo
				enddo
                !...calculate h(t)
                if ((dsqrt(1.0d0+t).le.SMALL_NUM).or.(t.le.SMALL_NUM)) then
                  	print*,'kernel1.f90: error, either dsqrt(1+t) or t is too small'
                    stop
                endif
                do i1 = 1,2
                  	do i2 = i1,2
                    	do i3 = 1,NDOF
                        	do i4 = i3,NDOF
                            	h(i1,i2,i3,i4) = f(i1,i2,i3,i4)*dlog(t)/dsqrt(1.0d0+t)           !f*log(t)/sqrt(1+t)
                            enddo
                        enddo
                    enddo
                enddo
                !...calculate D4
                do i1 = 1,2
                  	do i2 = i1,2
                    	do i3 = 1,NDOF
                        	do i4 = i3,NDOF
                            	D4(i1,i2,i3,i4) = D4(i1,i2,i3,i4) + h(i1,i2,i3,i4)*weight(n)	!multiplying h with weights to obtain D4
                            enddo
                        enddo
                    enddo
                enddo
            enddo !...of n = 1,nit2 for calculating D4
            !------------------------------------------------------------
            !...use symmetry of D4
			do i1=1,2
				do i2=1,2
					if (i2.ge.i1) then
                    !...case: already have data for D4(i1,i2,i3,i4) with i4>=i3
						do i3=1,NDOF
							do i4=1,NDOF
								if (i4.lt.i3) then
						    		D4(i1,i2,i3,i4)=D4(i1,i2,i4,i3)
						  		endif
							enddo
						enddo
					else
                    !...case: not yet have data for D4(i1,i2,i3,i4) at all
						do i3=1,NDOF
							do i4=1,NDOF
								if (i4.lt.i3) then
									D4(i1,i2,i3,i4)=D4(i2,i1,i4,i3)
								else
							    	D4(i1,i2,i3,i4)=D4(i2,i1,i3,i4)								!symmetry of D4
								endif
							enddo
					 	enddo
					endif
				enddo
			enddo
        	!...multiply D4 with constant
        	D4 = 0.5d0*(1.0d0 - c)*dsqrt(2.0d0)*D4
        	!------------------------------------------------------------
            !...loop over regular Gauss points to calculate D7
            do n = 1,nit2
              	!...t in term of t'
                t = 0.5d0*(bbound-abound)*(point(n)+1.0d0) + abound
                !------------------------------------------------------------
                !...compute z and inverse of (z,z)
                temp = dsin(dacos(t))
                z(1) = t*dcos(theta) + temp*dsin(theta)
                z(2) = -temp*dcos(theta) + t*dsin(theta)
                zz = 0.0d0
				do i1 = 1,NDOF
					do i2 = i1,NDOF
						do i3 = 1,2
							do i4 = 1,2
								zz(i1,i2) = zz(i1,i2) + z(i3)*E(i3,i1,i2,i4)*z(i4)
							enddo
						enddo
					enddo
				enddo
            	!...use symmetry of tensor (z,z)
            	do i1 = 2,NDOF
              		do i2 = 1,(i1-1)
                		zz(i1,i2) = zz(i2,i1)
                	enddo
            	enddo
				!...calculate inverse of (z,z)
                if (NDOF.eq.1) then
                  	zzin = 1.0d0/zz
                elseif (NDOF.eq.2) then
                  	call matrix2_inverse(zz,zzin)
                elseif (NDOF.eq.3) then
              		call matrix3_inverse(zz,zzin)
            	elseif (NDOF.eq.4) then
            		call matrix4_inverse(zz,zzin)
            	elseif (NDOF.eq.5) then
            		call matrix5_inverse(zz,zzin)
            	else
              		print*,'current version cannot compute inverse of zz for NDOF=',NDOF
                	stop
            	endif
            	!------------------------------------------------------------
				!...calculate f
				do i1 = 1,2
					do i2 = i1,2
						do i3 = 1,NDOF
							do i4 = i3,NDOF
								f(i1,i2,i3,i4) = z(i1)*z(i2)*zzin(i3,i4)
							enddo
						enddo
					enddo
				enddo
                !...calculate h(t)
                if ((dsqrt(1.0d0-t*t).le.SMALL_NUM).or.(t.le.SMALL_NUM)) then
                  	print*,'kernel1.f90: error, dsqrt(1-t*t) or t is too small'
                    stop
                endif
                do i1 = 1,2
                  	do i2 = i1,2
                    	do i3 = 1,NDOF
                        	do i4 = i3,NDOF
                            	h(i1,i2,i3,i4) = f(i1,i2,i3,i4)*dlog(t)/dsqrt(1.0d0-t*t)
                            enddo
                        enddo
                    enddo
                enddo
                !...calculate D7
                do i1 = 1,2
                  	do i2 = i1,2
                    	do i3 = 1,NDOF
                        	do i4 = i3,NDOF
                            	D7(i1,i2,i3,i4) = D7(i1,i2,i3,i4) + h(i1,i2,i3,i4)*weight(n)
                            enddo
                        enddo
                    enddo
                enddo
			enddo !...of n=1,nit2
            !------------------------------------------------------------
            !...use symmetry of D7
			do i1=1,2
				do i2=1,2
					if (i2.ge.i1) then
                    !...case: already have data for D7(i1,i2,i3,i4) with i4>=i3
						do i3=1,NDOF
							do i4=1,NDOF
								if (i4.lt.i3) then
						    		D7(i1,i2,i3,i4)=D7(i1,i2,i4,i3)
						  		endif
							enddo
						enddo
					else
                    !...case: not yet have data for D7(i1,i2,i3,i4) at all
						do i3=1,NDOF
							do i4=1,NDOF
								if (i4.lt.i3) then
									D7(i1,i2,i3,i4)=D7(i2,i1,i4,i3)
								else
							    	D7(i1,i2,i3,i4)=D7(i2,i1,i3,i4)
								endif
							enddo
					 	enddo
					endif
				enddo
			enddo
        	!...multiply D7 with constant
        	D7 = 0.5d0*(bbound - abound)*D7
            !------------------------------------------------------------
            D2a = D3 + D4 + D7
            !------------------------------------------------------------
			!...calculate D2b
            D2b = 0.0d0
            D3 = 0.0d0
            D4 = 0.0d0
            D5 = 0.0d0
            D6 = 0.0d0
            D7 = 0.0d0
            !------------------------------------------------------------
            !...loop over regular Gauss points to calculate D5
			do n=1,nit2
            	!...t' in term of t''
                tprime = 0.5d0*(point(n)+1.0d0)
                !...t in term of t'
                t = abound*tprime
                !------------------------------------------------------------
                !...compute z and inverse of (z,z)
                temp = dsin(dacos(t))
                z(1) = t*dcos(theta) - temp*dsin(theta)
                z(2) = temp*dcos(theta) + t*dsin(theta)
                zz = 0.0d0
				do i1 = 1,NDOF
					do i2 = i1,NDOF
						do i3 = 1,2
							do i4 = 1,2
								zz(i1,i2) = zz(i1,i2) + z(i3)*E(i3,i1,i2,i4)*z(i4)
							enddo
						enddo
					enddo
				enddo
            	!...use symmetry of tensor (z,z)
            	do i1 = 2,NDOF
              		do i2 = 1,(i1-1)
                		zz(i1,i2) = zz(i2,i1)
                	enddo
            	enddo
				!...calculate inverse of (z,z)
                if (NDOF.eq.1) then
                  	zzin = 1.0d0/zz
                elseif (NDOF.eq.2) then
                  	call matrix2_inverse(zz,zzin)
                elseif (NDOF.eq.3) then
              		call matrix3_inverse(zz,zzin)
            	elseif (NDOF.eq.4) then
            		call matrix4_inverse(zz,zzin)
            	elseif (NDOF.eq.5) then
            		call matrix5_inverse(zz,zzin)
            	else
              		print*,'current version cannot compute inverse of zz for NDOF=',NDOF
                	stop
            	endif
            	!------------------------------------------------------------
				!...calculate f
				do i1 = 1,2
					do i2 = i1,2
						do i3 = 1,NDOF
							do i4 = i3,NDOF
								f(i1,i2,i3,i4) = z(i1)*z(i2)*zzin(i3,i4)
							enddo
						enddo
					enddo
				enddo
                !...calculate h(t)
                if ((dsqrt(1.0d0+t).le.SMALL_NUM).or.(dsqrt(1.0d0-t).le.SMALL_NUM)) then
                  	print*,'kernel1.f90: error, dsqrt(1-t*t) is too small'
                    stop
                endif
                do i1 = 1,2
                  	do i2 = i1,2
                    	do i3 = 1,NDOF
                        	do i4 = i3,NDOF
                            	h(i1,i2,i3,i4) = f(i1,i2,i3,i4)/dsqrt(1.0d0-t*t)
                            enddo
                        enddo
                    enddo
                enddo
                !...calculate D5
                do i1 = 1,2
                  	do i2 = i1,2
                    	do i3 = 1,NDOF
                        	do i4 = i3,NDOF
                            	D5(i1,i2,i3,i4) = D5(i1,i2,i3,i4) + h(i1,i2,i3,i4)*weight(n)
                            enddo
                        enddo
                    enddo
                enddo
			enddo !...of n=1,nit2
            !------------------------------------------------------------
        	!...use symmetry of D5
			do i1=1,2
				do i2=1,2
					if (i2.ge.i1) then
                    !...case: already have data for D5(i1,i2,i3,i4) with i4>=i3
						do i3=1,NDOF
							do i4=1,NDOF
								if (i4.lt.i3) then
						    		D5(i1,i2,i3,i4)=D5(i1,i2,i4,i3)
						  		endif
							enddo
						enddo
					else
                    !...case: not yet have data for D5(i1,i2,i3,i4) at all
						do i3=1,NDOF
							do i4=1,NDOF
								if (i4.lt.i3) then
									D5(i1,i2,i3,i4)=D5(i2,i1,i4,i3)
								else
							    	D5(i1,i2,i3,i4)=D5(i2,i1,i3,i4)
								endif
							enddo
					 	enddo
					endif
				enddo
			enddo
        	!...multiply D5 with constant
        	D5 = 0.5d0*abound*dlog(abound)*D5
            !------------------------------------------------------------
            !...loop over logarith points to calculate D6
            do n = 1,nln
              	!...t in term of t'
                t = abound*pointln(n)
                !------------------------------------------------------------
                !...compute z and inverse of (z,z)
                temp = dsin(dacos(t))
                z(1) = t*dcos(theta) - temp*dsin(theta)
                z(2) = temp*dcos(theta) + t*dsin(theta)
                zz = 0.0d0
				do i1 = 1,NDOF
					do i2 = i1,NDOF
						do i3 = 1,2
							do i4 = 1,2
								zz(i1,i2) = zz(i1,i2) + z(i3)*E(i3,i1,i2,i4)*z(i4)
							enddo
						enddo
					enddo
				enddo
            	!...use symmetry of tensor (z,z)
            	do i1 = 2,NDOF
              		do i2 = 1,(i1-1)
                		zz(i1,i2) = zz(i2,i1)
                	enddo
            	enddo
				!...calculate inverse of (z,z)
                if (NDOF.eq.1) then
                  	zzin = 1.0d0/zz
                elseif (NDOF.eq.2) then
                  	call matrix2_inverse(zz,zzin)
                elseif (NDOF.eq.3) then
              		call matrix3_inverse(zz,zzin)
            	elseif (NDOF.eq.4) then
            		call matrix4_inverse(zz,zzin)
            	elseif (NDOF.eq.5) then
            		call matrix5_inverse(zz,zzin)
            	else
              		print*,'current version cannot compute inverse of zz for NDOF=',NDOF
                	stop
            	endif
            	!--------------------------------------------------------------------------------
				!...calculate f
				do i1 = 1,2
					do i2 = i1,2
						do i3 = 1,NDOF
							do i4 = i3,NDOF
								f(i1,i2,i3,i4) = z(i1)*z(i2)*zzin(i3,i4)
							enddo
						enddo
					enddo
				enddo
                !...calculate h(t)
                if ((dsqrt(1.0d0+t).le.SMALL_NUM).or.(dsqrt(1.0d0-t).le.SMALL_NUM)) then
                  	print*,'kernel1.f90: error, dsqrt(1-t*t) is too small'
                    stop
                endif
                do i1 = 1,2
                  	do i2 = i1,2
                    	do i3 = 1,NDOF
                        	do i4 = i3,NDOF
                            	h(i1,i2,i3,i4) = f(i1,i2,i3,i4)/dsqrt(1.0d0-t*t)
                            enddo
                        enddo
                    enddo
                enddo
                !...calculate D6
                do i1 = 1,2
                  	do i2 = i1,2
                    	do i3 = 1,NDOF
                        	do i4 = i3,NDOF
                            	D6(i1,i2,i3,i4) = D6(i1,i2,i3,i4) + h(i1,i2,i3,i4)*weightln(n)
                            enddo
                        enddo
                    enddo
                enddo
			enddo !...of n=1,nln
            !------------------------------------------------------------
        	!...use symmetry of D6
			do i1=1,2
				do i2=1,2
					if (i2.ge.i1) then
                    !...case: already have data for D6(i1,i2,i3,i4) with i4>=i3
						do i3=1,NDOF
							do i4=1,NDOF
								if (i4.lt.i3) then
						    		D6(i1,i2,i3,i4) = D6(i1,i2,i4,i3)
						  		endif
							enddo
						enddo
					else
                    !...case: not yet have data for D6(i1,i2,i3,i4) at all
						do i3=1,NDOF
							do i4=1,NDOF
								if (i4.lt.i3) then
									D6(i1,i2,i3,i4) = D6(i2,i1,i4,i3)
								else
							    	D6(i1,i2,i3,i4) = D6(i2,i1,i3,i4)
								endif
							enddo
					 	enddo
					endif
				enddo
			enddo
        	!...multiply D6 with constant
        	D6 = -1.0d0*abound*D6
            !------------------------------------------------------------
            D3 = D5 + D6
            !------------------------------------------------------------
            !...loop over regular Gauss points to calculate D4
            do n = 1,nit2
              	!...t' in term of t''
                tprime = 0.5d0*((1.0d0 - c)*point(n) + 1.0d0 + c)
                !...t in term of t'
                t = 1.0d0 - 0.5d0*(1.0d0 - tprime)*(1.0d0 - tprime)
                !------------------------------------------------------------
                !...compute z and inverse of (z,z)
                temp = dsin(dacos(t))
                z(1) = t*dcos(theta) - temp*dsin(theta)
                z(2) = temp*dcos(theta) + t*dsin(theta)
                zz = 0.0d0
				do i1 = 1,NDOF
					do i2 = i1,NDOF
						do i3 = 1,2
							do i4 = 1,2
								zz(i1,i2) = zz(i1,i2) + z(i3)*E(i3,i1,i2,i4)*z(i4)
							enddo
						enddo
					enddo
				enddo
            	!...use symmetry of tensor (z,z)
            	do i1 = 2,NDOF
              		do i2 = 1,(i1-1)
                		zz(i1,i2) = zz(i2,i1)
                	enddo
            	enddo
				!...calculate inverse of (z,z)
                if (NDOF.eq.1) then
                  	zzin = 1.0d0/zz
                elseif (NDOF.eq.2) then
                  	call matrix2_inverse(zz,zzin)
                elseif (NDOF.eq.3) then
              		call matrix3_inverse(zz,zzin)
            	elseif (NDOF.eq.4) then
            		call matrix4_inverse(zz,zzin)
            	elseif (NDOF.eq.5) then
            		call matrix5_inverse(zz,zzin)
            	else
              		print*,'current version cannot compute inverse of zz for NDOF=',NDOF
                	stop
            	endif
            	!--------------------------------------------------------------------------------
				!...calculate f
				do i1 = 1,2
					do i2 = i1,2
						do i3 = 1,NDOF
							do i4 = i3,NDOF
								f(i1,i2,i3,i4) = z(i1)*z(i2)*zzin(i3,i4)
							enddo
						enddo
					enddo
				enddo
                !...calculate h(t)
                if ((dsqrt(1.0d0+t).le.SMALL_NUM).or.(t.le.SMALL_NUM)) then
                  	print*,'kernel1.f90: error, either dsqrt(1+t) or t is too small'
                    stop
                endif
                do i1 = 1,2
                  	do i2 = i1,2
                    	do i3 = 1,NDOF
                        	do i4 = i3,NDOF
                            	h(i1,i2,i3,i4) = f(i1,i2,i3,i4)*dlog(t)/dsqrt(1.0d0+t)
                            enddo
                        enddo
                    enddo
                enddo
                !...calculate D4
                do i1 = 1,2
                  	do i2 = i1,2
                    	do i3 = 1,NDOF
                        	do i4 = i3,NDOF
                            	D4(i1,i2,i3,i4) = D4(i1,i2,i3,i4) + h(i1,i2,i3,i4)*weight(n)
                            enddo
                        enddo
                    enddo
                enddo
            enddo !...of n = 1,nit2 for calculating D4
            !------------------------------------------------------------
            !...use symmetry of D4
			do i1=1,2
				do i2=1,2
					if (i2.ge.i1) then
                    !...case: already have data for D4(i1,i2,i3,i4) with i4>=i3
						do i3=1,NDOF
							do i4=1,NDOF
								if (i4.lt.i3) then
						    		D4(i1,i2,i3,i4)=D4(i1,i2,i4,i3)								    !symmetry of D4
						  		endif
							enddo
						enddo
					else
                    !...case: not yet have data for D4(i1,i2,i3,i4) at all
						do i3=1,NDOF
							do i4=1,NDOF
								if (i4.lt.i3) then
									D4(i1,i2,i3,i4) = D4(i2,i1,i4,i3)
								else
							    	D4(i1,i2,i3,i4) = D4(i2,i1,i3,i4)								!symmetry of D4
								endif
							enddo
					 	enddo
					endif
				enddo
			enddo
        	!...multiply D4 with constant
        	D4 = 0.5d0*(1.0d0 - c)*dsqrt(2.0d0)*D4
        	!------------------------------------------------------------
            !...loop over regular Gauss points to calculate D7
            do n = 1,nit2
              	!...t in term of t'
                t = 0.5d0*(bbound-abound)*(point(n)+1.0d0) + abound							
                !------------------------------------------------------------
                !...compute z and inverse of (z,z)
                temp = dsin(dacos(t))														
                z(1) = t*dcos(theta) - temp*dsin(theta)
                z(2) = temp*dcos(theta) + t*dsin(theta)
                zz = 0.0d0
				do i1 = 1,NDOF
					do i2 = i1,NDOF
						do i3 = 1,2
							do i4 = 1,2
								zz(i1,i2) = zz(i1,i2) + z(i3)*E(i3,i1,i2,i4)*z(i4)					 !zz is computed this way, for more, got to understand Radon Transform
							enddo
						enddo
					enddo
				enddo
            	!...use symmetry of tensor (z,z)
            	do i1 = 2,NDOF
              		do i2 = 1,(i1-1)
                		zz(i1,i2) = zz(i2,i1)														 !symmetry of zz
                	enddo
            	enddo
				!...calculate inverse of (z,z)
                if (NDOF.eq.1) then
                  	zzin = 1.0d0/zz
                elseif (NDOF.eq.2) then
                  	call matrix2_inverse(zz,zzin)
                elseif (NDOF.eq.3) then
              		call matrix3_inverse(zz,zzin)
            	elseif (NDOF.eq.4) then
            		call matrix4_inverse(zz,zzin)
            	elseif (NDOF.eq.5) then
            		call matrix5_inverse(zz,zzin)													 !zzin is the inverse of zz
            	else
              		print*,'current version cannot compute inverse of zz for NDOF=',NDOF
                	stop
            	endif
            	!------------------------------------------------------------
				!...calculate f
				do i1 = 1,2
					do i2 = i1,2
						do i3 = 1,NDOF
							do i4 = i3,NDOF
								f(i1,i2,i3,i4) = z(i1)*z(i2)*zzin(i3,i4)							 !f is equation (47) in computational paper
							enddo
						enddo
					enddo
				enddo
                !...calculate h(t)
                if ((dsqrt(1.0d0-t*t).le.SMALL_NUM).or.(t.le.SMALL_NUM)) then
                  	print*,'kernel1.f90: error, dsqrt(1-t*t) or t is too small'
                    stop
                endif
                do i1 = 1,2
                  	do i2 = i1,2
                    	do i3 = 1,NDOF
                        	do i4 = i3,NDOF
                            	h(i1,i2,i3,i4) = f(i1,i2,i3,i4)*dlog(t)/dsqrt(1.0d0-t*t)			 !f*log(t)/sqrt(1-t^2)
                            enddo
                        enddo
                    enddo
                enddo
                !...calculate D7
                do i1 = 1,2
                  	do i2 = i1,2
                    	do i3 = 1,NDOF
                        	do i4 = i3,NDOF
                            	D7(i1,i2,i3,i4) = D7(i1,i2,i3,i4) + h(i1,i2,i3,i4)*weight(n)		 !Computing D7 by multiplying with weights
                            enddo
                        enddo
                    enddo
                enddo
			enddo !...of n=1,nit2
            !------------------------------------------------------------
            !...use symmetry of D7
			do i1=1,2
				do i2=1,2
					if (i2.ge.i1) then
                    !...case: already have data for D7(i1,i2,i3,i4) with i4>=i3
						do i3=1,NDOF
							do i4=1,NDOF
								if (i4.lt.i3) then
						    		D7(i1,i2,i3,i4)=D7(i1,i2,i4,i3)									 !symmetry of D7
						  		endif
							enddo
						enddo
					else
                    !...case: not yet have data for D7(i1,i2,i3,i4) at all
						do i3=1,NDOF
							do i4=1,NDOF
								if (i4.lt.i3) then
									D7(i1,i2,i3,i4)=D7(i2,i1,i4,i3)
								else
							    	D7(i1,i2,i3,i4)=D7(i2,i1,i3,i4)									 !symmetry of D7
								endif
							enddo
					 	enddo
					endif
				enddo
			enddo
        	!...multiply D7 with constant															 !What is this constant?
        	D7 = 0.5d0*(bbound - abound)*D7
            !------------------------------------------------------------							 !first check what is D7
            D2b = D3 + D4 + D7
            !------------------------------------------------------------							 !first check what is D2b
            D2 = D2a + D2b
            D2 = -0.5d0/(PI*PI)*D2																	 !1/(2\pi^2)*D2
        	!--------------------------------------------------------------------------------
			!...calculate table CI2_table(k,j,i) where i is the index for value of theta
        	A = 0.0d0
			do j=1,NDOF
        		do k=1,NDOF
                	do beta=1,2
                  		do gamma=1,2																										  
                    		do p=1,NDOF
                        		do m=1,NDOF
                            		do eita=1,2
                                		A(j,k,beta,gamma,p,m) = A(j,k,beta,gamma,p,m) + &
                                    	E(eita,k,p,beta)*E(eita,j,m,gamma) - &						 !equation (17) in computational paper
                                    	E(eita,j,k,eita)*E(beta,p,m,gamma)/dble(NDOF)
                                	enddo
                            		CI2_table(j,k,i) = CI2_table(j,k,i) + A(j,k,beta,gamma,p,m)*&	 !equation (***) in Han's Notes
                                		D2(beta,gamma,p,m)
                            	enddo
                        	enddo
                    	enddo
                	enddo
            	enddo
        	enddo
            !--------------------------------------------------------------------------------
			!...calculate table GI2_table(k,j,i) where i is the index for value of theta
			do j=1,NDOF
        		do p=1,NDOF
                	do k=1,NDOF
                  		do rho=1,2
                       		GI2_table(j,p,i) = GI2_table(j,p,i)+&
                            E(1,j,k,rho)*D2(2,rho,p,k) -E(2,j,k,rho)*D2(1,rho,p,k)					 !equation (**) in Han's Notes
                    	enddo
                	enddo
            	enddo
        	enddo
            !--------------------------------------------------------------------------------
			!...calculate table CI2_table(k,j,i) where i is the index for value of theta
			do j=1,NDOF
        		do p=1,NDOF
                	UI2_table(j,p,i) = D2(1,1,p,j) + D2(2,2,p,j)									 !equation (*) in Han's Notes
            	enddo
        	enddo
            !--------------------------------------------------------------------------------
    	enddo !loop over i=1,n_theta
    endif !of matl_type
End Subroutine Kernel1
