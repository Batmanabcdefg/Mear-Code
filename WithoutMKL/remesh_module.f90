MODULE Remesh
!*************************************************
!***  Module file for remeshing
!*************************************************
USE DefinitionConstant
TYPE NodeData
  INTEGER*4, ALLOCATABLE :: sys2user(:)
  INTEGER, ALLOCATABLE   :: id(:)
  REAL(KIND=DBL), ALLOCATABLE :: coor(:,:)
  INTEGER, ALLOCATABLE   :: ntype(:,:)
  INTEGER, ALLOCATABLE   :: unknown(:,:)
  REAL(KIND=DBL), ALLOCATABLE :: disp_val(:,:)
END TYPE NodeData
TYPE ElemData
  INTEGER*4, ALLOCATABLE :: sys2user(:)
  INTEGER, ALLOCATABLE   :: id(:)
  INTEGER, ALLOCATABLE   :: region(:)
  INTEGER, ALLOCATABLE   :: connect(:,:)
  INTEGER, ALLOCATABLE   :: loadid(:,:)
  REAL(KIND=DBL), ALLOCATABLE :: node_val(:,:,:)
  INTEGER, ALLOCATABLE   :: equno(:,:)
  INTEGER, ALLOCATABLE   :: indicator(:)
END TYPE ElemData
TYPE RemeshControl
  REAL(KIND=DBL)  :: grow_factor
  REAL(KIND=DBL)  :: gamma_crack
  REAL(KIND=DBL)  :: size_tip_min
  REAL(KIND=DBL)  :: size_tip_max
  REAL(KIND=DBL)  :: dcheck1_max
  REAL(KIND=DBL)  :: dcheck1_factor
  REAL(KIND=DBL)  :: dcheck2_max
  REAL(KIND=DBL)  :: dcheck2_factor
  REAL(KIND=DBL)  :: divide_ring_factor
  INTEGER				:: tri_combine
END TYPE RemeshControl
TYPE GrowthLaw
  INTEGER				:: law !...law = 1 for Paris law, law = 2 for bilinear law
  REAL (KIND=DBL)	:: paris_c
  REAL (KIND=DBL)	:: paris_m
  REAL (KIND=DBL)	:: bilinear_alpha
  REAL (KIND=DBL)	:: bilinear_beta
END TYPE GrowthLaw
TYPE MatProperty
  INTEGER				:: mtype
  REAL(KIND=DBL)  :: property(21)
  !...toughness(1)=K11; toughness(2)=K22; toughness(3)=K33
  REAL(KIND=DBL)  :: toughness(3)
END TYPE MatProperty
TYPE SifData
  INTEGER*4				:: ntipnode
  INTEGER*4, ALLOCATABLE :: sys2user(:)
  REAL(KIND=DBL), ALLOCATABLE :: k1(:)
  REAL(KIND=DBL), ALLOCATABLE :: k2(:)
  REAL(KIND=DBL), ALLOCATABLE :: k3(:)
END TYPE SifData
!REAL(KIND=DBL), PARAMETER :: PI = 3.141592653589793
!----------------------------------------------------------------------------------
CONTAINS
	SUBROUTINE edge_node_propagate(x1,x2,x3,x4,x5,xnew,cracksize,dcheck1,dcheck2,node4i,node8i,node8adj,&
    inewelem,node4new,node8new)
        IMPLICIT NONE
        REAL(KIND=DBL),INTENT(IN) ::x1(3),x2(3),x3(3),x4(3),x5(3),xnew(3),cracksize,dcheck1,dcheck2
        REAL(KIND=DBL),INTENT(OUT)::node4i(3),node8i(3),node8adj(3),node4new(3),node8new(3)
        INTEGER,INTENT(OUT)             ::inewelem
        !...local variables
        REAL(KIND=DBL)::halfl,dist,distt,dis1,dis2,dis3,dis4,dis5,dis6
        REAL(KIND=DBL)::zeta,zeta1,zeta2,dtt,dttt1,dttt2,adv
        REAL(KIND=DBL)::utemp(3)
        INTEGER             ::k,idcase
        !----------------------------------------------------------------------------------
	    !...Initilize inewelem
        inewelem=0
	    !...OBTAIN NEW COORDINATES OF node4i and node8i
	    !...obtain half length of crack tip element
	    halfl=0.5d0*cracksize
        !----------------------------------------------------------------------------------
        utemp = 0.d0
        !----------------------------------------------------------------------------------
	    !...obtain distances needed for calculation
	    dist=distance(xnew,x2)
	    distt=distance(xnew,x3)
	    adv=distance(xnew,x1)
        !----------------------------------------------------------------------------------
        if (adv.gt.cracksize) then
		    !...node8i and node4i are in between old node1i(x1) and new node1i(xnew)
		    !...obtain coordinate of node8i by linear interpolation
		    zeta=halfl/adv
	        do k=1,3
			    node8i(k)=(1.0d0-zeta)*xnew(k)+zeta*x1(k)
		    end do
		    !...obtain coordinate of node4i
		    zeta=2.0d0*halfl/adv
		    do k=1,3
		        node4i(k)=(1.0d0-zeta)*xnew(k)+zeta*x1(k)
                
		    end do
            idcase=3   !...node4i between (x1,xnew)
        else if (adv.ge.halfl) then
		    !...only node8i is in between old node1i(x1) and new node1i(xnew)
		    !...obtain coordinate of new ck_elnode(8,i)
		    zeta=halfl/adv
		    do k=1,3
			    node8i(k)=(1.0d0-zeta)*xnew(k)+zeta*x1(k)
		    end do
		    !...obtain distance between node8i and old node8i
		    dis1=distance(node8i,x2)
		    !...obtain distance between node8i and old node4i
		    dis2=distance(node8i,x3)
		    !...obtain distance between node8i and old node8adj
		    dis3=distance(node8i,x4)
		    if (dis1.ge.halfl) then
			    !...node4i is in between old node1i(x1) and old node8i(x2)
			    call findpoint(x1,x2,node8i,halfl,node4i)
			    idcase=4   !...node4i between (x1,x2)
		    else if (dis2.gt.halfl) then
			    !...node4i is in between old node8i(x2) and old node4i(x3)
                call findpoint(x2,x3,node8i,halfl,node4i) 
			    idcase=1   !...node4i between (x2,x3)
		    else if (dis3.gt.halfl) then
			    !...node4i is in between old node4i(x3) and old node8adj(x4)
			    !...obtain coordinate of new ck_elnode(4,i)
                call findpoint(x3,x4,node8i,halfl,node4i) 
                
			    idcase=2   !...node4i between (x3,x4)
            else
			    !...node4i is outside the range from old node8adj(x4) to old node1i(x1)
		        write(*,*)"CASE IS NOT COVERED YET: code#1"
                stop
            end if
        else if (dist.gt.halfl) then
		    !...node8i is in between old node1i(x1) and old node8i(x2)
		    call findpoint(x1,x2,xnew,halfl,node8i)
		    !...obtain distances needed for calculation
		    dis1=distance(node8i,x2)
		    dis2=distance(node8i,x3)
		    dis3=distance(node8i,x4)    
		    if (dis1.ge.halfl) then
			    !...node4i is in between old node1i(x1) and old node8i(x2)
			    call findpoint(x1,x2,node8i,halfl,node4i)
                
			    idcase=4
		    else if (dis2.ge.halfl) then
			    !...node4i is in between old node8i(x2) and old node4i(x3)
			    call findpoint(x2,x3,node8i,halfl,node4i)
                
			    idcase=1
		    else if (dis3.gt.halfl) then
			    !...node4i is in between old node4i(x3) and old node8adj(x4)
                call findpoint(x3,x4,node8i,halfl,node4i)
                
			    idcase=2
            else
			    !...node4i is outside the range from old node8adj(x4) to old node1i(x1)
			    write(*,*)"CASE IS NOT COVERED YET: code#2"
                stop
		    end if
	    else if (distt.gt.halfl) then
		    !...node8i is in between old node4i(x3) and old node8i(x2)
            call findpoint(x2,x3,xnew,halfl,node8i)
            
		    !...obtain distances need for calculation
		    dis2=distance(node8i,x3)
		    dis3=distance(node8i,x4)    
		    if (dis2.gt.halfl) then
			    !...node4i is in between old node4i(x3) and old node8i(x2)
			    call findpoint(x2,x3,node8i,halfl,node4i)
                
			    idcase=1
		    else if (dis3.gt.halfl) then
			    !...node4i is in between old node8adj(x4) and old node4i(x3)
			    call findpoint(x3,x4,node8i,halfl,node4i)
                
			    idcase=2
            else
			    !...node4i is outside the range
			    write(*,*)"CASE IS NOT COVERED YET: code#3"
                stop
		    end if
	    else
		    !...node8i is outside the range of old node4i(x3) to new node1i(xnew)
		    write(*,*)"CASE IS NOT COVERED YET: code#4"
            stop
	    end if
        !----------------------------------------------------------------------------------
	    !...determine coordinates of node8adj (no new element occurs) or node8adj and node4new and node8new (if new element occurs)
	    !...check location of new ck_elnode(4,i)
	    if (idcase.eq.1) then
        !...node4i between (x2,x3)
		    !...obtain distances needed for calculation
		    dis1=distance(node4i,x3)
		    dis2=distance(x3,x4)
		    dis3=distance(x4,x5)
		    dis4=dis1+dis2+dis3 
            !----------------------------------------------------------------------------------
		    !...check total distance between node4i and x5
		    if (dis4.ge.dcheck1) then
			    !...new element occurs
			    !...compute natural coordinate of node4new
			    zeta=2.0d0*(dis4-dcheck2)/(dis2+dis3)-1.0d0
                if (zeta.gt.1.0d0) then
				    !...node4new between old node4i (x3) and node4i
				    write(*,*)"CASE IS NOT COVERED YET: code#5"
                    stop
			    else
				    !...node4new between old node4i(x3) and node4adj(x5)
                    !...quadratic interpolate for node4new
				    call quad_interpolate(x5,x4,x3,zeta,node4new)
                    
				    !...obtain natural coordinate for node8adj
                    zeta=0.5d0*(zeta-1.0d0)
				    call quad_interpolate(x5,x4,x3,zeta,node8adj)
                    
                    !----------------------------------------------------------------------------------
				    !...obtain coordinates of node8new
                    dis5=distance(x3,node4new)
                    if (dis5.ge.dis1) then
					    !...node8new between old node4adj(x5) and old node4i(x3)
                        zeta=1.0d0-(dis5-dis1)/(dis2+dis3) !correct
					    call quad_interpolate(x5,x4,x3,zeta,node8new)
                        
				    else
					    !...node8new between old node1i(x1) and old node4i(x3)
					    zeta=1.0d0-(dis1-dis5)/(distance(x1,x2)+distance(x2,x3))
					    call quad_interpolate(x1,x2,x3,zeta,node8new)
                        
                    end if
			    end if
			    !...new element occur on this side
                inewelem=1
		    else
			    !...no new element
			    zeta=dis1/(dis2+dis3)
                if (zeta.le.1.0d0) then
				    !...node8adj between old node4adj(x5) and old node4i(x3)
				    call quad_interpolate(x5,x4,x3,zeta,node8adj)
                    
			    else
				    !...node8adj between old node4i(x3) and old node1i(x1)
				    zeta=1.0d0-2.0d0*(0.50d0*dis4-dis2-dis3)/(distance(x1,x2)+distance(x2,x3))
				    call quad_interpolate(x1,x2,x3,zeta,node8adj)
                    
			    end if
			    inewelem=0
		    end if
        !----------------------------------------------------------------------------------
	    else if (idcase.eq.2) then
        !...node4i between (x4,x3)
        	dis2=distance(x3,x4)
		    dis3=distance(x4,x5)
		    dis4=dis3+distance(x4,node4i)
            !----------------------------------------------------------------------------------
            if (dis4.ge.dcheck1) then
			    zeta=2.0d0*(dis4-dcheck2)/(dis2+dis3)-1.0d0
                if (zeta.gt.1.0d0) then
				    !...node4new is not in between (x5,x3)
				    write(*,*)"CASE IS IMPOSSIBLE (1):"
                    stop
			    else
				    !...node4new between (x5,x3)
				    !...OBTAIN COORDINATE OF node4new and node8adj
				    call quad_interpolate(x5,x4,x3,zeta,node4new)
                    zeta=0.5d0*(zeta-1.0d0)
				    call quad_interpolate(x5,x4,x3,zeta,node8adj)
                    
                    !----------------------------------------------------------------------------------
				    !...OBTAIN COORDINATE OF node8new
                    zeta=2.0d0*(dis4-0.5d0*dcheck2)/(dis2+dis3)-1.0d0
				    call quad_interpolate(x5,x4,x3,zeta,node8new)
                    
			    end if
			    !...new element occur on this side
                inewelem=1
		    else
			    !...no new element
			    zeta=-1.0d0+dis4/(dis2+dis3)
                if (zeta.le.1.0d0) then
				    !...node8adj between (x5,x3)
				    call quad_interpolate(x5,x4,x3,zeta,node8adj)
                    
			    else
				    write(*,*)"IMPOSSIBLE CASE (2):"
                    stop
			    end if
		    end if
        !----------------------------------------------------------------------------------
	    else if (idcase.eq.3) then
        !...node4i between (x1,xnew)
		    dis2=distance(x3,x4)
		    dis3=distance(x4,x5)
		    dis4=distance(node4i,x1)+distance(x1,x2)+distance(x2,x3)+dis2+dis3
            !----------------------------------------------------------------------------------
            if (dis4.ge.dcheck1) then
			    !...new element occurs
			    !...compute natural coordinate of node4new
			    zeta=2.0d0*(dis4-dcheck2)/(dis2+dis3)-1.0d0
                if (zeta.gt.1.0d0) then
				    !...node4new outside (x5,x3)
				    dtt=dis4-dcheck2-dis2-dis3
				    zeta1=2.0d0*dtt/(distance(x1,x2)+distance(x2,x3))-1.0d0
                    if (zeta1.lt.-1.0d0) then
				        write(*,*)"CASE IS IMPOSSIBLE (3)"
                        stop
				    elseif (zeta1.gt.1.0d0) then
					    !...node4new between old node1i(x1) and node4i (which is between (x1,xnew)): linear interpolation
					    zeta2=(dtt-distance(x1,x2)-distance(x2,x3))/distance(x1,node4i)
                        do k=1,3
						    node4new(k)=x1(k)*(1.0d0-zeta2)+node4i(k)*zeta2
						    node8new(k)=0.50d0*(node4new(k)+node4i(k))
					    end do
                        
					    if ((dis4-dcheck2)/2.0d0.gt.(dis2+dis3)) then
						    !...node8adj is outside of x5-x3
				            write(*,*)"CASE IS NOT COVERED YET: code#6"
                            stop
					    else
						    !...node8adj is in between x5-x3
						    zeta=(dis4-dcheck2)/(dis2+dis3)-1.0d0
						    call quad_interpolate(x5,x4,x3,zeta,node8adj)
                            
					    end if					  
				    else
					    !...node4new between (x1,x3)
                        call quad_interpolate(x3,x2,x1,zeta1,node4new)
                        
					    !...obtain coordinate of node8new
                        if (zeta1.le.0.0d0) then
						    !...case of node4new between (x2,x3)
						    dttt1=distance(x2,node4new)+distance(x1,x2)
                            dttt2=distance(node4i,x1)
 				        else
						    !...case of node4new between (x1,x2)
						    dttt1=distance(x1,node4new)
                            dttt2=distance(node4i,x1)
					    end if
					    if (dttt1.ge.dttt2) then
						    !...node8new between (x1,x3)
						    zeta=(dttt1-dttt2)/(distance(x1,x2)+distance(x2,x3))-1.0d0
						    call quad_interpolate(x1,x2,x3,zeta,node8new)
                            
					    else
						    !...node8new between old node1i(x1) and node4i (which is between (x1,xnew)): linear interpolation
						    zeta=(dttt1+dttt2)/dttt2/2.0d0
                            do k=1,3
							    node8new(k)=x1(k)*zeta+node4i(k)*(1.0d0-zeta)
				            end do
                            
					    end if
					    !...relocate node8adj
					    if ((dis4-dcheck2)/2.0d0.gt.(dis2+dis3)) then
						    !...node8adj outside of (x5,x3)
				            write(*,*)"CASE IS NOT COVERED YET: code#7"
                            stop
					    else
						    !...node8adj between (x5,x3)
						    zeta=(dis4-dcheck2)/(dis2+dis3)-1.0d0
						    call quad_interpolate(x5,x4,x3,zeta,node8adj)
                            
					    end if
				    end if
			    else
				    !...node4new inside (x5,x3)
				    call quad_interpolate(x5,x4,x3,zeta,node4new)
                    
                    zeta=0.5d0*(zeta-1.0d0)
				    call quad_interpolate(x5,x4,x3,zeta,node8adj)
                    
				    !...OBTAIN COORDINATE OF node8new
				    !...compute distances needed for calculation
                    dis5=distance(x3,node4new)
                    dis6=distance(x1,node4i)+distance(x1,x2)+distance(x2,x3)
                    if (dis5.ge.dis6) then
					    !...node8new between (x5,x3)
					    zeta=1.0d0-(dis5-dis6)/(dis2+dis3)
					    call quad_interpolate(x5,x4,x3,zeta,node8new)
                        
				    else
					    !...node8new outside (x5,x3)
					    zeta=(-distance(x1,node4i)+distance(x3,node4new))/(distance(x1,x2)+distance(x2,x3))
					    if (zeta.ge.-1.0d0) then
						    !...node8new between (x3,x1)
						    call quad_interpolate(x1,x2,x3,zeta,node8new)
                            
					    else
						    !...node8new between old node1i(x1) and node4i: linear interpolation
						    zeta=(distance(x1,node4i)+distance(x1,x2)+distance(x2,x3)+distance(x3,node4new))/distance(x1,node4i)/2.0d0
                            do k=1,3
							    node8new(k)=x1(k)*zeta+node4i(k)*(1.0d0-zeta)
					        end do
                            
					    end if
                    end if
			    end if
			    !...new element occur on this side
                inewelem=1
		    else if (dis4.lt.dcheck1) then
			    !...no new element: relcate node8adj
			    zeta=-1.0d0+dis4/(dis2+dis3)
                if (zeta.le.1.0d0) then
				    !...node8adj between x5-x3
				    call quad_interpolate(x5,x4,x3,zeta,node8adj)
                    
			    else
				    !...node8adj outside x5-x3
				    zeta=1.0d0-2.0d0*(dis4/2.0d0-dis2-dis3)/(distance(x1,x2)+distance(x2,x3))
				    if (zeta.ge.-1.0d0) then
					    !...node8adj between x1-x3
					    call quad_interpolate(x1,x2,x3,zeta,node8adj)
                        
				    else
					    !...node8adj between x1 and node4i: linear interpolation
					    zeta=dis4/2.0d0/distance(x1,node4i)
                        do k=1,3
						    node8adj(k)=x1(k)*zeta+node4i(k)*(1.0d0-zeta)
					    end do
                        
				    end if
			    end if
		    else
			    write(*,*)"CASE IS IMPOSSIBLE (4)"
                stop
		    end if
        !----------------------------------------------------------------------------------
	    else if (idcase.eq.4) then
        !...node4i between (x1,x2)
		    dis2=distance(x3,x4)
		    dis3=distance(x4,x5)
		    dis4=distance(node4i,x2)+distance(x2,x3)+dis2+dis3
            if (dis4.ge.dcheck1) then
			    !...new element occurs
			    zeta=2.0d0*(dis4-dcheck2)/(dis2+dis3)-1.0d0
                if (zeta.gt.1.0d0) then
				    !...node4new between x3 and node4i
				    dtt=dis4-dcheck2-dis2-dis3
				    zeta1=2.0d0*dtt/(distance(x1,x2)+distance(x2,x3))-1.0d0
                    if (zeta1.lt.-1.0d0) then
				        write(*,*)"CASE IS IMPOSSIBLE (5):"
                        stop
				    elseif (zeta1.gt.1.0d0) then
				        write(*,*)"CASE IS IMPOSSIBLE (6):"
                        stop
				    else    
					    !...node4new between (x3,x1)
					    call quad_interpolate(x3,x2,x1,zeta1,node4new) !correct with zeta1 defined above
                        
					    !...obtain coordinate of node8new
                        !...natural coord of node4i
					    zeta2=2.0d0*(distance(x3,x2)+distance(x2,node4i))/(distance(x1,x2)+distance(x2,x3))-1.0d0
					    !...then, natural coord of node8new
					    zeta=0.50d0*(zeta1+zeta2)
					    call quad_interpolate(x3,x2,x1,zeta,node8new)
                        
					    !...relocate node8adj
					    if ((dis4-dcheck2)/2.0d0.gt.(dis2+dis3)) then
						    !...node8adj outside (x5,x3)
				            write(*,*)"CASE IS NOT COVERED YET: code#8"
                            stop
					    else
						    !...node8adj between (x5,x3)
						    zeta=(dis4-dcheck2)/(dis2+dis3)-1.0d0
						    call quad_interpolate(x5,x4,x3,zeta,node8adj)
                            
					    end if
				    end if
			    else
				    !...node4new between (x5,x3)
				    call quad_interpolate(x5,x4,x3,zeta,node4new)
                    
                    !...relocated node8adj
				    zeta=0.5d0*(zeta-1.0d0)
				    call quad_interpolate(x5,x4,x3,zeta,node8adj)
                    
				    !...OBTAIN COORDINATE OF node8new
				    !...compute distances needed for calculation
                    dis5=distance(x3,node4new) !here, if more exact, need to divide into 2 cases: node4new between x5-x4 or x4-x3
                    dis6=distance(x3,x2)+distance(x2,node4i)
                    if (dis5.ge.dis6) then
					    !...node8new between (x5,x3)
					    zeta=1.0d0-(dis5-dis6)/(dis2+dis3)
					    call quad_interpolate(x5,x4,x3,zeta,node8new)
                        
				    else
					    !...node8new between (x3,node4i)
					    zeta=(-distance(x2,node4i)+distance(x3,node4new)+distance(x1,x2))/(distance(x1,x2)+distance(x2,x3))
					    if (zeta.ge.-1.0d0) then
						    !...node8new between (x1,x3)
						    call quad_interpolate(x1,x2,x3,zeta,node8new)
                            
					    else
				            write(*,*)"CASE IS IMPOSSIBLE (7):"
                            stop
					    end if
                    end if
			    end if
			    !...new element occur on this side
                inewelem=1
		    else
			    !...no new element: just relocate node8adj
			    zeta=-1.0d0+dis4/(dis2+dis3)
                if (zeta.le.1.0d0) then
				    !...node8adj between (x5,x3)
				    call quad_interpolate(x5,x4,x3,zeta,node8adj)
                    
			    else
				    !...node8adj outside (x5,x3)
				    zeta=1.0d0-2.0d0*(dis4/2.0d0-dis2-dis3)/(distance(x1,x2)+distance(x2,x3))
				    if (zeta.ge.-1.0d0) then
					    !...node8adj between (x3,x1)
					    call quad_interpolate(x1,x2,x3,zeta,node8adj)
                        
				    else
				        write(*,*)"CASE IS IMPOSSIBLE (8):"
                        stop
				    end if
			    end if
		    end if
	    else
		    write(*,*)"CASE IS NOT COVERED YET: code#9 (not in 4 cases of idcase)"
            stop
	    end if
    END SUBROUTINE edge_node_propagate
    !----------------------------------------------------------------------------------
	FUNCTION distance(xx,yy)
	    IMPLICIT NONE
	    REAL(KIND=DBL)::distance,xx(3),yy(3)
      
		distance=dsqrt((xx(1)-yy(1))**2.0d0+(xx(2)-yy(2))**2.0d0+(xx(3)-yy(3))**2.0d0)
    END FUNCTION distance
    !----------------------------------------------------------------------------------
	SUBROUTINE modtransform(C,eb1,eb2,eb3,Cb)
	!...transform contracted moduli matrix C (6x6) from global coordinates to local coordinate
	!...which has basis (eb1,eb2,eb3)
	!...results are record in Cb(6x6)
		IMPLICIT NONE
		REAL(KIND=DBL),INTENT(IN)	::C(6,6),eb1(3),eb2(3),eb3(3)
		REAL(KIND=DBL),INTENT(OUT)	::Cb(6,6)
		!...local varibles
		REAL(KIND=DBL)		::E(3,3,3,3),Eb(3,3,3,3),a(3,3)
		INTEGER				::i,j,k,l,p,q,m,n

		!...determine directional cosines of local frame
		do i=1,3
			a(i,1)=eb1(i)
		enddo
		do i=1,3
			a(i,2)=eb2(i)
		enddo
		do i=1,3
			a(i,3)=eb3(i)
		enddo
		
		!...calculate corresponding Eijkl with Cij
		E(1,1,1,1)=C(1,1)
		E(1,1,2,2)=C(1,2)
		E(1,1,3,3)=C(1,3)
		E(1,1,2,3)=C(1,4)
		E(1,1,3,1)=C(1,5)
		E(1,1,1,2)=C(1,6)
		E(2,2,2,2)=C(2,2)
		E(2,2,3,3)=C(2,3)
		E(2,2,2,3)=C(2,4)
		E(2,2,3,1)=C(2,5)
		E(2,2,1,2)=C(2,6)
		E(3,3,3,3)=C(3,3)
		E(3,3,2,3)=C(3,4)
		E(3,3,3,1)=C(3,5)
		E(3,3,1,2)=C(3,6)
		E(2,3,2,3)=C(4,4)
		E(2,3,3,1)=C(4,5)
		E(2,3,1,2)=C(4,6)
		E(3,1,3,1)=C(5,5)
		E(3,1,1,2)=C(5,6)
		E(1,2,1,2)=C(6,6)
		!...obtain the rest of matrix Eijkl by using symmetric property
		E(2,2,1,1)=E(1,1,2,2)
		E(3,3,1,1)=E(1,1,3,3)
		E(3,3,2,2)=E(2,2,3,3)
		E(2,3,1,1)=E(1,1,2,3)
		E(2,3,2,2)=E(2,2,2,3)
		E(2,3,3,3)=E(3,3,2,3)
		E(3,1,1,1)=E(1,1,3,1)
		E(3,1,2,2)=E(2,2,3,1)
		E(3,1,3,3)=E(3,3,3,1)
		E(3,1,2,3)=E(2,3,3,1)
		E(1,2,1,1)=E(1,1,1,2)
		E(1,2,2,2)=E(2,2,1,2)
		E(1,2,3,3)=E(3,3,1,2)
		E(1,2,2,3)=E(2,3,1,2)
		E(1,2,3,1)=E(3,1,1,2)

		do i=1,3
			E(i,i,3,2)=E(i,i,2,3)
			E(i,i,1,3)=E(i,i,3,1)
			E(i,i,2,1)=E(i,i,1,2)
		enddo

		E(2,3,3,2)=E(2,3,2,3)
		E(2,3,1,3)=E(2,3,3,1)
		E(2,3,2,1)=E(2,3,1,2)

		E(3,1,3,2)=E(3,1,2,3)
		E(3,1,1,3)=E(3,1,3,1)
		E(3,1,2,1)=E(3,1,1,2)

		E(1,2,3,2)=E(1,2,2,3)
		E(1,2,1,3)=E(1,2,3,1)
		E(1,2,2,1)=E(1,2,1,2)

		E(3,2,1,1)=E(2,3,1,1)
		E(3,2,2,2)=E(2,3,2,2)
		E(3,2,3,3)=E(2,3,3,3)
		E(3,2,2,3)=E(2,3,2,3)
		E(3,2,3,1)=E(2,3,3,1)
		E(3,2,1,2)=E(2,3,1,2)
		E(3,2,3,2)=E(2,3,3,2)
		E(3,2,1,3)=E(2,3,1,3)
		E(3,2,2,1)=E(2,3,2,1)

		E(1,3,1,1)=E(3,1,1,1)
		E(1,3,2,2)=E(3,1,2,2)
		E(1,3,3,3)=E(3,1,3,3)
		E(1,3,2,3)=E(3,1,2,3)
		E(1,3,3,1)=E(3,1,3,1)
		E(1,3,1,2)=E(3,1,1,2)
		E(1,3,3,2)=E(3,1,3,2)
		E(1,3,1,3)=E(3,1,1,3)
		E(1,3,2,1)=E(3,1,2,1)

		E(2,1,1,1)=E(1,2,1,1)
		E(2,1,2,2)=E(1,2,2,2)
		E(2,1,3,3)=E(1,2,3,3)
		E(2,1,2,3)=E(1,2,2,3)
		E(2,1,3,1)=E(1,2,3,1)
		E(2,1,1,2)=E(1,2,1,2)
		E(2,1,3,2)=E(1,2,3,2)
		E(2,1,1,3)=E(1,2,1,3)
		E(2,1,2,1)=E(1,2,2,1)
   
		!...compute the moduli in new coordinate system
		Eb=0.0d0
		do i=1,3
			do j=1,3
				do k=1,3
					do l=1,3
						do p=1,3
							do q=1,3
								do m=1,3
									do n=1,3
										Eb(p,q,m,n)=Eb(p,q,m,n)+a(i,p)*a(j,q)*a(k,m)*a(l,n)*E(i,j,k,l)
									enddo
								enddo
							enddo
						enddo
					enddo
				enddo
			enddo
		enddo

		!...obtain the reduced matrix (6x6) of moduli
		Cb(1,1)=Eb(1,1,1,1)
		Cb(1,2)=Eb(1,1,2,2)
		Cb(1,3)=Eb(1,1,3,3)
		Cb(1,4)=Eb(1,1,2,3)
		Cb(1,5)=Eb(1,1,3,1)
		Cb(1,6)=Eb(1,1,1,2)
		Cb(2,2)=Eb(2,2,2,2)
		Cb(2,3)=Eb(2,2,3,3)
		Cb(2,4)=Eb(2,2,2,3)
		Cb(2,5)=Eb(2,2,3,1)
		Cb(2,6)=Eb(2,2,1,2)
		Cb(3,3)=Eb(3,3,3,3)
		Cb(3,4)=Eb(3,3,2,3)
		Cb(3,5)=Eb(3,3,3,1)
		Cb(3,6)=Eb(3,3,1,2)
		Cb(4,4)=Eb(2,3,2,3)
		Cb(4,5)=Eb(2,3,3,1)
		Cb(4,6)=Eb(2,3,1,2)
		Cb(5,5)=Eb(3,1,3,1)
		Cb(5,6)=Eb(3,1,1,2)
		Cb(6,6)=Eb(1,2,1,2)
		!...symmetry
		do i=2,6
			do j=1,i-1
				Cb(i,j)=Cb(j,i)
			enddo
		enddo
	END SUBROUTINE modtransform
!--------------------------------------------------------------------------------------------	
	SUBROUTINE comtransform(Sc,eb1,eb2,eb3,Scb)
	!...transform contracted compliance matrix Sc (6x6) from global coordinates to local coordinate
	!...which has basis (eb1,eb2,eb3)
	!...results are record in Scb (6x6)
		IMPLICIT NONE
		REAL(KIND=DBL),INTENT(IN)	::Sc(6,6),eb1(3),eb2(3),eb3(3)
		REAL(KIND=DBL),INTENT(OUT)::Scb(6,6)
		!...local varibles
		REAL(KIND=DBL)			::S(3,3,3,3),Sb(3,3,3,3),a(3,3)
		INTEGER							::i,j,k,l,p,q,m,n

		!...determine directional cosines of local frame
		do i=1,3
			a(i,1)=eb1(i)
		enddo
		do i=1,3
			a(i,2)=eb2(i)
		enddo
		do i=1,3
			a(i,3)=eb3(i)
		enddo
		
		!...calculate corresponding Sijkl with Sij
		S(1,1,1,1)=Sc(1,1)
        S(1,1,2,2)=Sc(1,2)
        S(1,1,3,3)=Sc(1,3)
        S(1,1,2,3)=0.5d0*Sc(1,4)
        S(1,1,3,1)=0.5d0*Sc(1,5)
        S(1,1,1,2)=0.5d0*Sc(1,6)
        S(2,2,1,1)=Sc(2,1)
        S(2,2,2,2)=Sc(2,2)
        S(2,2,3,3)=Sc(2,3)
        S(2,2,2,3)=0.5d0*Sc(2,4)
        S(2,2,3,1)=0.5d0*Sc(2,5)
        S(2,2,1,2)=0.5d0*Sc(2,6)
        S(3,3,1,1)=Sc(3,1)
        S(3,3,2,2)=Sc(3,2)
        S(3,3,3,3)=Sc(3,3)
        S(3,3,2,3)=0.5d0*Sc(3,4)
        S(3,3,3,1)=0.5d0*Sc(3,5)
        S(3,3,1,2)=0.5d0*Sc(3,6)
        S(2,3,1,1)=0.5d0*Sc(4,1)
        S(2,3,2,2)=0.5d0*Sc(4,2)
        S(2,3,3,3)=0.5d0*Sc(4,3)
        S(2,3,2,3)=0.25d0*Sc(4,4)
        S(2,3,3,1)=0.25d0*Sc(4,5)
        S(2,1,1,2)=0.25d0*Sc(4,6)
        S(3,1,1,1)=0.5d0*Sc(5,1)
        S(3,1,2,2)=0.5d0*Sc(5,2)
        S(3,1,3,3)=0.5d0*Sc(5,3)
        S(3,1,2,3)=0.25d0*Sc(5,4)
        S(3,1,3,1)=0.25d0*Sc(5,5)
        S(3,1,1,2)=0.25d0*Sc(5,6)
        S(1,2,1,1)=0.5d0*Sc(6,1)
        S(1,2,2,2)=0.5d0*Sc(6,2)
        S(1,2,3,3)=0.5d0*Sc(6,3)
        S(1,2,2,3)=0.25d0*Sc(6,4)
        S(1,2,3,1)=0.25d0*Sc(6,5)
        S(1,2,1,2)=0.25d0*Sc(6,6)
   
        do i=1,3
            S(i,i,3,2)=S(i,i,2,3)
            S(i,i,1,3)=S(i,i,3,1)
            S(i,i,2,1)=S(i,i,1,2)
        enddo

        S(2,3,3,2)=S(2,3,2,3)
        S(2,3,1,3)=S(2,3,3,1)
        S(2,3,2,1)=S(2,3,1,2)

        S(3,1,3,2)=S(3,1,2,3)
        S(3,1,1,3)=S(3,1,3,1)
        S(3,1,2,1)=S(3,1,1,2)

        S(1,2,3,2)=S(1,2,2,3)
        S(1,2,1,3)=S(1,2,3,1)
        S(1,2,2,1)=S(1,2,1,2)
        
        S(3,2,1,1)=S(2,3,1,1)
        S(3,2,2,2)=S(2,3,2,2)
        S(3,2,3,3)=S(2,3,3,3)
        S(3,2,2,3)=S(2,3,2,3)
        S(3,2,3,1)=S(2,3,3,1)
        S(3,2,1,2)=S(2,3,1,2)
        S(3,2,3,2)=S(2,3,3,2)
        S(3,2,1,3)=S(2,3,1,3)
        S(3,2,2,1)=S(2,3,2,1)

        S(1,3,1,1)=S(3,1,1,1)
        S(1,3,2,2)=S(3,1,2,2)
        S(1,3,3,3)=S(3,1,3,3)
        S(1,3,2,3)=S(3,1,2,3)
        S(1,3,3,1)=S(3,1,3,1)
        S(1,3,1,2)=S(3,1,1,2)
        S(1,3,3,2)=S(3,1,3,2)
        S(1,3,1,3)=S(3,1,1,3)
        S(1,3,2,1)=S(3,1,2,1)

        S(2,1,1,1)=S(1,2,1,1)
        S(2,1,2,2)=S(1,2,2,2)
        S(2,1,3,3)=S(1,2,3,3)
        S(2,1,2,3)=S(1,2,2,3)
        S(2,1,3,1)=S(1,2,3,1)
        S(2,1,1,2)=S(1,2,1,2)
        S(2,1,3,2)=S(1,2,3,2)
        S(2,1,1,3)=S(1,2,1,3)
        S(2,1,2,1)=S(1,2,2,1)
   
        !Now do transformation on S
        Sb=0.0d0
        do i=1,3
            do j=1,3
                do k=1,3
                    do l=1,3
                        do p=1,3
                            do q=1,3
                                do m=1,3
                                    do n=1,3
                                        Sb(p,q,m,n)=Sb(p,q,m,n)+a(i,p)*a(j,q)*a(k,m)*a(l,n)*S(i,j,k,l);
                                    enddo
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo

        !obtain the contracted matrix (6x6) of moduli
        Scb(1,1)=Sb(1,1,1,1)
        Scb(1,2)=Sb(1,1,2,2)
        Scb(1,3)=Sb(1,1,3,3)
        Scb(1,4)=2*Sb(1,1,2,3)
        Scb(1,5)=2*Sb(1,1,3,1)
        Scb(1,6)=2*Sb(1,1,1,2)
        Scb(2,2)=Sb(2,2,2,2)
        Scb(2,3)=Sb(2,2,3,3)
        Scb(2,4)=2*Sb(2,2,2,3)
        Scb(2,5)=2*Sb(2,2,3,1)
        Scb(2,6)=2*Sb(2,2,1,2)
        Scb(3,3)=Sb(3,3,3,3)
        Scb(3,4)=2*Sb(3,3,2,3)
        Scb(3,5)=2*Sb(3,3,3,1)
        Scb(3,6)=2*Sb(3,3,1,2)
        Scb(4,4)=4*Sb(2,3,2,3)
        Scb(4,5)=4*Sb(2,3,3,1)
        Scb(4,6)=4*Sb(2,3,1,2)
        Scb(5,5)=4*Sb(3,1,3,1)
        Scb(5,6)=4*Sb(3,1,1,2)
        Scb(6,6)=4*Sb(1,2,1,2)
        do i=2,6
            do j=1,i-1
                Scb(i,j)=Scb(j,i)
            enddo
        enddo
	END SUBROUTINE comtransform
!--------------------------------------------------------------------------------------------
    SUBROUTINE matrix6_inverse(m6,invm6)
	!...subroutine to determine the inverse of matrix m (6x6), result is put on invm
	!...the algorithm is from http://en.wikipedia.org/wiki/Invertible_matrix
		IMPLICIT NONE
		REAL(KIND=DBL),INTENT(IN)	::m6(6,6)
		REAL(KIND=DBL),INTENT(OUT)::invm6(6,6)
		!...local variables
		REAL(KIND=DBL)	::a(3,3),b(3,3),c(3,3),d(3,3)
		REAL(KIND=DBL)	::inva(3,3),e(3,3),f(3,3),g(3,3),h(3,3),invh(3,3),p(3,3),q(3,3),r(3,3)
		INTEGER					::i,j
	
		!...partition
		do i=1,3
			do j=1,3
				a(i,j)=m6(i,j)
				d(i,j)=m6(i+3,j+3)
				b(i,j)=m6(i,j+3)
				c(i,j)=m6(i+3,j)
			enddo
		enddo
		!...calculate sub-defined matrices
		call matrix3_inverse(a,inva)
		call matrix3_multiply(c,inva,e)
		call matrix3_multiply(inva,b,f)
		call matrix3_multiply(e,b,g)
		do i=1,3
			do j=1,3
				h(i,j)=d(i,j)-g(i,j)
			enddo
		enddo
		call matrix3_inverse(h,invh)
		call matrix3_multiply(invh,e,p)
		call matrix3_multiply(f,p,q)
		call matrix3_multiply(f,invh,r)
		!...inverse matrix
		do i=1,3
			do j=1,3
				invm6(i,j)=inva(i,j)+q(i,j)
				invm6(i,j+3)=-r(i,j)
				invm6(i+3,j)=-p(i,j)
				invm6(i+3,j+3)=invh(i,j)
			enddo
		enddo
	END SUBROUTINE matrix6_inverse
!--------------------------------------------------------------------------------------------
    SUBROUTINE matrix5_inverse(m5,invm5)
	!...subroutine to determine the inverse of matrix m (5x5), result is put on invm
	!...the algorithm is from http://en.wikipedia.org/wiki/Invertible_matrix
		IMPLICIT NONE
		REAL(KIND=DBL),INTENT(IN)	::m5(5,5)
		REAL(KIND=DBL),INTENT(OUT)	::invm5(5,5)
		!...local variables
		REAL(KIND=DBL)	::a(2,2),b(2,3),c(3,2),d(3,3)
		REAL(KIND=DBL)	::inva(2,2),e(3,2),f(2,3),g(3,3),h(3,3),invh(3,3),p(3,2),q(2,2),r(2,3)
		INTEGER			::i,j,k
	
		!...partition
        do i = 1,2
          	do j = 1,2
            	a(i,j) = m5(i,j)
            enddo
        enddo
        do i = 1,3
          	do j = 1,3
            	d(i,j) = m5(i+2,j+2)
            enddo
        enddo
        do i = 1,2
          	do j = 1,3
            	b(i,j) = m5(i,j+2)
                c(j,i) = m5(j+2,i)
            enddo
        enddo
		!...calculate sub-defined matrices
		call matrix2_inverse(a,inva)
        !...calculate e = c*inv(a)
        e = 0.0d0
        do i = 1,3
          	do j = 1,2
            	do k = 1,2
                	e(i,j) = e(i,j) + c(i,k)*inva(k,j)
                enddo
            enddo
        enddo
		!...calculate f = inv(a)*b
        f = 0.0d0
        do i = 1,2
          	do j = 1,3
            	do k = 1,2
                	f(i,j) = f(i,j) + inva(i,k)*b(k,j)
                enddo
            enddo
        enddo
		!...calculate g = e*b
        g = 0.0d0
        do i = 1,3
          	do j = 1,3
            	do k = 1,2
                	g(i,j) = g(i,j) + e(i,k)*b(k,j)
                enddo
            enddo
        enddo
		!...calculate h = d - g
		do i=1,3
			do j=1,3
				h(i,j)=d(i,j)-g(i,j)
			enddo
		enddo
		call matrix3_inverse(h,invh)
        !...calculate p = inv(h)*e
        p = 0.0d0
        do i = 1,3
          	do j = 1,2
            	do k = 1,3
                	p(i,j) = p(i,j) + invh(i,k)*e(k,j)
                enddo
            enddo
        enddo
		!...calculate q = f*p
        q = 0.0d0
        do i = 1,2
          	do j = 1,2
            	do k = 1,3
                	q(i,j) = q(i,j) + f(i,k)*p(k,j)
                enddo
            enddo
        enddo
		!...calculate r = f*inv(h)
        r = 0.0d0
        do i = 1,2
          	do j = 1,3
            	do k = 1,3
                	r(i,j) = r(i,j) + f(i,k)*invh(k,j)
                enddo
            enddo
        enddo
		!...inverse matrix
        do i = 1,2
          	do j = 1,2
            	invm5(i,j) = inva(i,j) + q(i,j)
            enddo
        enddo
        do i = 1,2
          	do j = 1,3
            	invm5(i,j+2) = -1.0d0*r(i,j)
                invm5(j+2,i) = -1.0d0*p(j,i)
            enddo
        enddo
		do i=1,3
			do j=1,3
				invm5(i+2,j+2)=invh(i,j)
			enddo
		enddo
	END SUBROUTINE matrix5_inverse
!--------------------------------------------------------------------------------------------
    SUBROUTINE matrix4_inverse(m4,invm4)
	!...subroutine to determine the inverse of matrix m (4x4), result is put on invm
	!...the algorithm is from http://en.wikipedia.org/wiki/Invertible_matrix
		IMPLICIT NONE
		REAL(KIND=DBL),INTENT(IN)	::m4(4,4)
		REAL(KIND=DBL),INTENT(OUT)	::invm4(4,4)
		!...local variables
		REAL(KIND=DBL)	::a(2,2),b(2,2),c(2,2),d(2,2)
		REAL(KIND=DBL)	::inva(2,2),e(2,2),f(2,2),g(2,2),h(2,2),invh(2,2),p(2,2),q(2,2),r(2,2)
		INTEGER			::i,j
	
		!...partition
		do i=1,2
			do j=1,2
				a(i,j)=m4(i,j)
				d(i,j)=m4(i+2,j+2)
				b(i,j)=m4(i,j+2)
				c(i,j)=m4(i+2,j)
			enddo
		enddo
		!...calculate sub-defined matrices
		call matrix2_inverse(a,inva)
		call matrix2_multiply(c,inva,e)
		call matrix2_multiply(inva,b,f)
		call matrix2_multiply(e,b,g)
		do i=1,2
			do j=1,2
				h(i,j)=d(i,j)-g(i,j)
			enddo
		enddo
		call matrix2_inverse(h,invh)
		call matrix2_multiply(invh,e,p)
		call matrix2_multiply(f,p,q)
		call matrix2_multiply(f,invh,r)
		!...inverse matrix
		do i=1,2
			do j=1,2
				invm4(i,j)=inva(i,j)+q(i,j)
				invm4(i,j+2)=-1.0d0*r(i,j)
				invm4(i+2,j)=-1.0d0*p(i,j)
				invm4(i+2,j+2)=invh(i,j)
			enddo
		enddo
	END SUBROUTINE matrix4_inverse
!--------------------------------------------------------------------------------------------
	SUBROUTINE matrix3_inverse(m,invm)
	!...subroutine to calculate inverse of 3x3 matrix
		IMPLICIT NONE
		REAL(KIND=DBL), INTENT(IN)::m(3,3)
		REAL(KIND=DBL), INTENT(OUT)::invm(3,3)
		!...local variable
		REAL(KIND=DBL)::a(2,2),inva(2,2)
		REAL(KIND=DBL)::e(2),f(2),g,h
		
		a(1,1)=m(1,1)
		a(1,2)=m(1,2)
		a(2,1)=m(2,1)
		a(2,2)=m(2,2)
		call matrix2_inverse(a,inva)
		e(1)=m(3,1)*inva(1,1)+m(3,2)*inva(2,1)
		e(2)=m(3,1)*inva(1,2)+m(3,2)*inva(2,2)
		f(1)=inva(1,1)*m(1,3)+inva(1,2)*m(2,3)
		f(2)=inva(2,1)*m(1,3)+inva(2,2)*m(2,3)
		g=e(1)*m(1,3)+e(2)*m(2,3)
		h=m(3,3)-g
		if (abs(h).lt.(1/LARGE_NUM)) then
			write(*,*)'error in subroutine inverse33 (1), h=',h
			stop
		endif
		invm(1,1)=inva(1,1)+f(1)*e(1)/h
		invm(1,2)=inva(1,2)+f(1)*e(2)/h
		invm(2,1)=inva(2,1)+f(2)*e(1)/h
		invm(2,2)=inva(2,2)+f(2)*e(2)/h
		invm(3,3)=1/h
		invm(1,3)=-f(1)/h
		invm(2,3)=-f(2)/h
		invm(3,1)=-e(1)/h
		invm(3,2)=-e(2)/h
	END SUBROUTINE matrix3_inverse
!--------------------------------------------------------------------------------------------
	SUBROUTINE matrix2_inverse(m,invm)
	!...subroutine to calculate inverse of 2x2 matrix
		IMPLICIT NONE
		REAL(KIND=DBL), INTENT(IN)::m(2,2)
		REAL(KIND=DBL), INTENT(OUT)::invm(2,2)
		!...local variable
		REAL(KIND=DBL)::det
		
		det=m(1,1)*m(2,2)-m(1,2)*m(2,1)
		if (abs(det).lt.(1/LARGE_NUM)) then
			write(*,*)'error in subroutine inverse22 (1), det=',det
			stop
		endif
		invm(1,1)=m(2,2)/det
		invm(1,2)=-m(1,2)/det
		invm(2,1)=-m(2,1)/det
		invm(2,2)=m(1,1)/det
	END SUBROUTINE matrix2_inverse
!--------------------------------------------------------------------------------------------
	SUBROUTINE matrix3_multiply(a,b,c)
	!...subroutine to calculate multiplication of 3x3 matrices: c=a*b
		IMPLICIT NONE
		REAL(KIND=DBL), INTENT(IN)::a(3,3),b(3,3)
		REAL(KIND=DBL), INTENT(OUT)::c(3,3)
		!...local variable
		INTEGER		::i,j,k
		
		c=0.0d0
		do i=1,3
			do j=1,3
				do k=1,3
					c(i,j)=c(i,j)+a(i,k)*b(k,j)
				enddo
			enddo
		enddo
	END SUBROUTINE matrix3_multiply
!--------------------------------------------------------------------------------------------
	SUBROUTINE matrix2_multiply(a,b,c)
	!...subroutine to calculate multiplication of 2x2 matrices: c=a*b
		IMPLICIT NONE
		REAL(KIND=DBL), INTENT(IN)::a(2,2),b(2,2)
		REAL(KIND=DBL), INTENT(OUT)::c(2,2)
		!...local variable
		INTEGER		::i,j,k
		
		c=0.0d0
		do i=1,2
			do j=1,2
				do k=1,2
					c(i,j)=c(i,j)+a(i,k)*b(k,j)
				enddo
			enddo
		enddo
	END SUBROUTINE matrix2_multiply
!--------------------------------------------------------------------------------------------
	!...5/18/09: do not need vector zi (local 3-axis)
	SUBROUTINE growth_direction(sif1,sif2,sif3,ti,fa,S,Kij,prop_angle,ratioKKc)
	!...subroutine calculate the angle of propagation prop_angle following method in Hoenig's paper
	!...and the ratio of Ke/toughness (ratioKKc) at that angle (prop_angle)
    !...in detail, prop_angle is the direction where sigma_theta/toughness is maximized
    !...toughness is calculated by projection of tensor value Kij (K=K11n1+K22n2+K33n3)
        IMPLICIT NONE
	    REAL(KIND=DBL),INTENT(IN) ::sif1,sif2,sif3
	    REAL(KIND=DBL),INTENT(IN) ::ti(3),fa(3),S(6,6)
        REAL(KIND=DBL),INTENT(IN) ::Kij(3)
	    REAL(KIND=DBL),INTENT(OUT)::prop_angle,ratioKKc
	    !--------------------------------------------------
	    !...local variables
	    REAL(KIND=DBL)            ::theta
	    REAL(KIND=DBL)            ::normalvector(3),Sp(6,6)
	    REAL(KIND=DBL)            ::toughness
	    REAL(KIND=DBL)            ::coeff(7),sol_re(6),sol_im(6)
	    REAL(KIND=DBL)            ::sigma_xx,sigma_yy,sigma_xy,sigma_theta
	    REAL(KIND=DBL)            ::ratio_stk,ratio_stk_max
	    REAL(KIND=DBL),PARAMETER::upper_bound=85.0,lower_bound=-85.0,step=0.1 !bounds to find crack prop. angle
	    COMPLEX(KIND=DBL)         ::p(3),Q(3),lambda(3),detN,invN(3,3),A(3),temp
	    INTEGER                         ::i,j,k,icheck
        REAL(KIND=DBL), PARAMETER :: PI = 3.141592653589793d0
        !--------------------------------------------------
	    ratio_stk_max=0.0d0
	    theta=lower_bound !inital guess of theta
	    !...loop to find maximum of ratio sigma_theta/toughness
	    do i=1,int((upper_bound-lower_bound)/step)
	        theta=theta+step
	        !...calculate normal vector of predicted crack plane
	        do j=1,3
	            normalvector(j)=-dsin(theta*PI/180.0d0)*ti(j)+dcos(theta*PI/180.0d0)*fa(j)
	        enddo
	        !...calculate toughness in that plane
	        toughness=Kij(1)*normalvector(1)*normalvector(1)+Kij(2)*normalvector(2)*normalvector(2)+Kij(3)*normalvector(3)*normalvector(3)
	        !...calculate reduced compliance (plane strain state)
	        do j=1,6
	            if (j.eq.3) cycle
	            do k=1,6
	                if (k.eq.3) cycle
	                Sp(j,k)=S(j,k)-S(j,3)*S(3,k)/S(3,3)
	            enddo
	        enddo
	        !...calculate coefficients of characteristic equation
	        coeff(7)=Sp(1,1)*Sp(5,5)-Sp(1,5)*Sp(1,5)
		    coeff(6)=-2.0d0*Sp(1,1)*Sp(4,5)-2.0d0*Sp(1,6)*Sp(5,5)+2.0d0*Sp(1,5)*Sp(1,4)+2.0d0*Sp(1,5)*Sp(5,6)
		    coeff(5)=Sp(1,1)*Sp(4,4)+4.0d0*Sp(1,6)*Sp(4,5)+2.0d0*Sp(1,2)*Sp(5,5)+Sp(6,6)*Sp(5,5)-Sp(1,4)*Sp(1,4)-Sp(5,6)*Sp(5,6)-2.0d0*Sp(1,5)*Sp(2,5)-2.0d0*Sp(1,5)*Sp(4,6)-2.0d0*Sp(1,4)*Sp(5,6)
		    coeff(4)=-2.0d0*Sp(6,6)*Sp(4,5)-2.0d0*Sp(2,6)*Sp(5,5)-2.0d0*Sp(1,6)*Sp(4,4)-4.0d0*Sp(1,2)*Sp(4,5)+2.0d0*Sp(1,4)*Sp(2,5)+2.0d0*Sp(1,4)*Sp(4,6)+2.0d0*Sp(5,6)*Sp(2,5)+2.0d0*Sp(5,6)*Sp(4,6)+2.0d0*Sp(1,5)*Sp(2,4)
		    coeff(3)=4.0d0*Sp(2,6)*Sp(4,5)+2.0d0*Sp(1,2)*Sp(4,4)+Sp(6,6)*Sp(4,4)+Sp(2,2)*Sp(5,5)-Sp(2,5)*Sp(2,5)-Sp(4,6)*Sp(4,6)-2.0d0*Sp(1,4)*Sp(2,4)-2.0d0*Sp(5,6)*Sp(2,4)-2.0d0*Sp(2,5)*Sp(4,6)
		    coeff(2)=-2.0d0*Sp(2,6)*Sp(4,4)-2.0d0*Sp(2,2)*Sp(4,5)+2.0d0*Sp(2,5)*Sp(2,4)+2.0d0*Sp(4,6)*Sp(2,4)
		    coeff(1)=Sp(2,2)*Sp(4,4)-Sp(2,4)*Sp(2,4)
		    !...solve characteristic equation
		    call zrhqr(6,coeff,sol_re,sol_im)
		    !...get the solution with positive imaginary part
		    icheck=0
		    do j=1,6
		        if (sol_im(j).gt.0.0d0) then
		            icheck=icheck+1
		            if (icheck.gt.3) then
		                write(*,*)"number of solution which has positive imaginary part must be <=3"
		                stop
		            endif
		            p(icheck)=cmplx(sol_re(j),sol_im(j))
		        endif
		    enddo
		    !...calculate values of Q
		    do j=1,3
!...8/9/09: try to use function sqrt here, will double check later
		        Q(j)=sqrt(dcos(theta*PI/180.0d0)+p(j)*dsin(theta*PI/180.0d0))
		    enddo
		    !...calculate values of lambda
		    do j=1,3
		        lambda(j)=-(Sp(1,5)*p(j)**3.0-(Sp(1,4)+Sp(5,6))*p(j)**2.0d0+(Sp(2,5)+Sp(4,6))*p(j)-Sp(2,4))/(Sp(5,5)*p(j)**2.0d0-2.0d0*Sp(4,5)*p(j)+Sp(4,4))
		    enddo
		    !...calculate determinant of N
		    detN=p(2)*lambda(3)+p(3)*lambda(1)+p(1)*lambda(2)-p(2)*lambda(1)-p(1)*lambda(3)-p(3)*lambda(2)
		    !...calculate inverse of N
		    invN(1,1)=(p(2)*lambda(3)-p(3)*lambda(2))/detN
		    invN(1,2)=(lambda(3)-lambda(2))/detN
		    invN(1,3)=(p(2)-p(3))/detN
		    invN(2,1)=(p(3)*lambda(1)-p(1)*lambda(3))/detN
		    invN(2,2)=(lambda(1)-lambda(3))/detN
		    invN(2,3)=(p(3)-p(1))/detN
		    invN(3,1)=(p(1)*lambda(2)-p(2)*lambda(1))/detN
		    invN(3,2)=(lambda(2)-lambda(1))/detN
		    invN(3,3)=(p(1)-p(2))/detN
		    !...calculate values of A(i)=invN(i,j)*K(j)
		    do j=1,3
		        A(j)=invN(j,1)*sif1+invN(j,2)*sif2+invN(j,3)*sif3
		    enddo
		    !...calculate stresses
		    temp=(0.0d0,0.0d0)
		    do j=1,3
		        temp=temp+(p(j)**2.0d0)*A(j)/Q(j)
		    enddo
		    sigma_xx=real(temp,DBL)
		    temp=(0.0d0,0.0d0)
		    do j=1,3
		        temp=temp+A(j)/Q(j)
		    enddo
		    sigma_yy=real(temp,DBL)
		    temp=(0.0d0,0.0d0)
		    do j=1,3
		        temp=temp-p(j)*A(j)/Q(j)
		    enddo
		    sigma_xy=real(temp,DBL)
		    !...hoop stress in the predicted crack plane
		    sigma_theta=sigma_xx*dsin(theta*PI/180.0d0)*dsin(theta*PI/180.0d0)+ &
						sigma_yy*dcos(theta*PI/180.0d0)*dcos(theta*PI/180.0d0)-sigma_xy*dsin(2.0d0*theta*PI/180.0d0)
		    !...ratio between hoop stress and toughness
		    ratio_stk=sigma_theta/toughness
		    if (ratio_stk_max.lt.ratio_stk) then
		        ratio_stk_max=ratio_stk
		        prop_angle=theta
		    endif
		enddo
		!...get the ratio of Ke/Kec for output
		ratioKKc=ratio_stk_max
	END SUBROUTINE growth_direction
!--------------------------------------------------------------------------------------------
    SUBROUTINE growth_direction_sym(sif1,sif2,ti,fa,S,Kij,prop_angle,ratioKKc)
	!...subroutine calculate the angle of propagation prop_angle following method in Saouma's paper (symmetric plane)
	!...and the ratio of Ke/toughness (ratioKKc) at that angle (prop_angle)
    !...in detail, prop_angle is the direction where sigma_theta/toughness is maximized
    !...toughness is calculated by projection of tensor value Kij (K=K11n1+K22n2+K33n3)
        IMPLICIT NONE
	    REAL(KIND=DBL),INTENT(IN) ::sif1,sif2
	    REAL(KIND=DBL),INTENT(IN) ::ti(3),fa(3),S(6,6)
        REAL(KIND=DBL),INTENT(IN) ::Kij(3)
	    REAL(KIND=DBL),INTENT(OUT)::prop_angle,ratioKKc
	    
	    !...local variables
	    REAL(KIND=DBL)            ::theta
	    REAL(KIND=DBL)            ::normalvector(3),Sp(6,6)
	    REAL(KIND=DBL)            ::toughness
	    REAL(KIND=DBL)            ::coeff(5),sol_re(4),sol_im(4)
	    REAL(KIND=DBL)            ::sigma_xx,sigma_yy,sigma_xy,sigma_theta
	    REAL(KIND=DBL)            ::ratio_stk,ratio_stk_max
	    REAL(KIND=DBL),PARAMETER::upper_bound=85.0,lower_bound=-85.0,step=0.1 !bounds to find crack prop. angle
	    COMPLEX(KIND=DBL)         ::p(2),Q(2)
	    INTEGER                         ::i,j,k,icheck
        REAL(KIND=DBL), PARAMETER :: PI = 3.141592653589793
	    ratio_stk_max=0.0d0
	    theta=lower_bound !inital guess of theta
	    !...loop to find maximum of ratio sigma_theta/toughness
	    do i=1,int((upper_bound-lower_bound)/step)
	        theta=theta+step
	        !...calculate normal vector of predicted crack plane
	        do j=1,3
	            normalvector(j)=-dsin(theta*PI/180.0d0)*ti(j)+dcos(theta*PI/180.0d0)*fa(j)
	        enddo
	        !...calculate toughness in that plane
	        toughness=Kij(1)*normalvector(1)*normalvector(1)+Kij(2)*normalvector(2)*normalvector(2)+Kij(3)*normalvector(3)*normalvector(3)
	        !...calculate reduced compliance (plane strain state)
	        do j=1,6
	            if (j.eq.3) cycle
	            do k=1,6
	                if (k.eq.3) cycle
	                Sp(j,k)=S(j,k)-S(j,3)*S(3,k)/S(3,3)
	            enddo
	        enddo
	        !...calculate coefficients of characteristic equation
		    coeff(5)=Sp(1,1)
		    coeff(4)=-2.0d0*Sp(1,6)
		    coeff(3)=2.0d0*Sp(1,2)+Sp(6,6)
		    coeff(2)=-2.0d0*Sp(2,6)
		    coeff(1)=Sp(2,2)
		    !...solve characteristic equation
		    call zrhqr(4,coeff,sol_re,sol_im)
		    !...get the solution with positive imaginary part
		    icheck=0
		    do j=1,4
		        if (sol_im(j).gt.0.0d0) then
		            icheck=icheck+1
		            if (icheck.gt.2) then
		                write(*,*)"number of solution which has positive imaginary part must be <=2"
		                stop
		            endif
		            p(icheck)=cmplx(sol_re(j),sol_im(j))
		        endif
		    enddo
		    !...calculate values of Q
		    do j=1,2
		        Q(j)=sqrt(dcos(theta*PI/180.0d0)+p(j)*dsin(theta*PI/180.0d0))
		    enddo

		    !...calculate stresses
		    sigma_xx=sif1*real(p(1)*p(2)/(p(1)-p(2))*(p(2)/Q(2)-p(1)/Q(1)),DBL)+sif2*real(1.0d0/(p(1)-p(2))*(p(2)**2.0d0/Q(2)-p(1)**2.0d0/Q(1)),DBL)
            sigma_yy=sif1*real(1.0d0/(p(1)-p(2))*(p(1)/Q(2)-p(2)/Q(1)),DBL)+sif2*real(1.0d0/(p(1)-p(2))*(1.0d0/Q(2)-1.0d0/Q(1)),DBL)
            sigma_xy=sif1*real(p(1)*p(2)/(p(1)-p(2))*(1.0d0/Q(1)-1.0d0/Q(2)),DBL)+sif2*real(1.0d0/(p(1)-p(2))*(p(1)/Q(1)-p(2)/Q(2)),DBL)

		    !...hoop stress in the predicted crack plane
		    sigma_theta=sigma_xx*dsin(theta*PI/180.0d0)*dsin(theta*PI/180.0d0)+ &
						sigma_yy*dcos(theta*PI/180.0d0)*dcos(theta*PI/180.0d0)-sigma_xy*dsin(2.0d0*theta*PI/180.0d0)
		    !...ratio between hoop stress and toughness
		    ratio_stk=sigma_theta/toughness
		    if (ratio_stk_max.lt.ratio_stk) then
		        ratio_stk_max=ratio_stk
		        prop_angle=theta
		    endif
		enddo
		!...get the ratio of Ke/Kec for output
		ratioKKc=ratio_stk_max
	END SUBROUTINE growth_direction_sym
!--------------------------------------------------------------------------------------------
    SUBROUTINE findpoint(x1,x2,x3,dist,x4)
    !.....obtain coordinates of a point x4 located on a straight line connecting 
    !.....points x1 and x2 with a specified distance dist from point x3
        IMPLICIT NONE
	    REAL(KIND=DBL),INTENT(IN) ::x1(3),x2(3),x3(3),dist
	    REAL(KIND=DBL),INTENT(OUT)::x4(3)
	    !...local variables
	    REAL(KIND=DBL)            ::d1,a2,b,c2,ze1,ze2,ze
	    INTEGER                         ::k

	    a2=(x1(1)-x3(1))**2.0d0+(x1(2)-x3(2))**2.0d0+(x1(3)-x3(3))**2.0d0
	    b =(x1(1)-x3(1))*(x2(1)-x1(1))+(x1(2)-x3(2))*(x2(2)-x1(2))+(x1(3)-x3(3))*(x2(3)-x1(3))
	    c2=(x2(1)-x1(1))**2.0d0+(x2(2)-x1(2))**2.0d0+(x2(3)-x1(3))**2.0d0
	    d1=(b/c2)**2.0d0-(a2-dist*dist)/c2
	    if   (d1.ge.0.0d0) then
	        ze1=-b/c2+dsqrt(d1)
	        ze2=-b/c2-dsqrt(d1)
		    if (ze1.ge.0.0d0.and.ze1.le.1.0d0) then
		        if (ze2.ge.0.0d0.and.ze2.le.1.0d0) then
	                write(*,*)"CASE IS NOT COVERED IN SUBROUTINE: findpoint(1)"
			        write(*,*)(x1(k),k=1,3)
			        write(*,*)(x2(k),k=1,3)
			        write(*,*)(x3(k),k=1,3)
			        write(*,*)ze1,ze2,dist
		            stop
		        else
		            ze=ze1
		        end if
		    else
		        if (ze2.ge.0.0d0.and.ze2.le.1.0d0) then
		            ze=ze2
		        else
	                write(*,*)"CASE IS NOT COVERED IN SUBROUTINE: findpoint(2)"
			        write(*,*)(x1(k),k=1,3)
			        write(*,*)(x2(k),k=1,3)
			        write(*,*)(x3(k),k=1,3)
			        write(*,*)ze1,ze2,dist
		            stop
		        end if		  
		    end if
		    do k=1,3
		        x4(k)=(1.0d0-ze)*x1(k)+ze*x2(k)
		    end do
	    else
	        write(*,*)"CASE CANNOT OCCUR IN SUBROUTINE: findpoint(3)"
		    stop
	    end if
	END SUBROUTINE findpoint
!--------------------------------------------------------------------------------------------
    SUBROUTINE quad_interpolate(p1,p2,p3,xi,p)
	!...this subroutine calculate position of a node base on positions of 3 nodes
	!...p1, p2 and p3 by using quadratic interpolation. Result is put onto p
	    IMPLICIT NONE
	    REAL(KIND=DBL),INTENT(IN) ::p1(3),p2(3),p3(3)
	    REAL(KIND=DBL),INTENT(IN) ::xi
	    REAL(KIND=DBL),INTENT(OUT)::p(3)
	    !...local variables
	    REAL(KIND=DBL)::phi1,phi2,phi3
	    INTEGER         ::k
	
	    !...evaluate shape functions
		    phi1=0.5d0*(xi*xi-xi)
		    phi2=1.0d0-xi*xi
		    phi3=0.5d0*(xi*xi+xi)

		    do k=1,3
			    p(k)=p1(k)*phi1+p2(k)*phi2+p3(k)*phi3
		    end do
	END SUBROUTINE quad_interpolate
!--------------------------------------------------------------------------------------------
END MODULE Remesh
