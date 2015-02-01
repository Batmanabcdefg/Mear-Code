Subroutine ek_c(isrc,coor_xyz,ni,xic,wic,ekc)

!...this is for evaluating the C_p(X,Y) in equation (3.19) for all the equations
!...Subroutine to calculate ekc for all kinds of (outer) elements
!...use include in the module ElementMatrix, thus don't need GlobaData/DefinitionConstant

    IMPLICIT NONE
    
    INTEGER, INTENT(IN)		:: isrc,ni
    REAL(KIND=DBL), INTENT(IN)	:: coor_xyz(:,:),xic(:),wic(:)
    REAL(KIND=DBL), INTENT(OUT)   :: ekc(:,:)
	!...local variables
    INTEGER					:: i,j,k,ii,jj,kk,nodey
    REAL(KIND=DBL)          :: y,psiy_i(3),psiy_j(3),dpsiy(3),&
                                dyds(2),jacoby,wiy,C
    LOGICAL					:: CTIP_ELEMENT

    
	!...Set value of constant C depending on type of (outer) element
    if (ELTYPE(elemid(isrc)).eq.CTIP .or. ELTYPE(elemid(isrc)).eq.CREGULAR) then
      	!...outer element is on Sc, then C = 1
        C = 4.0d0*ht/3.0d0			
		!...changed here to account for the pseudo 3D profile
		!...the change is made only for source element being crack tip or regular crack element
    else
      	!...outer element is on St or Su (or interface element, need to verify!), then c = 1/2
        C = 0.5d0
		!...check out equation (3.4) in Han's thesis
    endif
    !...determine if (outer) element is a tip element so that stretching transformation will be applied
    if (ELTYPE(elemid(isrc)).eq.CTIP) then
      	CTIP_ELEMENT = .true.
    else
      	CTIP_ELEMENT = .false.
    endif
    !...get number of nodes of (outer) element
    nodey = NODE(elemid(isrc))
    !...initialize
    ekc = 0.0d0

    !...loop over integration points
    do k = 1,ni
      	!...stretching transformation if CTIP element
		!...check out section 3.3.2 in Han's thesis where a special crack tip element is explained
        if (CTIP_ELEMENT) then
          	!...(outer) element is tip, need more detail about type of tip element
          	if (elemid(isrc).eq.CTIP1) then
            	y = 0.5d0*(xic(k)+1.0d0)*(xic(k)+1.0d0) - 1.0d0
                wiy = wic(k)*(xic(k)+1.0d0)   !...weight includes jacobian of stretching transform
            elseif (elemid(isrc).eq.CTIP2) then
            	y = 1.0d0 - 0.5d0*(1.0d0-xic(k))*(1.0d0-xic(k))
                wiy = wic(k)*(1.0d0-xic(k))   !...weight includes jacobian of stretching transform
            else
              	print*,'wrong type of tip element: ',elem_sys2user(isrc)
                stop
            endif
        else
          	!...(outer) element is not a tip element
			!...we need not use special elements for not crack tip
            y = xic(k)
            wiy = wic(k)
        endif
        !...compute regular shape functions at integration points
        call reg_shape(elemid(isrc),y,psiy_j,dpsiy)
        !...compute derivative of position vector
        dyds = 0.0d0
        do i = 1,nodey
          	do kk = 1,2
            	dyds(kk) = dyds(kk) + coor_xyz(kk,i)*dpsiy(i)
            enddo
        enddo
        !...jacobian of (outer) element
        jacoby = dsqrt(dyds(1)*dyds(1)+dyds(2)*dyds(2))
        !...need tip shape functions if (outer) element is tip element
        if (CTIP_ELEMENT) then
          	call tip_shape(elemid(isrc),y,psiy_i)
        else
          	psiy_i = psiy_j
        endif
        !...compute ekc
        do i = 1,nodey
          	!ii = 3*i - 2
            ii = NDOF*i - NDOF + 1
            do j = 1,nodey
              	!jj = 3*j - 2
                jj = NDOF*j - NDOF + 1
                ekc(ii,jj) = ekc(ii,jj) + C*jacoby*psiy_i(i)*psiy_j(j)*wiy
            enddo
        enddo
    enddo   !...of loop over integration points
    !...continue to fill up ekc
    do i = 1,nodey
      	!ii = 3*i - 2
        ii = NDOF*i - NDOF + 1
        do j = 1,nodey
          	!jj = 3*j - 2
            jj = NDOF*j - NDOF + 1
            do k = 1,NDOF-1
            	ekc(ii+k,jj+k) = ekc(ii,jj)
            enddo
        enddo
    enddo
End Subroutine ek_c
