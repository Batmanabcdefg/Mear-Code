SUBROUTINE reg_shape(elid,xi,psi,dpsi)
  ! compute the shape function and it's derivatives
  !
  ! Usage: call (elid, xi, psi, dpsi)
  !
  !        xi     -- the coordinate(IN), in master element coordinate
  !        elid   -- the ID of element
  !        psi    -- the shape function at xi(OUT)
  !        dpsi   -- the derivative od shape function at xi(OUT)
  !
  USE DefinitionConstant
  IMPLICIT NONE  

  REAL(KIND=DBL), INTENT(in)    :: xi
  INTEGER, INTENT(in)           :: elid
  REAL(KIND=DBL), INTENT(out)   :: psi(3), dpsi(3)

  SELECT CASE (elid)
    !...4/28/09: add ILIN
    CASE (CLIN,BLIN,ILIN)
    
      ! 2-node linear element
      psi(1) = 0.5d0*(1.0d0 - xi)
      psi(2) = 0.5d0*(1.0d0 + xi)

      ! derivative w.r.t natural coordinate
      dpsi(1) = -0.5d0
      dpsi(2) = 0.5d0
	!...4/28/09: add IQUAD
    CASE (CQUAD,CTIP1,CTIP2,BQUAD,IQUAD)
      ! 3-node quadratic element 

      psi(1) = 0.5d0*xi*(xi - 1.0d0)
      psi(2) = 0.5d0*xi*(xi + 1.0d0)
      psi(3) = (1.0d0 - xi)*(1.0d0 + xi)

      ! compute partial derivatives w.r.t. natural coordinates 
      dpsi(1) = xi - 0.5d0
      dpsi(2) = xi + 0.5d0
      dpsi(3) = -2.0d0*xi

    CASE default
      PRINT *, elid, ' reg_shape: type element is not included in element library'
      STOP
  END SELECT
END SUBROUTINE reg_shape
!--------------------------------------------------------------------------------
SUBROUTINE tip_shape(elid, xi, psi)
  ! compute the shape function of tip element
  !
  ! Usage: call tip_shape(elid, xi, psi)
  !
  !     xi     -- the coordinate(IN), in master element coordinate
  !     elid   -- the ID of element
  !     psi    -- the shape function at x(OUT)
  !
  USE DefinitionConstant
  IMPLICIT NONE  

  REAL(KIND=DBL),INTENT(in)	:: xi
  INTEGER, INTENT(in)		:: elid
  REAL(KIND=DBL),INTENT(out):: psi(3)

  SELECT CASE (elid)
    CASE (CTIP1)  ! 3-node tip element with local node 1 is tip node
      	psi(1) = dsqrt(1.0d0+xi)*xi*(xi-1.0d0)
      	psi(2) = 0.5d0*dsqrt(0.5d0*(1.0d0+xi))*xi*(xi+1.0d0)
      	psi(3) = dsqrt(1.0d0+xi)*(1.0d0-xi)*(1.0d0+xi)
    CASE (CTIP2)  ! 3-node tip element with local node 2 is tip node
      	psi(1) = 0.5d0*dsqrt(0.5d0*(1.0d0-xi))*xi*(xi-1.0d0)
		psi(2) = dsqrt(1.0d0-xi)*xi*(xi+1.0d0)
      	psi(3) = dsqrt(1.0d0-xi)*(1.0d0-xi)*(1.0d0+xi)
    CASE default
      PRINT *, elid, ' type element is not included in tip-element library'
      STOP
  END SELECT
END SUBROUTINE tip_shape
!--------------------------------------------------------------------------------
SUBROUTINE tip_dshape(elid, xi, dpsi)
	USE DefinitionConstant
	IMPLICIT NONE

	REAL(KIND=DBL), INTENT(in) :: xi
	INTEGER, INTENT(in)        :: elid
	REAL(KIND=DBL), INTENT(out):: dpsi(3)

	!...internal variable
	REAL(KIND=DBL)	:: temp,temp2

	temp = dsqrt(1.0d0 + xi)
    temp2 = dsqrt(1.0d0 - xi)

    !...defensive programming
    if ((temp.lt.1.0E-15).or.(temp2.lt.1.0E-15)) then
      	print*,'error: Gauss point is too close to the -1 or 1 and makes derivative of tip shape singular'
        print*,'Gauss point=',xi
        stop
    endif

	SELECT CASE (elid)
    	CASE (CTIP1)  ! 3-node tip element with local node 1 is tip node
      		dpsi(1) = 0.5d0*xi*(xi-1.0d0)/temp + temp*(2.0d0*xi-1.0d0)
      		dpsi(2) = 0.5d0/dsqrt(2.0d0)*(0.5d0*xi*(xi+1.0d0)/temp + temp*(2.0d0*xi+1.0d0))
      		dpsi(3) = 0.5d0*(1.0d0+xi)*(1.0d0-xi)/temp - 2.0d0*temp*xi
    	CASE (CTIP2)  ! 3-node tip element with local node 2 is tip node
			dpsi(2) = -0.5d0*xi*(xi+1.0d0)/temp2 + temp2*(2.0d0*xi+1.0d0)
      		dpsi(1) = 0.5d0/dsqrt(2.0d0)*(0.5d0*xi*(1.0d0-xi)/temp2 + temp2*(2.0d0*xi-1.0d0))
      		dpsi(3) = -0.5d0*(1.0d0+xi)*(1.0d0-xi)/temp2 - 2.0d0*temp2*xi
    CASE default
      PRINT *, elid, ' type element is not included in tip-element library'
      STOP
	END SELECT

END SUBROUTINE tip_dshape