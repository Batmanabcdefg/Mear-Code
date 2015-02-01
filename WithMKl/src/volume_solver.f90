SUBROUTINE Volume_Solver(n,SIF_trac,SIF_pres,vol_trac,vol_pressure,vol_ratio,pressure_scale_temp)

  USE DefinitionConstant
  USE GlobalData
  IMPLICIT NONE
    
    INTEGER , INTENT(IN) :: n
    REAL(KIND=DBL),INTENT(IN), DIMENSION(NDOF,total_tip_node) :: SIF_trac  !(NDOF,total_tip_node)
    REAL(KIND=DBL),INTENT(IN), DIMENSION(total_cracks,NDOF,total_tip_node) :: SIF_pres  !(total_cracks,NDOF,total_tip_node)
    REAL(KIND=DBL),INTENT(IN), DIMENSION(total_cracks) :: vol_trac    !(total_cracks)
    REAL(KIND=DBL),INTENT(IN), DIMENSION(total_cracks,total_cracks) :: vol_pressure  !(total_cracks,total_cracks)
    REAL(KIND=DBL),INTENT(IN), DIMENSION(total_cracks) :: vol_ratio   !(total_cracks)
    REAL(KIND=DBL),INTENT(OUT), DIMENSION(total_cracks,total_tip_node):: pressure_scale_temp   !(total_cracks,total_tip_node)

    ! local variables
    INTEGER ::  i , j , k, nminus1
    REAL(KIND=DBL), DIMENSION(n-1,n-1)  :: A_L
    REAL(KIND=DBL), DIMENSION(n-1,n-1)  :: A_M
    REAL(KIND=DBL), DIMENSION(n-1) :: L
    REAL(KIND=DBL), DIMENSION(n-1) :: M
    REAL(KIND=DBL), DIMENSION(total_tip_node) :: pressure_scale_n_temp
    REAL(KIND=DBL), DIMENSION(n-1) :: pressure_scale_L
    REAL(KIND=DBL), DIMENSION(n-1) :: pressure_scale_M
    INTEGER, DIMENSION(n-1) :: Indx_L
    INTEGER, DIMENSION(n-1) :: Indx_M
    REAL(KIND=DBL) :: dot_L,dot_M


    nminus1=n-1

 !......set up reduced matrix A and two right-hand sides L and M
    do i = 1,nminus1

       do j = 1,nminus1

          A_L(i,j) = vol_pressure(i,j)-vol_ratio(i)*vol_pressure(n,j)

          A_M(i,j) = A_L(i,j)

       enddo

       L(i) = vol_ratio(i)*vol_trac(n)-vol_trac(i)

       M(i) = vol_ratio(i)*vol_pressure(n,n)-vol_pressure(i,n)

    enddo

    print*
    print*,'The reduced volume matrix A is: '
    print*, A_L
    print*,'The reduced RHS due to remote tractions is: '
    print*, L
    print*,'The reduced RHS due to the pressures are: '
    print*, M
    print*

!....pass A and L and A and R to Gauss Solver

    call Gauss_Solver(A_L,n-1,L,pressure_scale_L,Indx_L)
    
    call Gauss_Solver(A_M,n-1,M,pressure_scale_M,Indx_M)

    print*,'The pressure scaling solution due to remote tractions is: '
    print*,pressure_scale_L
    print*,'The pressure scaling solution due to pressure is: '
    print*,pressure_scale_M




    do k = 1,total_tip_node

       dot_L = 0.0d0
       dot_M = 0.0d0

       do j = 1,n-1
          
          dot_L = dot_L + SIF_pres(j,1,k)*pressure_scale_L(j)
          dot_M = dot_M + SIF_pres(j,1,k)*pressure_scale_M(j)

       enddo
              
          pressure_scale_n_temp(k)=(globalK(1)-SIF_trac(1,k)-dot_L)/(dot_M+SIF_pres(n,1,k))
              

    enddo



    do k = 1,total_tip_node
       do j = 1,n-1
              
          pressure_scale_temp(j,k)= pressure_scale_L(j)+pressure_scale_M(j)*pressure_scale_n_temp(k)
              
       enddo

       pressure_scale_temp(n,k)=pressure_scale_n_temp(k)

    enddo

    print*
    print*,'The temporary pressure scaling values for the last crack are: ', pressure_scale_n_temp
    print*


CONTAINS

SUBROUTINE Gauss_Solver(A,n,b,x,Indx)
!
! Subroutine to solve the equation A(n,n)*X(n) = B(n) with the
! partial-pivoting Gaussian elimination scheme.
! Copyright (c) Tao Pang 2001.
!

  INTEGER, INTENT(IN) :: n
  REAL(KIND=DBL), INTENT(INOUT), DIMENSION(n,n) :: A
  REAL(KIND=DBL), INTENT(INOUT), DIMENSION(n) :: b
  REAL(KIND=DBL), INTENT(OUT), DIMENSION(n):: x
  INTEGER, INTENT(OUT), DIMENSION(n) :: Indx


!local variables
  INTEGER :: i,j

!
  CALL Gauss_Elim (A,n,Indx)
!
  do i = 1, n-1
    do j = i+1, n

      b(Indx(j)) = b(Indx(j))-A(Indx(j),i)*b(Indx(i))

    enddo
  enddo
!
  x(n) = b(Indx(n))/A(Indx(n),n)

  do i = n-1, 1, -1
    x(i) = b(Indx(i))

    do j = i+1, n
      x(i) = x(i)-A(Indx(i),j)*x(j)
    enddo

    x(i) =  x(i)/A(Indx(i),i)

  enddo
!
END SUBROUTINE Gauss_Solver

!
SUBROUTINE Gauss_Elim (A,n,Indx)
!
! Subroutine to perform the partial-pivoting Gaussian elimination.
! A(n,n) is the original matrix in the input and transformed matrix
! plus the pivoting element ratios below the diagonal in the output.
! Indx(n) records the pivoting order.  Copyright (c) Tao Pang 2001.
!

  INTEGER, INTENT (IN) :: n
  REAL(KIND=DBL), INTENT (INOUT), DIMENSION(n,n) :: A
  INTEGER, INTENT (OUT), DIMENSION(n) :: Indx

! local variables
  INTEGER :: i,j,k,i_temp
  REAL(KIND=DBL) :: C1,PI,PI1,PJ
  REAL(KIND=DBL), DIMENSION(n) :: C
!
! Initialize the index
!
  do i = 1, n
    Indx(i) = i
  enddo
!
! Find the rescaling factors, one from each row
!
  do i = 1, n
    C1= 0.0d0
    do j = 1, n

      C1 = DMAX1(C1,ABS(A(i,j)))

   enddo

    C(i) = C1

 enddo
!
! Search the pivoting (largest) element from each column
!
  do j = 1, n-1

    PI1 = 0.0d0

    do i = j, n

      PI = ABS(A(Indx(i),j))/C(Indx(i))

      if (PI.GT.PI1) then

        PI1 = PI
        k   = i

     endif
  enddo
!
! Interchange the rows via Indx(n) to record pivoting order
!
    i_temp    = Indx(j)
    Indx(j) = Indx(k)
    Indx(k) = i_temp
    do i = j+1, n

      PJ  = A(Indx(i),j)/A(Indx(j),j)
!
! Record pivoting ratios below the diagonal
!
      A(Indx(i),j) = PJ
!
! Modify other elements accordingly
!
      do k = j+1, n
        A(Indx(i),k) = A(Indx(i),k)-PJ*A(Indx(j),k)
     enddo
  enddo
enddo
!
END SUBROUTINE Gauss_Elim

END SUBROUTINE VOLUME_SOLVER
