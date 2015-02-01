MODULE LinearSolver
  USE DefinitionConstant  ! SINGLE_KIND is definied in it
  IMPLICIT NONE
  PUBLIC :: dPCGSolver     ! PCG iterative solver(double)
!  REAL(KIND=DBL), PARAMETER :: tol = 1.0e-14

CONTAINS
      
  function dPCGSolver(N, A, b)
    !
    ! interface for the PCG linear solver(double precision)
    !
    ! Argument
    ! 
    ! N       (in) integer, 
    !         On entry, the dimension of the matrix
    ! 
    ! A       (in), DBL
    !         On entry, the coefficient matrix
    ! 
    ! b       (inout), DBL
    !         On entry, the right hand side vector
    !         On exit, the solution of the linear system
    !

    logical dPCGSolver
    integer, intent(in) :: N
    REAL(KIND=DBL), DIMENSION(:), INTENT(inout)    :: A
    REAL(KIND=DBL), DIMENSION(:), INTENT(inout) :: b

    ! local variables
    REAL(KIND=DBL), DIMENSION(N) :: x
    !real*8::b1(1,N)
    real(KIND=DBL) :: tol
    integer :: MAX_ITER, i
    ! tolerence
    tol = 1.0d-6 !this tolerance is from original Xiao's code
    
    ! maximum iteration number
    !...modified by Han to increase max number of iteration (for problems with small size)
    !MAX_ITER = 1.5 * N 
    MAX_ITER = 100 * N 
    !b1=0.0d0
    ! initial guess
    do i = 1, N
      x(i) = b(i) / A(i*(i-1)/2 + i)
    !  b1(1,i)=b(i)
    end do
    
    !call dtptrs('U','N','N',N,1,A,b1,1,info)
    ! calling CG
    dPCGSolver = dPCG(N, A, b, x, MAX_ITER, tol);
    !dPCGSolver=.true.
    print *
    print *, "=== Precoditioning(Jacobi) Conjugate Gradient Iterative Solver ==="
    print *, "Iteration number =", MAX_ITER
    print *, "Residual = ", tol
    !do i = 1,N
    !    b(i)=b1(1,i)
    !end do
    b=x;
  end function dPCGSolver

  FUNCTION dPCG(N, A, b, x, ITER, Resid)
    !
    ! (double precision)
    ! PCG solves the linear system Ax = b using 
    ! Preconditioner Conjugate Gradient iterative method
    ! (Algorithm see "Iterative Methods for Linear and Nonlinear Equations"
    ! By C.T.Kelley, SIAM
    !
    ! Convergence test : |b-A*x|/|b| < TOL
    !
    ! Argument:    
    !
    ! N        (inpupt) Integer
    !          On entry, the dimension of the matrix.
    !          On exit, Unchanged.
    !
    ! A        (input) DBL
    !          On entry, the coefficient matrix
    !          On exit, unchanged.
    !
    ! b        (input) DBL, dimension N
    !          On entry, right hand side vector B
    !          On exit, Unchanged.
    !
    ! x        (input/output) DBL, dimension N
    !          On entry, the initial guess. This is commonly set to 
    !          the zero vector.
    !          On exit, the iterated approximate solution
    !
    ! ITER     (input/output) Integer
    !          On entry, the maximum iterations to be performed.
    !          On exit, the actual number of iteration performed.
    !
    ! Resid    (input/output) DBL
    !          On extry, the allowable convergence measure for |b-A*x|/|b|.
    !          On exit, the final value of the measure.
    ! 
    ! 
    !=======================================================================

    LOGICAL dPCG
    INTEGER, INTENT(in)    :: N
    INTEGER, INTENT(inout) :: ITER
    REAL(KIND=DBL), DIMENSION(:), INTENT(in)    :: A
    REAL(KIND=DBL), DIMENSION(:), INTENT(in)    :: b
    REAL(KIND=DBL), DIMENSION(:), INTENT(inout) :: x
    REAL(KIND=DBL), INTENT(inout)               :: Resid

    ! local variables
    REAL(KIND=DBL), DIMENSION(N) :: r, p, w, z
    REAL(KIND=DBL) :: alpha,beta,rho,tau1,tau2
    REAL(KIND=DBL) :: TOL, bnorm, temp1
    INTEGER :: i,k
    real*8::ddot
    TOL = Resid*Resid
    ! compute the norm of b
    bnorm = 0.0d0
    !do i = 1, N
    !  bnorm = bnorm + b(i)*b(i)
    !end do
    bnorm=ddot(N,b,1,b,1) 
    ! the first step : compute r = b - Ax
    ! compute Ax, at exit : r = Ax
    !call D_mult_Ax(N, A, x, r)
    call dspmv('U',N,1.0d0,A,x,1,0.0d0,r,1)
    ! r = b - Ax
    do i = 1, N
      r(i) = b(i) - r(i)
    end do
    
    ! rho1 = |r|^2
    rho = 0.0d0
    !do i = 1, N
    !  rho = rho + r(i)*r(i)
    !end do
    rho=ddot(N,r,1,r,1)
	!...added by Han to initialize tau2
    tau2 = 0.d0

    ! begin the iteration
    do k = 1, ITER
      
      ! check the convergence
      if (rho < TOL * bnorm) then
        ! convergent
        ITER = k - 1
        Resid = dsqrt(rho)
        dPCG = .true.
        return
      end if
      
      ! call preconditioner
      call D_preconditioner(N, A, r, z)

      ! assign tau2 to tau1
      tau1 = tau2

      ! compute tau
      tau2 = 0.0d0
      !do i = 1, N
      !  tau2 = tau2 + r(i) * z(i)
      !end do
      tau2=ddot(N,r,1,z,1)
      if (k == 1) then
        do i = 1, N
          beta = 0.0d0
          p(i)  = z(i)
        end do
      else
        beta = tau2 / tau1
        do i = 1, N
          p(i) = z(i) + beta * p(i)
        end do
      end if

      ! w = Ap
      !call D_mult_Ax(N, A, p, w)
call dspmv('U',N,1.0d0,A,p,1,0.0d0,w,1)
      ! alpha = tau2/p^T * w
      temp1 = 0.0d0
      !do i = 1, N
      !  temp1 = temp1 + p(i)*w(i)
      !end do
      temp1=ddot(N,p,1,w,1)
      alpha = tau2 / temp1

      do i = 1, N
        x(i) = x(i) + alpha * p(i)
      end do
      do i = 1, N
        r(i) = r(i) - alpha * w(i)
      end do


      rho = 0.0d0
      !do i = 1, N
      !  rho = rho + r(i)*r(i)
      !end do
      rho=ddot(N,r,1,r,1)
    end do  ! loop k 

    ! failed to convergence
    dPCG = .false.
    ! return the current residual
    Resid = dsqrt(rho)

  END FUNCTION dPCG

  subroutine D_mult_Ax(N, A, x, y)
    !
    ! (double precision)
    ! Compulte the multiplicatiopn of matrix A and vector x
    ! This subroutine should be provide by the user according
    ! to the storage of matrix A
    ! In our case, the A is symmetry and the upper half is stored 
    ! in a one dimension array.
    !
    ! Argument:    
    !
    ! N        (inpupt) Integer
    !          On entry, the dimension of the matrix.
    !          On exit, Unchanged.
    !
    ! A        (input) DBL
    !          On entry, the coefficient matrix
    !          On exit, unchanged.
    !
    ! x        (input) DBL, dimension N
    !          On entry, right hand side vector B
    !          On exit, Unchanged.
    !
    ! y        (out) DBL, dimension N
    !          On exit, y = Ax
    !
    
    INTEGER, INTENT(in) :: N
    REAL(KIND=DBL), DIMENSION(:), INTENT(in)  :: A, x
    REAL(KIND=DBL), DIMENSION(:), INTENT(out) :: y

    ! local variables
    INTEGER :: i, j, k, kk
    REAL(KIND=DBL) :: temp1, temp2

    do i = 1, N
      y(i) = 0.0d0
    end do

    kk = 1
    do j = 1, N
      temp1 = x(j)
      temp2 = 0.0d0
      k = kk
      do i = 1, j-1
        y(i) = y(i) + temp1 * A(k)
        temp2 = temp2 + A(k) * x(i)
        k = k + 1
      end do
      y(j) = y(j) + temp1 * A(kk + j -1) + temp2
      kk = kk + j
    end do
  end subroutine D_mult_Ax

  subroutine D_preconditioner(N, A, x, y)
    !
    ! (double preicsion)
    ! Compulte the multiplicatiopn of preconditioner M and vector x
    ! This subroutine should be provide by the user according
    ! what kind of preconditioner they what to use, here we just 
    ! try a very simple preconditioner, Jacobi preconditioning, where
    ! M is the inverse of the diagonal part of of A.
    ! Note: In our case, the A is symmetry and the upper half is stored 
    ! in a one dimension array.
    !
    ! Argument:    
    !
    ! N        (inpupt) Integer
    !          On entry, the dimension of the matrix.
    !          On exit, Unchanged.
    !
    ! A        (input) DBL
    !          On entry, the coefficient matrix
    !          On exit, unchanged.
    !
    ! x        (input) DBL, dimension N
    !          On entry, right hand side vector B
    !          On exit, Unchanged.
    !
    ! y        (out) DBL, dimension N
    !          On exit, y = Mx
    !
    
    INTEGER, INTENT(in) :: N
    REAL(KIND=DBL), DIMENSION(:), INTENT(in)  :: A, x
    REAL(KIND=DBL), DIMENSION(:), INTENT(out) :: y

    ! local variables
    INTEGER :: i

    do i = 1, N
      y(i) =  x(i) / A(i*(i-1)/2+i)
    end do

  end subroutine D_preconditioner

END MODULE LinearSolver




