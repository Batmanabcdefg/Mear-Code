!...This subroutine is completely new by Han
!...Set integration points and weights for 1D: 1->max_nint=60 points
!...This subroutine is similar to setint.f90 of fadd3d
!...Example: xi(:,2),wi(:,2) are 2 points and weights of 2-point Gauss
Subroutine InitGaussIntegration
	USE DefinitionConstant
	USE GlobalData
	IMPLICIT NONE

	INTEGER		:: i,ierr

    !...allocate space for 1D integration points and weights
    allocate(xi(max_nint,max_nint),wi(max_nint,max_nint),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'xi,wi: allocation request denied!'
        stop
    endif
    allocate(xi_log(max_logint,max_logint),wi_log(max_logint,max_logint),STAT=ierr)
    if (ierr.ne.0) then
      	print*,'xi_log,wi_log: allocation request denied'
        stop
    endif

    !...initialize integration points and weights (regular Gauss integral)
    do i=1,max_nint
      	call GaussIntPoint(i,xi(:,i),wi(:,i))
    enddo
    
    !...initialize integration points and weights for logarith integral
    do i=2,max_logint
      	call LogIntPoint(i,xi_log(:,i),wi_log(:,i))
    enddo
    
End Subroutine InitGaussIntegration