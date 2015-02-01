SUBROUTINE RearrangeAKF
! Rearranging the global matrix AK and load vector F to remove 
! the rigid body motion
!
! Called by main
! 
! Calls  none
!
USE GlobalData
USE DefinitionConstant
IMPLICIT NONE
  
INTEGER :: i, j, k, m, ck
integer*4 :: ig1, mk
  
! do i=1,total_dof
!  print "(f10.7)",dF(i)
! end do
    !print'(f10.3)',ht
	IF (total_rigid_constrain_point == 0) RETURN
    ! remove rigid body motion
	!print'(f10.3)',ht
    DO i = 1, total_rigid_constrain_point
        DO j = 1,NDOF
            IF (rigid_constrain_dir(i,j) == 1) THEN  
                ! constrained in this direction, determine th global index of this constrain
                k = rigid_constrain_point_eqn(i) + j -1
                !...position in column-vector of the term right before the column k
                ig1 = k*(k-1)/2
                !...set all terms in column k equal to zero (here, only upper-half is set to zero)
                do m = 1, k-1
                    dAK(ig1+m) = 0.0d0  ! AK(:,k)
                    dAK_traction(ig1+m) = 0.0d0  ! AK_traction(:,k)
                    do ck = 1,total_cracks
                       dAK_pressure(ck,ig1+m) = 0.0d0  ! AK_pressure1(:,k)
                    enddo

                end do
                !...set all terms in that row equal to zero
                do m = k+1, total_dof
                    mk = m*(m-1)/2 +k
                    dAK(mk) = 0.0d0  ! Ak(k,:)
                    dAK_traction(mk) = 0.0d0  ! Ak_traction(k,:)
                    do ck = 1,total_cracks
                       dAK_pressure(ck,mk) = 0.0d0  ! Ak_pressure(ck,k,:)
                    enddo


                end do
                !...set the corresponding term in load vector equal to zero

                dF(k) = 0.0d0
                dF_traction(k) = 0.0d0
                
                do ck = 1,total_cracks
                   dF_pressure(ck,k) = 0.0d0
                enddo

                !...set the diagonal term to 1 (this may be left as it is, will consider carefully later)
                dAK(k*(k-1)/2+k) = 1.0d0   ! AK(k,k)
                dAK_traction(k*(k-1)/2+k) = 1.0d0   ! AK_traction(k,k)

                do ck = 1,total_cracks
                   dAK_pressure(ck,k*(k-1)/2+k) = 1.0d0   ! AK_pressure(ck,k,k)
                enddo

            END IF
        END DO
    END DO


END SUBROUTINE RearrangeAKF

