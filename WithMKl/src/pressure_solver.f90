SUBROUTINE Pressure_Solver(n,SIF_trac,SIF_pres,pressure_scale_temp)

  USE DefinitionConstant
  USE GlobalData
  IMPLICIT NONE
    
    INTEGER , INTENT(IN) :: n !total_cracks
    REAL(KIND=DBL),INTENT(IN), DIMENSION(NDOF,total_tip_node) :: SIF_trac  !(NDOF,total_tip_node)
    REAL(KIND=DBL),INTENT(IN), DIMENSION(total_cracks,NDOF,total_tip_node) :: SIF_pres  !(total_cracks,NDOF,total_tip_node)

    REAL(KIND=DBL),INTENT(OUT), DIMENSION(total_cracks,total_tip_node):: pressure_scale_temp   !(total_cracks,total_tip_node)

    ! local variables
    INTEGER ::  j , k
    REAL(KIND=DBL), DIMENSION(total_tip_node) :: SIF_pres_total

    do k = 1,total_tip_node

       SIF_pres_total(k) = 0

       do j = 1,n   ! total_cracks

          SIF_pres_total(k) = SIF_pres_total(k)+SIF_pres(j,1,k)

       enddo

       do j=1,n   ! total_cracks

          pressure_scale_temp(j,k) = (globalK(1)-SIF_trac(1,k))/SIF_pres_total(k)
          !...globalK(1) is fed as 1 in input file
		  !...SIF_trac(1,k) is K_I on crack tip k due to remote tractions
		  !...SIF_pres_total(k) is total K_I due to pressure on kth crack tip   
       enddo
    enddo

END SUBROUTINE PRESSURE_SOLVER
