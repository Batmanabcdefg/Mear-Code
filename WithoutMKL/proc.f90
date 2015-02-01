SUBROUTINE Processing
    ! Solve the linear equations
    USE GlobalData
    USE LinearSolver
    IMPLICIT NONE

    INTEGER :: ck

    LOGICAL :: successful_traction
    
	LOGICAL :: successful_pressure
	
	!-------------------------------------------------
	
	IF (sigma_x == 0.0 .AND. sigma_y == 0.0) THEN

	    PRINT *,'NO REMOTE STRESS!'

	ELSE
    
	successful_traction = dPCGSolver(total_dof, dAK_traction, dF_traction)

    IF (successful_traction) THEN
        
		PRINT *, 'Successfully solved the linear system for tractions!'
       
    ELSE
        
		PRINT *, 'Failed to solve the linear system for tractions, the system may be singular!'
        
		stop
    
	END IF

	END IF

    !---------------------------------------------------
	
	IF (pressure_value == 0.d0) THEN

	    PRINT *,'NO PRESSURE!'

	ELSE
	
        	do ck = 1,total_cracks

       
	   
	   successful_pressure = dPCGSolver(total_dof, dAK_pressure(ck,:), dF_pressure(ck,:))

       !.....this needs to be generalized to n cracks, not just 2....4/24/12...aje567

       !.....should be able to improve efficiency by solving multiple RHS at the same time

       IF (successful_pressure) THEN
          
		  PRINT *, 'Successfully solved the linear system for pressure for crack', ck, '!'
       
	   ELSE
          
		  PRINT *, 'Failed to solve the linear system for pressure for crack', ck, 'the system may be singular!'
          
		  stop
       
	   END IF
	
	enddo

	END IF

	!-----------------------------------------------------

END SUBROUTINE Processing







