Program main			 
!...Hydraulic fracture - crack(s) in unbounded domain
!...In this compilation, we solve crack(s) in unbounded domain using remote stress and pressure in crack(s). A pressure flag with two possible
!...non-zero values (1 for volume constraint and 2 for pressure constraint) is invoked. 1 refers to a strategy to simulate crack propagation
!...due to linearly dependent injection pressures (Pg.7 in Andrew's Report). 2 refers to a strategy to simulate crack propagation due to
!...linearly dependent injection volumes (Pg.8 in Andrew's Report).  
	USE GlobalData		 
	!...defined in global.f90,module globaldata uses module definitionconstant also
	IMPLICIT NONE
	INTEGER  	:: t0,t1,t2,t3,t4,t5,t6,mclock,i
	nstep = 1
	t0 = mclock()
	call Datapreparation 
	!...defined in prep.f90,module remesh defined in remesh_module.f90
	!...read the input file, allocate space for variables
    t1 = mclock()
    call FormAKF	                             
	!...defined in formakf.f90,uses definitionconstant,globaldata and element matrix defined in elem.f90
	!...forms the K matrix and force vector
    t2 = mclock()
    call RearrangeAKF	                         
	!...defined in rearge.f90, uses globaldata,definitionconstant 
    !...rearranges K matrix, force vector based on rigid body constraints
    t3 = mclock()
    call Processing	                             
	!...defined in proc.f90, uses globaldata and linearsolver defined in solver.f90
	!...solves the AX=b system
    t4 = mclock()
    call Post									 
	!...defined in post_remesh.f90,uses definitionconstant, globaldata and remesh
    !...outputs the results
	t5 = mclock()
    if (inner_flag.eq.1) then
!      	 call RigidInner
!        call RigidInner1
!        call RigidInner2
!        call RigidInner3(9)
		call RigidInner4
    endif
    if (tstress_flag.eq.YES) then
      	!...T-stress calculation by 1st approach, integrate on tip element only
    	call Tstress
        !...T-stress calculation by 1st approach, integrate on whole crack
        call Tstress1
        !...T-stress calculation by 2nd approach, integrate on whole crack, sum_uk at SBBL computed from integral equations
        !call Tstress2
        !...T-stress calculation by 2nd approach, integrate on whole crack, sum_uk at SBBL directly from solution vector
        call Tstress2a
        !...T-stress calculation by 2nd approach, integrate on tip element only (still have a bug when evaluate sum_uk at the other node of tip element)
        !call Tstress3
    endif
    t6 = mclock()
    !...display CPU time for each step of analysis
    print*,'---------------------------------------------------'
    print*,'          CPU TIME FOR EACH ANALYSIS STEP          '
    print*,'---------------------------------------------------'
    print'(3x,a30,f9.2,a7)','Data preparation             :',(t1-t0)/1000.0,'sec'
    print'(3x,a30,f9.2,a7)','Forming AK, F                :',(t2-t1)/1000.0,'sec'
    print'(3x,a30,f9.2,a7)','Constraint disp              :',(t3-t2)/1000.0,'sec'
    print'(3x,a30,f9.2,a7)','Solving linear eqns          :',(t4-t3)/1000.0,'sec'
    print'(3x,a30,f9.2,a7)','Post processing              :',(t5-t4)/1000.0,'sec'
    print'(3x,a30,f9.2,a7)','Tstress calculation          :',(t6-t5)/1000.0,'sec'
    print'(3x,a30,f9.2,a7)','Total CPU time               :',(t6-t0)/1000.0,'sec'
    if (growth_flag.eq.1) then
      	do i = 1,total_step
        	nstep = nstep + 1
        	call FormAKF
            call RearrangeAKF
            call Processing
            call Post
        enddo       
    endif 	
End program main
