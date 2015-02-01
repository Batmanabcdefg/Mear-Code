!...This subroutine is copied from Xiao's code, adding comments for clear
SUBROUTINE D_assemble(region, isrc, ifld, sequ, fequ, node_src, node_fld, ec, ea, &
                      eat, eb, ed, gk, af, af_traction, af_pressure)
        ! double precison
        USE DefinitionConstant
        USE GlobalData
        IMPLICIT NONE

        ! Arguments
        INTEGER, INTENT(in)   :: region, isrc, ifld, sequ(:), fequ(:), node_src, node_fld 
        REAL(KIND=DBL), INTENT(in)  :: ec(:,:), ea(:,:), eat(:,:), eb(:,:), ed(:,:)
        REAL(KIND=DBL), INTENT(out) :: gk(:), af(:), af_traction(:), af_pressure(:,:)


        ! local variables
        REAL(kind=DBL) :: eg(3*NDOF,3*NDOF), ef(3*NDOF),ef_traction(3*NDOF),ef_pressure(total_cracks,3*NDOF), eg_u(3*NDOF,3*NDOF),eg_t(3*NDOF,3*NDOF),&
                          ef_u(3*NDOF),ef_traction_u(3*NDOF),ef_pressure_u(total_cracks,3*NDOF),ef_t(3*NDOF),ef_traction_t(3*NDOF),&
                          ef_pressure_t(total_cracks,3*NDOF)
        
        INTEGER   :: i, j, ki, kj, ig, jg, I_sign, ck
        INTEGER*4 :: ij
        logical   :: SRC_INTERFACE, FLD_INTERFACE

        SRC_INTERFACE = .false.
        FLD_INTERFACE = .false.


        if (ELTYPE(elemid(isrc)) == INTERFE) SRC_INTERFACE = .true.
        if (ELTYPE(elemid(ifld)) == INTERFE) FLD_INTERFACE = .true.

  
        !...construct the element stiffness matrix eg and force vector ef
        call Compute_Egf(eg, ef, ef_traction, ef_pressure, eg_u, eg_t, ef_u,ef_traction_u,ef_pressure_u,&
                         ef_t,ef_traction_t,ef_pressure_t)

 !       print*,'-----------------------------------'
 !       print*,ef
 !       print*,'-----------------------------------'
 !       print*,ef_traction
 !       print*,'-----------------------------------'
 !       print*,'The elemental RHS is: '
 !       print*,ef_pressure
 !       print*,'-----------------------------------'

!....................................................................................
!print*, " "
!print*, "Sequ:", isrc
!do i = 1, node_src
!    print*, sequ(i)
!enddo
!print*, "Fequ:", ifld
!do i = 1, node_fld
!    print*, sequ(i)
!enddo
!print*, " "

!print*, "node_src: ", node_src
!print*, "node_fld: ", node_fld

!....................................................................................



        if ((.not.SRC_INTERFACE) .and. (.not.FLD_INTERFACE)) then
    ! neither source nor field are INTERFACE element
            DO i = 1, node_src               ! loop over nodes of source element(row)
                 DO ki = 1, NDOF                    ! loop over components of each node(row)
                      ! just for "corner problem" roller constrain case
                      ig = sequ(i) + ki -1             ! global index corresponding ki(row)
                      ! added for "corner problem" constrain case
                      !if (elnode_switch(ki,i,isrc) == 0) cycle  ! no contribution
                      af(ig) = af(ig) + ef(NDOF*(i-1)+ki) ! assemble "load" vector F
                      !...added traction and pressure specific load vectors ....3/5/12...aje567
                      af_traction(ig) = af_traction(ig) + ef_traction(NDOF*(i-1)+ki) !assemble "load" vector F for tractions only

                      do ck = 1,total_cracks
                            af_pressure(ck,ig) = af_pressure(ck,ig) + ef_pressure(ck,NDOF*(i-1)+ki) !assemble "load" vector F for pressure only
                      enddo

                      DO j = 1, node_fld               ! loop over nodes of field element(col)
                           DO kj = 1, NDOF              ! loop over components of each node(col)
                                 jg = fequ(j) + kj -1    ! global index corresponding kj(col)
!...........................................................................
!print*, isrc, ifld, i, j
!        if( isrc .ne. 1 .or. ifld .ne. 1) then
!            if ( isrc .eq. ifld ) then
!                if( i .eq. 1 ) then
!                   ig = ig - 3
!                endif
!                if( j .eq. 1 ) then 
!                   jg = jg - 3
!                endif
!            elseif ( ifld .eq. isrc - 1) then
!                if( i .eq. 1 ) then
!                    ig = ig - 3
!                    !jg = jg + 2
!                endif    
!!            elseif ( isrc .eq. ifld - 1) then
!!                if( j .eq. 1 ) then
!!                    ig = ig + 3
!!                    jg = jg - 3
!!                endif
!            else
!                if (j .eq. 1) then
!                    jg = jg - 3
!                endif
!            endif
!        endif
!...........................................................................
                                 ! added for "corner problem" constrain case
                                 !if (elnode_switch(kj,j,ifld) == 0 ) cycle ! no contribution
                                 if (ig > jg) CYCLE 
                                 ! global index in one-dimensional array(upper half)
                                 ij = jg*(jg-1)/2 + ig
                                 ! assemble global matrix AK
                                 gk(ij) = gk(ij) + eg(NDOF*(i-1)+ki,NDOF*(j-1)+kj)
!print*, "gk", ig, jg, ij, "eg", NDOF*(i-1)+ki, NDOF*(j-1)+kj
!...........................................................................
!        if( isrc .ne. 1 .or. ifld .ne. 1) then
!            if ( isrc .eq. ifld ) then
!                if( i .eq. 1 ) then
!                   ig = ig + 3
!                endif
!                if( j .eq. 1 ) then 
!                   jg = jg + 3
!                endif
!            elseif ( ifld .eq. isrc - 1) then
!                if( i .eq. 1 ) then
!                    ig = ig + 3
!                    !jg = jg - 2
!                endif    
!!            elseif ( isrc .eq. ifld - 1) then
!!                if( j .eq. 1 ) then
!!                    ig = ig - 3
!!                    jg = jg + 3
!!                endif
!            else
!                if (j .eq. 1) then
!                    jg = jg + 3
!                endif
!            endif
!        endif
!!...........................................................................
                           END DO
                           IF (node_id(elnode(j,ifld)) == SBLP .and. &
                                 (ELTYPE(elemid(ifld)) == CREGULAR .or. &
                                  ELTYPE(elemid(ifld)) == CTIP)) THEN
                               ! assemble to C-
                               DO kj = 1, NDOF            ! loop over components of each node(col)
                                     jg = fequ(j) + kj + NDOF - 1   ! global index corresponding kj(col)
!...........................................................................
!        if( isrc .ne. 1 .or. ifld .ne. 1) then
!            if ( isrc .eq. ifld ) then
!                if( i .eq. 1 ) then
!                   ig = ig - 3
!                endif
!                if( j .eq. 1 ) then 
!                   jg = jg - 3
!                endif
!            elseif ( ifld .eq. isrc - 1) then
!                if( i .eq. 1 ) then
!                   ig = ig - 3
!                   !jg = jg + 2
!                endif    
!!            elseif ( isrc .eq. ifld - 1) then
!                if( j .eq. 1 ) then
!                    ig = ig + 3
!                    jg = jg - 3
!                endif
!            else
!                if (j .eq. 1) then
!                    jg = jg - 3
!                endif
!            endif
!        endif
!...........................................................................
                                     if (ig > jg) CYCLE
                                     ! global index in one-dimensional array(upper half)
                                     ij = jg*(jg-1)/2 + ig
                                     ! assemble global matrix AK
                                     gk(ij) = gk(ij) - eg(NDOF*(i-1)+ki,NDOF*(j-1)+kj)
!...........................................................................
!        if( isrc .ne. 1 .or. ifld .ne. 1) then
!            if ( isrc .eq. ifld ) then
!                if( i .eq. 1 ) then
!                   ig = ig + 3
!                endif
!                if( j .eq. 1 ) then 
!                   jg = jg + 3
!                endif
!1            elseif ( ifld .eq. isrc - 1) then
!1                if( i .eq. 1 ) then
!                   ig = ig + 3
!                   !jg = jg - 2
!                endif    
!!            elseif ( isrc .eq. ifld - 1) then
!!                if( j .eq. 1 ) then
!!                    ig = ig - 3
!!                    jg = jg + 3
!!                endif
!            else
!                if (j .eq. 1) then
!                    jg = jg + 3
!                endif
!1            endif
!        endif
!!...........................................................................
                               END DO
                           END IF
                      END DO! loop j
                 END DO! loop ki

                 IF (node_id(elnode(i,isrc)) == SBLP .and.  &
                          (ELTYPE(elemid(isrc)) == CREGULAR .or.   &
                           ELTYPE(elemid(isrc)) == CTIP)) THEN
                     DO ki = 1, NDOF
                           ig = sequ(i) + ki + NDOF - 1            ! global index for C- 
                           af(ig) = af(ig) - ef(NDOF*(i-1)+ki) ! assemble "load" vector F
                           af_traction(ig) = af_traction(ig) - ef_traction(NDOF*(i-1)+ki) ! assemble "load" vector F for traction
                                
                           do ck = 1,total_cracks
                                 af_pressure(ck,ig) = af_pressure(ck,ig) - ef_pressure(ck,NDOF*(i-1)+ki) ! assemble "load" vector F for pressure
                           enddo

                           DO j = 1, node_fld               ! loop over nodes of field element(col)
                                DO kj = 1, NDOF             ! loop over component of each node(col)
                                      jg = fequ(j) + kj -1      ! global index corresponding kj(col)
                                      if (ig > jg) CYCLE
                                      ! global index in one-dimensional array
                                      ij = jg*(jg-1)/2 + ig
                                      ! assemble global matrix AK
                                      gk(ij) = gk(ij) - eg(NDOF*(i-1)+ki,NDOF*(j-1)+kj)
                                END DO
                                IF (node_id(elnode(j,ifld)) == SBLP .and. &
                                        (ELTYPE(elemid(ifld)) == CREGULAR .or.  &
                                         ELTYPE(elemid(ifld)) == CTIP)) THEN
                                    ! assemble to C-
                                    DO kj = 1, NDOF           ! loop over component of each node(col)
                                          jg = fequ(j) + kj + NDOF - 1   ! global index corresponding kj(col)
                                          if (ig > jg) CYCLE
                                          ! global index in one-dimensional array(upper half)
                                          ij = jg*(jg-1)/2 + ig
                                          ! assemble global matrix AK
                                          gk(ij) = gk(ij) + eg(NDOF*(i-1)+ki,NDOF*(j-1)+kj)
                                     END DO
                                 END IF
                           END DO! loop j
                     END DO! loop ki
                 END IF
            END DO! loop i
        else if((.not.SRC_INTERFACE) .and. FLD_INTERFACE) then
            ! just the field element is INTERFACE element
            ! no load vector for this case
            do i = 1, node_src
                 DO ki = 1, NDOF                    ! loop over component of each node(row)
                       ig = sequ(i) + ki -1             ! global index corresponding ki(row)
                       DO j = 1, node_fld               ! loop over nodes of field element(col)
                            DO kj = 1, NDOF              ! loop over component of each node(col)
                                  ! displacement unknown on INTERFACE
                                  jg = fequ(j) + kj - 1       ! global index corresponding kj(col)
                                  if (ig <= jg) then 
                                      ! global index in one-dimensional array(upper half)
                                      ij = jg*(jg-1)/2 + ig
                                      ! assemble global matrix AK
                                      gk(ij) = gk(ij) + eg_u(NDOF*(i-1)+ki,NDOF*(j-1)+kj)
                                  end if
                                  ! traction unknown on INTERFACE
                                  jg = fequ(j) + kj + NDOF - 1       ! global index corresponding kj(col)
                                  if (ig <= jg) then 
                                      ! global index in one-dimensional array(upper half)
                                      ij = jg*(jg-1)/2 + ig
                                      ! assemble global matrix AK
                                      gk(ij) = gk(ij) + eg_t(NDOF*(i-1)+ki,NDOF*(j-1)+kj)
                                  end if
                            END DO  ! loop kj
                       end DO  ! loop j
                 end DO  ! loop ki

                 ! for SBLP nodes
                 IF (node_id(elnode(i,isrc)) == SBLP .and.  &
                          (ELTYPE(elemid(isrc)) == CREGULAR .or.   &
                           ELTYPE(elemid(isrc)) == CTIP)) THEN
                     DO ki = 1, NDOF
                           ig = sequ(i) + ki + NDOF - 1            ! global index for C- 
                           DO j = 1, node_fld               ! loop over nodes of field element(col)
                                DO kj = 1, NDOF              ! loop over component of each node(col)
                                      ! displacement unknown on INTERFACE
                                      jg = fequ(j) + kj - 1       ! global index corresponding kj(col)
                                      if (ig <= jg) then 
                                          ! global index in one-dimensional array(upper half)
                                          ij = jg*(jg-1)/2 + ig
                                          ! assemble global matrix AK
                                          !...5/4/09: change the sign here
                                          gk(ij) = gk(ij) - eg_u(NDOF*(i-1)+ki,NDOF*(j-1)+kj)
                                      end if
                                      ! traction unknown on INTERFACE
                                      jg = fequ(j) + kj + NDOF -1       ! global index corresponding kj(col)
                                      if (ig <= jg) then 
                                          ! global index in one-dimensional array(upper half)
                                          ij = jg*(jg-1)/2 + ig
                                          ! assemble global matrix AK
                                          !...5/4/09: change the sign here
                                          gk(ij) = gk(ij) - eg_t(NDOF*(i-1)+ki,NDOF*(j-1)+kj)
                                      end if
                                END DO! loop kj
                           end DO! loop j
                     end DO! loop ki
                 end IF
            end do
        else if (SRC_INTERFACE .and. (.not. FLD_INTERFACE)) then
            ! source element is INTERFACE element
            do i = 1, node_src
                 DO ki = 1, NDOF                    ! loop over component of each node(row)
                      ! traction. eqn.
                      ig = sequ(i) + ki -1             ! global index corresponding ki(row)
                      af(ig) = af(ig) + ef_t(NDOF*(i-1)+ki) ! assemble "load" vector F
                      af_traction(ig) = af_traction(ig) + ef_traction_t(NDOF*(i-1)+ki)

                      do ck = 1,total_cracks
                           af_pressure(ck,ig) = af_pressure(ck,ig) + ef_pressure_t(ck,NDOF*(i-1)+ki)
                      enddo

                      DO j = 1, node_fld               ! loop over nodes of field element(col)
                           DO kj = 1, NDOF              ! loop over component of each node(col)
                               jg = fequ(j) + kj -1       ! global index corresponding kj(col)
                               if (ig > jg) CYCLE 
                               ! global index in one-dimensional array(upper half)
                               ij = jg*(jg-1)/2 + ig
                               ! assemble global matrix AK
                               gk(ij) = gk(ij) + eg_t(NDOF*(i-1)+ki,NDOF*(j-1)+kj)
                           END DO
                           IF (node_id(elnode(j,ifld)) == SBLP .and. &
                                    (ELTYPE(elemid(ifld)) == CREGULAR .or. &
                                     ELTYPE(elemid(ifld)) == CTIP)) THEN
                                ! assemble to C-
                                DO kj = 1, NDOF            ! loop over component of each node(col)
                                     jg = fequ(j) + kj + NDOF - 1    ! global index corresponding kj(col)
                                     if (ig > jg) CYCLE
                                     ! global index in one-dimensional array(upper half)
                                     ij = jg*(jg-1)/2 + ig
                                     ! assemble global matrix AK
                                     gk(ij) = gk(ij) - eg_t(NDOF*(i-1)+ki,NDOF*(j-1)+kj)
                                END DO
                           END IF
                      END DO! loop j
                      ! displacement eqn.
                      ig = sequ(i) + ki + NDOF - 1              ! global index corresponding ki(row)
                      af(ig) = af(ig) + ef_u(NDOF*(i-1)+ki) ! assemble "load" vector F
                      af_traction(ig) = af_traction(ig) + ef_traction_u(NDOF*(i-1)+ki) ! assemble "load" vector F

                      do ck = 1, total_cracks
                           af_pressure(ck,ig) = af_pressure(ck,ig) + ef_pressure_u(ck,NDOF*(i-1)+ki) ! assemble "load" vector F
                      enddo
                

                      DO j = 1, node_fld                 ! loop over nodes of field element(col)
                           DO kj = 1, NDOF              ! loop over component of each node(col)
                                jg = fequ(j) + kj -1       ! global index corresponding kj(col)
                                if (ig > jg) CYCLE 
                                    ! global index in one-dimensional array(upper half)
                                    ij = jg*(jg-1)/2 + ig
                                    ! assemble global matrix AK
                                    gk(ij) = gk(ij) + eg_u(NDOF*(i-1)+ki,NDOF*(j-1)+kj)
                           END DO
                           IF (node_id(elnode(j,ifld)) == SBLP .and. &
                                     (ELTYPE(elemid(ifld)) == CREGULAR .or. &
                                     ELTYPE(elemid(ifld)) == CTIP)) THEN
                               ! assemble to C-
                               DO kj = 1, NDOF            ! loop over component of each node(col)
                                     jg = fequ(j) + kj + NDOF - 1    ! global index corresponding kj(col)
                                     if (ig > jg) CYCLE
                                     ! global index in one-dimensional array(upper half)
                                     ij = jg*(jg-1)/2 + ig
                                     ! assemble global matrix AK
                                     gk(ij) = gk(ij) - eg_u(NDOF*(i-1)+ki,NDOF*(j-1)+kj)
                               END DO
                           END IF
                      END DO  ! loop j
                 END DO  ! loop ki
            end do  ! loop i
        else if (SRC_INTERFACE .and. FLD_INTERFACE) then
            ! both source and field element are INTERFACE element
            do i = 1, node_src
                 DO ki = 1, NDOF                    ! loop over component of each node(row)
                      ! traction eqn.
                      ig = sequ(i) + ki -1             ! global index corresponding ki(row)
                      DO j = 1, node_fld               ! loop over nodes of field element(col)
                           DO kj = 1, NDOF              ! loop over component of each node(col)
                                ! displacement unknown on INTERFACE
                                jg = fequ(j) + kj -1       ! global index corresponding kj(col)
                                if (ig <= jg) then
                                    ! global index in one-dimensional array(upper half)
                                    ij = jg*(jg-1)/2 + ig
                                    ! assemble global matrix AK
                                    gk(ij) = gk(ij) + I_sign * ed(NDOF*(i-1)+ki,NDOF*(j-1)+kj)
                                end if
                                ! traction unknown in INTERFACE
                                jg = fequ(j) + kj + NDOF - 1
                                if (ig <= jg) then
                                    ! global index in one-dimensional array(upper half)
                                    ij = jg*(jg-1)/2 + ig
                                    ! assemble global matrix AK
                                    gk(ij) = gk(ij) + I_sign * ea(NDOF*(i-1)+ki,NDOF*(j-1)+kj)
                                end if
                           end DO
                      END DO
                      ! displacement eqn.
                      ig = sequ(i) + ki + NDOF - 1            ! global index corresponding ki(row)
                      DO j = 1, node_fld               ! loop over nodes of field element(col)
                           DO kj = 1, NDOF                ! loop over component of each node(col)
                                ! displacement unknown on INTERFACE
                                jg = fequ(j) + kj -1       ! global index corresponding kj(col)
                                if (ig <= jg) then
                                    ! global index in one-dimensional array(upper half)
                                    ij = jg*(jg-1)/2 + ig
                                    ! assemble global matrix AK
                                    gk(ij) = gk(ij) + I_sign * eat(NDOF*(i-1)+ki,NDOF*(j-1)+kj)
                                end if
                                ! traction unknown in INTERFACE
                                jg = fequ(j) + kj + NDOF - 1
                                if (ig <= jg) then
                                    ! global index in one-dimensional array(upper half)
                                    ij = jg*(jg-1)/2 + ig
                                    ! assemble global matrix AK
                                    gk(ij) = gk(ij) - I_sign * eb(NDOF*(i-1)+ki,NDOF*(j-1)+kj)
                                end if
                           end DO
                      END DO
                 end DO  ! loop ki
            end do  ! loop i
     end if


!     print*,'-----------------------------------------------'
!     print*,'The total right hand side is: '
!     print*, af
!     print*,'-----------------------------------------------'
!     print*,'The right hand side due to tractions is: '
!     print*, af_traction
!     print*,'-----------------------------------------------'
!     print*,'The right hand side due to pressure (from assemb) is: '
!     print*, af_pressure
!     print*,'-----------------------------------------------'

     


CONTAINS
        SUBROUTINE Compute_Egf(eg,ef,ef_traction,ef_pressure,eg_u,eg_t,ef_u,ef_traction_u,ef_pressure_u,&
                               ef_t,ef_traction_t,ef_pressure_t)
        !...construct element stiffness eg, and load vector ef
        USE DefinitionConstant
        USE GlobalData
        IMPLICIT NONE

    REAL(KIND=DBL), INTENT(out) :: eg(3*NDOF,3*NDOF),ef(3*NDOF),ef_traction(3*NDOF),ef_pressure(total_cracks,3*NDOF), eg_u(3*NDOF,3*NDOF),&
                                   eg_t(3*NDOF,3*NDOF),ef_u(3*NDOF),ef_traction_u(3*NDOF),ef_pressure_u(total_cracks,3*NDOF),ef_t(3*NDOF),&
                                   ef_traction_t(3*NDOF),ef_pressure_t(total_cracks,3*NDOF)
  
    ! local variables
    !...ew is the matrix from which multiply with nodal value give force vector
    real(kind=DBL) :: ew(3*NDOF,3*NDOF),rb(3*NDOF),rb_traction(3*NDOF),rb_pressure(total_cracks,3*NDOF),ew_u(3*NDOF,3*NDOF),ew_t(3*NDOF,3*NDOF)
    integer :: i, j, ii, jj, alpha, beta, ck
    logical :: src_master = .false., fld_master = .false.
 
    ew = 0.0d0
    ew_u = 0.0d0
    ew_t = 0.0d0

    src_master = .false.
    fld_master = .false.

    if (region == elem_region(isrc)) src_master = .true.
    if (region == elem_region(ifld)) fld_master = .true.

    ! determin the I_sign
    if ((src_master .and. (.not. fld_master)) .or. &
        (fld_master .and. (.not. src_master))) then
      ! put "-" sign
      I_sign = -1
    else
      ! put "+" sign
      I_sign = 1
    end if

    ii = 0
    DO i = 1, node_src    ! over nodes of source element(equations)
         DO alpha = 1,NDOF  ! over NDOF "directions"
             ii = ii + 1
             SELECT CASE (elldid(alpha,isrc))
                CASE (BDISP)
                    ! using "disp. eqn"
                    jj = 0
                    DO j = 1, node_fld  ! over nodes of field element(unknowns)
                         DO beta = 1,NDOF ! over 3 directions
                             jj = jj + 1
                             SELECT CASE (elldid(beta,ifld))
                                 CASE (BDISP)
                                     ! unknownss are "TRACTION"
                                     eg(ii,jj) = -eb(ii,jj)
                                     IF (isrc == ifld) THEN
                                         ! source and field element are the same
                                         ew(ii,jj) = -ec(ii,jj) - eat(ii,jj)
                                     ELSE
                                         ! source and field element are different
                                         ew(ii,jj) = -eat(ii,jj)
                                     END IF
                                 !----------
                                 CASE(BTRAC, BTRFREE)
                                     IF (node_type(beta,elnode(j,ifld)) == EDGE_NODE) THEN
                                         ! for edge nodes(unknowns are "TRACTION")
                                         ! no contribations to EG
                                         eg(ii,jj) = 0.0d0
                                         ew(ii,jj) = -eat(ii,jj) 
                                     ELSE
                                         ! unknowns are "DISP"
                                         eg(ii,jj) = eat(ii,jj)
                                         ew(ii,jj) = eb(ii,jj)
                                     END IF
                                 !----------
                                 case (CTRAC, CTRFREE)
                                       IF (node_type(beta,elnode(j,ifld)) == EDGE_NODE) THEN
                                           ! for edge nodes(unknowns are "TRACTION")
                                           ! no contribations to EG
                                           eg(ii,jj) = 0.0d0
                                           ew(ii,jj) = -eat(ii,jj) 
                                        ELSE
                                           ! unknowns are "DISP"
                                           eg(ii,jj) = eat(ii,jj)
                                           ew(ii,jj) = 0.0d0 ! =eb(ii,jj), eb=0, for field element is on crack
                                        END IF
                                 !----------
                                 case (NOLOAD)  ! unknowns on INTERFACE surface, disp. eqn.
                                       !...not checked for 2D yet
                                       if (fld_master) then
                                           ! coefficients of u and t
                                           eg_u(ii,jj) = eat(ii,jj)
                                           eg_t(ii,jj) = -eb(ii,jj)
                                       else  ! put "-" sign
                                           ! coefficients of u and t
                                           eg_u(ii,jj) = -eat(ii,jj)
                                           eg_t(ii,jj) = eb(ii,jj)
                                       end if
                             END SELECT
                        END DO! loop beta
                    END DO! loop j
                !--------------------
                CASE (BTRAC, CTRAC, BTRFREE, CTRFREE)
                    ! using "traction eqn"
                    jj = 0
                    DO j = 1, node_fld  ! over nodes of field element(unknowns)
                        DO beta = 1,NDOF ! over 3 directions
                             jj = jj + 1
                             SELECT CASE (elldid(beta,ifld))
                                 CASE (BDISP)
                                     IF (node_type(alpha,elnode(i,isrc)) == EDGE_NODE) THEN
                                         ! unknowns are "TRACTION"
                                         ! "EDGE_NODE" has no contributions to the equation
                                         !...Su-Su
                                         eg(ii,jj) = 0.0d0
                                         ew(ii,jj) = 0.0d0
                                     ELSE
                                         ! unknowns are "TRACTION", using "traction eqn"
                                         !...St/Sc-Su
                                         eg(ii,jj) = ea(ii,jj)
                                         ew(ii,jj) = -ed(ii,jj)
                                     END IF
                                 !----------
                                 CASE (BTRAC)
                                     ! both source and field elements are traction prescribed
                                     IF (node_type(alpha,elnode(i,isrc)) == EDGE_NODE) THEN
                                         ! "EDGE_NODE" has no contribution to the equation
                                         !...Su-whatever
                                         eg(ii,jj) = 0.0d0
                                         ew(ii,jj) = 0.0d0
                                     ELSE
                                         IF (isrc == ifld) THEN
                                             IF (node_type(beta,elnode(j,ifld)) == EDGE_NODE) THEN
                                                 ! no contribution to "EDGE_NODE"
                                                 !...St/Sc-Su (same element)
                                                 eg(ii,jj) = 0.0d0   !...since we are on same element
                                                 ew(ii,jj) = -ed(ii,jj)   !...same element, but since node j is edge node ->using ed
                                             ELSE
                                                 !...St/Sc-St (same element)
                                                 eg(ii,jj) = ed(ii,jj)
                                                 ew(ii,jj) = ec(ii,jj) - ea(ii,jj)   !...this should be -ec ???
                                             END IF
                                         ELSE
                                             IF (node_type(beta,elnode(j,ifld)) == EDGE_NODE) THEN
                                                 ! no contribution to "EDGE_NODE" 
                                                 !...St/Sc-Su (different element)
                                                 eg(ii,jj) = 0.0d0
                                                 ew(ii,jj) = -ed(ii,jj)
                                             ELSE
                                                 !...St/Sc-St (different element)
                                                 eg(ii,jj) = ed(ii,jj) 
                                                 ew(ii,jj) = -ea(ii,jj)
                                             END IF
                                         END IF
                                     END IF
                                 !----------
                                 CASE (CTRAC)
                                     ! both source and field elements are traction prescribed
                                     IF (node_type(alpha,elnode(i,isrc)) == EDGE_NODE) THEN
                                         ! "EDGE_NODE" has no contribution to the equation
                                         !...Su-whatever
                                         eg(ii,jj) = 0.0d0
                                         ew(ii,jj) = 0.0d0
                                     ELSE
                                         IF (isrc == ifld) THEN
                                             IF (node_type(beta,elnode(j,ifld)) == EDGE_NODE) THEN
                                                 ! no contribution to "EDGE_NODE"
                                                 !...St/Sc-Su (same element)
                                                 eg(ii,jj) = 0.0d0
                                                 ew(ii,jj) = -ed(ii,jj)                     
                                             ELSE
                                                 !...St/Sc-Sc (same element)
                                                 eg(ii,jj) = ed(ii,jj)
                                                 ew(ii,jj) = ec(ii,jj) ! -ea(ii,jj), ea=0 for field element is on crack
                                             END IF
                                         ELSE
                                             IF (node_type(beta,elnode(j,ifld)) == EDGE_NODE) THEN
                                                 ! no contribution to "EDGE_NODE"
                                                 !...St/Sc-Su (different elements)
                                                 eg(ii,jj) = 0.0d0
                                                 ew(ii,jj) = -ed(ii,jj)
                                             ELSE
                                                 !...St/Sc-Sc (different elements)
                                                 eg(ii,jj) = ed(ii,jj) 
                                                 ew(ii,jj) = 0.0d0  !-ea(ii,jj), ea=0 for field element is on crack
                                             END IF
                                         END IF
                                     END IF
                                 !----------
                                 CASE (BTRFREE, CTRFREE)
                                     ! both source and field elements are traction prescribed
                                     IF (node_type(alpha,elnode(i,isrc)) == EDGE_NODE) THEN
                                         ! "EDGE_NODE" has no contribution to the equation
                                         !...Su-whatever
                                         eg(ii,jj) = 0.0d0
                                         ew(ii,jj) = 0.0d0
                                      ELSE
                                          !...St/Sc-
                                          IF (isrc == ifld) THEN
                                              IF (node_type(beta,elnode(j,ifld)) == EDGE_NODE) THEN
                                                  ! no contribution to "EDGE_NODE"
                                                  !...St/Sc-Su
                                                  eg(ii,jj) = 0.0d0
                                                  ew(ii,jj) = -ed(ii,jj)                     
                                              ELSE
                                                  !...St/Sc-St/Sc(traction free)
                                                  eg(ii,jj) = ed(ii,jj)
                                                  ew(ii,jj) = 0.0d0
                                              END IF
                                          ELSE
                                              IF (node_type(beta,elnode(j,ifld)) == EDGE_NODE) THEN
                                                  ! no contribution to "EDGE_NODE"
                                                  !...St/Sc-Su
                                                  eg(ii,jj) = 0.0d0
                                                  ew(ii,jj) = -ed(ii,jj)
                                              ELSE
                                                  !...St/Sc-St/Sc(traction free)
                                                  eg(ii,jj) = ed(ii,jj) 
                                                  ew(ii,jj) = 0.0d0
                                              END IF
                                          END IF
                                      END IF
                                 !----------
                                 case (NOLOAD)  ! unknowns on INTERFACE surface, traction eqn.
                                     IF (node_type(alpha,elnode(i,isrc)) == EDGE_NODE) THEN
                                         ! "EDGE_NODE" has no contributions to the traction equations
                                         eg_u(ii,jj) = 0.0d0
                                         eg_t(ii,jj) = 0.0d0
                                     else
                                         if (fld_master) then
                                             ! coefficients of u and t
                                             eg_u(ii,jj) = ed(ii,jj)
                                             eg_t(ii,jj) = ea(ii,jj)
                                         else  ! put "-" sign
                                             ! coefficients of u and t
                                             eg_u(ii,jj) = -ed(ii,jj)
                                             eg_t(ii,jj) = -ea(ii,jj)
                                         end if
                                     end IF
                             END SELECT
                        END DO  ! loop beta
                    END DO  ! loop j
                !--------------------
                CASE (NOLOAD)  ! source element is INTERFACE element
                    jj = 0
                    DO j = 1, node_fld  ! over nodes of field element(unknowns)
                         DO beta = 1,NDOF ! over 3 directions
                             jj = jj + 1
                             SELECT CASE (elldid(beta,ifld))
                                 CASE (BDISP)
                                     ! unknowns are "TRACTION"
                                     if (src_master) then 
                                         ! disp. eqn.
                                         eg_u(ii,jj) = -eb(ii,jj)
                                         ew_u(ii,jj) = -eat(ii,jj)
                                         ! traction eqn.
                                         eg_t(ii,jj) = ea(ii,jj)
                                         ew_t(ii,jj) = -ed(ii,jj)
                                     else ! put "-" sign
                                         ! disp. eqn.
                                         eg_u(ii,jj) = eb(ii,jj)
                                         ew_u(ii,jj) = eat(ii,jj)
                                         ! traction eqn.
                                         eg_t(ii,jj) = -ea(ii,jj)
                                         ew_t(ii,jj) = ed(ii,jj)
                                     end if

                                 CASE (BTRAC, BTRFREE)
                                     IF (node_type(beta,elnode(j,ifld)) == EDGE_NODE) THEN
                                         ! for edge nodes(unknowns are "TRACTION")
                                         ! no contribations to EG
                                         if (src_master) then
                                             ! disp. eqn.
                                             eg_u(ii,jj) = 0.0d0
                                             ew_u(ii,jj) = -eat(ii,jj) 
                                             ! traction eqn.
                                             eg_t(ii,jj) = 0.0d0
                                             ew_t(ii,jj) = -ed(ii,jj) 
                                         else  ! put "-" sign
                                             ! disp. eqn.
                                             eg_u(ii,jj) = 0.0d0
                                             ew_u(ii,jj) = eat(ii,jj) 
                                             ! traction eqn.
                                             eg_t(ii,jj) = 0.0d0
                                             ew_t(ii,jj) = ed(ii,jj) 
                                         end if
                                     ELSE  
                                         ! unknowns are "DISP"
                                         if (src_master) then
                                             ! disp. eqn.
                                             eg_u(ii,jj) = eat(ii,jj)
                                             ew_u(ii,jj) = eb(ii,jj)
                                             ! traction eqn.
                                             eg_t(ii,jj) = ed(ii,jj)
                                             ew_t(ii,jj) = -ea(ii,jj)
                                         else ! put "-" sign
                                             ! disp. eqn.
                                             eg_u(ii,jj) = -eat(ii,jj)
                                             ew_u(ii,jj) = -eb(ii,jj)
                                             ! traction eqn.
                                             eg_t(ii,jj) = -ed(ii,jj)
                                             ew_t(ii,jj) = ea(ii,jj)
                                         end if
                                     END IF

                                 case (CTRAC, CTRFREE)
                                     IF (node_type(beta,elnode(j,ifld)) == EDGE_NODE) THEN
                                         ! for edge nodes(unknowns are "TRACTION")
                                         ! no contribations to EG
                                         if (src_master) then
                                             ! disp. eqn.
                                             eg_u(ii,jj) = 0.0d0
                                             ew_u(ii,jj) = -eat(ii,jj) 
                                             ! traction eqn.
                                             eg_t(ii,jj) = 0.0d0
                                             ew_t(ii,jj) = -ed(ii,jj) 
                                         else  ! put "-" sign
                                             ! disp. eqn.
                                             eg_u(ii,jj) = 0.0d0
                                             ew_u(ii,jj) = eat(ii,jj) 
                                             ! traction eqn.
                                             eg_t(ii,jj) = 0.0d0
                                             ew_t(ii,jj) = ed(ii,jj) 
                                          end if
                                     ELSE  
                                         ! unknowns are "DISP"
                                         if (src_master) then
                                             ! disp. eqn.
                                             eg_u(ii,jj) = eat(ii,jj)
                                             ew_u(ii,jj) = 0.0d0  ! eb(ii,jj), eb=0 for field element is on crack(disp. eqn)
                                             ! traction eqn.
                                             eg_t(ii,jj) = ed(ii,jj)
                                             ew_t(ii,jj) = 0.0d0  ! -ea(ii,jj),ea=0 for field element is on crack(trac. eqn)
                                         else ! put "-" sign
                                             ! disp. eqn.
                                             eg_u(ii,jj) = -eat(ii,jj)
                                             ew_u(ii,jj) = 0.0d0  !-eb(ii,jj)
                                             ! traction eqn.
                                             eg_t(ii,jj) = -ed(ii,jj)
                                             ew_t(ii,jj) = 0.0d0  ! ea(ii,jj)
                                         end if
                                     END IF
                
                                 case (NOLOAD)  ! unknowns on INTERFACE surface, disp. eqn.
                                       ! put proper sign

                             end SELECT
                         end DO !...of beta
                    end DO !...of j
             END SELECT
         END DO! loop alpha
    END DO! loop i
    !--------------------------------------------------
    ! begin to form the right-hand side terms
    jj = 0
    DO j = 1, node_fld
         !...Load vector involves only on field element, since for the case of ec (source element)
         !->source element = field element
         DO beta = 1,NDOF
              jj = jj + 1
              SELECT CASE (elldid(beta,ifld))
                  CASE (BDISP)
                      ! known is displaceemnt
                      rb(jj) = elnode_val(beta,j,ifld)
                      !----------
                  CASE (BTRAC, CTRAC)
                      ! traction prescribed(field element)
                      !         IF (unknown(beta, elnode(j,ifld)) == DISPLACEMENT) THEN
                      !.........Modified by Jaroon on 01/17/2006
                      IF (unknown(beta, elnode(j,ifld)) == DISPLACEMENT.or.unknown(beta,elnode(j,ifld))==DISP_TRAC) THEN
                          ! known value is traction
                          !...This is the case of field element has some node with known traction, the other with known displacement (edge_node)
                          rb(jj) = elnode_val(beta,j,ifld)
                          rb_traction(jj) = elnode_val_trac(beta,j,ifld)
                        
                          do ck = 1, total_cracks
                              rb_pressure(ck,jj) = elnode_val_pressure(ck,beta,j,ifld)
                          enddo

                        

!                            print*,'---------------------------------------------------'
!                            print*,rb
!                            print*,'---------------------------------------------------'
!                            print*,rb_traction
!                            print*,'---------------------------------------------------'
!                            print*,'rb_pressure is: '
!                            print*,rb_pressure
!                            print*,'---------------------------------------------------'


                      ELSE
                          ! unknown is traction
                          rb(jj) = node_disp_val(beta,elnode(j,ifld))
                      END IF
                        !----------
                  CASE (BTRFREE, CTRFREE) 
                      IF (unknown(beta, elnode(j,ifld)) == DISPLACEMENT) THEN
                          rb(jj) = 0.0d0
                      ELSE
                          ! unknown is traction
                          rb(jj) = node_disp_val(beta,elnode(j,ifld))
                      END IF
                  !----------
                  case (NOLOAD)  ! INTERFACE element
                      ! do nothing
                      rb(jj) = 0.0d0

                  CASE default
                      WRITE(*,*) "ERROR : Unknow boundary type!"
                      WRITE(*,*) '(',beta,',',ifld,')=',elldid(beta,ifld)
              END SELECT
         END DO! loop over beta = 1,NDOF
    END DO! loop over j = 1, nodey
    !...multiply the "pre" load vector matrix ew with nodal value of field element (also source element for ec)
    DO i = 1, NDOF*node_src
         ef(i) = 0.0d0
         ef_u(i) = 0.0d0
         ef_t(i) = 0.0d0

         ef_traction(i) = 0.0d0
         ef_traction_u(i) = 0.0d0
         ef_traction_t(i) = 0.0d0

         do ck = 1, total_cracks
               ef_pressure(ck,i) = 0.0d0
               ef_pressure_u(ck,i) = 0.0d0
               ef_pressure_t(ck,i) = 0.0d0
         enddo

         DO j = 1, NDOF*node_fld
              ef(i) = ef(i) + ew(i,j) * rb(j)
              ef_u(i) = ef_u(i) + ew_u(i,j) * rb(j)
              ef_t(i) = ef_t(i) + ew_t(i,j) * rb(j)

              ef_traction(i) = ef_traction(i) + ew(i,j) * rb_traction(j)
              ef_traction_u(i) = ef_traction_u(i) + ew_u(i,j) * rb_traction(j)
              ef_traction_t(i) = ef_traction_t(i) + ew_t(i,j) * rb_traction(j)

              do ck = 1, total_cracks
                    ef_pressure(ck,i) = ef_pressure(ck,i) + ew(i,j) * rb_pressure(ck,j)
                    ef_pressure_u(ck,i) = ef_pressure_u(ck,i) + ew_u(i,j) * rb_pressure(ck,j)
                    ef_pressure_t(ck,i) = ef_pressure_t(ck,i) + ew_t(i,j) * rb_pressure(ck,j)
              enddo

         END DO
    END DO

!    print*,'---------------------------------------------------'
!    print*,ef
!    print*,'---------------------------------------------------'
!    print*,ef_traction
!    print*,'---------------------------------------------------'
!    print*,'The element RHS from Comp_Egf is: '
!    print*,ef_pressure
!    print*,'---------------------------------------------------'




!    !...debugging
!    print*,isrc,ifld
!    do i = 1,3
!              print'(3f15.8)',ef(3*i-2),ef(3*i-1),ef(3*i)
!    enddo

! if (isrc.eq.1.or.ifld.eq.1) then
!  print *,isrc,ifld
!  do i=1,9
!    print "(6f15.8)", ef(3*i-2),ef(3*i-1),ef(3*i),rb(3*i-2),rb(3*i-1),rb(3*i)
!  end do
!  print *
!  do i=1,9
!    print "(6f15.8)", ef_u(3*i-2),ef_u(3*i-1),ef_u(3*i),ef_t(3*i-2),ef_t(3*i-1),ef_t(3*i)
!  end do        
!  end if
!  do i=1,27
!   do j=1,27
!    print "(2i7,2f15.8)", i,j,ew(i,j)
!        end do
!  end do
! end if

!print*, "eg:"
!do i = 1, 9
!print "(ES10.3, ES10.3, ES10.3, ES10.3, ES10.3, ES10.3, ES10.3, ES10.3, ES10.3, ES10.3)",&
!        eg(1, i), eg(2,i), eg(3,i), eg(4,i), eg(5,i), eg(6,i), eg(7,i), eg(8,i), eg(9,i)
!enddo
        




        END SUBROUTINE Compute_Egf

END SUBROUTINE D_assemble
