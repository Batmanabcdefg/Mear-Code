MODULE ElementMatrix
USE GlobalData
USE DefinitionConstant
CONTAINS
SUBROUTINE element(region,isrc, ifld, sequ, fequ, ec, eaij, eatij, ebij, edij,         &
                                                        eaji, eatji, ebji, edji)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: isrc, ifld, region, sequ(3), fequ(3)
        REAL(KIND=DBL), INTENT(out)  :: ec(:,:), eaij(:,:), eatij(:,:),          &
                                          ebij(:,:), edij(:,:), eaji(:,:),         &
                                          eatji(:,:), ebji(:,:), edji(:,:)
    ! local variables
    logical                :: debug
    real(kind=DBL) :: temp(3*NDOF,3*NDOF)

    ! node coordinates of source and field element
    REAL(KIND=DBL), SAVE :: ncoor_src(2,3), ncoor_fld(2,3)

    INTEGER :: i, alpha, j, nodex, nodey, beta
    LOGICAL :: eaij_flag, eatij_flag, ebij_flag, edij_flag,  &
               eaji_flag, eatji_flag, ebji_flag, edji_flag, upper_half, lower_half,      &
               edgeij_elem, edgeji_elem, src_interface, fld_interface

    !debug = .true.
    debug = .false.

    ! material constant, shear modulus
    !mu = 0.5d0 * E / (1.d0 + nu)

    if (ELTYPE(elemid(isrc)) == INTERFE) then
            src_interface = .true.
    else
            src_interface = .false.
    end if

    if (ELTYPE(elemid(ifld)) == INTERFE) then
            fld_interface = .true.
    else
            fld_interface = .false.
    end if

    ! initialize matrices
    ec  = 0.0d0
    eaij  = 0.0d0
    eatij = 0.0d0
    ebij  = 0.0d0
    edij  = 0.0d0
    eaji  = 0.0d0
    eatji = 0.0d0
    ebji  = 0.0d0
    edji  = 0.0d0

    ! initialize flags
    eaij_flag  = .false. !.true.
    eatij_flag = .false. !.true.
    ebij_flag  = .false. !.true.
    edij_flag  = .true.
    eaji_flag  = .false. !.true.
    eatji_flag = .false. !.true.
    ebji_flag  = .false. !.true.
    edji_flag  = .true.
    upper_half = .false.
    lower_half = .false.
    edgeij_elem  = .false.  
    edgeji_elem  = .false.  


    nodex = NODE(elemid(ifld))  ! number of nodes in field element
    nodey = NODE(elemid(isrc))  ! number of nodes in source element

    ! to see if the field element is "EDGE_ELEMENT"
    do beta = 1, NDOF     ! over 3 directions 
        do j = 1, nodex  ! over the nodes of field element
            if (node_type(beta,elnode(j,ifld)) == EDGE_NODE) then
                edgeij_elem = .true.
                exit
            end if
        end do
    end do
    !...to see if the source element is "EDGE_ELEMENT"
    do alpha = 1, NDOF     ! over 3 directions 
         do i = 1, nodey  ! over the nodes of field element
            if (node_type(alpha,elnode(i,isrc)) == EDGE_NODE) then
                edgeji_elem = .true.
                exit
            end if
         end do
    end do

    ! the position of element matrices in GLOBAL matrix
    !...modified by Han: maxval(sequ(1:nodey))+3 and maxval(fequ(1:nodex))+3 for not interface
    if (src_interface) then
            if (maxval(sequ(1:nodey))+NDOF < minval(fequ(1:nodex))) upper_half = .true.
    else
              !if (maxval(sequ(1:nodey))< minval(fequ(1:nodex))) upper_half = .true.
            if (maxval(sequ(1:nodey))+NDOF< minval(fequ(1:nodex))) upper_half = .true.
    end if
    if (fld_interface) then
            if (minval(sequ(1:nodey)) > maxval(fequ(1:nodex))+NDOF) lower_half = .true.
    else
              !if (minval(sequ(1:nodey)) > maxval(fequ(1:nodex))) lower_half = .true.
            if (minval(sequ(1:nodey)) > maxval(fequ(1:nodex))+NDOF) lower_half = .true.
    end if


    ! node coordinates of source element
    ncoor_src = node_coor(:, elnode(:,isrc)) 
    ! node coordinates of field element
    ncoor_fld = node_coor(:, elnode(:,ifld))

    ! compute the distance between source element and field element
    !if (isrc /= ifld .and. adjacent(isrc,ifld) <= 0) then
    !  call compute_elem_dist(elem_dist, isrc, ifld, nodey, nodex, &
    !                         ncoor_src, ncoor_fld)
    !end if
         
    ! calculate ec
    IF (isrc == ifld) THEN  
            CALL ek_c(isrc, ncoor_src, ng_ef, xi(:,ng_ef), wi(:,ng_ef), ec)
    END IF

!............................................................................
!  lower_half = .false.
!print*, isrc, equno(1,isrc), equno(2,isrc), equno(3,isrc)
!print*, ifld, equno(1,ifld), equno(2,ifld), equno(3,ifld) 
!............................................................................

    if (lower_half) then
            ! for lower_half case , we still need  to calculate some matrices
            ! for calculating the right-hand side(load vector).
            ! determine which matrix need to calculate in this case
            !...even if this case happens, we still have make enough assembling for global
            !stiffness matrix since lower_half is the "else" of upper_half later, and then
            !all matrices ji... will be constructed, and in this case ji-matrix is upper half -> assembly to K
        !...debug
        if (debug) then
                print*,'lower_half'
        endif
        
        DO alpha = 1,NDOF  ! over NDOF directions
             SELECT CASE (elldid(alpha, isrc))
                !--------------------
                CASE (BDISP)
                  ! using "disp. eqn"
                     DO beta = 1,NDOF ! over NDOF directions
                        SELECT CASE (elldid(beta, ifld))
                          !----------
                          CASE (BDISP) 
                              if (eatij_flag) then
                                CALL ek_at(region,isrc, ifld, ncoor_src, ncoor_fld, eatij,.false.)
                                eatij_flag = .false.
                                !...debug
                                if (debug) then
                                    print*,'eatij from 1-a'
                                    do i =1,9
                                        print'(9(f12.6,2x))',(eatij(i,j),j=1,9)
                                    enddo
                                endif
                              end if
                          !----------
                          case (BTRAC)
                               if (ebij_flag) then
                                 CALL ek_b(region,isrc, ifld, ncoor_src, ncoor_fld, ebij,.false.)
                                 ebij_flag = .false.
                                 !...debug
                                 if (debug) then
                                    print*,'ebij from 1-b'
                                    do i =1,9
                                      print'(9(f12.6,2x))',(ebij(i,j),j=1,9)
                                    enddo
                                 endif
                               end if
                               if (edgeij_elem .and. eatij_flag) then
                                 CALL ek_at(region,isrc, ifld, ncoor_src, ncoor_fld, eatij,.false.)
                                 eatij_flag = .false.
                                 !...debug
                                 if (debug) then
                                    print*,'eatij from 1-b'
                                    do i =1,9
                                       print'(9(f12.6,2x))',(eatij(i,j),j=1,9)
                                    enddo
                                  endif
                                end if
                          !----------
                          CASE(BTRFREE)
                             ! for "edge_elem", we need EAT
                               if (edgeij_elem .and. eatij_flag) then
                                 CALL ek_at(region,isrc, ifld, ncoor_src, ncoor_fld, eatij,.false.)
                                 eatij_flag = .false.
                                 !...debug
                                 if (debug) then
                                    print*,'eatij from 1-c'
                                    do i =1,9
                                       print'(9(f12.6,2x))',(eatij(i,j),j=1,9)
                                    enddo
                                  endif
                                end if
                          !----------
                          case (CTRAC, CTRFREE)
                            ! do nothing
                          !----------
                          case (NOLOAD) ! INTERFACE element
                            ! do nothing
                       END SELECT
                     END DO   ! loop beta
                !--------------------
                CASE (BTRAC, CTRAC, BTRFREE, CTRFREE)
                    ! using "traction" equation
                    DO beta = 1,NDOF ! over 3 directions
                       SELECT CASE (elldid(beta, ifld))
                         !----------
                         CASE (BDISP)
                             if (edij_flag) then
                               CALL ek_d(region,isrc, ifld, ncoor_src, ncoor_fld, edij)
                               edij_flag = .false.
                               !...debug
                               if (debug) then
                                  print*,'edij from 2-a'
                                  do i =1,9
                                    print'(9(f12.6,2x))',(edij(i,j),j=1,9)
                                  enddo
                               endif
                             end if
                         !----------
                         case (BTRAC)
                             if (eaij_flag) then
                               CALL ek_at(region,ifld, isrc, ncoor_fld, ncoor_src, temp,.false.)
                               !...ea is transpose of eat
                               eaij = TRANSPOSE(temp)
                               eaij_flag = .false.
                               !...debug
                               if (debug) then
                                    print*,'eaij from 2-b'
                                    do i =1,9
                                       print'(9(f12.6,2x))',(eaij(i,j),j=1,9)
                                    enddo
                               endif
                             end if
                             ! for edge_elem 
                             if (edgeij_elem .and. edij_flag) then
                               CALL ek_d(region,isrc, ifld, ncoor_src, ncoor_fld, edij)
                               edij_flag = .false.
                               !...debug
                               if (debug) then
                                    print*,'edij from 2-b'
                                    do i =1,9
                                       print'(9(f12.6,2x))',(edij(i,j),j=1,9)
                                    enddo
                               endif
                             end if
                         !----------
                         CASE (BTRFREE)
                             ! for edge_elem 
                             if (edgeij_elem .and. edij_flag) then
                               CALL ek_d(region,isrc, ifld, ncoor_src, ncoor_fld, edij)
                               edij_flag = .false.
                               !...debug
                               if (debug) then
                                  print*,'edij from 2-c'
                                  do i =1,9
                                     print'(9(f12.6,2x))',(edij(i,j),j=1,9)
                                  enddo
                               endif
                             end if
                         !----------
                         case(CTRAC, CTRFREE)
                             ! do nothing
                         !----------
                         case (NOLOAD) ! INTERFACE 
                         ! do nothing
                       END SELECT
                    END DO ! loop beta
                !--------------------
                !...Lower-half: case 3
                case (NOLOAD) ! source element is on INTERFACE surface
                    DO beta = 1,NDOF ! over 3 directions
                       SELECT CASE (elldid(beta, ifld))
                         !----------
                         !...3a
                         CASE (BDISP)  ! displacement prescribed
                             ! traction equation
                             if (edij_flag) then
                                CALL ek_d(region,isrc, ifld, ncoor_src, ncoor_fld, edij)
                                edij_flag = .false.
                                !...debug
                                if (debug) then
                                    print*,'edij from 3-a, traction eq'
                                    do i =1,9
                                       print'(9(f12.6,2x))',(edij(i,j),j=1,9)
                                    enddo
                                endif
                             end if
                             ! disp. equaiton
                             if (eatij_flag) then
                                CALL ek_at(region,isrc, ifld, ncoor_src, ncoor_fld, eatij,.false.)
                                eatij_flag = .false.
                                !...debug
                                if (debug) then
                                    print*,'eatij from 3-a, disp eq'
                                    do i =1,9
                                        print'(9(f12.6,2x))',(eatij(i,j),j=1,9)
                                    enddo
                                endif
                             end if
                         !----------
                         !...3b
                         case (BTRAC, CTRAC)  ! traction prescribed
                             if (edgeij_elem) then  ! EDGE_NODE element
                                ! disp. equation
                                if (eatij_flag) then
                                   CALL ek_at(region,isrc, ifld, ncoor_src, ncoor_fld, eatij,.false.)
                                   eatij_flag = .false.
                                   !...debug
                                   if (debug) then
                                      print*,'eatij from 3-b, field elem has edge node, disp eq'
                                      do i =1,9
                                          print'(9(f12.6,2x))',(eatij(i,j),j=1,9)
                                      enddo
                                   endif
                                end if
                                ! traction equation
                                if (edij_flag) then
                                   CALL ek_d(region,isrc, ifld, ncoor_src, ncoor_fld, edij)
                                   edij_flag = .false.
                                   !...debug
                                   if (debug) then        
                                      print*,'edij from 3-b, field elem has edge node, traction eq'
                                      do i =1,9
                                          print'(9(f12.6,2x))',(edij(i,j),j=1,9)
                                      enddo
                                   endif
                                end if
                             end if
                             ! for other nodes of the element
                             ! disp. equation
                             if (eaij_flag) then
                                CALL ek_at(region,ifld, isrc, ncoor_fld, ncoor_src, temp,.false.)
                                eaij = TRANSPOSE(temp)
                                eaij_flag = .false.
                                !...debug
                                if (debug) then
                                    print*,'eaij from 3-b, disp eq'
                                    do i =1,9
                                       print'(9(f12.6,2x))',(eaij(i,j),j=1,9)
                                    enddo
                                endif
                             end if
                             ! traction equation
                             if (ebij_flag) then
                                CALL ek_b(region,isrc, ifld, ncoor_src, ncoor_fld, ebij,.false.)
                                ebij_flag = .false.
                                !...debug
                                if (debug) then
                                    print*,'ebij from 3-b, traction eq'
                                    do i =1,9
                                       print'(9(f12.6,2x))',(ebij(i,j),j=1,9)
                                    enddo
                                endif
                             end if
                         !----------
                         !...3c
                         case (BTRFREE, CTRFREE)
                             if (edgeij_elem) then  ! EDGE_NODE element
                                ! disp. equation
                                if (eatij_flag) then
                                   CALL ek_at(region,isrc, ifld, ncoor_src, ncoor_fld, eatij,.false.)
                                   eatij_flag = .false.
                                   !...debug
                                   if (debug) then
                                      print*,'eatij from 3-c, field elem has edge nodes, disp eq'
                                      do i =1,9
                                          print'(9(f12.6,2x))',(eatij(i,j),j=1,9)
                                      enddo
                                   endif
                                end if
                                ! traction equation
                                if (edij_flag) then
                                   CALL ek_d(region,isrc, ifld, ncoor_src, ncoor_fld, edij)
                                   edij_flag = .false.
                                   !...debug
                                   if (debug) then
                                      print*,'edij from 3-c, field elem has edge nodes, traction eq'
                                      do i =1,9
                                          print'(9(f12.6,2x))',(edij(i,j),j=1,9)
                                      enddo
                                   endif
                                end if
                             end if
                         END SELECT
                    END DO   !...of beta
             END SELECT !...of elldid(alpha,isrc)
        END DO !...of alpha
    else 
        ! for not lower_half case
        ! calculate matrices for element matrix and load vector
        !...this is the case that we want to assembly to upper half of K
        !...however, for lowerhalf, then "ELSE" of upperhalf will include, then ji-matrix is computed
        !...for assembly to K
        !...debug
        if (debug) then
           print*,'not lower half'
        endif
        
        DO alpha = 1,NDOF  ! over 3 directions
           SELECT CASE (elldid(alpha, isrc))
              !--------------------
              !...not lower-half, case 1
              case (BDISP)
                  ! using "disp. eqn"
                  if (eatij_flag) then
                     CALL ek_at(region,isrc, ifld, ncoor_src, ncoor_fld, eatij,.false.)
                     eatij_flag = .false.
                     !...debug
                     if (debug) then
                        print*,'eatij from 1'
                        do i =1,9
                           print'(9(f12.6,2x))',(eatij(i,j),j=1,9)
                        enddo
                  endif
                          end if
                
                  DO beta = 1,NDOF ! over 3 directions
                     SELECT CASE (elldid(beta, ifld))
                        !----------
                        !...1a
                        CASE (BDISP, BTRAC)
                            if (ebij_flag) then
                               CALL ek_b(region,isrc, ifld, ncoor_src, ncoor_fld, ebij,.false.)
                               ebij_flag = .false.
                               !...debug
                               if (debug) then
                                  print*,'ebij from 1-a'
                                  do i =1,9
                                       print'(9(f12.6,2x))',(ebij(i,j),j=1,9)
                                  enddo
                               endif
                            end if
                        !----------
                        !...1b
                        case (BTRFREE)
                          ! do nothing
                        !----------
                        !...1c
                        case (CTRAC, CTRFREE)
                          ! do nothing
                        !----------
                        !...1d
                        case (NOLOAD)  ! field element is INTERFACE element
                             if (ebij_flag) then
                                CALL ek_b(region,isrc, ifld, ncoor_src, ncoor_fld, ebij,.false.)
                                ebij_flag = .false.
                                !...debug
                                if (debug) then
                                    print*,'ebij from 1-d'
                                    do i =1,9
                                        print'(9(f12.6,2x))',(ebij(i,j),j=1,9)
                                    enddo
                                endif
                             end if
                     END SELECT
                  END DO! of beta
              !--------------------
              !...not lower-half, case 2
              CASE (BTRAC, BTRFREE, CTRAC, CTRFREE)
                  ! using traction equation
                  if (edij_flag) then
                     CALL ek_d(region,isrc, ifld, ncoor_src, ncoor_fld, edij)
                     edij_flag = .false.
                     !...debug
                     if (debug) then
!............................................................
                       if (lower_half .eq. .false.) then
                        print*,'edij from 2:'
                        do i =1,9
                             print'(9(f12.6,2x))',(edij(i,j),j=1,9)
                        enddo
                       endif
!............................................................
                     endif
                  end if
                  DO beta = 1,NDOF ! over 3 directions
                     SELECT CASE (elldid(beta, ifld))
                        !----------
                        !...2a
                        CASE (BDISP, BTRAC)
                             if (eaij_flag) then
                                CALL ek_at(region,ifld, isrc, ncoor_fld, ncoor_src, temp,.false.)
                                eaij = TRANSPOSE(temp)
                                eaij_flag = .false.
                                !...debug
                                if (debug) then
                                    print*,'eaij from 2-a'
                                    do i =1,9
                                         print'(9(f12.6,2x))',(eaij(i,j),j=1,9)
                                    enddo
                                endif
                             end if
                        !----------
                        !...2b
                        case (CTRAC, CTRFREE, BTRFREE)
                        !----------
                        !...2c
                        case (NOLOAD)  ! field element is INTERFACE eleement
                             if (eaij_flag) then
                                CALL ek_at(region,ifld, isrc, ncoor_fld, ncoor_src, temp,.false.)
                                eaij = TRANSPOSE(temp)
                                eaij_flag = .false.
                                !...debug
                                if (debug) then
                                    print*,'eaij from 2-c'
                                    do i =1,9
                                       print'(9(f12.6,2x))',(eaij(i,j),j=1,9)
                                    enddo
                                endif
                             end if
                     END SELECT
                  END DO! loop beta
              !--------------------
              !...not lower-half, case 3
              case (NOLOAD)  ! source element is INTERFACE element
                  ! traction equation
                  if (eaij_flag) then
                      !...what about the case ea(gamma_i,gamma_c)? because this case will never happen!
                      CALL ek_at(region,ifld, isrc, ncoor_fld, ncoor_src, temp,.false.)
                      eaij = TRANSPOSE(temp)
                      eaij_flag = .false.
                      if (debug) then
                         print*,'eaij from 3, traction'
                         do i =1,9
                              print'(9(f12.6,2x))',(eaij(i,j),j=1,9)
                         enddo
                      endif
                  end if
                  if (edij_flag) then
                     CALL ek_d(region,isrc, ifld, ncoor_src, ncoor_fld, edij)
                     edij_flag = .false.
                     if (debug) then
                        print*,'edij from 3, traction'
                        do i =1,9
                             print'(9(f12.6,2x))',(edij(i,j),j=1,9)
                        enddo
                     endif
                  end if

                  ! disp. equaiton
                  if (eatij_flag) then
                     CALL ek_at(region,isrc, ifld, ncoor_src, ncoor_fld, eatij,.false.)
                     eatij_flag = .false.
                     if (debug) then
                        print*,'eatij from 3, disp equation'
                        do i =1,9
                             print'(9(f12.6,2x))',(eatij(i,j),j=1,9)
                        enddo
                     endif
                  end if
                  if (ebij_flag) then
                     CALL ek_b(region,isrc, ifld, ncoor_src, ncoor_fld, ebij,.false.)
                     ebij_flag = .false.
                     if (debug) then
                        print*,'ebij from 3, disp equation'
                        do i =1,9
                             print'(9(f12.6,2x))',(ebij(i,j),j=1,9)
                        enddo
                     endif
                  end if
           END SELECT
        END DO! loop alpha
    end if  ! if (lower_half)
        
    if (isrc == ifld) return
        !...yes, correct, since if isrc=ifld, we don't need to repeat the same thing for ji-matrix
    
    if (upper_half) then
        !...this is also included of "not lowerhalf" above, so all ij-matrix were computed
        ! calculate matrices for load vector
        !...debug
        if (debug) then
           print*,'upper half'
        endif
        
        DO beta = 1,NDOF  ! over 3 directions
           SELECT CASE (elldid(beta, ifld))
               !--------------------
               CASE (BDISP)
                    ! using "disp. eqn"
                    DO alpha = 1,NDOF ! over 3 directions
                       SELECT CASE (elldid(alpha, isrc))
                         !----------
                         CASE (BDISP) 
                             if (eatji_flag) then
                                if (.not. eaij_flag) then
                                   eatji = transpose(eaij)
                                else
                                   CALL ek_at(region,ifld, isrc, ncoor_fld, ncoor_src, eatji,.false.)
                                end if
                                eatji_flag = .false.
                                !...debug
                                if (debug) then
                                    print*,'eatji from 1-a'
                                    do i =1,9
                                         print'(9(f12.6,2x))',(eatji(i,j),j=1,9)
                                    enddo
                                endif
                             end if
                         !----------
                         case (BTRAC)
                              if (ebji_flag) then
                                 if (.not. ebij_flag) then
                                    ebji = transpose(ebij)
                                 else
                                    CALL ek_b(region,ifld, isrc, ncoor_fld, ncoor_src, ebji,.false.)
                                 end if
                                 ebji_flag = .false.
                                 !...debug
                                 if (debug) then
                                    print*,'ebji from 1-b'
                                    do i =1,9
                                         print'(9(f12.6,2x))',(ebji(i,j),j=1,9)
                                    enddo
                                 endif
                              end if
                              if (edgeji_elem .and. eatji_flag) then 
                                 if (.not. eaij_flag) then
                                     eatji = transpose(eaij)
                                 else
                                     CALL ek_at(region,ifld, isrc, ncoor_fld, ncoor_src, eatji,.false.)
                                 end if
                                 eatji_flag = .false.
                                 !...debug
                                 if (debug) then
                                    print*,'eatji from 1-b'
                                    do i =1,9
                                         print'(9(f12.6,2x))',(eatji(i,j),j=1,9)
                                    enddo
                                 endif
                              end if
                         !----------
                         CASE(BTRFREE)
                             if (edgeji_elem .and. eatji_flag) then 
                                if (.not. eaij_flag) then
                                   eatji = transpose(eaij)
                                else
                                CALL ek_at(region,ifld, isrc, ncoor_fld, ncoor_src, eatji,.false.)
                                end if
                                eatji_flag = .false.
                                !...debug
                                if (debug) then
                                    print*,'eatji from 1-c'
                                    do i =1,9
                                         print'(9(f12.6,2x))',(eatji(i,j),j=1,9)
                                    enddo
                                endif
                             end if
                         !----------
                         case (CTRAC, CTRFREE)
                              ! do nothing
                         !----------
                         case (NOLOAD) ! INTERFACE 
                              ! do nothing
                       END SELECT
                    END DO! loop alpha
               !--------------------
               CASE (BTRAC, CTRAC, BTRFREE, CTRFREE)
                   ! using "traction eqn".
                   DO alpha = 1,NDOF ! over 3 directions
                      SELECT CASE (elldid(alpha,isrc))
                        !----------
                        CASE (BDISP)
                            if (edji_flag) then
                               if (.not. edij_flag) then
                                  edji = transpose(edij) 
                               else
                                  CALL ek_d(region,ifld, isrc, ncoor_fld, ncoor_src, edji)
                               end if
                               edji_flag = .false.
                               !...debug
                               if (debug) then
                                  print*,'edji from 2-a'
                                  do i =1,9
                                       print'(9(f12.6,2x))',(edji(i,j),j=1,9)
                                  enddo
                               endif
                            end if
                        !----------
                        case (BTRAC)
                            if (eaji_flag) then
                               if (.not. eatij_flag) then
                                  eaji = transpose(eatij) 
                               else
                                  CALL ek_at(region,isrc, ifld, ncoor_src, ncoor_fld, temp,.false.)
                                  eaji = TRANSPOSE(temp)
                               end if
                               eaji_flag = .false.
                               !...debug
                               if (debug) then
                                    print*,'eaji from 2-b'
                                    do i =1,9
                                         print'(9(f12.6,2x))',(eaji(i,j),j=1,9)
                                    enddo
                               endif
                            end if
                            if (edgeji_elem .and. edji_flag) then
                               if (.not. edij_flag) then
                                  edji = transpose(edij) 
                               else
                                  CALL ek_d(region,ifld, isrc, ncoor_fld, ncoor_src, edji)
                               end if
                               edji_flag = .false.
                               !...debug
                               if (debug) then
                                    print*,'edji from 2-b'
                                    do i =1,9
                                         print'(9(f12.6,2x))',(edji(i,j),j=1,9)
                                    enddo
                               endif
                            end if
                        !----------
                        CASE (BTRFREE)
                             if (edgeji_elem .and. edji_flag) then
                                if (.not. edij_flag) then
                                   edji = transpose(edij) 
                                else
                                   CALL ek_d(region,ifld, isrc, ncoor_fld, ncoor_src, edji)
                                end if
                                edji_flag = .false.
                                !...debug
                                if (debug) then
                                    print*,'edji from 2-c'
                                    do i =1,9
                                         print'(9(f12.6,2x))',(edji(i,j),j=1,9)
                                    enddo
                                endif
                             end if
                        !----------
                        case (CTRAC, CTRFREE)
                             ! do nothing
                      END SELECT
                   END DO! loop alpha
               !--------------------
               !...upper half, case 3
               case (NOLOAD)  ! source element is INTERFACE element
                    DO alpha = 1,NDOF ! over 3 directions
                       SELECT CASE (elldid(alpha, isrc))
                         !----------
                         CASE (BDISP)  ! displacement prescribed on field element
                             ! traction eqn.
                             if (edji_flag) then
                                if (.not. edij_flag) then
                                   edji = transpose(edij)
                                else
                                   CALL ek_d(region,ifld, isrc, ncoor_fld, ncoor_src, edji)
                                end if
                                edji_flag = .false.
                               if (debug) then
                                   print*,'edji from 3a, traction eq'
                                   do i =1,9
                                         print'(9(f12.6,2x))',(edji(i,j),j=1,9)
                                   enddo
                               endif
                             end if
                             ! displacement eqn.
                             if (eatji_flag) then
                                if (.not. eaij_flag) then
                                   eatji = transpose(eaij)
                                else
                                   CALL ek_at(region,ifld, isrc, ncoor_fld, ncoor_src, eatji,.false.)
                                end if
                                eatji_flag = .false.
                                if (debug) then
                                    print*,'eatji from 3a, disp eq'
                                    do i =1,9
                                         print'(9(f12.6,2x))',(eatji(i,j),j=1,9)
                                    enddo
                                endif
                             end if
                         !----------
                         case (BTRAC, CTRAC)  ! traction prescribed on field element
                              if (edgeji_elem) then ! EDGE_NODE element
                                  ! disp. eqn.
                                  if (eatji_flag) then 
                                     if (.not. eaij_flag) then
                                        eatji = transpose(eaij)
                                     else
                                        CALL ek_at(region,ifld, isrc, ncoor_fld, ncoor_src, eatji,.false.)
                                     end if
                                     eatji_flag = .false.
                                     if (debug) then
                                         print*,'eatji from 3b, source elem has edge nodes, disp eq'
                                         do i =1,9
                                               print'(9(f12.6,2x))',(eatji(i,j),j=1,9)
                                         enddo
                                     endif
                                  end if
                                  ! traction eqn.
                                  if (edji_flag) then
                                     if (.not. edij_flag) then
                                        edji = transpose(edij)
                                     else
                                        call ek_d(region,ifld, isrc, ncoor_fld, ncoor_src, edji)
                                     end if
                                     edji_flag = .false.
                                     if (debug) then
                                         print*,'edji from 3b, source elem has edge nodes, traction eq'
                                         do i =1,9
                                              print'(9(f12.6,2x))',(edji(i,j),j=1,9)
                                         enddo
                                     endif
                                  end if
                              end if
                              ! for other nodes of the element
                              ! disp. eqn.
                              if (eaji_flag) then 
                                  if (.not. eatij_flag) then
                                      eaji = transpose(eatij)
                                  else
                                      CALL ek_at(region,isrc, ifld, ncoor_src, ncoor_fld, temp,.false.)
                                      eaji = transpose(temp)
                                  end if
                                  eaji_flag = .false.
                                  if (debug) then
                                      print*,'eaji from 3b, disp eq'
                                      do i =1,9
                                           print'(9(f12.6,2x))',(eaji(i,j),j=1,9)
                                      enddo
                                  endif
                              end if
                              ! traction eqn.
                              if (ebji_flag) then
                                  if (.not. ebij_flag) then
                                      ebji = transpose(ebij)
                                  else
                                      CALL ek_b(region,ifld, isrc, ncoor_fld, ncoor_src, ebji,.false.)
                                  end if
                                  ebji_flag = .false.
                                  if (debug) then
                                      print*,'ebji from 3b, traction eq'
                                      do i =1,9
                                           print'(9(f12.6,2x))',(ebji(i,j),j=1,9)
                                      enddo
                                  endif
                              end if
                         !----------
                         case (BTRFREE, CTRFREE)
                              if (edgeji_elem) then ! EDGE_NODE element
                                 ! disp. eqn.
                                 if (eatji_flag) then 
                                    if (.not. eaij_flag) then
                                        eatji = transpose(eaij)
                                    else
                                        CALL ek_at(region,ifld, isrc, ncoor_fld, ncoor_src, eatji,.false.)
                                    end if
                                    eatji_flag = .false.
                                    if (debug) then
                                        print*,'eatji from 3c, source elem has edge nodes, disp eq'
                                        do i =1,9
                                             print'(9(f12.6,2x))',(eatji(i,j),j=1,9)
                                        enddo
                                    endif
                                 end if
                                 ! traction eqn.
                                 if (edji_flag) then
                                     if (.not. edij_flag) then
                                         edji = transpose(edij)
                                     else
                                         call ek_d(region,ifld, isrc, ncoor_fld, ncoor_src, edji)
                                     end if
                                     edji_flag = .false.
                                     if (debug) then
                                         print*,'edji from 3c, source elem has edge nodes, trac eq'
                                         do i =1,9
                                              print'(9(f12.6,2x))',(edji(i,j),j=1,9)
                                         enddo
                                     endif
                                 end if
                              end if
                       end SELECT
                 end DO !...of alpha
            END SELECT
        END DO! loop beta
    else   !...not upperhalf case
           !...this is also included lowerhalf case (above), thus ji-matrix is computed and
        !...this ji-matrix is now upperhalf, and necessary to be assembled to global K
        !to assure that all upper half of global K is fullfilled
        ! determine which matrix need to calculate
        !...debug
        if (debug) then
                print*,'not upper half'
        endif
        
        DO beta = 1,NDOF  ! over 3 directions
            SELECT CASE (elldid(beta, ifld))
              !--------------------
              CASE (BDISP)
                  ! using "disp. eqn"
                  if (eatji_flag) then
                     if (.not. eaij_flag) then
                         eatji = transpose(eaij)
                     else
                         CALL ek_at(region,ifld, isrc, ncoor_fld, ncoor_src, eatji,.false.)
                     end if
                     eatji_flag = .false.
                     !...debug
                     if (debug) then
                         print*,'eatji from 1'
                         do i =1,9
                              print'(9(f12.6,2x))',(eatji(i,j),j=1,9)
                         enddo
                     endif
                  end if
                  DO alpha = 1,NDOF ! over 3 directions
                     SELECT CASE (elldid(alpha, isrc))
                        !----------
                        CASE (BDISP, BTRAC)
                             if (ebji_flag) then
                                if (.not. ebij_flag) then
                                    ebji = transpose(ebij)
                                else
                                    CALL ek_b(region,ifld, isrc, ncoor_fld, ncoor_src, ebji,.false.)
                                end if
                                ebji_flag = .false.
                                !...debug
                                if (debug) then
                                    print*,'ebji from 1-a'
                                    do i =1,9
                                         print'(9(f12.6,2x))',(ebji(i,j),j=1,9)
                                    enddo
                                endif
                             end if
                        !----------
                        !...1a
                        CASE(BTRFREE, CTRAC, CTRFREE)
                            ! do nothing
                        !----------
                        !...1c
                        case (NOLOAD)  ! INTERFACE element
                            if (ebji_flag) then
                               if (.not. ebij_flag) then
                                   ebji = transpose(ebij)
                               else
                                   CALL ek_b(region,ifld, isrc, ncoor_fld, ncoor_src, ebji,.false.)
                               end if
                               ebji_flag = .false.
                               !...debug
                               if (debug) then
                                   print*,'ebji from 1c'
                                   do i =1,9
                                        print'(9(f12.6,2x))',(ebji(i,j),j=1,9)
                                   enddo
                               endif
                            end if
                     END SELECT
                  END DO! loop alpha
              !--------------------
              CASE (BTRAC, CTRAC, BTRFREE, CTRFREE)
                   ! using traction equation
                   if (edji_flag) then
                      if (.not. edij_flag) then
                          edji = transpose(edij)
                      else
                          CALL ek_d(region,ifld, isrc, ncoor_fld, ncoor_src, edji)
                      end if
                      edji_flag = .false.
                      !...debug
                      if (debug) then
!...............................................
                        if(lower_half .eq. .false.) then
                         print*,'edji from 2'
                         do i =1,9
                              print'(9(f12.6,2x))',(edji(i,j),j=1,9)
                         enddo
                        endif
!......................................................
                      endif
                   end if
                   DO alpha = 1,NDOF ! over 3 directions
                      SELECT CASE (elldid(alpha, isrc))
                         !----------
                         CASE (BDISP, BTRAC)
                             if (eaji_flag) then
                                if (.not. eatij_flag) then
                                    eaji = transpose(eatij)
                                else
                                    CALL ek_at(region,isrc, ifld, ncoor_src, ncoor_fld, temp,.false.)
                                    eaji = TRANSPOSE(temp)
                                end if
                                eaji_flag = .false.
                                !...debug
                                if (debug) then
                                    print*,'eaji from 2-a'
                                    do i =1,9
                                         print'(9(f12.6,2x))',(eaji(i,j),j=1,9)
                                    enddo
                                endif
                             end if
                         !----------
                         !...2b
                         CASE (BTRFREE, CTRAC, CTRFREE)
                             ! do nothing
                         !----------
                         !...2c
                         case (NOLOAD)  ! filed element is INTERFACE element
                             if (eaji_flag) then
                                 if (.not. eatij_flag) then
                                     eaji = transpose(eatij)
                                 else
                                     CALL ek_at(region,isrc, ifld, ncoor_src, ncoor_fld, temp,.false.)
                                     eaji = TRANSPOSE(temp)
                                 end if
                                 eaji_flag = .false.
                                 !...debug
                                 if (debug) then
                                     print*,'eaji from 2c'
                                     do i =1,9
                                          print'(9(f12.6,2x))',(eaji(i,j),j=1,9)
                                     enddo
                                 endif
                             end if
                      END SELECT
                   END DO! loop alpha
              !--------------------
              !...case 3
              case (NOLOAD)  ! source element is INTERFACE element
                   ! traction equation 
                   if (eaji_flag) then
                       if (.not. eatij_flag) then
                           eaji = transpose(eatij)
                       else
                           CALL ek_at(region,isrc, ifld, ncoor_src, ncoor_fld, temp,.false.)
                           eaji = TRANSPOSE(temp)
                       end if
                       eaji_flag = .false.
                       !...debug
                       if (debug) then
                           print*,'eaji from 3, traction eq'
                           do i =1,9
                                print'(9(f12.6,2x))',(eaji(i,j),j=1,9)
                           enddo
                       endif
                   end if
                   if (edji_flag) then
                       if (.not. edij_flag) then
                           edji = transpose(edij)
                       else
                           CALL ek_d(region,ifld, isrc, ncoor_fld, ncoor_src, edji)
                       end if
                       edji_flag = .false.
                       !...debug
                       if (debug) then
                           print*,'edji from 3, traction eq'
                           do i =1,9
                                print'(9(f12.6,2x))',(edji(i,j),j=1,9)
                           enddo
                       endif
                   end if

                   ! disp. equation
                   if (eatji_flag) then
                       if (.not. eaij_flag) then
                           eatji = transpose(eaij)
                       else
                           CALL ek_at(region,ifld, isrc, ncoor_fld, ncoor_src, eatji,.false.)
                       end if
                       eatji_flag = .false.
                       !...debug
                       if (debug) then
                            print*,'eatji from 3, disp eq'
                            do i =1,9
                                 print'(9(f12.6,2x))',(eatji(i,j),j=1,9)
                            enddo
                       endif
                   end if
                   if (ebji_flag) then
                       if (.not. ebij_flag) then
                           ebji = transpose(ebij)
                       else
                           CALL ek_b(region,ifld, isrc, ncoor_fld, ncoor_src, ebji,.false.)
                       end if
                       ebji_flag = .false.
                       !...debug
                       if (debug) then
                           print*,'ebji from 3, disp eq'
                           do i =1,9
                                print'(9(f12.6,2x))',(ebji(i,j),j=1,9)
                           enddo
                       endif
                   end if
            END SELECT
         END DO   !...loop beta
    end if   !...of upperhalf

END SUBROUTINE element

INCLUDE "ekc_2.f90"

!INCLUDE "ek_d.f90"
!INCLUDE "ek_d1.f90"
!INCLUDE "ek_d2.f90"
!...10/10/2010
INCLUDE "ekd_2.f90"

!INCLUDE "ek_b.f90"
INCLUDE "ek_b2.f90"

!INCLUDE "ek_at.f90"
!INCLUDE "ek_at1.f90"
INCLUDE "ek_at2.f90"

!...for T-stress calculation
!include "eatts.f90"
include "eatts2.f90"
include "eattip.f90"
!include "ebts.f90"
include "ebts2.f90"
include "ebtip.f90"
include "eat3.f90"
include "eb3.f90"

!...for rigid displacement calculation
include "eat_rigid.f90"
include "eb_rigid.f90"

END MODULE ElementMatrix
