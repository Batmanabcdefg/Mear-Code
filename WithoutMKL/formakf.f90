Subroutine FormAKF
        USE DefinitionConstant
    USE GlobalData
    USE ElementMatrix
    IMPLICIT NONE
    INTERFACE
        Subroutine D_assemble(region,isrc,ifld,sequ,fequ,node_src,node_fld,ec,ea,&
                                                  eat,eb,ed,gk,af,af_traction,af_pressure)
                    INTEGER, INTENT(IN)                        :: region,isrc,ifld,sequ(:),fequ(:),node_src,node_fld
                    REAL(SELECTED_REAL_KIND(15)), INTENT(IN)        :: ec(:,:),ea(:,:),eat(:,:),eb(:,:),ed(:,:)
            REAL(SELECTED_REAL_KIND(15)), INTENT(OUT)        :: gk(:),af(:),af_traction(:),af_pressure(:,:)

        End Subroutine D_assemble
    END INTERFACE
        !--------------------------------------------------------------------------------
        !...local variables: all local variables are size-declared by using global NDOF
    REAL(KIND=DBL)        :: ec(3*NDOF,3*NDOF),eaij(3*NDOF,3*NDOF),eatij(3*NDOF,3*NDOF),&
                                                       ebij(3*NDOF,3*NDOF),edij(3*NDOF,3*NDOF),eaji(3*NDOF,3*NDOF),&
                               eatji(3*NDOF,3*NDOF),ebji(3*NDOF,3*NDOF),edji(3*NDOF,3*NDOF)
        !--------------------------------------------------------------------------------
    INTEGER        :: isrc,jfld,IGsrc,JGfld,k,k0, ck
    !LOGICAL        :: debug = .true.
    LOGICAL        :: debug = .false.
    !--------------------------------------------------------------------------------

    k0 = 0
    do k = 1,total_region
              !...compute kernels CI1,CI2_table (and GI,UI)

        call Kernel1(k)

                !...loop over all elements of region k
        do isrc = 1,total_region_elem(k)
                   IGsrc = isrc + k0   !...global source element number
            do jfld = 1,isrc
                       JGfld = jfld + k0   !...global field element number
                !...debug
                            if (debug) then
                                    print*,'elements (user #): ',elem_sys2user(IGsrc),elem_sys2user(JGfld)
                                print*,'elements (system #): ',IGsrc,JGfld
                            endif
                !...compute element matrix and load vector
                call element(k,IGsrc,JGfld,equno(:,IGsrc),equno(:,JGfld),ec,eaij,eatij,ebij,edij,&
                                          eaji,eatji,ebji,edji)
                !...assemble eg and ef to global matrix AK and F: all ij-matrix
                call D_assemble(k,IGsrc,JGfld,equno(:,IGsrc),equno(:,JGfld),NODE(elemid(IGsrc)),&
                                                  NODE(elemid(JGfld)),ec,eaij,eatij,ebij,edij,dAK,dF,dF_traction,dF_pressure)
                
                if (IGsrc.ne.JGfld) then
                           !...assemble all ji-matrix
                    call D_assemble(k,JGfld,IGsrc,equno(:,JGfld),equno(:,IGsrc),NODE(elemid(JGfld)),&
                                                  NODE(elemid(IGsrc)),ec,eaji,eatji,ebji,edji,dAK,dF,dF_traction,dF_pressure)
                endif
            enddo   !...of loop over field elements
        enddo   !...of loop over source elements
        dAK_traction=dAK
        
        do ck=1,total_cracks
           dAk_pressure(ck,:)=dAK
        enddo
        !...update number of element if the problem has more than 1 region
        !...this algorith just work in the case of element # are in right order in each region
        k0 = k0 + total_region_elem(k)
     enddo   !...of loop over regions



!.....used for debugging....3/6/12...aje567
!     print*,'-----------------------------------------------'
!     print*,'The total right hand side is: '
!     print*, dF
!     print*,'-----------------------------------------------'
!     print*,'The right hand side due to tractions is: '
!     print*, dF_traction
!     print*,'-----------------------------------------------'
!     print*,'The right hand side due to pressure (from formakf) is: '
!     print*, dF_pressure
!     print*,'-----------------------------------------------'


End Subroutine FormAKF
