Subroutine assts_gl(neo,nei,sequ,fequ,el,gl)
!...subroutine to assemble element matrix el to global matrix gl
	USE DefinitionConstant
    USE GlobalData
    IMPLICIT NONE
    INTEGER, INTENT(IN)			:: neo,nei
    INTEGER, INTENT(IN)			:: sequ(:),fequ(:)
    REAL(KIND=DBL), INTENT(IN)	:: el(:,:)
    REAL(KIND=DBL), INTENT(OUT)	:: gl(:)
    !...local variables
    INTEGER			:: i,j,ki,kj,ig,jg,ij

	do i = 1,neo
       	do ki = 1,NDOF
           	!...global index associated with row
            ig = sequ(i) + ki - 1
            do j = 1,nei
               	do kj = 1,NDOF
               	!...global index associated with column
                    jg = fequ(j) + kj - 1
                    !...only assemble the upper half of gl
                    if (ig.gt.jg) cycle
                    !...transform to one-dimensional array
                    ij = jg*(jg-1)/2 + ig
                    !...assemble to global gl
                    gl(ij) = gl(ij) + el(NDOF*(i-1)+ki,NDOF*(j-1)+kj)
                enddo
            enddo
        enddo
    enddo
End Subroutine assts_gl
!------------------------------------------------------------
Subroutine assts_grhs(neo,sequ,rhs,grhs)
!...subroutine to assemble element matrix el(9,9) to global matrix gl
	USE DefinitionConstant
    USE GlobalData
 	IMPLICIT NONE
   	INTEGER, INTENT(IN)			:: neo
   	INTEGER, INTENT(IN)			:: sequ(:)
   	REAL(KIND=DBL), INTENT(IN)	:: rhs(:)
   	REAL(KIND=DBL), INTENT(OUT)	:: grhs(:)
    !...local variables
    INTEGER			:: i,ki,ig

	do i = 1,neo
       	do ki = 1,NDOF
           	!...global index associated with row
            ig = sequ(i) + ki - 1
            !...assemble to global gl
            grhs(ig) = grhs(ig) + rhs(NDOF*(i-1)+ki)
        enddo
    enddo
End Subroutine assts_grhs