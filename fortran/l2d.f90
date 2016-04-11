!
!
 SUBROUTINE L2D(d,nb,l,pl)
! **********************************************************************
! *                                                                    *
! * Returns elements of D and N matrices from L element and equil pop  *
! *                                                                    *
! **********************************************************************
! *                                                                    *
! * Author:     Boris Fackovec                                         *
! * Date:       13/02/2016                                             *
! * Version:    1.0                                                    *
! *                                                                    *
! **********************************************************************
!
   IMPLICIT NONE
   DOUBLE PRECISION, INTENT(IN):: l, pl
   DOUBLE PRECISION, INTENT(OUT) :: d
   LOGICAL, INTENT(OUT) :: nb
!
! -------------------------------------------------------------------
!
!
! -------------------------------------------------------------------
!
   IF (l.GT.pl) THEN
      nb = .TRUE.
      CALL LogDiffExp(l-pl,0.0d0,d)
   ELSE
      nb = .FALSE.
      CALL LogDiffExp(0.0d0,l-pl,d)
   ENDIF
!
   RETURN
!
! -------------------------------------------------------------------
!
 END SUBROUTINE L2D
!
!
