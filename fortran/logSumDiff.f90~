!
!
 SUBROUTINE LogSumDiff(n,vL,vbPM,lSD,bPM)
! **********************************************************************
! *                                                                    *
! * Returns sum of (minus if vbPM is .false.) exponentials of          *
! * terms of vL                                                        *
! *                                                                    *
! **********************************************************************
! *                                                                    *
! * Author:     Boris Fackovec                                         *
! * Date:       03/02/2016                                             *
! * Version:    1.1                                                    *
! *                                                                    *
! **********************************************************************
!
   USE M_MRGRNK, ONLY: MRGRNK
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n
   DOUBLE PRECISION, INTENT(IN), DIMENSION(n) :: vL
   LOGICAL, INTENT(IN), DIMENSION(n) :: vbPM
   DOUBLE PRECISION, INTENT(IN) :: lSD
   LOGICAL, INTENT(IN) :: bPM
!
! -------------------------------------------------------------------
!
   INTEGER :: i
   INTEGER, DIMENSION(n) :: vInd, vI
   DOUBLE PRECISION :: l, l0
   DOUBLE PRECISION, DIMENSION(2) :: vUjo
   LOGICAL :: b, b0
!
! -------------------------------------------------------------------
!
   IF (n.EQ.1) THEN
      lSD = vL(1)
      bPM = vbPM(1)
      RETURN
   ENDIF
!
   CALL MRGRNK(vL,vInd)
   cycInds: DO i=1,n
      vI(vInd(i)) = i
   ENDDO cycSetMasses
!
!  now start from smallest term and add/subtract
!
   lSD = vL(vI(1))
   bPM = vbPM(vI(1))
   cycAddSub: DO i=2:n
      l0 = vL(vI(i))
      b0 = vbPM(vI(i))
      IF (b0.EQV.bPM) THEN
         vUjo(1) = lSD
         vUjo(2) = l0
         CALL logSumExp(2,vUjo,lSD)
      ELSE IF (lSD.GT.l0) THEN
         CALL LogDiffExp(lSD,l0,l)
         lSD = l
      ELSE IF (l0.GT.lSD) THEN
         CALL LogDiffExp(l0,lSD,l)
         lSD = l
         bPM = b0
      ELSE
         lSD = -1e99
         bPM = .TRUE.
      ENDIF
   ENDDO cycAddSub
!
   RETURN
!
! -------------------------------------------------------------------
!
 END SUBROUTINE LogSumDiff
!
!
!
 SUBROUTINE LogSumExp(n,vL,lSE)
! **********************************************************************
! *                                                                    *
! * Returns log of Sum of Exp of terms                                 *
! *                                                                    *
! **********************************************************************
! *                                                                    *
! * Author:     Boris Fackovec                                         *
! * Date:       21/01/2016                                             *
! * Version:    1.0                                                    *
! *                                                                    *
! **********************************************************************
!
   USE M_REFSOR, ONLY: REFSOR
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n
   DOUBLE PRECISION, INTENT(IN), DIMENSION(n) :: vL
   DOUBLE PRECISION, INTENT(IN) :: lSE
!
! -------------------------------------------------------------------
!
   INTEGER :: i
   DOUBLE PRECISION :: l
   DOUBLE PRECISION, DIMENSION(n) :: vLT
!
! -------------------------------------------------------------------
!
   IF (n.EQ.1) THEN
      lSE = vL(1)
      RETURN
   ELSE IF (n.EQ.2) THEN
      lSE = MAX(vL(1),vL(2))
      l = MIN(vL(1),vL(2))-lSE
      lSE = lSE + LOG(1+EXP(l))
      RETURN
   ENDIF
!
   vLT = vL
   CALL REFSOR(vLT)
!
!  now start adding from the smallest term to the largest
!
   lSE = vLT(n)
   cycInds: DO i=1,n-1
      l = vLT(n-i)
      lSE = l + LOG(1+EXP(lSE-l))
   ENDDO cycInds
!  a few orders expression for small values of lSE-l could be used
!  in order to speed up the calculation
!
   RETURN
!
! -------------------------------------------------------------------
!
 END SUBROUTINE LogSumExp
!
!
!
 SUBROUTINE LogDiffExp(l1,l2,ldif)
! **********************************************************************
! *                                                                    *
! * Returns log of Diff of Exp of l1 and l2                            *
! *                                                                    *
! **********************************************************************
! *                                                                    *
! * Author:     Boris Fackovec                                         *
! * Date:       21/01/2016                                             *
! * Version:    1.0                                                    *
! *                                                                    *
! **********************************************************************
!
   IMPLICIT NONE
   DOUBLE PRECISION, INTENT(IN) :: l1
   DOUBLE PRECISION, INTENT(IN) :: l2
   DOUBLE PRECISION, INTENT(IN) :: ldif
!
! -------------------------------------------------------------------
!
   DOUBLE PRECISION :: l
!
! -------------------------------------------------------------------
!
   IF (l1.LT.l2) THEN
      WRITE(*,*) "warning: switching elements in LogDiffExp!"
      l = l1
      l1 = l2
      l2 = l
   ELSE IF (l1.EQ.l2) THEN
      ldif = -1e99
      RETURN
   ENDIF
   l = l2-l1
   lSE = l1 + LOG(1+EXP(l))
!
   RETURN
!
! -------------------------------------------------------------------
!
 END SUBROUTINE LogSumExp
!
