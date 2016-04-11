!
!
 SUBROUTINE LPopul(n,L,D,Nb,vpl,va,xl)
! **********************************************************************
! *                                                                    *
! * Initialize geometry, velocities and parameters for MD simulation   *
! *                                                                    *
! **********************************************************************
! *                                                                    *
! * Author:     Boris Fackovec                                         *
! * Date:       30/01/2016                                             *
! * Version:    1.0                                                    *
! *                                                                    *
! **********************************************************************
!
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n
   DOUBLE PRECISION, INTENT(IN), DIMENSION(n,n) :: L
   DOUBLE PRECISION, INTENT(IN), DIMENSION(n,n) :: D
   LOGICAL, INTENT(IN), DIMENSION(n,n) :: Nb
   DOUBLE PRECISION, INTENT(IN), DIMENSION(n) :: vpl
   LOGICAL, INTENT(IN), DIMENSION(n) :: va
   DOUBLE PRECISION, INTENT(OUT) :: x
!
! -------------------------------------------------------------------
!
   INTEGER :: na, i, j, k
   DOUBLE PRECISION :: pl, pl1, pl2, pltemp
   INTEGER, DIMENSION(n) :: nind
   DOUBLE PRECISION, DIMENSION(n) :: vpla, vplnow, vpltemp
   DOUBLE PRECISION, DIMENSION(n**2) :: vf1
   LOGICAL :: bf1
   LOGICAL, DIMENSION(n**2) :: vbf1
!
! -------------------------------------------------------------------
!
!
!  find the size of state A and calculate its population
!
   na = 0
   cycCountA: DO i=1,n
      IF (va(i)) THEN 
         na = na + 1 
         vpla(na) = vpl(i)
         vind(na) = i
      ENDIF
   ENDDO cycCountA
   CALL LogSumExp(na,vpl,pl1)
   CALL LogDiffExp(0.0d0,pl1,pl2)
!
!  get the population from log transition matrix
!
   cycSumCols: DO i=1,na
      cycSumRows: DO j=1,na
         vpltemp(j) = L(vind(i),vind(j))+vpl(vind(j))-pl1
      ENDDO cycSumRows
      CALL LogSumExp(na,vpltemp,pltemp)
      vplnow(i) = pltemp
   ENDDO cycSumCols
   CALL LogSumExp(na,vplnow,pl)
   xl = logDiffExp(pl,pl1)-pl2
   IF (xl.LT.-3d0) THEN
      RETURN
   ENDIF
!
!   get the population from log(1-t)
!
   cycAddCols: DO i=1,na
      cycAddRows: DO j=1,na
         k = (i-1)*na+j
         vf1(k) = vpl(vind(i)) + vpl(vind(j)) + D(vind(i),vind(j))
         vb1f(k) = N(vind(i),vind(j))
      ENDDO cycAddRows
   ENDDO cycAddCols
   CALL LogSumDiff(na**2,vf1,vbf1,pl,bf1)
   CALL LogDiffExp(0.0d0,pl,xl)
   xl = xl - pl1 - pl2
!
   RETURN
!
! -------------------------------------------------------------------
!
 END SUBROUTINE LPopul
!
