!
!
 SUBROUTINE LPopul(n,L,D,Nb,vpl,vba,xl)
! **********************************************************************
! *                                                                    *
! * Calculates progress of relaxation from current transition matrix   *
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
   LOGICAL, INTENT(IN), DIMENSION(n) :: vba
   DOUBLE PRECISION, INTENT(OUT) :: xl
!
! -------------------------------------------------------------------
!
   INTEGER :: na, i, j, k
   DOUBLE PRECISION :: pl, pl1, pl2, pltemp
   INTEGER, DIMENSION(n) :: vind
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
!   WRITE(*,*) "Subroutine LPopul"
   na = 0
   cycCountA: DO i=1,n
      IF (vba(i)) THEN 
         na = na + 1 
         vpla(na) = vpl(i)
         vind(na) = i
      ENDIF
   ENDDO cycCountA
   CALL LogSumExp(na,vpla,pl1)
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
!   WRITE(*,*) "individual vpls", vplnow(1:na)
   CALL LogSumExp(na,vplnow,pl)
   WRITE(*,'(A20,F12.7)') "total state a", pl
   CALL LogDiffExp(pl,pl1,xl)
   xl = xl - pl2
   WRITE(*,'(A20,F12.7)') "xl from L", xl
!   IF (xl.GT.-3d0) THEN
 !     RETURN
  ! ENDIF
!
!   get the population from log(1-t)
!
   cycAddCols: DO i=1,na
      cycAddRows: DO j=1,na
         k = (i-1)*na+j
         vf1(k) = vpl(vind(i)) + vpl(vind(j)) + D(vind(i),vind(j))
         vbf1(k) = Nb(vind(i),vind(j))
      ENDDO cycAddRows
   ENDDO cycAddCols
   CALL LogSumDiff(na**2,vf1,vbf1,pl,bf1)
!   WRITE(*,*) "from 1 subtract", pl
!   CALL LogDiffExp(0.0d0,pl,xl)
   xl = pl - pl1 - pl2
   WRITE(*,'(A20,F12.7)') "xl from D", xl
!
!
!
   IF (xl.GT.1e-8) THEN
      WRITE(*,*) "Error: xl can only be negative"
      OPEN(UNIT=999,STATUS='NEW',FILE='Dlast')
      OPEN(UNIT=998,STATUS='NEW',FILE='Llast')
      OPEN(UNIT=997,STATUS='NEW',FILE='Nblast')
      DO i=1,n
         WRITE(999,*) D(i,:)
         WRITE(998,*) L(i,:)
         WRITE(997,*) Nb(i,:)
      ENDDO
      CLOSE(999)
      CLOSE(998)
      CLOSE(997)
      STOP
   ENDIF
!
   RETURN
!
! -------------------------------------------------------------------
!
 END SUBROUTINE LPopul
!
