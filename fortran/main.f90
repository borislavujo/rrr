!
!
 PROGRAM CalcRXN
! **********************************************************************
! *                                                                    *
! * Parametres for all functions used by relaxation rate calculation   *
! *                                                                    *
! **********************************************************************
! *                                                                    *
! * Author:     Boris Fackovec                                         *
! * Date:       04/02/2016                                             *
! * Version:    1.0                                                    *
! *                                                                    *
! **********************************************************************
!
   IMPLICIT NONE
   INTEGER :: ns, nts, its, i, j
   DOUBLE PRECISION xl, klab, klba
   LOGICAL, DIMENSION(:), ALLOCATABLE :: vba
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: vpl
   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Kl
!
! -------------------------------------------------------------------
!
   OPEN(UNIT=1,STATUS='OLD',FILE='N')
   READ(1,*) ns, nts
   CLOSE(1)
   ALLOCATE(vpl(ns))
   ALLOCATE(vba(ns))
   ALLOCATE(Kl(ns,ns))
   cycDoRows: DO i=1,ns
      cycDoCols: DO j=1,ns
         Kl(i,j) = -99e9
      ENDDO cycDoCols
   ENDDO cycDoRows
!
!  read vpl
!
   OPEN(UNIT=2,STATUS='OLD',FILE='vpl')
   cycREADvpl: DO i=1,ns
      READ(2,*) vpl(i)
   ENDDO cycREADvpl
   CLOSE(2)
!
!  read Kl  
!
   OPEN(UNIT=3,STATUS='OLD',FILE='Kl3')
   cycREADKl: DO its=1,nts
      READ(3,*) i, j, xl
      Kl(i,j) = xl
      Kl(j,i) = xl - vpl(i) + vpl(j)
   ENDDO cycREADKl
   CLOSE(3)
!
!  read vb
!
   OPEN(UNIT=4,STATUS='OLD',FILE='vb')
   cycREADvb: DO i=1,ns
      READ(4,*) j
      IF (j.EQ.0) THEN
         vba(i) = .FALSE.
      ELSE
         vba(i) = .TRUE.
      ENDIF
   ENDDO cycREADvb
   CLOSE(4)
!
!  do relaxation
!
   CALL RXN(ns,Kl,vpl,vba,klab,klba)
!
   STOP
!
! -------------------------------------------------------------------
!
 END PROGRAM CalcRXN
!
