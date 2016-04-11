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
! * Date:       03/02/2016                                             *
! * Version:    1.0                                                    *
! *                                                                    *
! **********************************************************************
!
   IMPLICIT NONE
   INTEGER :: ns, nts, i
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: vpl
   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Kl3
   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Kl
!
! -------------------------------------------------------------------
!
   OPEN(UNIT=1,STATUS='OLD',FILE='N')
   READ(1,*) ns, nts
   CLOSE(1)
   ALLOCATE(vpl(ns))
   ALLOCATE(Kl3(nts,3))
   ALLOCATE(Kl(ns,ns))
!
! -------------------------------------------------------------------
!
 END PROGRAM CalcRXN
!
