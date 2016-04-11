!
!
 MODULE ParamsRXN
! **********************************************************************
! *                                                                    *
! * Parametres for all functions used by relaxation rate calculation   *
! *                                                                    *
! **********************************************************************
! *                                                                    *
! * Author:     Boris Fackovec                                         *
! * Date:       11/02/2016                                             *
! * Version:    1.1                                                    *
! *                                                                    *
! **********************************************************************
!
   INTEGER :: howFine = 2
   DOUBLE PRECISION :: lstartdt = -6.0d0
   DOUBLE PRECISION :: thresh = -6.0d0
   DOUBLE PRECISION :: minxl = -5.0d0
   DOUBLE PRECISION :: mindoldd = -3.0d0
   DOUBLE PRECISION :: ln2 = LOG(2.0d0)
   DOUBLE PRECISION :: h = -3.0d0
!
! -------------------------------------------------------------------
!
 END MODULE ParamsRXN
!
