!
!
 SUBROUTINE RXN(n,Kl0,vpl0,vba,lkab,lkba)
! **********************************************************************
! *                                                                    *
! * Calculates log rate constants between 2 states                     *
! *   using enhanced exponential lumping                               *
! **********************************************************************
! *                                                                    *
! * Author:     Boris Fackovec                                         *
! * Date:       11/02/2016                                             *
! * Version:    1.2                                                    *
! *                                                                    *
! **********************************************************************
!
   USE ParamsRXN, ONLY: howFine, lstartdt, thresh, ln2, minxl
!, mindoldd
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n
   DOUBLE PRECISION, INTENT(IN), DIMENSION(n) :: vpl0
   DOUBLE PRECISION, INTENT(IN), DIMENSION(n,n) :: Kl0
   LOGICAL, INTENT(IN), DIMENSION(n) :: vba
   DOUBLE PRECISION, INTENT(OUT) :: lkab, lkba
!
! -------------------------------------------------------------------
!
   INTEGER :: i, j, k, wcol, wrow
   DOUBLE PRECISION :: pl, ltemp, xl, dmax, pl1, pl2 ! doldd
   DOUBLE PRECISION :: lds, ldt, ltnow, ltau
   INTEGER :: nFine, na, nb
   DOUBLE PRECISION, DIMENSION(n,n) :: Kl, L1, D1, Ls, Ds, Lp, Dp, &
        Ltem, Dtem    ! Dold
   LOGICAL, DIMENSION(n,n) :: Nb1, Nbs, Nbp, Nbtem
   DOUBLE PRECISION, DIMENSION(n) :: vlnow, vpla, vplb, vpl, vtemp
   DOUBLE PRECISION, DIMENSION(2) :: vUjo
   DOUBLE PRECISION, DIMENSION(3*n) :: vtr
   LOGICAL :: btemp, bl
   LOGICAL, DIMENSION(3*n) :: vbtr
   DOUBLE PRECISION, DIMENSION(9999) :: vXdata1, vXdata2 
! max ratio of timescales = 2^(9999/(2^3))
   INTEGER :: sizeX, nind
!
! -------------------------------------------------------------------
!
   WRITE(*,*) "Subroutine RXN"
!
!  make sure populations sum to 1 and rate mat satisfies detailed bal
!
   CALL LogSumExp(n,vpl0,pl)
   cycNormvpl: DO i=1,n
      vpl(i) = vpl0(i) - pl
   ENDDO cycNormvpl
   Kl = Kl0
   CALL SymmetriseRateMat(n,Kl,vpl)
!
!  calculate quilibrium populations
!
   WRITE(*,*) "Calc equil populs"
   na = 0
   nb = 0
   cycCountA: DO i=1,n
      IF (vba(i)) THEN 
         na = na + 1 
         vpla(na) = vpl(i)
      ELSE
         nb = nb + 1
         vplb(nb) = vpl(i)
      ENDIF
   ENDDO cycCountA
   CALL LogSumExp(na,vpla,pl1)
   CALL LogSumExp(nb,vplb,pl2) ! more robust than CALL LogDiffExp(0.0d0,pl1,pl2)
!
!  calculate optimal minimum (log) time step lds
!
   nFine = 2**howFine
   cycDiagRates: DO j=1,n
      cycRows: DO i=1,n
         vlnow(i) = Kl(i,j)
      ENDDO cycRows
      vlnow(j) = vlnow(n)
      CALL LogSumExp(n-1,vlnow,pl)
      vtemp(j) = pl
   ENDDO cycDiagRates
   lds = -MAXVAL(vtemp) + lstartdt - REAL(howFine)*ln2
   WRITE(*,'(A8,E20.7)') "lds", lds
!
!  produce log transition matrix L = log(t)
!
   WRITE(*,*) "Make the log trans mat"
   cycTMCol: DO j=1,n
      cycTMRow: DO i=1,n
         Ls(i,j) = Kl(i,j) + lds
         vlnow(i) = Ls(i,j)
      ENDDO cycTMRow
      vlnow(j) = vlnow(n)
      CALL LogSumExp(n-1,vlnow,ltemp)
      CALL LogDiffExp(0.0d0,ltemp,pl)
      Ls(j,j) = pl
   ENDDO cycTMCol
!
!  produce D and N matrices for log(1-t) formalism from L
!
   WRITE(*,*) "produce D and N matrices"
   cycDoRows: DO i=1,n
      cycDoCols: DO j=1,n
         CALL L2D(Ds(i,j),Nbs(i,j),Ls(i,j),vpl(i))
      ENDDO cycDoCols
   ENDDO cycDoRows
!
!  log transition matrices for dt = nFine*exp(lds)
!
   WRITE(*,*) "get the matrices for nFine*lds"
   L1 = Ls
   D1 = Ds
   Nb1 = Nbs
   IF (howFine.GT.0) THEN
      cycDoubling: DO i=1,howFine
         CALL MultLogMat(n,vpl,L1,D1,Nb1,L1,D1,Nb1,Ltem,Dtem,Nbtem)
         L1 = Ltem
         D1 = Dtem
         Nb1 = Nbtem
      ENDDO cycDoubling
   ENDIF
   xl = 0
   nind = 1
   vXdata1(nind) = -99e9
   vXdata2(nind) = xl
   ldt = lds + REAL(howFine)*ln2
!   doldd = 1.0d1
   dmax = 1.0d1
!
!  main cycle
!
   WRITE(*,*) "Main cycle"
!   cycMain: DO WHILE( ((xl.GT.minxl).OR.(doldd.GT.mindoldd))&
!        .AND.(dmax.GT.thresh))
   cycMain: DO WHILE((xl.GT.minxl).OR.(dmax.GT.thresh))
!      Dold = D1
      CALL LPopul(n,L1,D1,Nb1,vpl,vba,xl)
      vXdata1(nind) = ldt
      vXdata2(nind) = xl
      WRITE(*,'(A10,I6,A10,E20.7,A10,E20.7)') "nind", nind, "ltnow", ltnow, "xl", xl
      nind = nind + 1
      Lp = L1
      Dp = D1
      Nbp = Nb1
      ltnow = ldt
      IF (nFine.GT.1) THEn
      cycLinProp: DO i=1,nFine-1
         CALL MultLogMat(n,vpl,Ls,Ds,Nbs,Lp,Dp,Nbp,Ltem,Dtem,Nbtem)
         Lp = Ltem
         Dp = Dtem
         Nbp = Nbtem
         CALL LPopul(n,Lp,Dp,Nbp,vpl,vba,xl)
         vUjo(1) = ltnow
         vUjo(2) = lds
         CALL LogSumExp(2,vUjo,ltnow)
         vXdata1(nind) = ltnow
         vXdata2(nind) = xl
         WRITE(*,'(A10,I6,A10,E20.7,A10,E20.7)') "nind", nind, "ltnow", ltnow, "xl", xl
         nind = nind + 1
      ENDDO cycLinProp
      ENDIF
      CALL MultLogMat(n,vpl,L1,D1,Nb1,L1,D1,Nb1,Ltem,Dtem,Nbtem)
      L1 = Ltem
      D1 = Dtem
      Nb1 = Nbtem
      CALL MultLogMat(n,vpl,Ls,Ds,Nbs,Ls,Ds,Nbs,Ltem,Dtem,Nbtem)
      Ls = Ltem
      Ds = Dtem
      Nbs = Nbtem
      ldt = ldt + ln2
      lds = lds + ln2
!      doldd = LOG(SUM(ABS(D1-Dold)))
      dmax = -99e9
      wcol = 0
      wrow = 0
      DO i=1,n-1
        DO j=i+1,n
          IF (D1(i,j).GT.dmax) THEN
            dmax = D1(i,j)
            wrow = i
            wcol = j
          ENDIF
          IF (D1(j,i).GT.dmax) THEN
            dmax = D1(j,i)
            wrow = j
            wcol = i
          ENDIF
        ENDDO
      ENDDO
      WRITE(*,*) "xl", xl, "dmax", dmax, "which", wrow, wcol
!      WRITE(*,*) "xl", xl, "doldd", doldd, "dmax", dmax, "which", MAXLOC(D1)
   ENDDO cycMain
   CALL TrapzLog(nind-1,vXdata1,vXdata2,ltau)
   WRITE(*,*) "ltau", ltau
   WRITE(*,*) "pl1", pl1, "pl2", pl2, "p1", EXP(pl1), "p2", EXP(pl2)
   lkab = pl1 - ltau
   lkba = pl2 - ltau
   WRITE(*,*) "lkab", lkab, "lkba", lkba
!
   RETURN
!
! -------------------------------------------------------------------
!
 END SUBROUTINE RXN
!
!
!
 SUBROUTINE SymmetriseRateMat(n,Kl,vpl)
! **********************************************************************
! *                                                                    *
! * Ensures the rate matrix satisfies detailed balance                 *
! *   Kl must have -99e9 for zero rates and -1e99 at its diagonal      *
! *                                                                    *
! **********************************************************************
! *                                                                    *
! * Author:     Boris Fackovec                                         *
! * Date:       02/02/2016                                             *
! * Version:    1.0                                                    *
! *                                                                    *
! **********************************************************************
!
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n
   DOUBLE PRECISION, INTENT(IN), DIMENSION(n) :: vpl
   DOUBLE PRECISION, INTENT(INOUT), DIMENSION(n,n) :: Kl
!
! -------------------------------------------------------------------
!
   INTEGER :: i, j
   DOUBLE PRECISION :: recTau, pl
   DOUBLE PRECISION, DIMENSION(2) :: vUjo
   LOGICAL, DIMENSION(3*n) :: vbtr
!
! -------------------------------------------------------------------
!
!   cycDoDiags: DO i=1,n
!      Kl(i,i) = -1e99
!   ENDDO cycDoDiags
   cycDoRows: DO i=1,n-1
      cycDoCols: DO j=i+1,n
         vUjo(1) = Kl(i,j)
         vUjo(2) = Kl(j,i)
         CALL LogSumExp(2,vUjo,recTau)
         vUjo(1) = vpl(i)
         vUjo(2) = vpl(j)
         CALL LogSumExp(2,vUjo,pl)
         Kl(i,j) = vpl(i) + recTau - pl
         Kl(j,i) = vpl(j) + recTau - pl
      ENDDO cycDoCols
   ENDDO cycDoRows
!
   RETURN
!
! -------------------------------------------------------------------
!
 END SUBROUTINE SymmetriseRateMat
!
!
!
 SUBROUTINE TrapzLog(n,vX1,vX2,ltau)
! **********************************************************************
! *                                                                    *
! * Integrates a vector in log formalism using trapezoidal ruel        *
! *                                                                    *
! **********************************************************************
! *                                                                    *
! * Author:     Boris Fackovec                                         *
! * Date:       11/02/2016                                             *
! * Version:    1.1                                                    *
! *                                                                    *
! **********************************************************************
!
   USE ParamsRXN, ONLY: ln2
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n
   DOUBLE PRECISION, INTENT(IN), DIMENSION(n) :: vX1, vX2
   DOUBLE PRECISION, INTENT(OUT) :: ltau
!
! -------------------------------------------------------------------
!
   INTEGER :: i
   DOUBLE PRECISION :: ldt, pl
   DOUBLE PRECISION, DIMENSION(2) :: vUjo
!
! -------------------------------------------------------------------
!
   ltau = -99d9
   cycTrapez: DO i=1,n-1
      CALL LogDiffExp(vX1(i+1),vX1(i),ldt)
      vUjo(1) = vX2(i)
      vUjo(2) = vX2(i+1)
      CALL LogSumExp(2,vUjo,pl)
      pl = pl - ln2
      vUjo(1) = ltau
      vUjo(2) = ldt + pl
      CALL LogSumExp(2,vUjo,ltau)
   ENDDO cycTrapez
!
   RETURN
!
! -------------------------------------------------------------------
!
 END SUBROUTINE TrapzLog
!
