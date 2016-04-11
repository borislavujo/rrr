!
!
 SUBROUTINE MultLogMat(n,vpl,L1,D1,Nb1,L2,D2,Nb2,L3,D3,Nb3)
! **********************************************************************
! *                                                                    *
! * Multiplies transition matrices in the log and log(1-t) formalisms  *
! *                                                                    *
! **********************************************************************
! *                                                                    *
! * Author:     Boris Fackovec                                         *
! * Date:       11/02/2016                                             *
! * Version:    1.2                                                    *
! *                                                                    *
! **********************************************************************
!
   USE ParamsRXN, ONLY: ln2, h
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n
   DOUBLE PRECISION, INTENT(IN), DIMENSION(n) :: vpl
   DOUBLE PRECISION, INTENT(IN), DIMENSION(n,n) :: L1, L2
   DOUBLE PRECISION, INTENT(IN), DIMENSION(n,n) :: D1, D2
   LOGICAL, INTENT(IN), DIMENSION(n,n) :: Nb1, Nb2
   DOUBLE PRECISION, INTENT(OUT), DIMENSION(n,n) :: L3
   DOUBLE PRECISION, INTENT(OUT), DIMENSION(n,n) :: D3
   LOGICAL, INTENT(OUT), DIMENSION(n,n) :: Nb3
!
! -------------------------------------------------------------------
!
   INTEGER :: i, j, k
   DOUBLE PRECISION :: pl, ltemp, pl2, pplus, pminus, popdif, tempdif
   INTEGER :: nind, nplus, nminus, wrow, wcol
   DOUBLE PRECISION, DIMENSION(n) :: vlnow, vplus, vminus
   DOUBLE PRECISION, DIMENSION(2) :: vUjo
   DOUBLE PRECISION, DIMENSION(3*n) :: vtr
   LOGICAL :: btemp, bl
   LOGICAL, DIMENSION(3*n) :: vbtr
   LOGICAL, DIMENSION(2) :: vbUjo
   LOGICAL, DIMENSION(n,n) :: Nbted
   CHARACTER(LEN=32) :: fmt_f
!
! -------------------------------------------------------------------
!
   fmt_f = '(16E14.5)'
!
!  L3 = L1 * L2
!
   cycCalcRowL3: DO j=1,n
      cycCalcColL3: DO i=1,n
         cycCorrL12: DO k=1,n
            vlnow(k) = L1(i,k)+L2(k,j)
         ENDDO cycCorrL12
         CALL LogSumExp(n,vlnow,ltemp)
         L3(i,j) = ltemp
      ENDDO cycCalcColL3
   ENDDO cycCalcRowL3
!
!  symmetrise L3
!
   cycDoRows: DO i=1,n-1
      cycDoCols: DO j=i+1,n
         vUjo(1) = L3(i,j)
         vUjo(2) = L3(j,i)
         CALL LogSumExp(2,vUjo,ltemp)
         vUjo(1) = vpl(i)
         vUjo(2) = vpl(j)
         CALL LogSumExp(2,vUjo,pl)
         L3(i,j) = vpl(i) + ltemp - pl
         L3(j,i) = vpl(j) + ltemp - pl
      ENDDO cycDoCols
   ENDDO cycDoRows
!
!  normalise L3
!
   cycL3norm: DO j=1,n
      cycGetL3Trms: DO i=1,n
         vlnow(i) = L3(i,j)
      ENDDO cycGetL3Trms
      CALL LogSumExp(n,vlnow,ltemp)
      cycNormL3: DO i=1,n
         L3(i,j) = L3(i,j) - ltemp
      ENDDO cycNormL3
   ENDDO cycL3norm
!   WRITE(*,*) "L3"
!   DO i=1,n
!      WRITE(*,fmt_f) L3(i,:)
!   ENDDO
!
!  calculate D3 AND Nb3
!
   cycCalcRowD3: DO i=1,n
      cycCalcColD3: DO j=1,n
         IF (L3(i,j)-vpl(i).GT.h) THEN
            Nbted(i,j) = .TRUE.
            cycFillvtr: DO k=1,n
               vtr((k-1)*3+1)  = vpl(k) + D1(i,k)
               vtr((k-1)*3+2)  = vpl(k)           + D2(k,j)
               vtr(k*3)        = vpl(k) + D1(i,k) + D2(k,j)
               vbtr((k-1)*3+1) = Nb1(i,k)
               vbtr((k-1)*3+2) = Nb2(k,j)
               vbtr(k*3)       = Nb1(i,k).EQV.Nb2(k,j)
            ENDDO cycFillvtr
            CALL LogSumDiff(3*n,vtr,vbtr,ltemp,btemp)
            Nb3(i,j) = btemp
            D3(i,j) = ltemp
            IF ((.NOT.Nb3(i,j)).AND.(D3(i,j).GT.0)) THEN
               WRITE(*,*) "warning: negative popul was calculated for i, j", i, &
                    j, "D(i,j):", D3(i,j), "D(i,j)+vpl(i)", D3(i,j)+vpl(i)
               D3(i,j) = 0.0d0
!            WRITE(*,*) "i, j", i, j, "Nb3(i,j)", Nb3(i,j), "D3(i,j)", D3(i,j)
!            WRITE(*,*) "vtr", vtr
!            WRITE(*,*) "vbtr", vbtr
            ENDIF
         ELSE
            Nbted(i,j) = .FALSE.
            CALL L2D(D3(i,j),Nb3(i,j),L3(i,j),vpl(i))
         ENDIF
      ENDDO cycCalcColD3
   ENDDO cycCalcRowD3
!   WRITE(*,*) "how produced"
!   DO i=1,n
!      WRITE(*,*) Nbted(i,:)
!   ENDDO
!   WRITE(*,*) "before symmetris"
!   DO i=1,n
!      WRITE(*,fmt_f) D3(i,:)
!   ENDDO
!   WRITE(*,*) "N"
!   DO i=1,n
!      WRITE(*,*) Nb3(i,:)
!   ENDDO
!
!  symmetrise D3
!
   popdif = -99e9
   wrow = 0
   wcol = 0
   cycSymmD3Row: DO i=1,n-1
      cycSymmD3Col: DO j=i+1,n
         IF (D3(i,j).GT.D3(j,i)) THEN
            Nb3(j,i) = Nb3(i,j)
         ELSE
            Nb3(i,j) = Nb3(j,i)
         ENDIF
         vUjo(1) = D3(i,j)+vpl(i)
         vUjo(2) = D3(j,i)+vpl(j)
         CALL LogSumExp(2,vUjo,ltemp)
         vUjo(1) = vpl(i)
         vUjo(2) = vpl(j)
         CALL LogSumExp(2,vUjo,pl)
         ltemp = ltemp - pl
         CALL LogDiffExp(MAX(ltemp,D3(i,j)),MIN(ltemp,D3(i,j)),tempdif)
         IF(vpl(i)+tempdif.GT.popdif) THEN
           popdif = vpl(i) + tempdif
           wrow = i
           wcol = j
         ENDIF
         CALL LogDiffExp(MAX(ltemp,D3(j,i)),MIN(ltemp,D3(j,i)),tempdif)
         IF(vpl(j)+tempdif.GT.popdif) THEN
           popdif = vpl(j) + tempdif
           wrow = j
           wcol = i
         ENDIF
         D3(i,j) = ltemp
         D3(j,i) = ltemp
      ENDDO cycSymmD3Col
   ENDDO cycSymmD3Row
   WRITE(*,*) "symmetris difference", popdif, wrow, wcol
!   DO i=1,n
!      WRITE(*,fmt_f) D3(i,:)-Ddif(i,:)
!   ENDDO
!   WRITE(*,*) "after symmetris"
!   DO i=1,n
!      WRITE(*,fmt_f) D3(i,:)
!   ENDDO
!   WRITE(*,*) "N"
!   DO i=1,n
!      WRITE(*,*) Nb3(i,:)
!   ENDDO
!
!  normalise d3
!
   popdif = -99e9
   wcol = 0
   cycNormCols: DO j=1,n
      cycGetTrms: DO i=1,n
         vtr(i)  = D3(i,j)+vpl(i)
         vbtr(i) = Nb3(i,j)
      ENDDO cycGetTrms
      vtr(j) = vtr(n)
      vbtr(j) = vbtr(n)
      CALL LogSumDiff(n-1,vtr,vbtr,pl,bl)
      Nb3(j,j) = .NOT.bl
      CALL LogDiffExp(MAX(D3(j,j)+vpl(j),pl),MIN(D3(j,j)+vpl(j),pl),tempdif)
      IF (tempdif.GT.popdif) THEN
        wcol = j
        popdif = tempdif
      ENDIF
      D3(j,j) = pl - vpl(j)
      IF ((D3(j,j).GT.0).AND.bl) THEN
         WRITE(*,*) "Warning: column ", j, " sums to x = 1 + exp", D3(j,j)+vpl(j)
!         WRITE(*,*) "The problematic column", j, ":", D3(1:n,j)
         D3(j,j) = 0.0d0
! here perhaps this excess population should be distributed to remainind D(i,i)
      ENDIF
   ENDDO cycNormCols
   WRITE(*,*) "normalisation difference", popdif, wcol
!   DO i=1,n
!      WRITE(*,fmt_f) D3(i,:)-Ddif(i,:)
!   ENDDO
!   WRITE(*,*) "after normalis"
!   DO i=1,n
!      WRITE(*,fmt_f) D3(i,:)
!   ENDDO
!   WRITE(*,*) "N"
!   DO i=1,n
!      WRITE(*,*) Nb3(i,:)
!   ENDDO
!   WRITE(*,*) "final D*P"
!   DO i=1,n
!      WRITE(*,fmt_f) D3(i,:)+vpl
!   ENDDO
!
   RETURN
!
! -------------------------------------------------------------------
!
 END SUBROUTINE MultLogMat
!
