!
!
 SUBROUTINE BXD2D_Initialize
! **********************************************************************
! *                                                                    *
! * Initialize geometry, velocities and parameters for MD simulation   *
! *                                                                    *
! **********************************************************************
! *                                                                    *
! * Author:     Boris Fackovec                                         *
! * Date:       19/07/2013                                             *
! * Version:    1.0                                                    *
! *                                                                    *
! **********************************************************************
!
   USE ModBXD2D
   IMPLICIT NONE
   INTEGER i, j
   WRITE(*,*) "Subroutine BXD2D_Initialize"
!
! -------------------------------------------------------------------
!
!
!   allocate memory
!
   ALLOCATE (vMasses(nAtoms))
   ALLOCATE (Coords(nAtoms,nDim))
   ALLOCATE (vSigma(nAtoms))
   ALLOCATE (vEps(nAtoms))
   ALLOCATE (Vat(nAtoms,nDim))
   ALLOCATE (Accel(nAtoms,nDim))
   ALLOCATE (Forces(nAtoms,nDim))
   ALLOCATE (vRecord(nSteps))
   ALLOCATE (vEpot0(nPoints))
   ALLOCATE (PointsCoords(nPoints,nAtoms,nDim))
   ALLOCATE (RemCoords1(nRemember,nAtoms,nDim))
   ALLOCATE (RemCoords2(nRemember,nAtoms,nDim))
   WRITE(*,*) "Memory succesfully allocated"
   iRemStep = 1
   iRefl = 1
   dt0 = dt
!
!   tooFar scaled by square of number of atoms -> then user does not need
!     to increase if number of atoms is increased
!
   tooFarScaled = tooFar*SQRT(REAL(nAtoms))
!
!   filling the allocated points
!
   cycLJParams: do i=1,nAtoms
      vSigma(i) = sigmaLJ
      vEps(i) = epsilonLJ
   enddo cycLJParams
   cycSetMasses: DO i=1,nAtoms
      vMasses(i) = 1.0d0
   ENDDO cycSetMasses
!
!   if restarting trajectory from a file, read
!
   WRITE(myUnit,*) "Reading coordinates"
   IF (b_restart) THEN
      CALL BXD2D_ReadRestartData
   ELSE
      CALL BXD2D_ReadCoords
      CALL BXD2D_SeedVels
   ENDIF
!
!   move center of mass into (0,0)
!     set drift and angular velocities to 0
!
   CALL BXD2D_RPoints2D
   WRITE(myUnit,*) "Centering the structure"
   CALL BXD2D_Center(.TRUE.)
!
!   compute forces just to calculate potential energy
!
   WRITE(myUnit,*) "Computing gradient of potential energy"
   CALL BXD2D_ComputeForces
   RETURN
!
! -------------------------------------------------------------------
!
   END SUBROUTINE BXD2D_Initialize
!
!
   SUBROUTINE BXD2D_Propagator
! **********************************************************************
! *                                                                    *
! * Velocity Verlet integrator                                         *
! *   for both microcanonical and canonical simulations                *
! *                                                                    *
! **********************************************************************
! *                                                                    *
! * Author:     Boris Fackovec                                         *
! * Date:       19/07/2013                                             *
! * Version:    1.0                                                    *
! *                                                                    *
! **********************************************************************
! *                                                                    *
! * Last changes:                                                      *
! *                                                                    *
! **********************************************************************
!
   USE ModBXD2D
   IMPLICIT NONE
   DOUBLE PRECISION t_taut, dtSmall, rAcc, eTot
   DOUBLE PRECISION temp, qq, pi, shat, eta, pistar, pi2star, g
   REAL cpuTime
   INTEGER i
   INTEGER, DIMENSION(1) :: rVSeed
!
! -------------------------------------------------------------------
!
   WRITE(*,*) "Starting MD evolution"
   IF (thermType.EQ.0) THEN
      WRITE(myUnit,*) "Using Velocity Verlet integrator"
   ELSE IF (thermType.EQ.2) THEN
!
!   g - number of degrees of freedom
!
      g = REAL(nAtoms*nDim - 3)
!
!   artificial mass is a global variable
!
      qq = 1.0d0 / artifMass
!
!   pi - something arbitrary 
!      - will be changed during equilibr period
!
      pi = 1.0d0 
      WRITE(myUnit,*) "Using Nose-Hoover scheme"
   ENDIF
   time = 0
!
!   simulation starts here
!
   CALL BXD2D_Output(1)
   cycForEachStep: DO iStep=1,nSteps
      IF (reflecting.GT.0) THEN
         dt = dt * (1.0e6**(1/REAL(nEquil)))
!         WRITE(*,*) "dt", dt
         reflecting = reflecting - 1
         CALL BXD2D_AdjustEtot
      ENDIF
!      WRITE(*,*) "Step", iStep
      IF ((thermType.EQ.0).OR.(thermType.EQ.1)) THEN
!
!   use velocity verlet scheme
!
!
!   v(t+0.5dt)
!
         Vat = Vat + 5.0d-1 * Accel * dt
!
!   x(t+dt)
!
         Coords = Coords + Vat * dt
!
!   v(t+dt)
!
!         WRITE(*,*) "Computing forces"
         CALL BXD2D_ComputeForces
         cycUpdateVeloc: DO i=1,nAtoms
            Vat(i,1) = Vat(i,1) + 5.0d-1 * Accel(i,1) * dt
            Vat(i,2) = Vat(i,2) + 5.0d-1 * Accel(i,2) * dt
         ENDDO cycUpdateVeloc
!         Vat = Vat + 5.0d-1 * Accel * dt
!         WRITE(*,*) "New Vat"
!
!   center, but do not care about angular velocity
!
         time = time + dt
!         WRITE(*,*) "Calling thermostat"
         IF (thermType.EQ.1) THEN
            CALL BXD2D_Berendsen
         ENDIF
!         WRITE(*,*) "Step finished"
      ELSE
!
!   use Nose-Hoover propagator by Kleinerman et al. 2008
!
         CALL BXD2D_ComputeForces
!
!   equations are numbered in accord with Kleinerman et al. 2008
!
!
!   eqn (6):    pi* = pi + dt/4 * (2 E_k - g k T)
!
         pistar = pi + dt * 2.5d-1 * ( 2.0d0 * eKin - g * desiredT)
!
!   eqn (7):   shat = exp(- pistar dt / (2 Q) )
!
         shat = EXP(- 5.0d-1 * pistar * qq * dt )
!
!   eqn (8):   eta(n+1/2) = eta(n) + (pi* dt) / (2 Q)
!
         eta = eta + pistar * dt * 5.0d-1 * qq
!
!   calculation of p(n)*shat
!
         eKin = 0
         cycCalcEkshat1: DO i=1,nAtoms
            eKin = eKin + 5.0d-1 * shat * vMasses(i) * &
                 SUM(Vat(i,1:nDim)*Vat(i,1:nDim))
         ENDDO cycCalcEkshat1
!
!   eqn (9):   pi(n+1/2) = pi* + dt/4 ( 2 E_k(p(n)shat) - g T k )
!
         pi = pistar + 2.5d-1 * dt * ( 2.0d0 * eKin - g * desiredT )
!
!   eqn (10):  p(n+1/2) = p(n)*shat - dV/dq|(n) * dt/2
!
         Vat = Vat * shat + 5.0d-1 * Accel * dt
!
!   eqn (11):  q(n+1) = q(n) + p(n+1/2)/m * dt
!
         Coords = Coords + Vat * dt
!
!  eqn (12):  p* = p(n+1/2) - dV/dq|(n+1) * dt/2
!
         CALL BXD2D_ComputeForces
         Vat = Vat + 5.0d-1 * Accel * dt
!
!  eqn (13):  pi** = pi(n+1/2) + dt/4 ( 2 Ek(p*) - g T k )
!
         eKin = 0
         cycCalcEk2: DO i=1,nAtoms
            eKin = eKin + 5.0d-1 * vMasses(i) * &
                 SUM(Vat(i,1:nDim)*Vat(i,1:nDim))
         ENDDO cycCalcEk2
         pi2star = pi + 2.5d-1 * dt * ( 2.0d0 * eKin - g * desiredT )
!
!  eqn (14):  shat = exp( - (pi** dt) / (2 Q) )
!
         shat = EXP( - 5.0d-1 * pi2star * qq * dt )
!
!  eqn (15):  eta(n+1) = eta(n+1/2) + (pi** dt)/(2 Q)
!
         eta = eta + pi2star * dt * 5.0d-1 * qq
!
!  eqn (16): pi(n+1) = pi** + dt/4 ( 2 Ek(p*shat - g T k )
!
         eKin = 0
         cycCalcEkshat2: DO i=1,nAtoms
            eKin = eKin + 5.0d-1 * shat * vMasses(i) * &
                 SUM(Vat(i,1:nDim)*Vat(i,1:nDim))
         ENDDO cycCalcEkshat2
         pi = pi2star + 2.5d-1 * dt * ( 2.0d0 * eKin - g * desiredT )
!
!  eqn (17):  p(n+1) = p* shat
!
         Vat = Vat * shat
         time = time + dt
      ENDIF
!
!   check for evaporated atoms
!
!      WRITE(*,*) "Calling travellers"
!      CALL MD_Travellers
!
!   print output
!
!
!   zakomentovanie nasledujuceho prikazu ma za nasledok seg fault!!!???!!!
!
!      WRITE(myUnit,*)
      IF (MOD(iStep,nSkipOutput) .EQ. 0) THEN
         CALL BXD2D_Output(1)
!
!   for temperature  statistics
!
         CALL BXD2D_GetMe(1, temp)
         sumTemp = sumTemp + temp
         iTemp = iTemp + 1
      ENDIF
!
!   in microcanonical runs, numerical errors cause problems
!     with conservation of momenta
!
!      WRITE(*,*) "Centering"
      IF ((MOD(iStep,nStepsAdjust) .EQ.0).AND.(thermType.EQ.0)) THEN
         CALL BXD2D_Center(.TRUE.)
         WRITE(*,*) "Step", iStep
         CALL BXD2D_GetMe(2, eTot)
         WRITE(*,*) "Etot", eTot
         WRITE(*,*) "Adjusting"
         CALL BXD2D_AdjustEtot
!         CALL MD_Output(1)
      ENDIF
      vRecord(iStep) = ePot
!      CALL MD_GetMe(2, eTot)
!      WRITE(*,*) "Step really finished"
!      WRITE(*,*) "Etot", eTot
!
!   save the current geometry
!
      CALL BXD2D_RememberMe(0)
!
!   check box
!
      CALL BXD2D_CheckBox
   ENDDO cycForEachStep
   RETURN
!
! -------------------------------------------------------------------
!
   END SUBROUTINE BXD2D_Propagator
!
!
!
   SUBROUTINE BXD2D_Terminate
! **********************************************************************
! *                                                                    *
! * Print final coordinates                                            *
! *                                                                    *
! **********************************************************************
! *                                                                    *
! * Author:     Boris Fackovec                                         *
! * Date:       19/07/2013                                             *
! * Version:    1.0                                                    *
! *                                                                    *
! **********************************************************************
!
   USE ModBXD2D
   IMPLICIT NONE
   INTEGER i, iNow, iNeighbor
   DOUBLE PRECISION maxTime
   WRITE(myUnit,*) "Subroutine Terminate"
!
! -------------------------------------------------------------------
!
   WRITE(myUnit,*) "Final coordinates"
   cycWriteFinal: DO i=1,nAtoms
     WRITE(myUnit,*) Coords(i,1:nDim)
   ENDDO cycWriteFinal
!
   WRITE(myUnit,*) "Average temperature", sumTemp/iTemp
   WRITE(myUnit,*) "Total simulation time", time
   WRITE(myUnit,*) "recorded quantity"
!   DO i=1,nSteps
!      WRITE(myUnit,*) i, vRecord(i)
!   ENDDO
   CLOSE(myUnit)

!
!   deallocate
!
   DEALLOCATE (Coords)
   DEALLOCATE (vMasses)
   DEALLOCATE (vSigma)
   DEALLOCATE (vEps)
   DEALLOCATE (Vat)
   DEALLOCATE (Accel)
   DEALLOCATE (Forces)
   DEALLOCATE (vRecord)
   DEALLOCATE (PointsCoords)
!
! -------------------------------------------------------------------
!
   END SUBROUTINE BXD2D_Terminate
!
