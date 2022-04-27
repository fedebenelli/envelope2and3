module parameters
   use iso_fortran_env
   implicit none
   real(8), save :: Tc, dc
   real(8), save :: a, b, del1
end module parameters


module converter
   use iso_fortran_env
   use constants
   implicit none

contains
   ! =============================================================================
   !  Subroutines to obtain EOS parameters from critical constants
   ! -----------------------------------------------------------------------------
   subroutine pr_params_from_crit(Tc, Pc, w, R, ac, b, m)
      real(8), intent(in) :: Tc ! Critical temperature
      real(8), intent(in) :: Pc ! Critical pressure
      real(8), intent(in) :: w ! Accentric factor
      real(8), intent(in) :: R ! Gas constant

      real(8), intent(out) :: ac ! ac parameter
      real(8), intent(out) :: b ! Covolume
      real(8), intent(out) :: m ! Constant to calculate a(T)

      ! Internal varibles
      real(8) :: del1 = 1.0d0 + sqrt(2.d0) ! delta_1 for PR since the parameters are calculated based on the RKPR method
      real(8) :: OMa, OMb, Zc, Vceos ! 
      real(8) :: RT ! R*T product

      RT = R*Tc

      call get_Zc_OMa_OMb(del1, Zc, OMa, OMb)

      ac = OMa*RT**2/Pc
      b = OMb*RT/Pc

      Vceos = Zc*R*Tc/Pc

      ! m constant to calculate a depending on temperature
      if (w <= 0.491) then
         ! m from PR
         m = 0.37464 + 1.54226*w - 0.26992*w**2
      end if
      if (w > 0.491) then
         ! PR78
         m = 0.379642 + 1.48503*w - 0.164423*w**2 + 0.016666*w**3
      end if
   end subroutine pr_params_from_crit

   subroutine srk_params_from_crit(Tc, Pc, w, R, ac, b, m)
      real(8), intent(in) :: Tc
      real(8), intent(in) :: Pc
      real(8), intent(out) :: w
      real(8), intent(in) :: R

      real(8), intent(out) :: ac
      real(8), intent(out) :: b
      real(8), intent(out) :: m

      real(8) :: del1 = 1.0D0
      real(8) :: OMa, OMb, Zc, Vceos
      real(8) :: RT

      RT = R*Tc

      call get_Zc_OMa_OMb(del1, Zc, OMa, OMb)

      ac = OMa*RT**2/Pc
      b = OMb*RT/Pc
      Vceos = Zc*R*Tc/Pc

      m = 0.48 + 1.574*w - 0.175*w**2
   end subroutine srk_params_from_crit

   subroutine rkpr_params_from_crit( &
      Tc, Pc, w, R, ac, b, k, del1, Vceos, T, Pvdat, RhoLsat &
      )
      ! Get the RKPR EOS parameters from the substance critical properties.

      ! Input
      real(8), intent(in) :: Tc ! Critical temperature
      real(8), intent(in) :: Pc ! Critical pressure
      real(8), intent(in) :: w  ! accentric factor
      real(8), intent(in) :: R  ! Gas constant

      ! Optional cases where extra specifications are made
      real(8), intent(inout), optional :: del1 ! delta_1
      real(8), intent(inout), optional :: Vceos ! Critical volume at specification
      real(8), intent(in), optional :: T       ! Temperature used to either estimate k or del1
      real(8), intent(inout), optional :: Pvdat    ! Vapor pressure used to estimate k
      real(8), intent(in), optional :: RhoLsat ! Saturation density used to estimate del1

      !------------------------------------------------------------------------
      real(8), intent(out) :: ac   ! ac parameter
      real(8), intent(out) :: b    ! covolume
      real(8), intent(out) :: k    ! k to calculate "a" with ac and T

      real(8) :: OMa, OMb, Zc, RT, del1ini, dc, Tr, &
                      a, Pv, RHOL, RHOV, phiL, delta_k, Pold, oldk, &
                      Trho, RHOld, del1_old, delta_del1

      logical :: del1_spec, Pv_spec, rhoL_spec
      del1_spec = .false.
      Pv_spec = .false.
      rhoL_spec = .false.

      if (present(del1)) then
         del1_spec = .true.
         !print *, "Specified del1"
      end if

      if (present(T)) then
         if (present(Pvdat)) then
            Pv_spec = .true.
            !print *, "Specified Pv"
         end if
      end if

      if (present(T)) then
         if (present(RhoLsat)) then
            rhoL_spec = .true.
            !print *, "Specified RhoLsat"
         end if
      end if

      ! Initialize delta_1 and get the value that statisfies the Zc
      ! condition
      ! -----------------------------------------------------------------------
      RT = R*Tc
      if (present(Vceos)) then
         ! Usual specification with Vceos
         Zc = Pc*Vceos/RT
         del1ini = D(1) + D(2)*(D(3) - Zc)**D(4) + D(5)*(D(3) - Zc)**D(6)
         call getdel1(Zc, del1ini, del1)
      end if

      if (rhoL_spec) then
         Trho = T/Tc
         del1 = 2.0  ! initial value
         RHOld = 0.d0
      end if

      ! Calculate the inbetween parameters and finally get ac and b
      ! -----------------------------------------------------------------------
17    call get_Zc_OMa_OMb(del1, Zc, OMa, OMb)

      ac = OMa*RT**2/Pc
      b = OMb*RT/Pc

      ! Obtain the k parameter
      ! -----------------------------------------------------------------------
      dc = Pc/Zc/RT
      Vceos = 1.0d0/dc
      ! initial guess for k parameter
      k = (A1*Zc + A0)*w**2 + (B1*Zc + B0)*w + (C1*Zc + C0)
      a = ac*(3/(2 + Tr))**k

      if (del1_spec) then
         Vceos = Zc*RT/Pc
      end if

      if (Pv_spec) then
         ! added 29/06/2013 in order to allow for better reproductions of Pv curves
         Tr = T/Tc
      else
         Tr = 0.7d0
         Pvdat = Pc*10**-(1.0d0 + w)
      end if

      ! Find the value of k that fits with the expected Pv
      call VaporPressure(a, b, del1, Tc, dc, Tr, Pvdat, Pv, RHOL, RHOV, phiL)
      if (Pv > Pvdat) then
         delta_k = 0.1
      else
         delta_k = -0.1
      end if

      do while (abs(Pv - Pvdat)/Pvdat > 0.005)
         Pold = Pv
         oldk = k
         k = k + delta_k
         a = ac*(3/(2 + Tr))**k
         call VaporPressure(a, b, del1, Tc, dc, Tr, Pvdat, Pv, RHOL, RHOV, phiL)
         delta_k = -(Pv - Pvdat)*(k - oldk)/(Pv - Pold)
      end do

      if (rhoL_spec) then
         ! November 2011 for RKPR specifying T, RHOLsat
         if (abs(Trho - 0.70) > 1.d-2) then
            ! get calculated RHOL when Trho is no 0.70
            Pvdat = Pc*10**-((1./Trho - 1d0)*7*(1.0D0 + w)/3)
            a = ac*(3/(2 + Trho))**k
            call VaporPressure(a, b, del1, Tc, dc, Trho, Pvdat, Pv, RHOL, RHOV, phiL)
         end if

         if (RHOld == 0.d0) then
            del1_old = del1
            ! condition for the strange case that del1=2 is solution
            if (abs(RHOL - RHOLSAT)/RHOLSAT > 1.d-4) del1 = 2.1
         else
            delta_del1 = -(RHOL - RhoLsat)*(del1 - del1_old)/(RHOL - RHOld)
            del1_old = del1
            del1 = del1 + delta_del1
         end if

         RHOld = RHOL

         if (abs(RHOL - RHOLSAT)/RHOLSAT > 1.d-4) go to 17

         call get_Zc_OMa_OMb(del1, Zc, OMa, OMb)
         ac = OMa*RT**2/Pc
         b = OMb*RT/Pc
      end if
   end subroutine rkpr_params_from_crit
   ! =============================================================================

   
   ! ==========================================================================
   !  Subroutines to obtain critical constants from EOS parameters
   ! --------------------------------------------------------------------------
   subroutine pr_critical_from_params(ac, b, m, R, Tc, Pc, w, Vceos)
      real(8), intent(in) :: ac ! a parameter
      real(8), intent(in) :: b  ! b parameter
      real(8), intent(in) :: m  ! m parameter
      real(8), intent(in) :: R  ! Gas constant

      real(8), intent(out) :: Tc ! Critical temperature
      real(8), intent(out) :: Pc ! Critical pressure
      real(8), intent(out) :: w ! accentric factor
      real(8), intent(out) :: Vceos ! Critical volume according to eos

      real(8) :: del1 = 1.0D0 + sqrt(2.0)
      real(8) :: OMa, OMb, Zc, al, be, ga

      call get_Zc_OMa_OMb(del1, Zc, OMa, OMb)

      Tc = OMb*ac/(OMa*R*b)
      Pc = OMb*R*Tc/b
      Vceos = Zc*R*Tc/Pc

      al = -0.26992
      be = 1.54226
      ga = 0.37464 - m

      w = 0.5*(-be + sqrt(be**2 - 4*al*ga))/al
   end subroutine pr_critical_from_params

   subroutine srk_critical_from_params(ac, b, m, Tc, Pc, w, Vceos)
      real(8), intent(in) :: ac ! a parameter
      real(8), intent(in) :: b  ! b parameter
      real(8), intent(in) :: m  ! m parameter

      real(8), intent(out) :: Tc ! Critical temperature
      real(8), intent(out) :: Pc ! Critical pressure
      real(8), intent(out) :: w ! accentric factor
      real(8), intent(out) :: Vceos ! Critical volume according to eos

      real(8) :: del1 = 1.0d0
      real(8) :: OMa, OMb, Zc, al, be, ga

      call get_Zc_OMa_OMb(del1, Zc, OMa, OMb)

      Tc = OMb*ac/(OMa*RGAS*b)
      Pc = OMb*RGAS*Tc/b
      Vceos = Zc*RGAS*Tc/Pc

      al = -0.175
      be = 1.574
      ga = 0.48 - m
      w = 0.5*(-be + sqrt(be**2 - 4*al*ga))/al
   end subroutine srk_critical_from_params

   subroutine rkpr_critical_from_params(ac, b, del1, k, Tc, Pc, w, Vceos)
      use constants
      real(8), intent(in) :: ac ! ac parameter
      real(8), intent(in) :: b  ! b parameter
      real(8), intent(in) :: del1 ! delta_1 parameter
      real(8), intent(in) :: k ! k parameter

      real(8), intent(out) :: Tc ! Critical temperature
      real(8), intent(out) :: Pc ! Critical pressure
      real(8), intent(out) :: w ! Accentric factor
      real(8), intent(out) :: Vceos ! Critical volume

      real(8) :: Zc, OMa, OMb, al, be, ga

      call get_Zc_OMa_OMb(del1, Zc, OMa, OMb)

      Tc = OMb*ac/(OMa*RGAS*b)
      Pc = OMb*RGAS*Tc/b
      Vceos = Zc*RGAS*Tc/Pc

      al = A1*Zc + A0
      be = B1*Zc + B0
      ga = C1*Zc + C0 - k
      w = 0.5*(-be + sqrt(be**2 - 4*al*ga))/al
   end subroutine rkpr_critical_from_params
   ! ==========================================================================

   ! ==========================================================================
   !  Extra subroutines
   ! --------------------------------------------------------------------------
   subroutine get_Zc_OMa_OMb(del1, Zc, OMa, OMb)
      real(8), intent(in) :: del1 ! RKPR delta_1 parameter
      real(8), intent(out) :: Zc ! Critical compressibility factor
      real(8), intent(out) :: OMa !
      real(8), intent(out) :: OMb !

      real(8) :: d1, y

      d1 = (1.d0 + del1**2.d0)/(1.d0 + del1)
      y = 1.d0 + (2.d0*(1.d0 + del1))**(1.0d0/3.d0) + (4.d0/(1.d0 + del1))**(1.0d0/3)
      OMa = (3.d0*y*y + 3.d0*y*d1 + d1**2.d0 + d1 - 1.0d0)/(3.d0*y + d1 - 1.0d0)**2.d0
      OMb = 1.d0/(3.d0*y + d1 - 1.0d0)
      Zc = y/(3.d0*y + d1 - 1.0d0)
   end subroutine get_Zc_OMa_OMb

   subroutine getdel1(Zc_in, del1_ini, del1)
      real(8), intent(in) ::  Zc_in
      real(8), intent(in) ::  del1_ini
      real(8), intent(out) ::  del1

      real(8) :: d1, y, del1_old, Zc, Z_old, aux, error=1.d0

      del1 = del1_ini
      d1 = (1 + del1**2)/(1 + del1)
      y = 1 + (2*(1 + del1))**(1.0d0/3) + (4/(1 + del1))**(1.0d0/3)
      Zc = y/(3*y + d1 - 1.0d0)

      del1_old = del1

      if (Zc .gt. Zc_in) then
         del1 = 1.01*del1
      else
         del1 = 0.99*del1
      end if

      do while (error >= 1.0d-6)
         d1 = (1 + del1**2)/(1 + del1)
         y = 1 + (2*(1 + del1))**(1.0d0/3) + (4/(1 + del1))**(1.0d0/3)
         Z_old = Zc
         Zc = y/(3*y + d1 - 1.0d0)
         aux = del1
         del1 = del1 - (Zc - Zc_in)*(del1 - del1_old)/(Zc - Z_old)
         del1_old = aux
         error = abs(Zc - Zc_in)
      end do
   end subroutine getdel1

   recursive subroutine VaporPressure(&
         a, b, del1, Tc, dc, Tr, PVini, Pv, RHOL, RHOV, phiL&
         )
      use constants
      real(8), intent(in) :: a
      real(8), intent(in) :: b
      real(8), intent(in) :: del1
      real(8), intent(in) :: Tc
      real(8), intent(in) :: dc
      real(8), intent(in) :: Tr
      real(8), intent(in) :: PVini

      real(8), intent(out) :: Pv
      real(8), intent(out) :: RHOL
      real(8), intent(out) :: RHOV
      real(8), intent(out) :: phiL

      real(8) :: dphi= 0.d0, P, T, V, phi, phiV, dphiold, Pold, Plast
      P = PVini
      T = Tr*Tc

      do while (RHOL < 0.9*dc .or. RHOV > dc)
         if (RHOL < 0.9*dc) then
            P = 1.01*P
         else if (RHOV > dc) then
            P = 0.99*P
         end if
         call VCALC(1, a, b, del1, T, P, V)
         RHOL = 1/V
         call VCALC(-1, a, b, del1, T, P, V) ! SOLVE for vapor density
         RHOV = 1/V
      end do

      call FUG_CALC(a, b, del1, T, P, 1/RHOL, phi)
      phiL = phi
      call FUG_CALC(a, b, del1, T, P, V, phi)
      phiV = phi
      dphiold = dphi
      dphi = phiV - phiL

      ! ASK: Is this really a recursion?
      Plast = P
      if (ABS(dphi) .gt. ERRMAX) then
         Pold = Plast
         Plast = P
         if (dphiold == 0.0D0 .or. Tr .gt. 0.975) then
            P = P*(phiL/phiV)
         else
            P = Plast - dphi*(Plast - Pold)/(dphi - dphiold)
         end if
         call VaporPressure(a, b, del1, Tc, dc, Tr, P, Pv, RHOL, RHOV, phiL)
      end if
      PV = P
      return
   end

   recursive subroutine VCALC(ITYP, a, b, del1, T, P, V)
      ! ROUTINE FOR CALCULATION OF VOLUME, GIVEN PRESSURE
      use constants, only: RGAS

      integer, intent(in) :: ITYP ! Type of root desired 1 for liquid -1 for vapor
      real(8), intent(in) :: a
      real(8), intent(in) :: b
      real(8), intent(in) :: del1
      real(8), intent(in) :: T
      real(8), intent(in) :: P

      real(8), intent(out) :: V

      ! Internal variables
      integer :: ITER
      logical :: FIRST_RUN
      real(8) :: ZETMIN, ZETMAX, ZETA, F, F_V, F_2V, F_N, &
         PCALC, del, AT, DER, VVAP, AVAP

      FIRST_RUN = .TRUE.
      ITER = 0
      ZETMIN = 0.D0
      ZETMAX = .99D0

      if (ITYP .GT. 0) then
         ! Liquid estimate
         ZETA = .5D0
      else
         ! Ideal gas estimate
         ZETA = MIN(.5D0, b*P/(RGAS*T))
      end if

      del = 1.d0
      do while (abs(del) > 1d-10)
         V = b/ZETA
         ITER = ITER + 1
         CALL vdWg_Derivs(a, b, del1, T, V, F, F_V, F_2V, F_N)
         PCALC = RGAS*T*(1/V - F_V)

         if (PCALC .GT. P) then
            ZETMAX = ZETA
         else
            ZETMIN = ZETA
         end if

         AT = F - LOG(V) + V*P/(T*RGAS)
         DER = RGAS*T*(F_2V + 1.d0)/b
         DEL = -(PCALC - P)/DER
         ZETA = ZETA + MAX(MIN(DEL, 0.1d0), -0.1d0)

         if (ZETA .GT. ZETMAX .OR. ZETA .LT. ZETMIN) then
            ZETA = .5D0*(ZETMAX + ZETMIN)
         end if
      end do

      if (ITYP == 0) then
         VVAP = V
         AVAP = AT
         ! Calculate Liquid volume and It's energy
         call VCALC(1, a, b, del1, T, P, V)
         call vdWg_Derivs(a, b, del1, T, V, F, F_V, F_2V, F_N)
         AT = F - LOG(V) + V*P/(T*RGAS)
         if (AT .GT. AVAP) V = VVAP
      end if
   end subroutine VCALC

   subroutine vdWg_Derivs(a, b, del1, T, V, F, F_V, F_2V, F_N)
      ! CALCULATES THE CONTRIBUTION TO THE RESIDUAL, REDUCED HELMHOLZ ENERGY (F)
      ! AND ITS FIRST AND SECOND DERIVATIVE WRT V
      use constants

      real(8), intent(in) :: a
      real(8), intent(in) :: b
      real(8), intent(in) :: del1
      real(8), intent(in) :: T ! Temperature [K]
      real(8), intent(in) :: V ! Volume: [mL/mol] or [mL] for checking n-derivatives

      real(8), intent(out) :: F    ! A^RES/RT CONTRIBUTION (DIMENSIONLESS) or (MOLES)
      real(8), intent(out) :: F_V  ! 1ST V-DERIVATIVE OF F
      real(8), intent(out) :: F_2V ! 1ST V-DERIVATIVE OF F_V  (*V**2)
      real(8), intent(out) :: F_N  ! 1ST N-DERIVATIVE OF F

      real(8) :: C, aRT, ETA, SUMC, SUMD, REP, ATT, ATTV, REPV, REP2V, ATT2V

      C = (1 - del1)/(1 + del1)
      aRT = a/(RGAS*T)
      ETA = 0.25*b/V
      SUMC = c*b + V
      SUMD = del1*b + V
      REP = -log(1 - 4*ETA)
      ATT = aRT*LOG(SUMD/SUMC)/(b*(C - del1))
      ATTV = aRT/SUMC/SUMD
      REPV = 1/(1 - 4*ETA) - 1
      REP2V = 1/(1 - 4*ETA)**2 - 1
      ATT2V = aRT*V**2*(1/SUMD**2 - 1/SUMC**2)/(b*(C - del1))
      F = REP + ATT
      F_V = (-REPV/V + ATTV)
      F_2V = REP2V - ATT2V
      F_N = REP + ATT - V*F_V
   end

   subroutine FUG_CALC(a, b, del1, T, P, V, phi)
      use constants
      real(8), intent(in) :: a
      real(8), intent(in) :: b
      real(8), intent(in) :: del1
      real(8), intent(in) :: T
      real(8), intent(in) :: P
      real(8), intent(in) :: V
      real(8), intent(out) :: phi

      real(8) :: RT, Z, F, F_V, F_2V, F_N
      
      RT = RGAS*T
      Z = P*V/RT
      call vdWg_Derivs(a, b, del1, T, V, F, F_V, F_2V, F_N)
      phi=exp(F_N)/Z
   end
   ! ==========================================================================
end module converter
