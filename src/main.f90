
! Listado de commons utilizados, que podrían pasar a un módulo:
!   common /CRIT/TC(nco),PC(nco),DCeos(nco),omg(nco)
!        common /COMPONENTS/ ac(nco),b(nco),delta1(nco),rk(nco),Kij_or_K0,NTDEP
!        common /MODEL/ NMODEL
!        common /rule/ncomb
!        common /bcross/bij(nco,nco)
!        common /Tdep/ Kinf,Tstar


program calc_envelope2and3
   implicit double precision(A - H, O - Z)
   logical Comp3ph
   common/writeComp/Comp3ph, i1, i2

   open (1, FILE='envelIN.txt')
   open (2, FILE='envelOUT.txt')
   read (1, *) N

   !write (6, *) 'write extra output with compositions for 2 compounds along 3-phase lines?'
   !write (6, *) 'Enter 1 for YES. Otherwise, any other number.'
   ! read (5, *) i

   i = 0

   if (i == 1) Comp3ph = .true.
   if (Comp3ph) then
      open (3, FILE='Comp3phOUT.txt')
      write (6, *) 'Enter order numbers for two selected compounds, separated by space.'
      read (5, *) i1, i2
      write (3, *) 'Molar fractions along three-phase boundaries are printed below for compounds with order:', i1, i2
   end if

   call readcase(n)

end program

subroutine readcase(n)
   use constants
   use dtypes, only: envelope, kfcross, cross
   use array_operations, only: find_cross

   implicit double precision(A - H, O - Z)

   parameter(nco=64, Pmax=700.0)
   dimension z(n), xx(n), w(n)

   double precision Kinf
   double precision, dimension(n) :: Kfact, KFsep
   double precision, dimension(nco) :: KFcr1, Kscr1, KFcr2, Kscr2, KlowT, PHILOGxlowT, KFsep1 ! go in commonS (cannot have n dimension)
   ! pure compound physical constants
   double precision, dimension(n) :: tcn
   double precision, dimension(n) :: pcn
   double precision, dimension(n) :: omgn

   ! eos parameters
   double precision, dimension(n) :: acn  ! in bar*(L/mol)**2
   double precision, dimension(n) :: bn   ! in L/mol
   double precision, dimension(n) :: delta1n  !only required for RKPR
   double precision, dimension(n) :: k_or_mn  ! k for RKPR ; m for SRK/PR

   ! interaction parameters matrices
   double precision, dimension(n, n) :: Kij_or_K0n, Lijn
   double precision, dimension(n, n) :: Tstarn

   ! interaction parameters matrices
   double precision, dimension(nco, nco) :: Kij_or_K0, Lij, Tstar

   ! T, P and Density of the calculated envelope
   double precision, dimension(800) :: Tv
   double precision, dimension(800) :: Pv
   double precision, dimension(800) :: Dv

   ! number of valid elements in To, Po and Do arrays
   integer :: n_points

   ! positions of the last saturation points before each critical point
   integer, dimension(4) :: icri
   ! T, P and Density of critical points
   double precision, dimension(4) :: Tcri
   double precision, dimension(4) :: Pcri
   double precision, dimension(4) :: Dcri

   ! number of valid elements in icri, Tcri, Pcri and Dcri arrays
   integer :: ncri

   double precision, dimension(n) :: y, PHILOGy, PHILOGx
   double precision, dimension(n) :: DLPHITx, DLPHIPx, DLPHITy, DLPHIPy
   double precision, dimension(n, n) :: FUGNx, FUGNy

   type(cross), allocatable :: crossings(:)

   character*4 :: spec
   logical FIRST

   common/MODEL/NMODEL
   common/CRIT/TC(nco), PC(nco), DCeos(nco), omg(nco)
   common/COMPONENTS/ac(nco), b(nco), delta1(nco), rk_or_m(nco), Kij_or_K0, NTDEP
   common/rule/ncomb
   common/bcross/bij(nco, nco)
   common/Tdep/Kinf, Tstar
   common/lforin/lij
   common/DewCurve/ilastDewC, TdewC(800), PdewC(800), dewK(800, nco)
   !common/CrossingPoints/Tcr1, Pcr1, Tcr2, Pcr2, KFcr1, Kscr1, KFcr2, Kscr2
   common/lowTbub/TlowT, PlowT, KlowT, PHILOGxlowT !shared with envelope2
   !common/lowTKsep/KFsep1     !shared with envelope3

   type(envelope) :: dew_envelope, low_t_envelope, high_p_envelope

   Tcr1 = 0.d0 ! T of 1st crossing point detected between different envelope segments
   Tcr2 = 0.d0

   read (1, *) (z(j), j=1, N)
   read (1, *) nmodel
   if (nmodel < 3) then
      call read2PcubicNC(N, 1, 2)
   else if (nmodel == 3) then
      call readRKPRNC(N, 1, 2)
   end if
   write (2, *)
   write (2, 4) (z(i), i=1, n)
4  format('Molar fractions: ', 20F7.4)

   ! Passing values from commons(nco) to input arguments (n)
   TCn = tc(:n)
   PCn = pc(:n)
   OMGn = omg(:n)
   acn = ac(:n)
   bn = b(:n)
   delta1n = delta1(:n)
   k_or_mn = rk_or_m(:n)
   Kij_or_K0n = Kij_or_K0(:n, :n)
   Tstarn = Tstar(:n, :n)
   lijn = lij(:n, :n)

   ! Start from dew curve at low T and P
   ! Although the algorithm is written to start from a bubble point, here we invert the Wilson K factors
   ! and y will actually mean x (ichoice = 2)

   ichoice = 2
   P = 1.0
   T = 315.0
   do while (P > 0.1)
      T = T - 5.D0
      P = 1.d0/sum(z/(PCn*exp(5.373*(1 + omgn)*(1 - TCn/T))))
   end do
   KFACT = PCn*exp(5.373*(1 + omgn)*(1 - TCn/T))/P  ! standard Wilson K factors
   KFACT = 1.d0/KFACT  ! inversion

   call envelope2(ichoice, nmodel, n, z, T, P, KFACT, tcn, pcn, omgn, acn, bn, k_or_mn, delta1n, &
                  Kij_or_K0n, Tstarn, Lijn, n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri, dew_envelope)
   call WriteEnvel(n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)
   TdewC = Tv
   PdewC = Pv
   ilastDewC = n_points

   if (Tcr1 > 0.d0) then  ! self crossing detected (usual in some asymmetric hc mixtures)

   else if (P > Pmax) then  ! now run from Low T Bubble point
      ichoice = 1
      P = 11.0   ! 11.0
      T = 205.0
      do while (P > 10) ! > 10
         T = T - 5.D0
         P = sum(z*PCn*exp(5.373*(1 + omgn)*(1 - TCn/T)))
      end do
      KFACT = PCn*exp(5.373*(1 + omgn)*(1 - TCn/T))/P
      call envelope2(ichoice, nmodel, n, z, T, P, KFACT, tcn, pcn, omgn, acn, bn, k_or_mn, delta1n, &
                     Kij_or_K0n, Tstarn, Lijn, n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri, low_t_envelope)
      call WriteEnvel(n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)

      ! ========================================================================
      !  Crossings finding
      !   To find the crossings between envelopes a generic cross finding
      !   subroutine is called. After the... crosses? have been found,
      !   the K-Factors of each envelope (and each cross) are calculated
      !   with an interpolation based on temperature.
      ! ------------------------------------------------------------------------

      ! Find the cross between two lines (in this case, dew envelope and
      ! (starting from) low temperature bubble envelope
      call find_cross(dew_envelope%t, low_t_envelope%t, &
                      dew_envelope%p, low_t_envelope%p, crossings)

      if (size(crossings) > 1) then
         ! At least two crosses were found

         ! Left Cross (since the search is made from the dew line, the right cross
         !             will be found first, but keeping the 1 and 2 naming for
         !             compatibility, should be fixed after standarizing)
         Tcr1 = crossings(2)%x
         Pcr1 = crossings(2)%y
         icross = crossings(2)%i
         jcross = crossings(2)%j

         ! New Kfactors interpolated for the Left cross
         kfcr1 = kfcross(jcross, low_t_envelope%t, low_t_envelope%logk, Tcr1)
         kscr1 = kfcross(icross, dew_envelope%t, dew_envelope%logk, Tcr1)

         ! Right cross
         Tcr2 = crossings(1)%x
         Pcr2 = crossings(1)%y
         icross = crossings(1)%i
         jcross = crossings(1)%j

         ! New Kfactors interpolated for the right cross
         kfcr2 = kfcross(jcross, low_t_envelope%t, low_t_envelope%logk, Tcr2)
         kscr2 = kfcross(icross, dew_envelope%t, dew_envelope%logk, Tcr2)
      else
         Tcr1 = 0
         Pcr1 = 0
      end if
      ! ========================================================================

   else  ! now run from High P L-L saturation (incipient phase rich in last comp.)
      ichoice = 3
      P = Pmax
      T = 710.0
      iy = 1
      ix = 1
      y = 0.d0
      y(n) = 1.d0
      difw = -1.d0
      ! WRITE (2,*) '  P = ', P
      ! WRITE (2,*) '     T    dif LogFug'
      do while (difW < 0.d0)
         T = T - 10.D0
         call TERMO(n, ix, 1, T, P, z, Vx, PHILOGx, DLPHIPx, DLPHITx, FUGNx) ! for fluid
         call TERMO(n, iy, 1, T, P, y, Vy, PHILOGy, DLPHIPy, DLPHITy, FUGNy) ! for pure Water/Asphaltene
         difW = log(z(n)) + PHILOGx(n) - log(y(n)) - PHILOGy(n)
         !            WRITE (2,*) T, difW
      end do
      KFACT = 1.D-3
      KFACT(n) = 1/z(n)
      call envelope2(ichoice, nmodel, n, z, T, P, KFACT, tcn, pcn, omgn, acn, bn, k_or_mn, delta1n, &
                     Kij_or_K0n, Tstarn, Lijn, n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri, high_p_envelope)
      call WriteEnvel(n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)

      if (P < 1.0) then ! 5/6/22
         ichoice = 1   ! now run from Low T Bubble point, after isolated LL saturation curve
         P = 11.0
         T = 205.0
         do while (P > 10) ! > 10
            T = T - 5.D0
            P = sum(z*PCn*exp(5.373*(1 + omgn)*(1 - TCn/T)))
         end do
         KFACT = PCn*exp(5.373*(1 + omgn)*(1 - TCn/T))/P
         call envelope2(ichoice, nmodel, n, z, T, P, KFACT, tcn, pcn, omgn, acn, bn, k_or_mn, delta1n, &
                        Kij_or_K0n, Tstarn, Lijn, n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri, low_t_envelope)
         call WriteEnvel(n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)

         ! Find cross between high P LL and bubble
         call find_cross(high_p_envelope%t, low_t_envelope%t, &
                         high_p_envelope%p, low_t_envelope%p, crossings)
         Tcr1 = crossings(1)%x
         Pcr1 = crossings(1)%y
         icross = crossings(1)%i
         jcross = crossings(1)%j

         ! New Kfactors interpolated for the Left cross
         kfcr1 = kfcross(jcross, low_t_envelope%t, low_t_envelope%logk, Tcr1)
         kscr1 = kfcross(icross, high_p_envelope%t, high_p_envelope%logk, Tcr1)

         ! Find cross between dew and bubble
         call find_cross(dew_envelope%t, low_t_envelope%t, &
                         dew_envelope%p, low_t_envelope%p, crossings)
         Tcr2 = crossings(1)%x
         Pcr2 = crossings(1)%y
         icross = crossings(1)%i
         jcross = crossings(1)%j

         ! New Kfactors interpolated for the right cross
         kfcr2 = kfcross(jcross, low_t_envelope%t, low_t_envelope%logk, Tcr2)
         kscr2 = kfcross(icross, dew_envelope%t, dew_envelope%logk, Tcr2)
      end if
   end if

   print *, Tcr1, Pcr1, "cross1"
   print *, Tcr2, Pcr2, "cross2"

   if (abs(Tcr1) < 1d-5) then ! no crossings
      print *, "no cross found"
      ichoice = 1   ! first inner curve: 3-phase bubble curve (incipient: V)
      T = TlowT
      P = PlowT
      beta = z(n) ! asphaltenes (or last, separating compound) molar fraction
      KFACT = KlowT(1:n)
      y = 0.d0
      y(n) = 1.d0
      call TERMO(n, 1, 1, T, P, y, Vy, PHILOGy, DLPHIPy, DLPHITy, FUGNy)
      KFsep = exp(PHILOGxlowT(1:n) - PHILOGy)
      call envelope3(ichoice, nmodel, n, z, T, P, beta, KFACT, KFsep, tcn, pcn, omgn, acn, bn, k_or_mn, delta1n, &
                     Kij_or_K0n, Tstarn, Lijn, n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)
      call WriteEnvel(n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)

      ichoice = 2   ! second inner curve: Lower AOP curve (incipient: Asph. or 2nd L)
      T = TlowT
!          beta = sum(z(1:9)) ! volatile compounds until C5  (GENERALIZE LATER as sum of all compounds with Tc<470)
!          s = 2
!          do while (s>1.2)
!              P = P/2
!              beta = 0.8 * beta
!              KFsep = PCn * EXP(5.373*(1+omgn)*(1-TCn/T)) / P  ! standard Wilson K factors
!              xx = z/(1-beta+beta*KFsep)
!              s = sum(xx)
!          end do
!          KFACT(n) = 1/xx(n)
!          P = P1 + (P2-P1) * (1-S1) / (S2-S1)
!          KFACT = KFsep1(1:n) ! from converged first point for 3-phase bubble curve
      P = PlowT/2
      FIRST = .true.
      spec = 'TP'
      call flash(spec, FIRST, nmodel, n, z, tcn, pcn, omgn, acn, bn, k_or_mn, delta1n, &
                 Kij_or_K0n, Tstarn, Lijn, t, p, v, xx, w, rho_x, rho_y, beta, iter)
      call TERMO(n, 1, 1, T, P, xx, Vx, PHILOGx, DLPHIPy, DLPHITy, FUGNy)
      y = 0.d0
      y(n) = 1.d0
      call TERMO(n, 1, 1, T, P, y, Vy, PHILOGy, DLPHIPy, DLPHITy, FUGNy)
      dif = PHILOGy(n) - (PHILOGx(n) + log(xx(n)))
      Pold = P
      if (dif < 0) then
         P = 0.9*P
      else
         P = 1.1*P
      end if
      
!          yn = xx(n)*KFACT(n)
!          do while (yn>2.0)
      do while (abs(dif) > 0.1 .and. P > 0.9)
         call flash(spec, FIRST, nmodel, n, z, tcn, pcn, omgn, acn, bn, k_or_mn, delta1n, &
                    Kij_or_K0n, Tstarn, Lijn, t, p, v, xx, w, rho_x, rho_y, beta, iter)
         call TERMO(n, 1, 1, T, P, xx, Vx, PHILOGx, DLPHIPy, DLPHITy, FUGNy)
         call TERMO(n, 1, 1, T, P, y, Vy, PHILOGy, DLPHIPy, DLPHITy, FUGNy)
         dold = dif
         dif = PHILOGy(n) - (PHILOGx(n) + log(xx(n)))
         aux = P
         P = max(P/10, Pold - dold*(P - Pold)/(dif - dold))
         Pold = aux
      end do
      if (abs(dif) > 0.1) then
         FIRST = .true.
         Told = T
         T = T + 10.0
      end if
      do while (abs(dif) > 0.1)
         call flash(spec, FIRST, nmodel, n, z, tcn, pcn, omgn, acn, bn, k_or_mn, delta1n, &
                    Kij_or_K0n, Tstarn, Lijn, t, p, v, xx, w, rho_x, rho_y, beta, iter)
         call TERMO(n, 1, 1, T, P, xx, Vx, PHILOGx, DLPHIPy, DLPHITy, FUGNy)
         call TERMO(n, 1, 1, T, P, y, Vy, PHILOGy, DLPHIPy, DLPHITy, FUGNy)
         dold = dif
         dif = PHILOGy(n) - (PHILOGx(n) + log(xx(n)))
         aux = T
         T = min(Told - dold*(T - Told)/(dif - dold), T + 20.0)
         Told = aux
      end do
      KFACT = exp(PHILOGx - PHILOGy)
      KFsep = w/xx
      call envelope3(ichoice, nmodel, n, z, T, P, beta, KFACT, KFsep, tcn, pcn, omgn, acn, bn, k_or_mn, delta1n, &
                     Kij_or_K0n, Tstarn, Lijn, n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)
      call WriteEnvel(n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)

   else  ! at least one crossing
      if (ichoice /= 3) ichoice = 1
      print *, "at least one cross", Tcr1, Pcr1
      T = Tcr1
      P = Pcr1
      beta = 0.0d0
      KFACT = exp(KFcr1(:n))
      KFsep = exp(Kscr1(:n))
      call envelope3(ichoice, nmodel, n, z, T, P, beta, KFACT, KFsep, tcn, pcn, omgn, acn, bn, k_or_mn, delta1n, &
                     Kij_or_K0n, Tstarn, Lijn, n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)
      call WriteEnvel(n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)

      if (ichoice /= 3) ichoice = 2
      T = Tcr1
      P = Pcr1
      beta = 0.0d0
      KFACT = exp(Kscr1(:n)) ! now y (incipient phase in envelope3) will be the second liquid
      KFsep = exp(KFcr1(:n)) ! w will be vapor, with phase fraction beta
      call envelope3(ichoice, nmodel, n, z, T, P, beta, KFACT, KFsep, tcn, pcn, omgn, acn, bn, k_or_mn, delta1n, &
                     Kij_or_K0n, Tstarn, Lijn, n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)
      call WriteEnvel(n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)
   end if

   ! THIS IF SHOULD BE INSIDE THE ELSE ABOVE!!
   if (abs(Tcr2) > 1d-5) then
      print *, "second cross", Tcr2, Pcr2
      ichoice = 1
      T = Tcr2
      P = Pcr2
      beta = 0.0d0
      KFACT = exp(KFcr2(:n))
      KFsep = exp(Kscr2(:n))
      call envelope3(ichoice, nmodel, n, z, T, P, beta, KFACT, KFsep, tcn, pcn, omgn, acn, bn, k_or_mn, delta1n, &
                     Kij_or_K0n, Tstarn, Lijn, n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)
      call WriteEnvel(n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)

      ichoice = 2
      T = Tcr2
      P = Pcr2
      beta = 0.0d0
      KFACT = exp(Kscr2(:n)) ! now y (incipient phase in envelope3) will be the second liquid
      KFsep = exp(KFcr2(:n)) ! w will be vapor, with phase fraction beta
      call envelope3(ichoice, nmodel, n, z, T, P, beta, KFACT, KFsep, tcn, pcn, omgn, acn, bn, k_or_mn, delta1n, &
                     Kij_or_K0n, Tstarn, Lijn, n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)
      call WriteEnvel(n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)
   end if

end subroutine readcase

subroutine WriteEnvel(n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)
   use file_operations, only: out_i, outfile

   ! T, P and Density of the calculated envelope
   double precision, dimension(800) :: Tv
   double precision, dimension(800) :: Pv
   double precision, dimension(800) :: Dv

   ! number of valid elements in To, Po and Do arrays
   integer :: n_points

   ! positions of the last saturation points before each critical point
   integer, dimension(4) :: icri

   ! T, P and Density of critical points
   double precision, dimension(4) :: Tcri
   double precision, dimension(4) :: Pcri
   double precision, dimension(4) :: Dcri

   ! number of valid elements in icri, Tcri, Pcri and Dcri arrays
   integer :: ncri

   character(len=200) :: filename

   out_i = out_i + 1
   filename = "envelout"
   filename = outfile(filename, out_i)

   open (unit=out_i, file=filename)

   !WRITE (out_i, *)
   write (out_i, *) '   T(K)        P(bar)        D(mol/L)'
   do i = 1, n_points
      write (out_i, 1) Tv(i), Pv(i), Dv(i)
   end do

1  format(F12.4, 2E14.4, x, I4)
   write (out_i, *)
   write (out_i, *) ' Number of critical points found: ', ncri
   write (out_i, *) '   T(K)        P(bar)        D(mol/L)'
   do i = 1, ncri
      write (out_i, 1) Tcri(i), Pcri(i), Dcri(i), icri(i)
   end do
   close (out_i)
end subroutine WriteEnvel

subroutine envelope2(ichoice, model, n, z, T, P, KFACT, tcn, pcn, omgn, acn, bn, k_or_mn, delta1n, &
                     Kij_or_K0n, Tstarn, Lijn, n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri, &
                     this_envelope)

   use dtypes, only: envelope

   implicit double precision(A - H, O - Z)
   parameter(nco=64)

   ! M&M means the book by Michelsen and Mollerup, 2nd Edition (2007)

   ! eos id, number of compounds in the system and starting point type
   integer, intent(in) :: model, n, ichoice

   double precision Kinf

   ! estimated T and P for first point (then used for every point)
   double precision :: T, P

   ! estimated K factors for first point (then used for every point)
   double precision, dimension(n) :: KFACT
   double precision, dimension(nco) :: KFcr1, Kscr1, KFcr2, Kscr2, KlowT, PHILOGxlowT ! go in commonS (cannot have n dimension)

   ! composition of the system
   double precision, dimension(n), intent(in) :: z

   ! pure compound physical constants
   double precision, dimension(n), intent(in) :: tcn
   double precision, dimension(n), intent(in) :: pcn
   double precision, dimension(n), intent(in) :: omgn

   ! eos parameters
   double precision, dimension(n), intent(in) :: acn  ! in bar*(L/mol)**2
   double precision, dimension(n), intent(in) :: bn   ! in L/mol
   double precision, dimension(n), intent(in) :: delta1n  !only required for RKPR
   double precision, dimension(n), intent(in) :: k_or_mn  ! k for RKPR ; m for SRK/PR

   ! interaction parameters matrices
   double precision, dimension(n, n), intent(in) :: Kij_or_K0n
   double precision, dimension(n, n), intent(in) :: Tstarn
   double precision, dimension(n, n), intent(in) :: Lijn

   ! T, P and Density of the calculated envelope
   double precision, dimension(800), intent(out) :: Tv
   double precision, dimension(800), intent(out) :: Pv
   double precision, dimension(800), intent(out) :: Dv

   ! number of valid elements in Tv, Pv and Dv arrays
   integer, intent(out) :: n_points

   ! positions of the last saturation points before each critical point
   integer, dimension(4), intent(out) :: icri
   ! T, P and Density of critical points
   double precision, dimension(4), intent(out) :: Tcri
   double precision, dimension(4), intent(out) :: Pcri
   double precision, dimension(4), intent(out) :: Dcri

   ! number of valid elements in icri, Tcri, Pcri and Dcri arrays
   integer, intent(out) :: ncri

   ! Intermediate variables during calculation process
   double precision, dimension(n) :: y, PHILOGy, PHILOGx
   double precision, dimension(n) :: DLPHITx, DLPHIPx, DLPHITy, DLPHIPy
   double precision, dimension(n, n) :: FUGNx, FUGNy
   integer, dimension(n + 2) :: ipiv
   double precision, dimension(n + 2) :: X, Xold, Xold2, delX, bd, F, dFdS, dXdS
   double precision, dimension(n + 2, n + 2) :: JAC, AJ
   double precision :: Vy, Vx
   double precision, dimension(2) :: TpairA, PpairA, TpairB, PpairB
   logical :: run, passingcri, Cross, minT, minmaxT

   double precision, dimension(nco, nco) :: Kij_or_K0, Tstar
   common/CRIT/TC(nco), PC(nco), DCeos(nco), omg(nco)
   common/COMPONENTS/ac(nco), b(nco), delta1(nco), rk_or_m(nco), Kij_or_K0, NTDEP
   common/MODEL/NMODEL
   common/rule/ncomb
   common/bcross/bij(nco, nco)
   common/Tdep/Kinf, Tstar
   common/DewCurve/ilastDewC, TdewC(800), PdewC(800), dewK(800, nco)
   common/CrossingPoints/Tcr1, Pcr1, Tcr2, Pcr2, KFcr1, Kscr1, KFcr2, Kscr2
   common/lowTbub/TlowT, PlowT, KlowT, PHILOGxlowT

   real(8) :: tmp_logk(800, n)
   type(envelope) :: this_envelope

   minT = .false.
   minmaxT = .false.
   Told2 = 0.0
   Told = 10.0
   maxP = 0.d0
   ! Charging the commons(nco) from input arguments (n)
   NMODEL = model
   TC(:n) = tcn
   PC(:n) = pcn
   OMG(:n) = omgn
   ac(:n) = acn
   b(:n) = bn
   delta1(:n) = delta1n
   rk_or_m(:n) = k_or_mn
   Kij_or_K0(:n, :n) = Kij_or_K0n
   Kinf = 0.0d0
   ncomb = 0  ! only  vdW combining rules and quadratic mixing rules by  the moment
   Tstar(:n, :n) = Tstarn

   ! b matrix for Classical or van der Waals combining rules:
   do i = 1, n
      do j = i, n
         bij(i, j) = (1 - lijn(i, j))*(b(i) + b(j))/2
         bij(j, i) = bij(i, j)
      end do
   end do
   !
   !-----------------------------------------------------------
   ! Continuation method for tracing the envelope starts here
   run = .true.
   i = 0
   ncri = 0
   JAC(n + 1, :) = 0.d0
   lda = n + 2
   ldb = n + 2
   X(:n) = log(KFACT)
   X(n + 1) = log(T)
   X(n + 2) = log(P)
   iy = 1
   ix = 1
   if (ichoice <= 2) then  ! low T bub (1) or dew (2)
      if (ichoice == 1) iy = -1
      if (ichoice == 2) ix = -1  ! x will be vapor phase during the first part, and liquid after a critical point is crossed
      ns = n + 1
      S = log(T)
      delS = 0.005
      y = KFACT*z ! Wilson estimate for vapor (or liquid) composition
   else    ! (ichoice==3) high P L-L sat
      PmaxDewC = maxval(PdewC(1:ilastDewC))
      ns = n + 2
      S = log(P)
      delS = -0.005
      y = 0.d0
      y(n) = 1.d0
   end if
   Xold = 0.d0
   dFdS = 0.d0
   dFdS(n + 2) = -1.d0

   do while (run)
      i = i + 1
      ! Newton starts here
      delX = 1.0
      iter = 0
      max_iter = 100
      do while (maxval(abs(delX)) > 1.d-5 .and. iter <= max_iter)
         iter = iter + 1
         !      nc,MTYP,INDIC,T,P,rn,V,PHILOG,DLPHI,DLPHIP,DLPHIT,FUGN
         call TERMO(n, iy, 4, T, P, y, Vy, PHILOGy, DLPHIPy, DLPHITy, FUGNy)
         call TERMO(n, ix, 2, T, P, z, Vx, PHILOGx, DLPHIPx, DLPHITx, FUGNx)
         F(:n) = X(:n) + PHILOGy - PHILOGx  ! X(:n) are LOG_K
         F(n + 1) = sum(y - z)
         F(n + 2) = X(ns) - S
         do j = 1, n
            JAC(1:n, j) = FUGNy(:, j)*y(j)  ! z*K=y
            JAC(j, j) = JAC(j, j) + 1.d0
         end do
         JAC(1:n, n + 1) = T*(DLPHITy - DLPHITx)
         JAC(1:n, n + 2) = P*(DLPHIPy - DLPHIPx)
         JAC(n + 1, 1:n) = y
         JAC(n + 2, :) = 0.d0
         JAC(n + 2, ns) = 1.d0
         ! call dgesv( n, nrhs, a, lda, ipiv, b, ldb, info )
         bd = -F
         AJ = JAC
         call dgesv(n + 2, 1, AJ, lda, ipiv, bd, ldb, info)
         if (info .ne. 0) then
            print *, "error with dgesv in parameter ", info, "540"
         end if
         delX = bd
         if (i == 1) then
            do while (maxval(abs(delX)) > 5.0)   ! Too large Newton step --> Reduce it
               delX = delX/2
            end do
         else
            do while (maxval(abs(delX)) > 0.08)   ! Too large Newton step --> Reduce it
               delX = delX/2
            end do
            if (iter > 8) delX = delX/2  ! too many iterations (sometimes due to oscillatory behavior near crit) --> Reduce it
         end if

         X = X + delX

         if (.not. passingcri .and. i /= 1 .and. iter > 20 .and. maxval(abs(delX)) > 0.001) then ! Too many iterations-->Reduce step to new point
            delS = delS*3.0/4.0
            S = S - delS
            X = Xold + dXdS*delS
            !iter = 0
         end if

         KFACT = exp(X(:n))
         y = z*KFACT
         T = exp(X(n + 1))
         P = exp(X(n + 2))
      end do
      ! Point converged (unless it jumped out because of high number of iterations)
      if (ichoice == 2) dewK(i, :n) = X(:n)
      if (ichoice == 1 .and. i == 1) then ! save for starting later a 3-envel
         TlowT = T
         PlowT = P
         KlowT(1:n) = KFACT
         PHILOGxlowT(1:n) = PHILOGx
      end if
      if (iter > max_iter) run = .false.
      if (P > maxP) maxP = P

      if (ichoice == 2 .and. i > 1) then
         ! TODO: If this is the way the low p dew line finishes, I think this should be more strict
         if (P < Pv(i - 1) .and. P < maxP/5 .and. T > 300) then
            run = .false.  ! to finish envelope going to low T bubble
         end if
      end if

      print *, T, P, ns, iter
      Tv(i) = T
      Pv(i) = P
      Dv(i) = 1/Vx    ! saturated phase density

      tmp_logk(i, :n) = X(:n)
      ! rho_y = 1/Vy     incipient phase density

      if (ichoice == 3 .and. P < 1.0) then
         run = .false.   ! isolated LL line detected. Stop and start a new one from low T false bubble point
      end if

      if (sum(X(:n)*Xold(:n)) < 0) then  ! critical point detected
         ncri = ncri + 1
         icri(ncri) = i - 1
         frac = -Xold(ns)/(X(ns) - Xold(ns))
         Tcri(ncri) = Tv(i - 1) + frac*(T - Tv(i - 1))
         Pcri(ncri) = Pv(i - 1) + frac*(P - Pv(i - 1))
         Dcri(ncri) = Dv(i - 1) + frac*(Dv(i) - Dv(i - 1))
         iy = -iy
         ix = -yx
      end if

      if (run) then
         ! Calculation of sensitivities (dXdS)
         ! dgesv( n, nrhs, a, lda, ipiv, b, ldb, info )
         bd = -dFdS
         AJ = JAC
         call dgesv(n + 2, 1, AJ, lda, ipiv, bd, ldb, info)
         if (info .ne. 0) then
            print *, "error with dgesv in parameter ", info, "760"
         end if
         dXdS = bd
         ! Selection of (the most changing) variable to be specified for the next point
         nsold = ns
         ns = maxloc(abs(dXdS), DIM=1)
         if (maxval(abs(X(:n))) < 0.2) then
            ns = maxloc(abs(dXdS(1:n)), DIM=1)  ! T and P not allowed to be chosen close to a critical point
         end if
         if (ns /= nsold) then
            delS = dXdS(ns)*delS    ! translation of delS to the  new specification variable
            dXdS = dXdS/dXdS(ns)  ! translation of sensitivities
            S = X(ns)               ! update of S
         end if
         ! Setting step in S for the next point to be calculated
         delmax = max(sqrt(abs(X(ns)))/10, 0.1)
         updel = delS*3/iter
         if (passingcri) updel = delS
         if (delS > 0) then
            delS = min(updel, delmax)
         else
            delS = max(updel, -delmax)
         end if
         S = S + delS
         ! Generation of estimates for the next point
         Told2 = Told
         Told = T
         Xold2 = Xold
         Xold = X
         X = Xold + dXdS*delS

         black_i = 0
         step_fixer = 1.d0
         if (passingcri) passingcri = .false.

         do while (maxval(abs(X(:n))) < 0.03)  ! approaching the black hole... get out of there! (0.03)
            black_i = black_i + 1
            if (black_i > 50) then
               print *, "Stuck on the black hole 1"
               step_fixer = 3.d0
            end if
            stepX = maxval(abs(X(:n) - Xold(:n))) ! the step given by the most changing logK to fall into the black hole
            passingcri = .true.
            if (stepX > 0.07) then
               S = S - delS/2
               X = X - dXdS*delS/2   !  half step back
            else
               S = S + delS
               X = X + step_fixer*dXdS*delS   ! one more step to jump over the critical point
            end if
         end do

         T = exp(X(n + 1))
         if (.not. passingcri .and. abs(T - Told) > 7) then ! Delta T estimations > 7K are not allowed
            delS = delS/2
            S = S - delS
            X = Xold + dXdS*delS
            T = exp(X(n + 1))
         end if

         P = exp(X(n + 2))
         KFACT = exp(X(:n))
         y = z*KFACT
         if ((dXdS(n + 1)*delS < 0 .and. P < 0.1 .or. T < 120.0) &  ! dew line stops when P<0.1 bar or T<150K
             .or. (P > 1.0 .and. T < 150.0) &   ! bubble line stops when T<150K
             .or. (P > 1500) &
             .or. (abs(dels) < 1.d-4)) then
            run = .false.
         end if
      end if
   end do

   n_points = i
   !-----------------------------------------------------------

   ! Define envelope values, omit the last point to avoid not really
   ! converged cases
   this_envelope%logk = tmp_logk(:n_points - 1, :)
   this_envelope%t = Tv(:n_points - 1)
   this_envelope%p = Pv(:n_points - 1)

   print *, y
   print *, rho_x
   print *, rho_y
   print *, beta
end subroutine envelope2

subroutine envelope3(ichoice, model, n, z, T, P, beta, KFACT, KFsep, tcn, pcn, omgn, acn, bn, k_or_mn, delta1n, &
                     Kij_or_K0n, Tstarn, Lijn, n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)

! Routine for tracing the boundaries of a three-phase region, starting from a known previously detected point.
! Developed as an extension and modification of "envelope2". March 2017.

!   In the main case (ichoice=1) y, xx and w denote molar fractions in the following phases:
!   y: vapor (incipient phase)
!  xx: hydrocarbon liquid
!   w: the other liquid (e.g. aqueous of asphaltenic phase)
!   beta is the fraction of w phase. Then we have (1-beta) of xx phase and 0 of y (incipient phase)
!   KFACT = y/xx
!   KFsep = w/xx

!   ichoice=2 is used for cases where the incipient phase "y" is the 2nd liquid.
!   w will correspond to vapor phase, with fraction beta

!   ichoice=3 is used for cases where the initial saturated phase "xx" is the vapor (e.g. OilB with water from Lindeloff-Michelsen).
!   y and w will correspond to the two liquid phases, with a beta fraction for w.

   implicit double precision(A - H, O - Z)
   parameter(nco=64)

   ! M&M means the book by Michelsen and Mollerup, 2nd Edition (2007)

   ! eos id, number of compounds in the system and starting point type
   integer, intent(in) :: model, n, ichoice

   double precision Kinf

   ! estimated T, P and "w" phase fraction for first point (then used for every point)
   double precision :: T, P, beta

   ! estimated K factors for first point (then used for every point)
   double precision, dimension(n) :: KFACT, KFsep

   ! composition of the system
   double precision, dimension(n), intent(in) :: z

   ! pure compound physical constants
   double precision, dimension(n), intent(in) :: tcn
   double precision, dimension(n), intent(in) :: pcn
   double precision, dimension(n), intent(in) :: omgn

   ! eos parameters
   double precision, dimension(n), intent(in) :: acn  ! in bar*(L/mol)**2
   double precision, dimension(n), intent(in) :: bn   ! in L/mol
   double precision, dimension(n), intent(in) :: delta1n  !only required for RKPR
   double precision, dimension(n), intent(in) :: k_or_mn  ! k for RKPR ; m for SRK/PR

   ! interaction parameters matrices
   double precision, dimension(n, n), intent(in) :: Kij_or_K0n
   double precision, dimension(n, n), intent(in) :: Tstarn
   double precision, dimension(n, n), intent(in) :: Lijn

   ! T, P and Density of the calculated envelope
   double precision, dimension(800), intent(out) :: Tv
   double precision, dimension(800), intent(out) :: Pv
   double precision, dimension(800), intent(out) :: Dv

   ! number of valid elements in Tv, Pv and Dv arrays
   integer, intent(out) :: n_points

   ! positions of the last saturation points before each critical point
   integer, dimension(4), intent(out) :: icri
   ! T, P and Density of critical points
   double precision, dimension(4), intent(out) :: Tcri
   double precision, dimension(4), intent(out) :: Pcri
   double precision, dimension(4), intent(out) :: Dcri

   ! number of valid elements in icri, Tcri, Pcri and Dcri arrays
   integer, intent(out) :: ncri

   ! Intermediate variables during calculation process
   double precision, dimension(n) :: y, xx, w, PHILOGy, PHILOGx, PHILOGw
   double precision, dimension(n) :: dxdB, dydB, dwdB, dxdKs, dydKs, dwdKs, aux
   double precision, dimension(n) :: DLPHITx, DLPHIPx, DLPHITy, DLPHIPy, DLPHITw, DLPHIPw
   double precision, dimension(n, n) :: FUGNx, FUGNy, FUGNw
   integer, dimension(2*n + 3) :: ipiv
   double precision, dimension(2*n + 3) :: X, Xold, delX, bd, F, dFdS, dXdS
   double precision, dimension(2*n + 3) :: Fp, JACnumK2, JACnumKn, JACnumKs2, JACnumKsn, JACnumlT, JACnumlP, JACnumB
   double precision, dimension(2*n + 3, 2*n + 3) :: JAC, AJ
   double precision :: Vy, Vx, Vw
!        DOUBLE PRECISION, dimension(2) :: TpairA, PpairA, TpairB, PpairB    ! crossing vars
   logical :: run, passingcri, Comp3ph  !, Cross                       ! crossing var (Cross)

   double precision, dimension(nco, nco) :: Kij_or_K0, Tstar
   double precision, dimension(nco) :: KFsep1
   common/CRIT/TC(nco), PC(nco), DCeos(nco), omg(nco)
   common/COMPONENTS/ac(nco), b(nco), delta1(nco), rk_or_m(nco), Kij_or_K0, NTDEP
   common/MODEL/NMODEL
   common/rule/ncomb
   common/bcross/bij(nco, nco)
   common/Tdep/Kinf, Tstar
   common/lowTKsep/KFsep1
   common/writeComp/Comp3ph, i1, i2
   ! common /DewCurve/ ilastDewC, TdewC(500), PdewC(500)     ! crossing vars
   ! common /CrossingPoints/ Tcr1,Pcr1,Tcr2,Pcr2,KFcr1,Kscr1,KFcr2,Kscr2

! Charging the commons(nco) from input arguments (n)
   NMODEL = model
   TC(:n) = tcn
   PC(:n) = pcn
   OMG(:n) = omgn
   ac(:n) = acn
   b(:n) = bn
   delta1(:n) = delta1n
   rk_or_m(:n) = k_or_mn
   Kij_or_K0(:n, :n) = Kij_or_K0n
   Kinf = 0.0d0
   ncomb = 0  ! only  vdW combining rules and quadratic mixing rules by  the moment
   Tstar(:n, :n) = Tstarn

   ! b matrix for Classical or van der Waals combining rules:
   do i = 1, n
      do j = i, n
         bij(i, j) = (1 - lijn(i, j))*(b(i) + b(j))/2
         bij(j, i) = bij(i, j)
      end do
   end do
!
   !-----------------------------------------------------------

   ! Continuation method for tracing the envelope starts here
   run = .true.
   i = 0
   ncri = 0
   JAC = 0.d0
   lda = 2*n + 3
   ldb = 2*n + 3
   X(:n) = log(KFACT)
   X(n + 1:2*n) = log(KFsep)
   X(2*n + 1) = log(T)
   X(2*n + 2) = log(P)
   X(2*n + 3) = beta

!
   ! To test Jacobian insert here block for numerical derivatives at the end of this code
!

   iy = 1
   ix = 1
   iw = 1

   if (ichoice == 1) iy = -1
   if (ichoice == 2) iw = -1  ! w will be vapor phase during the first part
   if (ichoice == 3) ix = -1  ! x will be vapor phase, saturated in the first point

   if (beta == 0.0d0) then
      ns = 2*n + 3
      S = 0.00   ! for beta
      delS = 0.001
   else ! start from low T "bubble" point or low T "Lower AOP"
      ns = 2*n + 1
      S = X(2*n + 1)   ! for log(T)
      delS = 0.001
   end if

   xx = z/(1 - beta + beta*KFsep)
   y = KFACT*xx
   w = KFsep*xx
   Xold = 0.d0
   dFdS = 0.d0
   dFdS(2*n + 3) = -1.d0
   write (2, *) '     T       P      beta     X(1)     X(n)     X(n+1)     X(2*n)    ns  iter'
   if (Comp3ph) write (3, *) '     T       P      beta     xa     xb     ya     yb    wa     wb'

   do while (run)
      i = i + 1  ! number of point to be calculated along the line
      ! Newton starts here
      delX = 1.0
      iter = 0
      max_iter = 100
      reps = 0

      do while (maxval(abs(delX)) > 1.d-5 .and. iter <= max_iter)
         iter = iter + 1
         ! nc,MTYP,INDIC,T,P,rn,V,PHILOG,DLPHI,DLPHIP,DLPHIT,FUGN
         call TERMO(n, iy, 4, T, P, y, Vy, PHILOGy, DLPHIPy, DLPHITy, FUGNy)
         call TERMO(n, ix, 4, T, P, xx, Vx, PHILOGx, DLPHIPx, DLPHITx, FUGNx)
         call TERMO(n, iw, 4, T, P, w, Vw, PHILOGw, DLPHIPw, DLPHITw, FUGNw)
         F(:n) = X(:n) + PHILOGy - PHILOGx  ! X(:n) are LOG_K
         F(n + 1:2*n) = X(n + 1:2*n) + PHILOGw - PHILOGx  ! X(:n) are LOG_K
         F(2*n + 1) = sum(y - xx)
         F(2*n + 2) = sum(w - xx)
         F(2*n + 3) = X(ns) - S
         dxdB = -(KFsep - 1.d0)*xx*xx/z
         dydB = KFACT*dxdB
         dwdB = KFsep*dxdB
         aux = -beta*xx/z
         dxdKs = aux*xx
         dydKs = aux*y
         dwdKs = xx*(1 + aux*KFsep)

         do j = 1, n
            JAC(1:n, j) = FUGNy(:, j)*y(j)  ! y=K*xx
            JAC(j, j) = JAC(j, j) + 1.d0
         end do

         do j = n + 1, 2*n    ! wrt Ks
            JAC(1:n, j) = KFsep(j - n)*(FUGNy(:, j - n)*dydKs(j - n) - FUGNx(:, j - n)*dxdKs(j - n))
         end do

         JAC(1:n, 2*n + 1) = T*(DLPHITy - DLPHITx)  ! wrt T
         JAC(1:n, 2*n + 2) = P*(DLPHIPy - DLPHIPx)  ! wrt P

         do l = 1, n    ! wrt beta
            JAC(l, 2*n + 3) = sum(FUGNy(l, :)*dydB - FUGNx(l, :)*dxdB)
         end do

         ! ders of F(n+1:2*n) wrt K = 0

         do j = n + 1, 2*n  ! wrt Ks
            JAC(n + 1:2*n, j) = KFsep(j - n)*(FUGNw(:, j - n)*dwdKs(j - n) - FUGNx(:, j - n)*dxdKs(j - n))
            JAC(j, j) = JAC(j, j) + 1.d0
         end do

         JAC(n + 1:2*n, 2*n + 1) = T*(DLPHITw - DLPHITx)  ! wrt T
         JAC(n + 1:2*n, 2*n + 2) = P*(DLPHIPw - DLPHIPx)  ! wrt P

         do l = n + 1, 2*n    ! wrt beta
            JAC(l, 2*n + 3) = sum(FUGNw(l - n, :)*dwdB - FUGNx(l - n, :)*dxdB)
         end do

         JAC(2*n + 1, 1:n) = y             ! sum(y-x) wrt K
         JAC(2*n + 1, n + 1:2*n) = KFsep*(dydKs - dxdKs)    ! sum(y-x) wrt Ks
         ! sum(y-x) wrt T or P = 0
         JAC(2*n + 1, 2*n + 3) = sum(dydB - dxdB)  ! wrt beta
         ! sum(w-x) wrt K = 0
         JAC(2*n + 2, n + 1:2*n) = KFsep*(dwdKs - dxdKs)    ! sum(w-x) wrt Ks
         ! sum(w-x) wrt T or P = 0
         JAC(2*n + 2, 2*n + 3) = sum(dwdB - dxdB)  ! wrt beta
         JAC(2*n + 3, :) = 0.d0
         JAC(2*n + 3, ns) = 1.d0

         bd = -F
         AJ = JAC

         call dgesv(2*n + 3, 1, AJ, lda, ipiv, bd, ldb, info)

         if (info .ne. 0) then
            print *, "error with dgesv in parameter ", info
            print *, "error at", T, P
         end if

         delX = bd

         if (i == 1) then
            do while (maxval(abs(delX)) > 1.0)   ! Too large Newton step --> Reduce it
               delX = delX/2
            end do
         else
            do while (maxval(abs(delX)) > 0.08)   ! Too large Newton step --> Reduce it
               delX = delX/2
            end do
            if (iter > 10) delX = delX/2  ! too many iterations (sometimes due to oscillatory behavior near crit) --> Reduce it
         end if

         X = X + delX

         if (.not. passingcri .and. i /= 1 .and. iter > 20 .and. maxval(abs(delX)) > 0.001) then ! Too many iterations--> Reduce step to new point
            delS = delS*3.d0/4.d0
            S = S - delS
            X = Xold + dXdS*delS
         end if

         KFACT = exp(X(:n))
         KFsep = exp(X(n + 1:2*n))
         T = exp(X(2*n + 1))
         P = exp(X(2*n + 2))
         beta = X(2*n + 3)
         xx = z/(1 - beta + beta*KFsep)
         y = KFACT*xx
         w = KFsep*xx
      end do

      ! Point converged (unless it jumped out because of high number of iterations)
      if (iter > max_iter) run = .false.
      if (beta < 0) run = .false.
      if (ichoice == 1 .and. i == 1) KFsep1(1:n) = KFsep

      print *, T, P, ns, iter
      write (2, 1) T, P, beta, X(1), X(n), X(n + 1), X(2*n), ns, iter
      if (Comp3ph) write (3, 3) T, P, beta, xx(i1), xx(i2), y(i1), y(i2), w(i1), w(i2)
      Tv(i) = T
      Pv(i) = P
      Dv(i) = 1/Vx    ! saturated phase density
      rho_x = 1/Vx
      rho_y = 1/Vy    ! incipient phase density
      if (sum(X(:n)*Xold(:n)) < 0) then  ! critical point detected between x and y phases
         ncri = ncri + 1
         icri(ncri) = i - 1
         frac = -Xold(ns)/(X(ns) - Xold(ns))
         Tcri(ncri) = Tv(i - 1) + frac*(T - Tv(i - 1))
         Pcri(ncri) = Pv(i - 1) + frac*(P - Pv(i - 1))
         Dcri(ncri) = Dv(i - 1) + frac*(Dv(i) - Dv(i - 1))
         iy = -iy
         ix = -yx
      end if
      if (sum(X(n + 1:2*n)*Xold(n + 1:2*n)) < 0) then  ! critical point detected between x and w phases
         ncri = ncri + 1
         icri(ncri) = i - 1
         frac = -Xold(ns)/(X(ns) - Xold(ns))
         Tcri(ncri) = Tv(i - 1) + frac*(T - Tv(i - 1))
         Pcri(ncri) = Pv(i - 1) + frac*(P - Pv(i - 1))
         Dcri(ncri) = Dv(i - 1) + frac*(Dv(i) - Dv(i - 1))
         iw = -iw
         ix = -yx
      end if
      if (run) then
         ! Calculation of sensitivities (dXdS)
         ! dgesv( n, nrhs, a, lda, ipiv, b, ldb, info )
         bd = -dFdS
         AJ = JAC
         call dgesv(2*n + 3, 1, AJ, lda, ipiv, bd, ldb, info)
         if (info .ne. 0) then
            print *, "error with dgesv in parameter ", info, "1125"
         end if
         dXdS = bd
         ! Selection of (the most changing) variable to be specified for the next point
         nsold = ns
         ns = maxloc(abs(dXdS), DIM=1)
         if (maxval(abs(X(:n))) < 0.2) then
            ns = maxloc(abs(dXdS(1:n)), DIM=1)  ! vars other than logK not allowed to be chosen close to a y-x critical point
         end if
         if (maxval(abs(X(n + 1:2*n))) < 0.2) then
            ns = maxloc(abs(dXdS(n + 1:2*n)), DIM=1)  ! vars other than logKs not allowed to be chosen close to a w-x critical point
         end if
         if (ns /= nsold) then
            delS = dXdS(ns)*delS    ! translation of delS to the  new specification variable
            dXdS = dXdS/dXdS(ns)    ! translation of sensitivities
            S = X(ns)               ! update of S
         end if
         ! Setting step in S for the next point to be calculated
         delmax = max(sqrt(abs(X(ns)))/10, 0.1)
         updel = delS*3/iter
         if (passingcri) updel = delS
         if (delS > 0) then
            delS = min(updel, delmax)
         else
            delS = max(updel, -delmax)
         end if
         S = S + delS
         ! Generation of estimates for the next point
         Told = T
         Xold = X
         X = Xold + dXdS*delS
         if (passingcri) passingcri = .false.

         ! do while (maxval(abs(X(:n)))<0.03)  ! approaching the black hole... get out of there! (0.03)
         !       stepX = maxval(abs(X(:n)-Xold(:n))) ! the step given by the most changing logK to fall into the black hole
         !       passingcri = .true.
         !       if (stepX > 0.07) then
         !           S = S - delS/2
         !           X = X - dXdS * delS/2   !  half step back
         !       else
         !           S = S + delS
         !           X = X + dXdS * delS   ! one more step to jump over the critical point
         !       end if
         ! end do

         if (maxval(abs(X(:n))) < 0.03) then  ! approaching the black hole... get out of there!
            if (delS > 0) delS = 0.04 - S
            if (delS < 0) delS = -0.04 - S
            S = S + delS
            X = X + dXdS*delS   ! one more step to jump over the critical point
            passingcri = .true.
         end if

         black_i = 0

         do while (maxval(abs(X(n + 1:2*n))) < 0.03)  ! approaching the black hole... get out of there! (0.03)
            black_i = black_i + 1
            if (black_i > 10) then
               print *, "Stuck on the black hole 2"
            end if
            stepX = maxval(abs(X(n + 1:2*n) - Xold(n + 1:2*n))) ! the step given by the most changing logKs to fall into the black hole
            passingcri = .true.
            if (stepX > 0.07) then
               S = S - delS/2
               X = X - dXdS*delS/2   !  half step back
            else
               S = S + delS
               X = X + dXdS*delS   ! one more step to jump over the critical point
            end if
         end do

         T = exp(X(2*n + 1))

         if (.not. passingcri .and. abs(T - Told) > 7) then ! Delta T estimations > 7K are not allowed
            delS = delS/2
            S = S - delS
            X = Xold + dXdS*delS
            T = exp(X(2*n + 1))
         end if
         P = exp(X(2*n + 2))
         KFACT = exp(X(:n))
         KFsep = exp(X(n + 1:2*n))
         beta = X(2*n + 3)
         xx = z/(1 - beta + beta*KFsep)
         y = KFACT*xx
         w = KFsep*xx
         if ((dXdS(2*n + 1)*delS < 0 .and. P < 0.1 .or. T < 120.0) &  ! this may need adjustment
             .or. (P > 1.0 .and. T < 150.0) &
             .or. (P > 1500) &
             .or. (abs(dels) < 1.d-4)) then
            run = .false.
         end if
      end if
   end do
   n_points = i
   !-----------------------------------------------------------

   print *, y
   print *, rho_x
   print *, rho_y
   print *, beta
1  format(7F10.4, 2I4)
3  format(3F10.4, 6E12.3)
end subroutine envelope3

subroutine EvalFEnvel3(n, z, X, F)

   implicit double precision(A - H, O - Z)
   double precision, dimension(n) :: KFACT, KFsep
   double precision, dimension(n) :: z, y, xx, w, PHILOGy, PHILOGx, PHILOGw
   double precision, dimension(n) :: DLPHITx, DLPHIPx, DLPHITy, DLPHIPy, DLPHITw, DLPHIPw
   double precision, dimension(n, n) :: FUGNx, FUGNy, FUGNw
   double precision, dimension(2*n + 3) :: X, F

   S = 0.001
   KFACT = exp(X(:n))
   KFsep = exp(X(n + 1:2*n))
   T = exp(X(2*n + 1))
   P = exp(X(2*n + 2))
   beta = X(2*n + 3)
   xx = z/(1 - beta + beta*KFsep)
   y = KFACT*xx
   w = KFsep*xx
   call TERMO(n, 1, 4, T, P, y, Vy, PHILOGy, DLPHIPy, DLPHITy, FUGNy)
   call TERMO(n, -1, 4, T, P, xx, Vx, PHILOGx, DLPHIPx, DLPHITx, FUGNx)
   call TERMO(n, 1, 4, T, P, w, Vw, PHILOGw, DLPHIPw, DLPHITw, FUGNw)
   F(:n) = X(:n) + PHILOGy - PHILOGx  ! X(:n) are LOG_K
   F(n + 1:2*n) = X(n + 1:2*n) + PHILOGw - PHILOGx  ! X(:n) are LOG_K
   F(2*n + 1) = sum(y - xx)
   F(2*n + 2) = sum(w - xx)
   F(2*n + 3) = X(25) - S
end subroutine EvalFEnvel3

!        eps = 1.d-6
!        call EvalFEnvel3(n,z,X,F)
!        x(2) = x(2)+eps
!        call EvalFEnvel3(n,z,X,Fp)
!        x(2) = x(2)-eps
!        JACnumK2 = (Fp-F)/eps
!        write(2,*)'JACnumK2: ',JACnumK2
!
!        x(n) = x(n)+eps
!        call EvalFEnvel3(n,z,X,Fp)
!        x(n) = x(n)-eps
!        JACnumKn = (Fp-F)/eps
!        write(2,*)'JACnumKn: ',JACnumKn
!
!        x(n+2) = x(n+2)+eps
!        call EvalFEnvel3(n,z,X,Fp)
!        x(n+2) = x(n+2)-eps
!        JACnumKs2 = (Fp-F)/eps
!        write(2,*)'JACnumKs2: ',JACnumKs2
!
!        x(2*n) = x(2*n)+eps
!        call EvalFEnvel3(n,z,X,Fp)
!        x(2*n) = x(2*n)-eps
!        JACnumKsn = (Fp-F)/eps
!        write(2,*)'JACnumKsn: ',JACnumKsn
!
!        x(2*n+1) = x(2*n+1)+eps
!        call EvalFEnvel3(n,z,X,Fp)
!        x(2*n+1) = x(2*n+1)-eps
!        JACnumlT = (Fp-F)/eps
!        write(2,*)'JACnumlT: ',JACnumlT
!
!        x(2*n+2) = x(2*n+2)+eps
!        call EvalFEnvel3(n,z,X,Fp)
!        x(2*n+2) = x(2*n+2)-eps
!        JACnumlP = (Fp-F)/eps
!        write(2,*)'JACnumlP: ',JACnumlP
!
!        x(2*n+3) = x(2*n+3)+eps
!        call EvalFEnvel3(n,z,X,Fp)
!        x(2*n+3) = x(2*n+3)-eps
!        JACnumB = (Fp-F)/eps
!        write(2,*)'JACnumB: ',JACnumB

