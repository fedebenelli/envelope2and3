module nfile
   implicit none
   integer :: xit_file = 0
   character(len=250) :: filename
end module nfile


program calc_envelope2and3
   use constants
   use io_nml, only: setup_input, read_system
   use system, only: z, nmodel => thermo_model, n => nc
   use io, only: str
   
   implicit none

   common/writeComp/Comp3ph, i1, i2

   logical :: Comp3ph = .false.

   ! cli args
   character(len=254) :: infile
   character(len=3) :: from_nml
   character(len=3) :: three_phase_arg
   
   logical :: three_phase = .false.

   ! Simple benchmark
   real :: start_time, end_time

   integer :: i, j
   integer :: i1, i2

   call execute_command_line("mkdir -p env23out && rm env23out/envelout*")
   from_nml=""
   call get_command_argument(1, infile)
   call get_command_argument(2, three_phase_arg)
   call get_command_argument(3, from_nml)

   open (2, FILE='envelOUT.txt')
   if (from_nml == "yes") then
      ! NML Based IO
      call read_system(infile)
   else
      ! Modified Legacy IO
      open (1, FILE=infile)! 'envelIN.txt')
      read (1, *) N

      allocate(z(n))
      read (1, *) (z(j), j=1, N)
      read (1, *) nmodel

      if (nmodel < 4) then
         call read2PcubicNC(N, 1, 2)
      else if (nmodel == 4) then
         call readRKPRNC(N, 1, 2)
      end if
   end if

   write (2, *)
   write (2, 4) (z(i), i=1, n)

   ! Normalize compositions
   z = z/sum(z)

   call cpu_time(start_time)
   call readcase(n, three_phase_arg)
   call cpu_time(end_time)

   print *, "Finished in " // str(end_time - start_time) // "seconds"

4  format('Molar fractions: ', 20F7.4)
end program

subroutine readcase(n, three_phase)
   use constants
   use system, only: z, nmodel => thermo_model, tc, pc, dceos => dc, omg => w, &
                     ac, b, delta1 => del1, rk_or_m => k, kij_or_k0 => kij, &
                     ntdep => tdep, ncomb => mixing_rule, bij, kinf, tstar, lij
   use dtypes, only: envelope, kfcross, point, print_header, env3, find_cross, find_self_cross
   use array_operations, only: diff
   use envelopes, only: max_points

   implicit none
   
   integer, intent(in) :: n
   character(len=3), intent(in) :: three_phase ! Calculate 3-phase lines?

   integer, parameter :: nco=64
   real(pr), parameter :: Pmax=700

   real(pr) :: xx(n), w(n), y(n)

   real(pr), dimension(n) :: Kfact, KFsep
   real(pr), dimension(nco) :: KFcr1, Kscr1, KFcr2, Kscr2, & 
                               KlowT, PHILOGxlowT ! go in commons (cannot have n dimension)
   
   ! T, P and Density of the calculated envelope
   real(pr) :: T, P, V
   real(pr) :: Tv(max_points)
   real(pr) :: Pv(max_points)
   real(pr) :: Dv(max_points)

   ! number of valid elements in To, Po and Do arrays
   integer :: n_points

   ! positions of the last saturation points before each critical point
   integer, dimension(4) :: icri

   ! T, P and Density of critical points
   real(pr), dimension(4) :: Tcri, Pcri, Dcri

   ! number of valid elements in icri, Tcri, Pcri and Dcri arrays
   integer :: ncri

   ! phases variables
   real(pr) :: Vx, PHILOGx(n), DLPHITx(n), DLPHIPx(n), FUGNx(n, n)
   real(pr) :: Vy, PHILOGy(n), DLPHITy(n), DLPHIPy(n), FUGNy(n, n)
   real(pr) :: rho_x
   real(pr) :: rho_y
   integer :: ix
   integer :: iy

   ! simple benchmark
   real(pr) :: start_time, end_time

   ! Flash
   character(len=4) :: spec
   logical :: first
   integer :: iter
   
   ! Envelopes
   real(pr) :: Tcr1, Tcr2, Pcr1, Pcr2
   integer :: ichoice
   type(point), allocatable :: crossings(:)
   integer :: icross, jcross
   type(envelope) :: dew_envelope, low_t_envelope, high_p_envelope
   type(env3) :: triphasic
   real(pr) :: beta

   logical :: highPLL_converged

   real(pr) :: aux, difw, dif, Pold, Told, dold, extrapolation

   ! Iterations variables
   integer :: i, j

   Tcr1 = 0.d0 ! T of 1st crossing point detected between different envelope segments
   Tcr2 = 0.d0

   ! Start from dew curve at low T and P
   ! Although the algorithm is written to start from a bubble point,
   ! here we invert the Wilson K factors
   ! and y will actually mean x (ichoice = 2)
   print *, "Running Dew Line"
   call cpu_time(start_time)
   ichoice = 2
   P = 1.0
   T = 315.0

   do while (P > 0.1)
      T = T - 5.D0
      P = 1.d0/sum(z/(pc*exp(5.373*(1 + omg)*(1 - tc/T))))
   end do

   KFACT = pc*exp(5.373*(1 + omg)*(1 - tc/T))/P  ! standard Wilson K factors
   KFACT = 1.d0/KFACT  ! inversion

   call envelope2(&
      ichoice, n, z, T, P, KFACT, &
      n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri, &
      dew_envelope)
   
   call cpu_time(end_time)
   call WriteEnvel(n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)
   call dew_envelope%write("env23out/envelout2-DEW")
   print *, "Done in: ", end_time-start_time


   call find_self_cross(dew_envelope%t, dew_envelope%p, crossings)
   if (allocated(crossings) .and. size(crossings) > 0) then
      ! self crossing detected (usual in some asymmetric hc mixtures)
      
      ! Self Cross found
      Tcr2 = crossings(1)%x
      Pcr2 = crossings(1)%y

      ! New Kfactors interpolated for the right cross
      kfcr2 = kfcross( &
         crossings(1)%i, dew_envelope%t, dew_envelope%logk, Tcr2 &
      )

      kscr2 = kfcross( &
         crossings(1)%j, dew_envelope%t, dew_envelope%logk, Tcr2 &
      )
      print *, "SELF CROSS", Tcr2, Pcr2

   end if

   if (P > Pmax) then  ! now run from Low T Bubble point
      print *, "Running LowTBub Line"
      call cpu_time(start_time)
      ichoice = 1
      P = 11.0   ! 11.0
      T = 205.0

      do while (P > 10) ! > 10
         T = T - 5.D0
         P = sum(z*pc*exp(5.373*(1 + omg)*(1 - tc/T)))
      end do

      KFACT = pc*exp(5.373*(1 + omg)*(1 - tc/T))/P

      call envelope2(&
         ichoice, n, z, T, P, KFACT, n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri, low_t_envelope &
      )
      call cpu_time(end_time)
      print *, "Done in: ", end_time-start_time
      call WriteEnvel(n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)
      call low_t_envelope%write("env23out/envelout2-LTBUB")

      ! ========================================================================
      !  Crossings finding
      !   To find the crossings between envelopes a generic cross finding
      !   subroutine is called. After the... crosses? have been found,
      !   the K-Factors of each envelope (and each cross) are calculated
      !   with an interpolation based on temperature.
      ! ------------------------------------------------------------------------

      ! Find the cross between two lines (in this case, dew envelope and
      ! (starting from) low temperature bubble envelope
      if (Tcr2 < 1d-5) then
         ! If there wasn't a self-cross in the dew line, find crossings between
         ! dew and bub
         call find_cross(dew_envelope%t, low_t_envelope%t, &
                        dew_envelope%p, low_t_envelope%p, crossings)
      end if

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

   else
      ! In the case the low pressure dew line went near negative pressure
      ! run from High P L-L saturation (incipient phase rich in last comp.)
      print *, "Running HighPLL Line"
      call cpu_time(start_time)
      ichoice = 3
      P = Pmax
      T = 710.0
      iy = 1   ! Liquid root
      ix = 1   ! Liquid root
      y = 0.d0
      y(n) = 1.d0
      difw = -1.d0

      do while (difW < 0.d0)
         T = T - 10.D0
         call TERMO(n, ix, 1, T, P, z, Vx, PHILOGx, DLPHIPx, DLPHITx, FUGNx) ! for fluid
         call TERMO(n, iy, 1, T, P, y, Vy, PHILOGy, DLPHIPy, DLPHITy, FUGNy) ! for pure Water/Asphaltene
         difW = log(z(n)) + PHILOGx(n) - log(y(n)) - PHILOGy(n)
      end do

      KFACT = 1.d-3
      KFACT(n) = 1/z(n)
      call envelope2(ichoice, n, z, T, P, KFACT, n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri, high_p_envelope)
      call cpu_time(end_time)
      print *, "Done in: ", end_time - start_time
      call high_p_envelope%write("env23out/envelout2-HPLL")
      call WriteEnvel(n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)

      if (P < 15.0) then ! 5/6/22
         ichoice = 1   ! now run from Low T Bubble point, after isolated LL saturation curve
         P = 11.0
         T = 205.0

         do while (P > 10) ! > 10, temporary was 7 for some reason
            T = T - 5.D0
            P = sum(z*pc*exp(5.373*(1 + omg)*(1 - TC/T)))
         end do

         KFACT = PC*exp(5.373*(1 + omg)*(1 - TC/T))/P
         call envelope2(ichoice, n, z, T, P, KFACT, n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri, low_t_envelope)
         call low_t_envelope%write("env23out/envelout2-LTBUB")
         
         call WriteEnvel(n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)

         ! Find cross between high P LL and bubble
         call find_cross(high_p_envelope%t, low_t_envelope%t, &
                         high_p_envelope%p, low_t_envelope%p, crossings)
         
         if (size(crossings) < 1) then
            call find_cross(high_p_envelope%t, dew_envelope%t, &
                            high_p_envelope%p, dew_envelope%p, crossings)
         if (size(crossings) > 1) then
            Tcr1 = crossings(1)%x
            Pcr1 = crossings(1)%y
            icross = crossings(1)%i
            jcross = crossings(1)%j

            ! New Kfactors interpolated for the Left cross
            kfcr1 = kfcross(jcross, dew_envelope%t, dew_envelope%logk, Tcr1)
            kscr1 = kfcross(icross, high_p_envelope%t, high_p_envelope%logk, Tcr1)
         endif
         else
            Tcr1 = crossings(1)%x
            Pcr1 = crossings(1)%y
            icross = crossings(1)%i
            jcross = crossings(1)%j

            ! New Kfactors interpolated for the Left cross
            kfcr1 = kfcross(jcross, low_t_envelope%t, low_t_envelope%logk, Tcr1)
            kscr1 = kfcross(icross, high_p_envelope%t, high_p_envelope%logk, Tcr1)
         end if

         if (Tcr2 < 1) then
            ! Find cross between dew and bubble
            call find_cross(dew_envelope%t, low_t_envelope%t, &
                           dew_envelope%p, low_t_envelope%p, crossings)
            if (size(crossings) > 1) then
               Tcr2 = crossings(1)%x
               Pcr2 = crossings(1)%y
               icross = crossings(1)%i
               jcross = crossings(1)%j

               ! New Kfactors interpolated for the right cross
               kfcr2 = kfcross(jcross, low_t_envelope%t, low_t_envelope%logk, Tcr2)
               kscr2 = kfcross(icross, dew_envelope%t, dew_envelope%logk, Tcr2)
            end if
         end if
      
      else
         low_t_envelope = dew_envelope
         low_t_envelope%logk = low_t_envelope%logk(size(low_t_envelope%t):1:-1, :)
         low_t_envelope%t = low_t_envelope%t(size(low_t_envelope%t):1:-1)
         low_t_envelope%p = low_t_envelope%p(size(low_t_envelope%t):1:-1)
      end if
   end if

   print *, Tcr1, Pcr1, "cross1"
   print *, Tcr2, Pcr2, "cross2"

   open(42, file="./env23out/DSPs")
   if (allocated(crossings)) then
      do i = 1, size(crossings)
         write(42, *) crossings(i)%x, crossings(i)%y
      end do
   end if
   close(42)

   print *, "============="
   print *, Tcr1, Pcr1
   print *, kfcr1(:4)
   print *, kscr1(:4)
   print *, ""
   print *, Tcr2, Pcr2
   print *, kfcr2(:4)
   print *, kscr2(:4)
   print *, "-------------"

   block_cross: block
      use envelopes, only: find_crossings
      call find_crossings(&
         dew_envelope, low_t_envelope, high_p_envelope, &
         Tcr1, Pcr1, Tcr2, Pcr2, &
         kfcr1, kscr1, kfcr2, kscr2)
         print *, "============="
         print *, Tcr1, Pcr1
         print *, kfcr1(:4)
         print *, kscr1(:4)
         print *, ""
         print *, Tcr2, Pcr2
         print *, kfcr2(:4)
         print *, kscr2(:4)
         print *, "-------------"
   end block block_cross

   if (three_phase /= "yes") then
      call exit(0)
   end if

   if (abs(Tcr1) < 1d-5) then ! no crossings
      print *, "no cross found"

      ! ========================================================================
      !  Inner envelopes calculation
      ! ------------------------------------------------------------------------
      print *, "Running 3phase with incipient vapor"

      ! First inner curve: 3-phase bubble curve (incipient: V)
      ichoice = 1

      ! Start from the low T bubble envelope first point
      if (allocated(low_t_envelope%t)) then
         T = low_t_envelope%t(1)
         P = low_t_envelope%p(1)
      else
         T = 150.d0
         P = 10.d0
      end if

      beta = z(n) ! asphaltenes (or last, separating compound) molar fraction WHY
      KFACT = exp(low_t_envelope%logk(1, :))
      
      ! Second liquid fugacity
      y = 0
      y(n) = 1
      call TERMO(n, 1, 1, T, P, y, Vy, PHILOGy, DLPHIPy, DLPHITy, FUGNy)

      ! K = bubble_fug/pure_asph_fug, for two liquid phases 
      call TERMO(n, 1, 1, T, P, z, Vx, PHILOGx, DLPHIPx, DLPHITx, FUGNx)
      KFsep = exp(PHILOGx - philogy)

      call envelope3(ichoice, n, z, T, P, beta, KFACT, KFsep, n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri, triphasic)
      call WriteEnvel(n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)
      call triphasic%write("env23out/envelout3-NOCROSS1")

      ! ------------------------------------------------------------------------
      ! second inner curve: Lower AOP curve (incipient: Asph. or 2nd L)
      ichoice = 2 
      T = low_t_envelope%t(1)
      P = low_t_envelope%p(1)/2
      
      ! ========================================================================
      !  Find starting point for second inner curve
      ! ------------------------------------------------------------------------
      FIRST = .true.
      spec = 'TP'

      call flash(spec, FIRST, nmodel, n, z, tc, pc, omg, ac, b, rk_or_m, delta1, &
                 Kij_or_K0, Tstar, Lij, t, p, v, xx, w, rho_x, rho_y, beta, iter)
      call TERMO(n, 1, 1, T, P, xx, Vx, PHILOGx, DLPHIPy, DLPHITy, FUGNy)

      y = 0 
      y(n) = 1

      call TERMO(n, 1, 1, T, P, y, Vy, PHILOGy, DLPHIPy, DLPHITy, FUGNy)
      dif = PHILOGy(n) - (PHILOGx(n) + log(xx(n)))
      Pold = P

      if (dif < 0) then
         P = 0.9*P
      else
         P = 1.1*P
      end if

      do while (abs(dif) > 0.1 .and. P > 0.9)
         call flash(spec, FIRST, nmodel, n, z, tc, pc, omg, ac, b, rk_or_m, delta1, &
                    Kij_or_K0, Tstar, Lij, t, p, v, xx, w, rho_x, rho_y, beta, iter)
         call TERMO(n, 1, 1, T, P, xx, Vx, PHILOGx, DLPHIPy, DLPHITy, FUGNy)
         call TERMO(n, 1, 1, T, P, y, Vy, PHILOGy, DLPHIPy, DLPHITy, FUGNy)
         dold = dif

         dif = PHILOGy(n) - (PHILOGx(n) + log(xx(n)))

         aux = P

         extrapolation = (dold*(P - Pold)/(dif - dold))

         P = max(P/10, Pold - extrapolation)
         Pold = aux
      end do

      if (abs(dif) > 0.1) then
         FIRST = .true.
         Told = T
         T = T + 10.0
      end if

      do while (abs(dif) > 0.1)
         call flash(spec, FIRST, &
            nmodel, n, z, tc, pc, omg, ac, b, rk_or_m, delta1, & ! spec
            Kij_or_K0, Tstar, Lij, & ! spec
            t, p, v, xx, w, rho_x, rho_y, beta, iter & ! out
         )
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
      ! ========================================================================

      call envelope3(&
         ichoice, n, z, T, P, beta, KFACT, KFsep, &
         n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri, triphasic &
      )
      call WriteEnvel(n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)
      call triphasic%write("env23out/envelout3-NOCROSS2")
   else  ! at least one crossing
      if (ichoice /= 3) ichoice = 1
      print *, "at least one cross", Tcr1, Pcr1
      T = Tcr1
      P = Pcr1
      beta = 0.0d0
      KFACT = exp(KFcr1(:n))
      KFsep = exp(Kscr1(:n))
      call envelope3(&
         ichoice, n, z, T, P, beta, KFACT, KFsep, &
         n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri, triphasic &
      )
      call WriteEnvel(n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)
      call triphasic%write("env23out/envelout3-CROSS11")

      if (ichoice /= 3) ichoice = 2
      T = Tcr1
      P = Pcr1
      beta = 0.0d0
      KFACT = exp(Kscr1(:n)) ! now y (incipient phase in envelope3) will be the second liquid
      KFsep = exp(KFcr1(:n)) ! w will be vapor, with phase fraction beta

      call envelope3(&
         ichoice, n, z, T, P, beta, KFACT, KFsep, &
         n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri, triphasic &
      )
      call WriteEnvel(n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)
      call triphasic%write("env23out/envelout3-CROSS12")
   end if

   ! Maybe this should be inside the else above?
   ! Not necesarily, there can be a single dsp. Still the whole logic
   ! should be revised
   if (abs(Tcr2) > 1d-5) then
      print *, "second cross", Tcr2, Pcr2

      ! ========================================================================
      ! Setup envelope
      ! ------------------------------------------------------------------------
      ! Check if the DSP is after the two-phase bubble line critical point
      ! if it is, initialize as a incipient liquid line
      if (size(dew_envelope%critical_points) > 0) then
         if (Tcr2 > dew_envelope%critical_points(1)%t) then
            print *, "CASE: DSP AFTER CP"
            ichoice = 3
         end if
      else
         print *, "CASE: DSP BEFORE CP"
         ichoice = 1
      end if
         
      T = Tcr2
      P = Pcr2
      beta = 0.0d0
      KFACT = exp(KFcr2(:n))
      KFsep = exp(Kscr2(:n))
      call envelope3(&
         ichoice, n, z, T, P, beta, KFACT, KFsep, &
         n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri, triphasic &
      )
      call WriteEnvel(n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)
      call triphasic%write("env23out/envelout3-CROSS21")
      ! ========================================================================

      ! ========================================================================
      ! Setup envelope
      ! ------------------------------------------------------------------------
      ichoice = 2
      T = Tcr2
      P = Pcr2
      beta = 0.0d0
      KFACT = exp(Kscr2(:n)) ! now y (incipient phase in envelope3) will be the second liquid
      KFsep = exp(KFcr2(:n)) ! w will be vapor, with phase fraction beta
      call envelope3(&
         ichoice, n, z, T, P, beta, KFACT, KFsep, &
         n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri, triphasic &
      )
      call WriteEnvel(n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)
      call triphasic%write("env23out/envelout3-CROSS22")
      ! ========================================================================
   end if
end subroutine readcase

subroutine WriteEnvel(n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)
   use constants
   use file_operations, only: out_i, outfile
   use envelopes, only: max_points

   ! T, P and Density of the calculated envelope
   real(pr), dimension(max_points) :: Tv
   real(pr), dimension(max_points) :: Pv
   real(pr), dimension(max_points) :: Dv

   ! number of valid elements in To, Po and Do arrays
   integer :: n_points

   ! positions of the last saturation points before each critical point
   integer, dimension(4) :: icri

   ! T, P and Density of critical points
   real(pr), dimension(4) :: Tcri
   real(pr), dimension(4) :: Pcri
   real(pr), dimension(4) :: Dcri

   character(len=*), parameter :: fmt = "(3E14.4)"

   ! number of valid elements in icri, Tcri, Pcri and Dcri arrays
   integer :: ncri

   character(len=200) :: filename
   integer :: file_unit

   out_i = out_i + 1
   filename = "envelout"
   filename = outfile(filename, out_i)

   open(newunit=file_unit, file=filename)

   write (file_unit, *) 'T(K)        P(bar)        D(mol/L)'
   do i = 1, n_points
      write (file_unit, *) Tv(i), Pv(i), Dv(i)
   end do
!
!1  format(F12.4, 2E14.4, x, I4)
!   write (file_unit, *)
!   write (file_unit, *) ' Number of critical points found: ', ncri
!   write (file_unit, *) '   T(K)        P(bar)        D(mol/L)'
!   do i = 1, ncri
!      write (file_unit, 1) Tcri(i), Pcri(i), Dcri(i), icri(i)
!   end do
   close (file_unit)
end subroutine WriteEnvel

subroutine envelope2(ichoice, n, z, T, P, KFACT, & ! This will probably always exist
                     n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri, & ! This shouldnt be here in the future
                     this_envelope) ! This output should encapsulate everything
   use constants
   use dtypes, only: envelope, find_cross, point, critical_point
   use system, only: nc, nmodel => thermo_model, &
                     tc, pc, dceos => dc, omg => w, &
                     ac, b, delta1 => del1, rk_or_m => k, &
                     kij_or_k0 => kij, ntdep => tdep, ncomb => mixing_rule, bij,&
                     kinf, tstar, lij
   use linalg, only: solve_system
   use envelopes, only: F2, X2, fix_delx, update_specification, max_points, &
                        env_number
   implicit none

   ! number of compounds in the system and starting point type
   integer, intent(in) :: n, ichoice

   ! estimated T and P for first point (then used for every point)
   real(pr) :: T, P

   ! Maximun pressure
   real(pr) :: maxP

   ! estimated K factors for first point (then used for every point)
   real(pr), intent(in out) :: KFACT(n)

   ! composition of the system
   real(pr), intent(in) :: z(n)

   ! T, P and Density of the calculated envelope
   real(pr), intent(out) :: Tv(max_points)
   real(pr), intent(out) :: Pv(max_points)
   real(pr), intent(out) :: Dv(max_points)

   ! number of valid elements in Tv, Pv and Dv arrays
   integer, intent(out) :: n_points

   ! positions of the last saturation points before each critical point
   integer, dimension(4), intent(out) :: icri

   ! T, P and Density of critical points
   real(pr), dimension(4), intent(out) :: Tcri(4), Pcri(4), Dcri(4)

   ! number of valid elements in icri, Tcri, Pcri and Dcri arrays
   integer, intent(out) :: ncri

   ! Intermediate variables during calculation process
   real(pr), dimension(n) :: y
   integer, dimension(n + 2) :: ipiv
   real(pr), dimension(n + 2) :: X, Xold, Xold2, delX, bd, F, dFdS, dXdS
   real(pr), dimension(n + 2, n + 2) :: JAC, AJ
   real(pr) :: Vy, Vx
   logical :: run, passingcri, minT, minmaxT

   character(len=:), allocatable :: incipient_phase

   type(envelope), intent(out) :: this_envelope
   real(pr) :: tmp_logk(max_points, n)
   real(pr) :: tmp_logphi(max_points, n)
   
   ! Extrapolation of variables to detect critical points
   real(pr) :: extra_slope(n + 2)
   real(pr) :: lnK_extrapolated(n)
   real(pr) :: delta_t

   integer :: i
   integer :: iy, ix ! Vapor or liquid selectors

   ! Specification value, delta and index
   real(pr) :: S, delS
   integer :: ns

   real(pr) :: Told2, Told
   real(pr) :: frac

   ! Netwon method
   integer :: iter ! Iteration
   integer :: max_iter

   ! Critical Points
   type(critical_point), allocatable :: critical_points(:)
   integer :: black_i ! Number of steps while trying to escape the CP
   real(pr) :: stepx

   integer :: funit_it
   integer :: funit_env
   character(len=20) :: fname_it
   character(len=20) :: fname_env

   ! =============================================================================
   !  OUTPUT file
   ! -----------------------------------------------------------------------------
   env_number = env_number + 1
   write(fname_env, *) env_number
   print *, fname_env
   fname_env = "ENV2_OUT" // "_" // trim(adjustl(fname_env))
   print *, fname_env
   open(newunit=funit_env, file=fname_env)
   ! =============================================================================

   ! Initialize with zero Tv and Pv
   allocate(this_envelope%vars(max_points-50, n+2))
   Tv = 0
   Pv = 0

   minT = .false.
   minmaxT = .false.
   passingcri = .false.
   Told2 = 0.0
   Told = 10.0
   maxP = 0.d0

   !-----------------------------------------------------------
   ! Continuation method for tracing the envelope starts here
   run = .true.
   i = 0
   ncri = 0
   JAC(n + 1, :) = 0.d0
   ! lda = n + 2
   ! ldb = n + 2
   X(:n) = log(KFACT)
   X(n + 1) = log(T)
   X(n + 2) = log(P)
   iy = 1
   ix = 1

   select case(ichoice)
   case (1)
      incipient_phase = "vapor"
      iy = -1
   case (2)
      incipient_phase = "liquid"
      ix = -1
   case (3)
      incipient_phase = "2ndliquid"
   end select
   write(funit_env, *) incipient_phase

   if (ichoice <= 2) then  
      ! low T bub (1) or dew (2)
      ! x will be vapor phase during the first part, 
      ! and liquid after a critical point is crossed
      if (ichoice == 1) iy = -1
      if (ichoice == 2) ix = -1  
      ns = n + 1
      S = log(T)
      delS = 0.005

      ! Wilson estimate for vapor (or liquid) composition
      y = KFACT*z 
   else    
      ! (ichoice==3) high P L-L sat
      ! PmaxDewC = maxval(PdewC(1:ilastDewC))
      ns = n + 2
      S = log(P)
      delS = -0.005
      y = 0.d0
      y(n) = 1.d0
   end if

   Xold = 0.d0
   dFdS = 0.d0
   dFdS(n + 2) = -1.d0

   write(fname_it, *) ichoice
   fname_it = trim(adjustl("X_it_2ph" // "_" //adjustl(trim(fname_it))))

   open(newunit=funit_it, file=fname_it)
   do while (run)
      i = i + 1
      if (i > max_points - 50) then
         exit
      end if
      ! Newton starts here
      delX = 1.0
      iter = 0
      max_iter = 500

      do while (maxval(abs(delX)) > 1.d-9 .and. iter <= max_iter)
         ! Solve point with full Newton method
         call F2(incipient_phase, z, y, X, S, ns, F, JAC)

         iter = iter + 1

         ! TODO: For some reason passing JAC and -F breaks the system
         !       update: it's because dgesv redefines A and B
         bd = -F
         AJ = JAC
         delX = solve_system(AJ, bd)
         call fix_delX(i, iter, 3, 10.0_pr, 0.08_pr, delX)

         X = X + delX

         if (.not. passingcri .and. i /= 1 &
             .and. iter > 10 &
             .and. maxval(abs(delX)) > 0.001) then 
            ! Too many iterations-->Reduce step to new point

            delS = delS*2.0/4.0
            S = S - delS
            X = Xold + dXdS*delS
         end if

         KFACT = exp(X(:n))
         y = z*KFACT
         T = exp(X(n + 1))
         P = exp(X(n + 2))

         write(funit_it, *) i, iter, T, P, KFACT
      end do
      write(funit_it, *) " "
      write(funit_it, *) " "

      ! Point converged (unless it jumped out because of high number of iterations)
      write(funit_env, *) T, P, exp(X(:n))
      if (iter > max_iter) run = .false.
      if (P > maxP) maxP = P

      if (incipient_phase == "liquid" .and. i > 1) then
         ! TODO: If this is the way the low p dew line finishes, 
         ! I think this could be better, like using dPdT
         if (P < Pv(i - 1) .and. P < maxP/5 .and. T > 300) then
            run = .true.  ! to finish envelope going to low T bubble
         end if
      end if

      Tv(i) = T
      Pv(i) = P
      Dv(i) = 1/Vx    ! saturated phase density

      tmp_logk(i, :n) = X(:n)
      ! tmp_logphi(i, :n) = philogx(:n)

      ! rho_y = 1/Vy     incipient phase density

      if (incipient_phase == "2ndliquid" .and. P < 1.0) then
         ! isolated LL line detected. 
         ! Stop and start a new one from low T false bubble point
         run = .false.
      end if
      
      print *, incipient_phase, i, T, P, ns, iter
      if (i > max_points - 50) exit

      if (sum(X(:n) * Xold(:n)) < 0) then  ! critical point detected
         print *, "Found critical!"
         ncri = ncri + 1
         icri(ncri) = i - 1
         frac = -Xold(ns)/(X(ns) - Xold(ns))
         Tcri(ncri) = Tv(i - 1) + frac*(T - Tv(i - 1))
         Pcri(ncri) = Pv(i - 1) + frac*(P - Pv(i - 1))
         Dcri(ncri) = Dv(i - 1) + frac*(Dv(i) - Dv(i - 1))


         select case (incipient_phase)
         case("liquid")
            incipient_phase = "vapor"
         case("vapor")
            incipient_phase = "liquid"
         end select

         write(funit_env, *) " "
         write(funit_env, *) " "
         write(funit_env, *) incipient_phase
      end if

      if (run) then
         ! Calculation of sensitivities (dXdS)
         ! dgesv( n, nrhs, a, lda, ipiv, b, ldb, info )

         call update_specification(iter, passingcri, X, JAC, ns, S, delS, dXdS)

         ! Generation of estimates for the next point
         Told2 = Told
         Told = T
         Xold2 = Xold
         Xold = X
         X = Xold + dXdS*delS

         black_i = 0
         step_fixer = 1.d0
         if (passingcri) passingcri = .false.

         do while (maxval(abs(X(:n))) < 0.03)
            print *, "Jumping critical"
            ! approaching the black hole... get out of there! (0.03)
            black_i = black_i + 1
            if (black_i > 50) then
                print *, "Stuck on the black hole"
                if (black_i > 100) stop
            end if

            stepX = maxval(abs(X(:n) - Xold(:n))) 

            ! the step given by the most changing logK to fall into the black hole
            passingcri = .true.
            if (stepX > 0.07) then
               S = S - delS/2
               X = X - dXdS*delS/2   !  half step back
            else
               S = S + delS
               X = X + dXdS*delS   ! one more step to jump over the critical point
            end if
         end do

         ! Extrapolate lnK ten degrees to detect CP
         if (i > 10) then
            extra_slope = (tmp_logk(i, :n) - tmp_logk(i-1, :n))/(Tv(i) - Tv(i-1))


         T = exp(X(n + 1))

         if (.not. passingcri .and. abs(T - Told) > 7) then 
            ! Delta T estimations > 7K are not allowed
            delS = delS/2
            S = S - delS
            X = Xold + dXdS*delS
            T = exp(X(n + 1))
         end if

         P = exp(X(n + 2))
         KFACT = exp(X(:n))
         y = z*KFACT

         ! Finish conditions
         if ((dXdS(n + 1)*delS < 0 .and. P < 0.1 .or. T < 120.0) &  ! dew line stops when P<0.1 bar or T<150K
             .or. (P > 1.0 .and. T < 150.0) &   ! bubble line stops when T<150K
             .or. (P > 1500) &
             .or. (abs(dels) < 1.d-6)) then
            run = .false.
         end if
      end if
   end do
   !-----------------------------------------------------------

   n_points = i

   write(funit_env, *) " "
   write(funit_env, *) " "
   write(funit_env, *) "critical"
   do i=1, ncri
      write(funit_env, *) Tcri(i), Pcri(i)
   end do

   ! Define envelope values, omit the last point to avoid not really
   ! converged cases
   close(funit_it)
   close(funit_env)
   this_envelope%logk = tmp_logk(:n_points - 1, :)
   this_envelope%logphi = tmp_logphi(:n_points - 1, :)
   this_envelope%t = Tv(:n_points - 1)
   this_envelope%p = Pv(:n_points - 1)
   this_envelope%z = z

   allocate(critical_points(ncri))
   critical_points%t = tcri(:ncri)
   critical_points%p = pcri(:ncri)
   this_envelope%critical_points = critical_points
end subroutine envelope2

subroutine envelope3(&
         ichoice, n, z, T, P, beta, KFACT, KFsep, & ! Inputs
         n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri, this_envelope & !Outputs
   )
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
   use nfile, only: xit_file, xit_filename => filename
   use constants
   use dtypes, only: env3, critical_point
   use system, only: & 
                     nmodel => thermo_model, ncomb => mixing_rule, ntdep => tdep, &
                     ac, b, delta1 => del1, rk_or_m => k, &
                     tc, pc, dceos => dc, omg => w, &
                     kij_or_k0 => kij, lij, bij, moles
   use linalg, only: solve_system
   use io, only: str
   use thermo, only: pressure
   use envelopes, only: max_points, env_number

   implicit none

   integer, parameter :: max_iter=500   !! Max number of iterations per Newton step
   integer, parameter :: isot_n=5000    !! Number of points at the final isotherm

   ! number of compounds in the system and starting point type
   integer, intent(in) :: n, ichoice

   ! estimated T, P and "w" phase fraction for first point (then used for every point)
   real(pr) :: T, P, beta

   ! estimated K factors for first point (then used for every point)
   real(pr) :: KFACT(n), KFsep(n)

   ! composition of the system
   real(pr), intent(in) :: z(n)

   ! T, P and Density of the calculated envelope
   real(pr), intent(out) :: Tv(max_points)
   real(pr), intent(out) :: Pv(max_points)
   real(pr), intent(out) :: Dv(max_points)

   ! number of valid elements in Tv, Pv and Dv arrays
   integer, intent(out) :: n_points

   ! positions of the last saturation points before each critical point
   integer, intent(out) :: icri(4)

   ! T, P and Density of critical points
   real(pr), intent(out) :: Tcri(4)
   real(pr), intent(out) :: Pcri(4)
   real(pr), intent(out) :: Dcri(4)

   ! number of valid elements in icri, Tcri, Pcri and Dcri arrays
   integer, intent(out) :: ncri

   ! Intermediate variables during calculation process
   real(pr), dimension(n) :: y, xx, w, PHILOGy, PHILOGx, PHILOGw
   real(pr), dimension(n) :: dxdB, dydB, dwdB, dxdKs, dydKs, dwdKs, aux
   real(pr), dimension(n) :: DLPHITx, DLPHIPx, DLPHITy, DLPHIPy, DLPHITw, DLPHIPw
   real(pr), dimension(n, n) :: FUGNx, FUGNy, FUGNw
   integer, dimension(2*n + 3) :: ipiv
   real(pr), dimension(2*n + 3) :: X, Xold, delX, F, dFdS, dXdS
   real(pr) :: bd(2*n + 3), AJ(2*n + 3, 2*n + 3)
   real(pr), dimension(2*n + 3, 2*n + 3) :: JAC
   real(pr) :: Vy, Vx, Vw
   logical :: run, passingcri, Comp3ph  !, Cross                 ! crossing var (Cross)

   integer  :: i, j, l      ! Iteration variables
   integer  :: i1, i2       ! components to save out info about during calculations
   integer  :: ix, iy, iw   ! Incipient phase?
   integer  :: lda, ldb     ! Dimension of the system of equations
   integer  :: iter         ! Number of iterations during full newton
   real(pr) :: rho_x, rho_y ! Phases densities
   real(pr) :: frac

   ! Newton method variables
   real(pr) :: Told ! Previous iteration temperature
   integer :: info  ! LAPACK's solve status
  
   ! Specification related variables 
   integer  :: ns, nsold ! Which specification is used
   real(pr) :: s, dels   ! Specification value and deltaS
   real(pr) :: delmax    ! Maximum delta S based on deltaX
   real(pr) :: updel     ! Maximum delta S based on iterations and previus delS

   ! Critical points related
   integer :: black_i ! Number of iterations trying to escape the black hole
   real(pr) :: stepx  ! Maximum changing logK
   real(pr) :: stepx_k, stepx_ks ! Maximun changing logK for each set
   logical :: K_critical, KS_critical

   ! Phase envelope points
   type(env3), intent(in out) :: this_envelope
   real(pr) :: incipient(n), saturated(n), minoritary(n)
   real(pr) :: tmp_logk(max_points, n)
   real(pr) :: tmp_logks(max_points, n)
   real(pr) :: tmp_x(max_points, n)
   real(pr) :: tmp_y(max_points, n)
   real(pr) :: tmp_w(max_points, n)
   real(pr) :: tmp_beta(max_points)
   type(critical_point), allocatable :: ll_critical_points(:), critical_points(:)

   integer :: unit_iterations, tmp_i

   character(len=:), allocatable :: incipient_phase
   

   ! Variables to draw an isotherm at the end of the algorithm
   real(pr) :: isot_v(isot_n), &
               isot_pz(isot_n), &
               isot_pxx(isot_n), &
               isot_py(isot_n), &
               isot_pw(isot_n)

   integer :: funit_it
   integer :: funit_env
   character(len=20) :: fname_it
   character(len=20) :: fname_env

   ! =============================================================================
   !  OUTPUT file
   ! -----------------------------------------------------------------------------
   env_number = env_number + 1
   write(fname_env, *) env_number
   print *, fname_env
   fname_env = "ENV3_OUT" // "_" // trim(adjustl(fname_env))
   print *, fname_env
   open(newunit=funit_env, file=fname_env)
   ! =============================================================================
   
   ix = 1
   iy = 1
   iw = 1
   allocate(this_envelope%vars(max_points - 50, 2*n+3), stat=i)

   select case(ichoice)
      case (1)
         incipient_phase = "Vapor"
         iy = -1
      case (2)
         incipient_phase = "2nd Liquid"
         iw = -1
      case (3)
         incipient_phase = "Saturated vapor"
         ix = -1
   end select
   print *, "Running: ", incipient_phase
   print *, "init: ", T, P
   !-----------------------------------------------------------

   ! Continuation method for tracing the envelope starts here
   passingcri = .false.
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

   if (beta .le. 1.0e-12) then
      ns = 2*n + 3
      S = 0.00   ! for beta
      delS = 0.001
   else 
      ! start from low T "bubble" point or low T "Lower AOP"
      ns = 2*n + 2 ! ln(P) specification
      S = X(ns)
      delS = 0.001
   end if

   xx = z/(1 - beta + beta*KFsep)
   y = KFACT*xx
   w = KFsep*xx

   Xold = 0.d0
   dFdS = 0.d0
   dFdS(2*n + 3) = -1.d0

   i = 0

   xit_file = xit_file + 1
   write(xit_filename, *) xit_file
   xit_filename = trim(adjustl("Xit" // adjustl(trim(xit_filename))))
   open(newunit=unit_iterations, file=xit_filename)
   ! Save all the newton iterations into a file
   write(unit_iterations, *) &
   "point ", ("K" // str(tmp_i) // " ", tmp_i=1,n), &
   ("KS" // str(tmp_i) // " ", tmp_i=1,n), &
   "T ", "P ", "beta ",  &
   ("del_K" // str(tmp_i) // " ", tmp_i=1,n), &
   ("del_KS" // str(tmp_i) // " ", tmp_i=1,n), &
   "del_T ", "del_P ", "del_beta ", "vx ", "vy ", "vw "

   write(funit_env, *) incipient_phase
   do while (run .and. i < max_points)
      i = i + 1  ! number of point to be calculated along the line

      ! ========================================================================
      ! Newton for i point starts here
      ! ------------------------------------------------------------------------
      delX = 1.0
      iter = 0
      do while (maxval(abs(delX)) > 1.d-6 .and. iter <= max_iter)
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
            JAC(n + 1:2*n, j) =   KFsep(j - n)*(FUGNw(:, j - n)*dwdKs(j - n) &
                                - FUGNx(:, j - n)*dxdKs(j - n))
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

         delX = bd

         if (i == 1) then
            do while (maxval(abs(delX)) > 0.1)   
               ! Too large Newton step --> Reduce it
               delX = delX/2
            end do
         else
            do while (maxval(abs(delX)) > 0.08)   
               ! Too large Newton step --> Reduce it
               delX = delX/2
            end do
            if (iter > 10) delX = delX/2  ! too many iterations (sometimes due to oscillatory behavior near crit) --> Reduce it
         end if

         write(unit_iterations, *) i, X, delX, vx, vy, vw

         X = X + delX

         if (.not. passingcri .and. i/= 1 .and. iter > max_iter/10 .and. maxval(abs(delX)) > 0.001) then ! Too many iterations--> Reduce step to new point
            delS = delS * 3.d0/4.d0
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
      !if (ichoice == 1) call exit()
      ! ========================================================================

      ! Point converged 
      ! unless it jumped out because of high number of iterations
      ! or the phase fraction went out of 0 < beta < 1
      if (iter > max_iter) run = .false.
      if (beta < 0 .or. beta > 1) run = .false.

      Tv(i) = T
      Pv(i) = P
      Dv(i) = 1/Vx    ! saturated phase density
      rho_x = 1/Vx
      rho_y = 1/Vy    ! incipient phase density
      this_envelope%vars(i, :) = X

      ! ========================================================================
      ! Reasign phases after a critical point
      if (passingcri) then
         ncri = ncri + 1
         icri(ncri) = i - 1
         frac = -Xold(ns)/(X(ns) - Xold(ns))
         Tcri(ncri) = Tv(i - 1) + frac*(T - Tv(i - 1))
         Pcri(ncri) = Pv(i - 1) + frac*(P - Pv(i - 1))
         Dcri(ncri) = Dv(i - 1) + frac*(Dv(i) - Dv(i - 1))
         iy = -iy
         !ix = -ix
      end if

      if (sum(X(n + 1:2*n)*Xold(n + 1:2*n)) < 0) then  
         ! critical point detected between x and w phases
         ncri = ncri + 1
         icri(ncri) = i - 1
         frac = -Xold(ns)/(X(ns) - Xold(ns))
         Tcri(ncri) = Tv(i - 1) + frac*(T - Tv(i - 1))
         Pcri(ncri) = Pv(i - 1) + frac*(P - Pv(i - 1))
         Dcri(ncri) = Dv(i - 1) + frac*(Dv(i) - Dv(i - 1))
         iw = -iw
         ix = -ix
      end if
      ! ========================================================================

      if (run) then
         ! Calculation of sensitivities (dXdS)
         ! dgesv( n, nrhs, a, lda, ipiv, b, ldb, info )
         bd = -dFdS
         AJ = JAC

         call dgesv(2*n + 3, 1, AJ, lda, ipiv, bd, ldb, info)
         dXdS = bd

         ! Selection of (the most changing) variable to be specified for the next point
         nsold = ns
         ns = maxloc(abs(dXdS), DIM=1)

         if (maxval(abs(X(:n))) < 0.2) then
            ! vars other than logK not allowed to be chosen close to a y-x critical point
            ns = maxloc(abs(dXdS(1:n)), DIM=1)  
         end if

         if (maxval(abs(X(n + 1:2*n))) < 0.2) then
            ! vars other than logKs not allowed to be chosen close to a w-x critical point
            ns = maxloc(abs(dXdS(n + 1:2*n)), DIM=1) 
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

         ! =====================================================================
         ! Critical point detection
         ! ---------------------------------------------------------------------

         if ( maxval(abs(X(:n))) < 0.03) then
            ! approaching the black hole... get out of there!
            if (delS > 0) then 
               delS = 0.04 - S
            else
               delS = -0.04 - S
            end if

            S = S + delS
            X = X + dXdS*delS ! one more step to jump over the critical point
            passingcri = .true.
            print *, "crit K: ", T, P
         end if

         black_i = 0

         do while ( &
               maxval(abs(KFact/KFsep)) < 1.01 & 
               ! .and. &
               ! maxval(abs(X(n+1:2*n)/X(:n))) < 1.05 &
            )
            black_i = black_i + 1
            if (black_i > 10) then
               print *, "Stuck on K/Ks hole"
            else if (black_i > 25) then 
               return
            end if
            ! approaching the black hole... get out of there!
            print *, "crit K: ", T, P, S, & 
               maxval(abs(X(:n)/X(n+1:2*n))), maxval(abs(X(n+1:2*n)/X(:n)))
            
            stepX = maxval(abs(X(n + 1:2*n) - Xold(n + 1:2*n))) 
            ! the step given by the most changing logKs to fall into the black hole

            if (stepX > 0.07) then
               S = S - delS/2
               X = X - dXdS*delS/2 ! half step back
            else
               S = S + 5*delS
               X = X + 5*dXdS*delS ! one more step to jump over the critical point
            end if

            passingcri = .true.
         end do
         

         do while (maxval(abs(X(n + 1:2*n))) < 0.03)
            ! approaching the black hole... get out of there! (0.03)
            black_i = black_i + 1

            if (black_i > 10) then
               print *, "Stuck on the black hole 2"
            end if

            stepX = maxval(abs(X(n + 1:2*n) - Xold(n + 1:2*n))) 
            ! the step given by the most changing logKs to fall into the black hole
            passingcri = .true.

            if (stepX > 0.07) then
               S = S - delS/2
               X = X - dXdS*delS/2 ! half step back
            else
               S = S + delS
               X = X + dXdS*delS ! one more step to jump over the critical point
            end if
            print *, "cirt KS: ", T, P
         end do

         T = exp(X(2*n + 1))

         if (.not. passingcri .and. abs(T - Told) > 7) then 
            ! Delta T estimations > 7K are not allowed
            delS = delS/2
            S = S - delS
            X = Xold + dXdS*delS
            T = exp(X(2*n + 1))
         end if

         ! =====================================================================

         P = exp(X(2*n + 2))
         KFACT = exp(X(:n))
         KFsep = exp(X(n + 1:2*n))
         beta = X(2*n + 3)
         xx = z/(1 - beta + beta*KFsep)
         y = KFACT*xx
         w = KFsep*xx

         if (&! this might need adjustment
             (dXdS(2*n + 1)*delS < 0 .and. P < 0.1 .or. T < 120.0) &  ! Finish on low P and T
             .or. (P > 1.0 .and. T < 150.0) & ! Finish on low T
             .or. (P > 1500) & ! Finish on high P
             .or. (abs(dels) < 1.d-8)) then ! Finish on low delS
            run = .false.
         end if
      end if

      tmp_logk(i, :) = kfact
      tmp_logks(i, :) = KFsep
      tmp_x(i, :) = xx
      tmp_y(i, :) = y
      tmp_w(i, :) = w
      tmp_beta(i) = beta
      write(funit_env, *) T, P, exp(X(:2*n))
      print *, ix, iy, iw, T, P, rho_x, rho_y, delS, ns, iter
   end do
   close(unit=unit_iterations)
   close(unit=funit_env)
   n_points = i - 1
   
   !-----------------------------------------------------------
   ! Save data into the envelope object
   this_envelope%z = z
   this_envelope%p = pv(:n_points)
   this_envelope%t = tv(:n_points)
   this_envelope%logk = tmp_logk(:n_points, :)
   this_envelope%logks = tmp_logks(:n_points, :)
   this_envelope%x = tmp_x(:n_points, :)
   this_envelope%y = tmp_y(:n_points, :)
   this_envelope%w = tmp_w(:n_points, :)
   this_envelope%beta = tmp_beta

   allocate(critical_points(ncri))
   critical_points%t = tcri(:ncri)
   critical_points%p = pcri(:ncri)
   this_envelope%critical_points = critical_points


   ! t = 653
   ! p =  180

   ! j = findloc(&
   !       abs(tv - t) - 1 < 1 .and. abs(pv - p) - 1 < 1, .true., dim=1 &
   !    )

   ! j = n_points

   ! print *, j, tv(j), pv(j)

   ! t = tv(j)

   ! open(newunit=nunit, file="isots")
   ! write(nunit, *) "v pz px py pw"

   ! do i=1,size(isot_v)
   !    isot_v(i) = real(i, pr)/500
   ! end do

   ! moles = z
   ! call pressure(isot_v, t, isot_pz)
   ! moles = this_envelope%x(j, :)
   ! call pressure(isot_v, t, isot_pxx)
   ! moles = this_envelope%y(j, :)
   ! call pressure(isot_v, t, isot_py)
   ! moles = this_envelope%w(j, :)
   ! call pressure(isot_v, t, isot_pw)

   ! do i=1,size(isot_v)
   !    write(nunit, *) isot_v(i), isot_pz(i), isot_pxx(i), isot_py(i), isot_pw(i)
   ! end do

   ! close(nunit)


   ! print *, y
   ! print *, rho_x
   ! print *, rho_y
   ! print *, beta
end subroutine envelope3

