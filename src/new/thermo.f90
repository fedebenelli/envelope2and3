module system
   ! Module for a cubic eos system,
   ! this should be later adapted into a simple oop system where an eos object
   ! stores the relevant parameters
   use constants, only: pr, R
   implicit none

   ! Model settings
   integer :: thermo_model !! Which thermodynamic model
   integer :: tdep !! Temperature dependance of kij
   integer :: mixing_rule !! What mixing rule to use
   integer :: nc !! Number of components

   ! Mole fracions
   real(pr), allocatable :: z(:)

   ! Critical constants
   real(pr), allocatable :: tc(:) !! Critical temperature
   real(pr), allocatable :: pc(:) !! Critical pressure
   real(pr), allocatable :: dc(:) !! Critical density
   real(pr), allocatable :: w(:)  !! Acentric factor

   ! Model parameters
   real(pr), allocatable :: ac(:) !! Critical attractive parameter
   real(pr), allocatable :: b(:)  !! repulsive parameter
   real(pr), allocatable :: del1(:) !! $$\delta_1$$ parameter
   real(pr), allocatable :: k(:) !! Attractive parameter constant
   real(pr), allocatable :: kij(:, :) !! Attractive BIP
   real(pr), allocatable :: lij(:, :) !! Repulsive BIP

   ! T dependant mixing rule parameters
   real(pr), allocatable :: kinf(:, :), tstar(:, :)
   real(pr), allocatable :: bij(:, :)

contains

   subroutine setup(n, nmodel, ntdep, ncomb)
      integer, intent(in) :: n
      integer, intent(in) :: nmodel
      integer, intent(in) :: ntdep
      integer, intent(in) :: ncomb

      thermo_model = nmodel
      tdep = ntdep
      mixing_rule = ncomb
      nc = n

      ! allocate(z(n))
      allocate(tc(n))
      allocate(pc(n))
      allocate(dc(n))
      allocate(w(n))
      allocate(ac(n))
      allocate(b(n))
      allocate(del1(n))
      allocate(k(n))
      allocate(kij(n, n))
      allocate(lij(n, n))
      allocate(kinf(n, n))
      allocate(tstar(n, n))
      ! allocate(aij(n, n))
      ! allocate(daijdt(n, n))
      ! allocate(daijdt2(n, n))
      allocate(bij(n, n))
   end subroutine setup

   subroutine PR_factory(moles_in, ac_in, b_in, tc_in, pc_in, w_in, k_in)
        !! PengRobinson factory
        real(pr), intent(in) :: moles_in(nc)
        real(pr), optional, intent(in) :: ac_in(nc)
        real(pr), optional, intent(in) :: b_in(nc)
        real(pr), optional, intent(in) :: tc_in(nc)
        real(pr), optional, intent(in) :: pc_in(nc)
        real(pr), optional, intent(in) :: w_in(nc)
        real(pr), optional, intent(in) :: k_in(nc)

        integer :: i

        logical :: params_spec, critical_spec
        real(pr) :: zc(nc), oma(nc), omb(nc)
        real(pr) :: vceos(nc), al, be, ga(nc)
        real(pr) :: RTc(nc)

        del1 = 1 + sqrt(2.0_pr)
        z = moles_in

        params_spec = (present(ac_in) .and. present(b_in) .and. present(k_in))
        critical_spec = (present(tc_in) .and. present(pc_in) .and. present(w_in))

        if (params_spec) then
            ac = ac_in
            b = b_in
            k = k_in

            call get_Zc_OMa_OMb(del1, zc, oma, omb)
            Tc = OMb * ac / (OMa * R* b)
            RTc = R * Tc
            Pc = OMb * RTc / b
            Vceos = Zc * RTc / Pc
            al = -0.26992
            be = 1.54226
            ga = 0.37464 - k
            w = 0.5 * (-be + sqrt(be**2 - 4 * al * ga)) / al
        else if (critical_spec) then
            tc = tc_in
            pc = pc_in
            w = w_in
            RTc = R*Tc

            call get_Zc_OMa_OMb(del1, Zc, OMa, OMb)

            ac = OMa * RTc**2 / Pc
            b = OMb * RTc / Pc
            Vceos = Zc * RTc / Pc
            ! k (or m) constant to calculate attractive parameter depending on temperature
            do i=1,nc
               if (w(i) <= 0.491) then
                  ! m from PR
                  k(i) = 0.37464 + 1.54226 * w(i) - 0.26992 * w(i)**2
               else
                  ! PR78
                  k(i) = 0.379642 + 1.48503 * w(i) - 0.164423 * w(i)**2 + 0.016666 * w(i)**3
               end if
            end do
        end if
   end subroutine

   subroutine SRK_factory(moles_in, ac_in, b_in, tc_in, pc_in, w_in, k_in)
        !! PengRobinson factory
        real(pr), intent(in) :: moles_in(nc)
        real(pr), optional, intent(in) :: ac_in(nc)
        real(pr), optional, intent(in) :: b_in(nc)
        real(pr), optional, intent(in) :: tc_in(nc)
        real(pr), optional, intent(in) :: pc_in(nc)
        real(pr), optional, intent(in) :: w_in(nc)
        real(pr), optional, intent(in) :: k_in(nc)

        logical :: params_spec, critical_spec
        real(pr) :: zc(nc), oma(nc), omb(nc)
        real(pr) :: vceos(nc), al, be, ga(nc)
        real(pr) :: RTc(nc)

        del1 = 1
        z = moles_in

        params_spec = (present(ac_in) .and. present(b_in) .and. present(k_in))
        critical_spec = (present(tc_in) .and. present(pc_in) .and. present(w_in))

        if (params_spec) then
            ac = ac_in
            b = b_in
            k = k_in

            call get_Zc_OMa_OMb(del1, zc, oma, omb)
            Tc = OMb * ac / (OMa * R* b)
            RTc = R * Tc
            Pc = OMb * RTc / b
            Vceos = Zc * RTc / Pc
            dc = 1/vceos
            al = -0.26992
            be = 1.54226
            ga = 0.37464 - k
            w = 0.5 * (-be + sqrt(be**2 - 4 * al * ga)) / al
        else if (critical_spec) then
            tc = tc_in
            pc = pc_in
            w = w_in
            RTc = R * Tc

            call get_Zc_OMa_OMb(del1, Zc, OMa, OMb)

            ac = OMa * RTc**2 / Pc
            b = OMb * RTc / Pc
            Vceos = Zc * RTc / Pc

            k = 0.48 + 1.574 * w - 0.175 * w**2
        end if
   end subroutine

   subroutine get_Zc_OMa_OMb(del1, Zc, OMa, OMb)
      !! Calculate Zc, OMa and OMb from the delta_1 parameter.
      real(pr), intent(in)  :: del1(:) !! delta_1 parameter
      real(pr), intent(out) :: Zc(:)   !! Critical compressibility factor
      real(pr), intent(out) :: OMa(:)  !! OMa
      real(pr), intent(out) :: OMb(:)  !! OMb

      real(pr) :: d1(size(del1)), y(size(del1))

      d1 = (1._pr + del1**2._pr)/(1._pr + del1)
      y = 1._pr + (2._pr*(1._pr + del1))**(1.0_pr/3._pr) + (4._pr/(1._pr + del1))**(1.0_pr/3)
      OMa = (3._pr*y*y + 3._pr*y*d1 + d1**2._pr + d1 - 1.0_pr)/(3._pr*y + d1 - 1.0_pr)**2._pr
      OMb = 1._pr/(3._pr*y + d1 - 1.0_pr)
      Zc = y/(3._pr*y + d1 - 1.0_pr)
   end subroutine get_Zc_OMa_OMb
end module system