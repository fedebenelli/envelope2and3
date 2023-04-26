module envelopes
   !! Functions to be used in the different continuation methods to trace
   !! phase envelopes 
   use constants, only: pr
   use linalg, only: solve_system
   use system, only: nc
   use dtypes, only: envelope
   implicit none

   integer, parameter :: max_points = 2000
   integer :: env_number = 0

   interface F
      module procedure :: F2
   end interface

contains

   ! ===========================================================================
   ! General routines
   ! ---------------------------------------------------------------------------
   subroutine update_specification(iter, passingcri, X, dF, ns, S, delS, dXdS)
      integer,  intent(in)     :: iter
      logical,  intent(in)     :: passingcri
      real(pr), intent(in)     :: X(nc + 2)
      real(pr), intent(in)     :: dF(nc + 2, nc + 2)
      integer,  intent(in out) :: ns
      real(pr), intent(in out) :: S
      real(pr), intent(in out) :: delS
      real(pr), intent(in out) :: dXdS(nc + 2)

      real(pr) :: dF_dS(nc + 2)
      real(pr) :: bd(nc + 2)
      real(pr) :: AJ(nc + 2, nc + 2)
      real(pr) :: delmax, updel
      integer  :: nsold

      dF_dS = 0
      call dFdS(dF_dS)

      bd = -dF_dS
      AJ = dF
      dXdS = solve_system(AJ, bd)

      ! Selection of (the most changing) variable to be specified for the next point
      nsold = ns

      ns = maxloc(abs(dXdS), dim=1)

      if (maxval(abs(X(:nc))) < 0.2) then
         ns = maxloc(abs(dXdS(:nc)), dim=1)  ! T and P not allowed to be chosen close to a critical point
      end if

      if (ns /= nsold) then
         delS = dXdS(ns)*delS  ! translation of delS to the  new specification variable
         dXdS = dXdS/dXdS(ns)  ! translation of sensitivities
         S = X(ns)             ! update of S
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
   end subroutine
   ! ===========================================================================

   ! ===========================================================================
   ! Specification function derivatives
   ! ---------------------------------------------------------------------------
   subroutine dFdS(dF_dS)
      use system, only: nc
      real(pr), intent(out) :: dF_dS(nc + 2)

      dF_dS = 0
      dF_dS(nc + 2) = -1
   end subroutine
   ! ---------------------------------------------------------------------------

   ! ===========================================================================
   ! Three Phase envelopes
   ! ---------------------------------------------------------------------------
   function X3(kx, kw, beta, p, t) result(X)
      ! Setup the X variable
      real(pr), intent(in) :: kx(nc)
      real(pr), intent(in) :: kw(nc)
      real(pr), intent(in) :: beta
      real(pr), intent(in) :: p
      real(pr), intent(in) :: t

      real(pr) :: X(2*nc)

      X(:nc) = log(kx)
      X(nc+1:2*nc) = log(kw)
      X(2*nc + 1) = beta
      X(2*nc + 2) = log(p)
      X(2*nc + 3) = log(t)
   end function
   ! ---------------------------------------------------------------------------
   
   ! ===========================================================================
   ! Two Phase envelopes
   ! ---------------------------------------------------------------------------
   function X2(kfact, P, T) result(X)
      real(pr), intent(in) :: kfact(nc)
      real(pr), intent(in) :: P
      real(pr), intent(in) :: T

      real(pr) :: X(nc + 2)

      integer :: n

      n = size(kfact)

      X(:n) = log(kfact)
      X(n + 1) = log(T)
      X(n + 2) = log(P)
   end function

   subroutine F2(incipient, z, y, X, S, ns, F, dF)
      character(len=:), allocatable, intent(in) :: incipient
      real(pr), intent(in) :: z(:)
      real(pr), intent(in) :: X(nc + 2)
      real(pr), intent(in) :: y(nc)
      real(pr), intent(in) :: S
      integer, intent(in) :: ns

      real(pr), intent(out) :: F(nc + 2)
      real(pr), intent(out) :: dF(nc + 2, nc + 2)

      real(pr) :: Vx, Vy, lnfug_x(nc), lnfug_y(nc)
      real(pr) :: dlnphi_dt_x(nc), dlnphi_dt_y(nc)
      real(pr) :: dlnphi_dp_x(nc), dlnphi_dp_y(nc)
      real(pr) :: dlnphi_dn_x(nc, nc), dlnphi_dn_y(nc, nc)

      real(pr) :: T, P

      integer :: ix, iy, n, j

      n = size(z)
      F = 0
      dF = 0

      T = exp(X(n+1))
      P = exp(X(n+2))

      select case(incipient)
      case ("liquid")
         ix = -1
         iy = 1
      case ("vapor")
         ix = 1
         iy = -1
      case ("2ndliquid")
         ix = 1
         iy = 1
      end select

      call TERMO(n, iy, 4, T, P, y, Vy, lnfug_y, dlnphi_dp_y, dlnphi_dt_y, dlnphi_dn_y)
      call TERMO(n, ix, 2, T, P, z, Vx, lnfug_x, dlnphi_dp_x, dlnphi_dt_x, dlnphi_dn_x)

      F(:n) = X(:n) + lnfug_y - lnfug_x  ! X(:n) are LOG_K
      F(n + 1) = sum(y - z)
      F(n + 2) = X(ns) - S

      ! Jacobian Matrix
      do j=1,n
         df(:n, j) = dlnphi_dn_y(:, j) * y(j)
         df(j, j) = dF(j, j) + 1
      end do

      df(:n, n + 1) = T * (dlnphi_dt_y - dlnphi_dt_x)
      df(:n, n + 2) = P * (dlnphi_dp_y - dlnphi_dp_x)

      df(n + 1, :n) = y

      df(n + 2, :) = 0
      df(n + 2, ns) = 1
   end subroutine F2

   subroutine fix_delx(&
         point, iterations, desired_iterations, first_tol, tol, delX &
      )
      integer, intent(in)  :: point
      integer, intent(in)  :: iterations
      integer, intent(in)  :: desired_iterations
      real(pr), intent(in) :: first_tol
      real(pr), intent(in) :: tol
      real(pr), intent(in out) :: delX(:)

      if (point == 1) then
         do while (maxval(abs(delX)) > first_tol) 
            ! Too large Newton step --> Reduce it
            delX = delX/2
         end do
      else
         do while (maxval(abs(delX)) > tol)   
            ! Too large Newton step --> Reduce it
            delX = delX/2
         end do
         if (iterations > desired_iterations)  then
            ! too many iterations (sometimes due to oscillatory behavior 
            ! near critical point) --> Reduce it
            delX = delX/2
         endif
      end if
   end subroutine
   ! ===========================================================================

   ! =============================================================================
   !  Crossing related
   ! -----------------------------------------------------------------------------
   subroutine find_crossings(&
         dew, bub, hpl, &
         Tcr1, Pcr1, Tcr2, Pcr2, &
         kfcr1, kscr1, kfcr2, kscr2 &
      )
      ! Find the crossings between the whole set of two-phase lines
      use dtypes, only: envelope, kfcross, point, find_cross, find_self_cross
      implicit none

      type(envelope),        intent(in out) :: dew !! Dew envelope (AOP)
      type(envelope),        intent(in out) :: bub !! Bubble envelope
      type(envelope),        intent(in out) :: hpl !! HPLL envelope
      real(pr),                 intent(out) :: Tcr1, Pcr1, Tcr2, Pcr2
      real(pr),                 intent(out) :: kfcr1(:), kscr1(:), kfcr2(:), kscr2(:)

      type(point), allocatable :: self_cross(:)
      type(point), allocatable :: dew_bub_cross(:)
      type(point), allocatable :: dew_hpl_cross(:)
      type(point), allocatable :: bub_hpl_cross(:)

      logical :: has_hpll_line

      logical :: crossed_dew_hpl
      logical :: crossed_bub_hpl
      logical :: crossed_dew_dew
      logical :: crossed_dew_bub
      logical :: crossed_self

      tcr1 = 0
      tcr2 = 0
      pcr1 = 0
      pcr2 = 0

      has_hpll_line = allocated(hpl%t)

      ! ========================================================================
      !  Find the crossings
      ! ------------------------------------------------------------------------

      ! First check if the dew_envelope or the low_t_envelope self_cross
      call find_self_cross(dew%t, dew%p, self_cross, crossed_self)

      ! Then:
      ! - Check if HPLL line has been traced
      ! - [Cross dew bub?]|[Cross bub HPLL andor Cross dew HPLL]
      if (has_hpll_line) then
         call find_cross(dew%t, hpl%t, dew%p, hpl%p, dew_hpl_cross, crossed_dew_hpl)
         call find_cross(bub%t, hpl%t, bub%p, hpl%p, bub_hpl_cross, crossed_bub_hpl)
         call find_cross(dew%t, bub%t, dew%p, bub%p, dew_bub_cross, crossed_dew_bub)
      else
         call find_cross(dew%t, bub%t, dew%p, bub%p, dew_bub_cross, crossed_dew_bub)
      end if
      ! ========================================================================

      if (has_hpll_line) then
         if (crossed_bub_hpl .and. crossed_dew_bub) then
            ! HPLL line crossed with bubble, and bubble crossed with dew
            call get_values(bub_hpl_cross, 1, bub, hpl, tcr1, pcr1, kfcr1, kscr1)
            call get_values(dew_bub_cross, 1, bub, dew, tcr2, pcr2, kfcr2, kscr2)
         else if (crossed_dew_hpl) then
            ! HPLL line crossed with dew line
            call get_values(dew_hpl_cross, 1, dew, hpl, tcr2, pcr2, kfcr2, kscr2)
         end if
      else
         if (crossed_dew_bub) then
            call get_values(dew_bub_cross, 1, bub, dew, tcr2, pcr2, kfcr2, kscr2)
            call get_values(dew_bub_cross, 2, bub, dew, tcr1, pcr1, kfcr1, kscr1)
         end if
      end if
      
      if (crossed_self) then
         call get_values(self_cross, 1, dew, dew, tcr2, pcr2, kfcr2, kscr2)
      end if
   contains
       subroutine get_values(cross, index, envelope1, envelope2, t, p, kf, ks)
          use dtypes, only: point
          type(point), allocatable, intent(in) :: cross(:)
          type(envelope), intent(in) :: envelope1, envelope2
          integer, intent(in) :: index
          real(pr), intent(out) :: t
          real(pr), intent(out) :: p
          real(pr), intent(out) :: kf(size(kfcr1))
          real(pr), intent(out) :: ks(size(kscr1))

          integer :: icross, jcross

          icross = cross(index)%i
          jcross = cross(index)%j

          t = cross(index)%x
          p = cross(index)%y

          kf = kfcross(jcross, envelope1%t, envelope1%logk, t)
          ks = kfcross(icross, envelope2%t, envelope2%logk, t)
       end subroutine
   end subroutine find_crossings
   ! ===========================================================================
end module envelopes