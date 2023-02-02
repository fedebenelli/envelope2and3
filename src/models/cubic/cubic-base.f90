module cubic_base
   use constants, only: pr
   use helmholtz, only: ArModel, size
   implicit none


   type, extends(ArModel) :: CubicEOS
      real(pr), allocatable :: ac(:)
      real(pr), allocatable :: b(:)
      real(pr), allocatable :: k(:)
      real(pr), allocatable :: del1(:)
      real(pr), allocatable :: del2(:)

      real(pr), allocatable :: pc(:)
      real(pr), allocatable :: tc(:)
      real(pr), allocatable :: zc(:)
      real(pr), allocatable :: vc(:)
      real(pr), allocatable :: w(:)
   end type


   interface residual_helmholtz
      module procedure :: residual_helmholtz_cubic
   end interface residual_helmholtz

   interface atractive_parameter
        module procedure :: atractive_generic_cubic
   end interface atractive_parameter

   interface repulsive_parameter
        module procedure :: repulsive_generic_cubic
   end interface

contains

   pure function atractive_generic_cubic(system, p, v, t) result(a)
         class(CubicEOS), intent(in) :: system
         real(pr), intent(in) :: p
         real(pr), intent(in) :: v
         real(pr), intent(in) :: t

         real(pr) :: a(size(system))

         real(pr) :: tr(size(system))

         associate (ac => system%ac, Tc => system%Tc, k => system%k)
            Tr = T/Tc
            a = ac*(1 + k*(1 - sqrt(Tr)))**2
         end associate
   end function

   pure function repulsive_generic_cubic(system, p, v, t) result(b)
         class(CubicEOS), intent(in) :: system
         real(pr), intent(in) :: p
         real(pr), intent(in) :: v
         real(pr), intent(in) :: t

         real(pr) :: b(size(system))

         b = system%b
   end function

   subroutine residual_helmholtz_cubic(system, v, t, ar)
      type(CubicEoS) :: system
      real(pr), intent(in) :: v
      real(pr), intent(in) :: t
      real(pr), intent(out) :: ar
   end subroutine

   subroutine get_Zc_OMa_OMb(system, OMa, OMb)
      !! Calculate Zc, OMa and OMb from the delta_1 parameter.
      class(CubicEOS), intent(in out) :: system
      real(pr), intent(out) :: OMa(size(system)) !! OMa
      real(pr), intent(out) :: OMb(size(system)) !! OMb

      real(pr) :: d1(size(system)), y(size(system))

      associate(Zc => system%Zc, del1 => system%del1)
      d1 = (1._pr + del1**2._pr)/(1._pr + del1)
      y = 1._pr + (2._pr*(1._pr + del1))**(1.0_pr/3._pr) + (4._pr/(1._pr + del1))**(1.0_pr/3)
      OMa = (3._pr*y*y + 3._pr*y*d1 + d1**2._pr + d1 - 1.0_pr)/(3._pr*y + d1 - 1.0_pr)**2._pr
      OMb = 1._pr/(3._pr*y + d1 - 1.0_pr)
      Zc = y/(3._pr*y + d1 - 1.0_pr)
      end associate
   end subroutine get_Zc_OMa_OMb
end module cubic_base
