module helmholtz
   use constants
   type, abstract :: ArModel
   contains
      procedure(residual_helmholtz), deferred :: residual_helmholtz
   end type

   type, extends(ArModel) :: CubicEOS
      real(pr), allocatable :: ac(:)
      real(pr), allocatable :: b(:)
      real(pr), allocatable :: k(:)
      real(pr), allocatable :: del1(:)
      real(pr), allocatable :: del2(:)

      real(pr), allocatable :: pc(:)
      real(pr), allocatable :: tc(:)
      real(pr), allocatable :: w(:)
   contains
      procedure :: residual_helmholtz => dummy_ar
   end type

   abstract interface
      pure subroutine residual_helmholtz(self, t, v, ar)
         use constants
         import ArModel
         class(ArModel), intent(in out) :: self
         real(pr), intent(in) :: t
         real(pr), intent(in) :: v
         real(pr), intent(out) :: ar
      end subroutine
   end interface

contains

   pure subroutine dummy_ar(self, t, v, ar)
      class(CubicEOS), intent(in out) :: self
      real(pr), intent(in) :: t
      real(pr), intent(in) :: v
      real(pr), intent(out) :: ar
   end subroutine
end module helmholtz

module peng_robinson
   use constants
   use helmholtz

   type, extends(CubicEOS) :: PengRobinson
   end type

   interface setup_cubic
      module procedure :: setup_pengrobinson
   end interface
contains

   pure subroutine setup_pengrobinson(self, n, tc, pc, w)
      class(PengRobinson), intent(in out) :: self
      integer, intent(in) :: n
      real(pr), intent(in) :: tc(n)
      real(pr), intent(in) :: pc(n)
      real(pr), intent(in) :: w(n)

      self%tc = tc
      self%pc = pc
      self%w = w
   end subroutine

end module peng_robinson
